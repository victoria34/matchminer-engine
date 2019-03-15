"""Copyright 2016 Dana-Farber Cancer Institute"""

import re
import os
import sys
import yaml
import json
import logging
import requests
import pandas as pd
import datetime as dt
from pymongo import MongoClient

import oncotreenx
from matchengine.settings import months, TUMOR_TREE, mmr_map, mmr_map_rev


def build_gquery(field, txt):
    """Builds the Mongo query from the genomic criteria"""

    # unless instructed otherwise, construct a positive query
    neg = False
    sv = False

    # Structural variants
    if field.lower() == 'variant_category' and txt in ['SV', '!SV']:
        sv = True

    # MMR and MS Status
    if field.lower() == 'mmr_status' or field.lower() == 'ms_status':
        key = '$eq'
        txt = mmr_map[txt]

    # Wildcard Protein Change
    elif field.lower() == 'wildcard_protein_change':

        # Negative queries will be run as positive queries and the matched sample ids will be subtracted from
        # the set of all sample ids in the database
        if txt.startswith('!'):
            neg = True
            txt = txt[1:]

        # By convention, all protein changes being with "p."
        if not txt.startswith('p.'):
            txt = 'p.' + txt

        key = '$regex'
        txt = '^%s[A-Z]' % txt

    # Match any variant category
    elif field.lower() == 'variant_category' and txt.lower() == 'any variation':
        txt = ['MUTATION', 'CNV']
        key = '$in'

    # Not equal to given value
    elif (isinstance(txt, str) and txt.startswith('!')) or (isinstance(txt, unicode) and txt.startswith('!')):
        key = '$eq'
        neg = True

        # Exon field is queried as an integer. All other fields as strings
        if field.lower() == 'exon':
            txt = int(txt.replace('!', ''))
        else:
            txt = txt.replace('!', '')

    # Otherwise set equal to
    else:
        key = '$eq'

    return key, txt, neg, sv


def build_cquery(c, norm_field, txt):
    """Builds the Mongo query from the clinical criteria"""

    if isinstance(txt, list):
        c[norm_field] = {}
        for i in txt:
            if i.startswith('!'):
                if '$nin' in c[norm_field]:
                    c[norm_field]['$nin'].append(i.replace('!', ''))
                else:
                    c[norm_field]['$nin'] = [i.replace('!', '')]
            else:
                if '$in' in c[norm_field]:
                    c[norm_field]['$in'].append(i)
                else:
                    c[norm_field]['$in'] = [i]

    else:
        if txt.startswith('!'):
            key = '$ne'
            txt = txt.replace('!', '')
        else:
            key = '$eq'

        c[norm_field] = {key: txt}

    return c


def build_oncotree():
    """Builds oncotree"""
    return oncotreenx.build_oncotree(file_path=TUMOR_TREE)


def normalize_fields(mapping, field):
    """Translates yaml field name into the database field name."""

    # parse map from db
    old_keys = [i['key_old'] for i in mapping]
    new_keys = [i['key_new'] for i in mapping]
    vals = [i['values'] for i in mapping]
    val_map = dict(zip(new_keys, vals))

    # translate keys
    field = field.upper()
    if field in old_keys:
        field = new_keys[old_keys.index(field)]
        return field, val_map[field]
    else:
        return field, val_map[field]


def normalize_values(mapping, field, val):
    """Translates yaml fields and values into database fields and values"""

    # first normalize keys
    field, mapping = normalize_fields(mapping, field)

    # exclude "!" from mapping
    ne = False
    map_by = val
    if (isinstance(val, str) and val[0] == "!") or (isinstance(val, unicode) and val[0] == "!"):
        map_by = val[1:]
        ne = True

    # return the translated keys and values
    if map_by in mapping:
        if ne:
            return field, '!%s' % str(mapping[map_by])
        else:
            return field, mapping[map_by]
    else:
        return field, val


def samples_from_mrns(db, mrns):
    """Returns a dictionary mapping each MRN to all of its associated SAMPLE_IDs"""

    clinical = list(db.clinical.find(
        {'MRN': {'$in': mrns}},
        {'MRN': 1, 'SAMPLE_ID': 1}
    ))

    mrn_map = {}
    for c in clinical:
        mrn_map[c['SAMPLE_ID']] = c['MRN']

    return mrn_map


def search_birth_date(c):
    """Converts query to filter by birth date based on the given age"""
    txt = c['BIRTH_DATE']['$eq'].lstrip()

    # translate to mongo query
    if txt.startswith('>='):
        key = '$lte'
    elif txt.startswith('<='):
        key = '$gte'
    elif txt.startswith('>'):
        key = '$lt'
    elif txt.startswith('<'):
        key = '$gt'
    else:
        raise ValueError

    # get age
    idx = 1
    if txt[1] == '=':
        idx = 2
    abs_age = str(txt[idx:])

    # date today
    today = dt.datetime.today()

    # calculate date to query
    if '.' in abs_age:
        month, year = get_months(abs_age, today)
        query_date = today.replace(month=month, year=(today.year + year))
    else:
        try:
            query_date = today.replace(year=today.year - int(abs_age))
        except ValueError:
            query_date = today + (dt.datetime(today.year + int(abs_age), 1, 1) - dt.datetime(today.year, 1, 1))

    return {key: query_date}


def get_months(abs_age, today):
    """Given a decimal, returns the number of months and number of years to subtract from today"""

    split_age = str(abs_age).split('.')

    # month
    month = split_age[1]
    month = int((float(month) * 12) / (10 ** len(month)))  # e.g. convert 5/10 to x/12

    # year
    if split_age[0] == '':
        year = 0
    else:
        year = int(split_age[0])

    # handle crossing over a year boundary
    if today.month - month <= 0:
        month = months.index(months[-(abs(today.month - month))])
        year = -(year + 1)
    else:
        month = today.month - month

    if month == 0:
        month = 1

    return month, year


def add_trials(trial_path, db):
    """Adds all ymls in the "trial_path" to the db"""

    inserted_ids = 0

    # search directory for ymls
    for yml in os.listdir(trial_path):
        ymlpath = os.path.join(trial_path, yml)

        # only add files of extension ".yml"
        if yml.split('.')[-1] != 'yml':
            continue

        # convert yml to json format
        with open(ymlpath) as f:
            t = yaml.load(f.read())

            # add trial to db
            result = db.trial.insert_one(t)
            if result.inserted_id:
                inserted_ids += 1

    return inserted_ids


def format_genomic_alteration(g, query):
    """Format the genomic alteration that matched a particular trial"""

    if g is None:
        return g

    # for clarity
    gene = 'TRUE_HUGO_SYMBOL'
    mut = 'TRUE_PROTEIN_CHANGE'
    cnv = 'CNV_CALL'
    var = 'TRUE_VARIANT_CLASSIFICATION'
    sv = 'VARIANT_CATEGORY'
    wt = 'WILDTYPE'
    mmr = 'MMR_STATUS'

    alteration = ''
    is_variant = 'gene'

    # Ignore wildtype when determining if match was gene- or variant-level
    if query.keys()[0] == '$and':
        query = query['$and'][0]

    # determine if match was gene- or variant-level
    if mut in query and query[mut] is not None:
        is_variant = 'variant'

    # add wildtype calls
    if wt in g and g[wt] is True:
        alteration += 'wt '

    # add gene
    if gene in g and g[gene] is not None:
        alteration += g[gene]

    # add mutation
    if mut in g and g[mut] is not None:
        alteration += ' %s' % g[mut]

    # add cnv call
    elif cnv in g and g[cnv] is not None:
        alteration += ' %s' % g[cnv]

    # add variant classification
    elif var in g and g[var] is not None:
        alteration += ' %s' % g[var]

    # add structural variation
    elif sv in g and g[sv] == 'SV':
        alteration += ' Structural Variation'

    # add mutational signtature
    elif sv in g and g[sv] == 'SIGNATURE' and mmr in g and g[mmr] is not None and g[mmr] in mmr_map_rev:
        alteration += mmr_map_rev[g[mmr]]

    return alteration, is_variant


def format_not_match(g):
    """Format the genomic alteration for genomic documents that matched a negative clause of a match tree"""

    alteration = ''
    is_variant = 'gene'

    # for clarity
    gene = 'TRUE_HUGO_SYMBOL'
    mut = 'TRUE_PROTEIN_CHANGE'
    cnv = 'CNV_CALL'
    var = 'TRUE_VARIANT_CLASSIFICATION'
    sv = 'VARIANT_CATEGORY'

    # Ignore wildtype when formatting genomic alteration
    if g.keys()[0] == '$and':
        g = g['$and'][0]

    # add gene
    if gene in g and g[gene] is not None:
        alteration = format_query(g[gene], gene=True)

    # add mutation
    if mut in g and g[mut] is not None:
        alteration += ' %s' % format_query(g[mut])
        is_variant = 'variant'

    # add cnv call
    elif cnv in g and g[cnv] is not None:
        alteration += ' %s' % format_query(g[cnv])

    # add variant classification
    elif var in g and g[var] is not None:
        alteration += ' %s' % format_query(g[var])

    # add structural variation
    elif sv in g and g[sv][g[sv].keys()[0]] == 'SV':
        alteration += ' Structural Variation'

    # if no gene is specified, the ! is added manually
    if gene not in g or g[gene] is None:
        alteration = '!' + alteration[1:]

    return alteration, is_variant


def format_query(g, gene=False):
    """Turns the mongo query into a formatted genomic alteration"""

    alteration = ''
    key = g.keys()[0]

    if key == '$regex':
        alteration += '!%s' % g[key].replace('^', '').replace('[A-Z]', '')
    elif key == '$in':
        for item in g[key][:-1]:
            alteration += '!%s, ' % item
        alteration += '!%s' % g[key][-1]
    else:
        alteration += '!%s' % g[key]

    if not gene:
        alteration = alteration.replace('!', '')

    return alteration


def add_matches(trial_matches_df, db):
    """Add the match table to the database or update what already exists theres"""

    if 'clinical_id' in trial_matches_df.columns:
        trial_matches_df['clinical_id'] = trial_matches_df['clinical_id'].apply(lambda x: str(x))

    if 'genomic_id' in trial_matches_df.columns:
        trial_matches_df['genomic_id'] = trial_matches_df['genomic_id'].apply(lambda x: str(x))

    if 'report_date' in trial_matches_df.columns:
        trial_matches_df['report_date'] = trial_matches_df['report_date'].apply(
            lambda x: dt.datetime.strftime(x, '%Y-%m-%d %X') if pd.notnull(x) else x)

    if len(trial_matches_df.index) > 0:
        db.trial_match.drop()
        for i in range(0, trial_matches_df.shape[0], 1000):
            records = json.loads(trial_matches_df[i:i + 1000].T.to_json()).values()
            db.trial_match.insert_many(records)


def get_db(uri):
    """Returns a Mongo connection"""

    if uri:
        MONGO_URI = uri
    else:

        # sanity check
        MONGO_URI = ""
        file_path = os.getenv("SECRETS_JSON", None)
        if file_path is None:
            uri = os.getenv("MONGO_URI")
            if uri:
                MONGO_URI = uri
            else:
                logging.error("ENVAR SECRETS_JSON not set")
                sys.exit(1)

        else:
            # pull values.
            with open(file_path) as fin:
                vars = json.load(fin)
                for name, value in vars.iteritems():
                    if name == "MONGO_URI":
                        MONGO_URI = value

    if not MONGO_URI:
        logging.error("MONGO_URI not set in SECRETS_JSON")
    else:
        os.environ["MONGO_URI"] = MONGO_URI
        connection = MongoClient(MONGO_URI)
        if 'mlab' in MONGO_URI:
            # heroku will assign a auto-generated mongodb that cannot be renamed
            return connection.get_default_database()
        else:
            return connection["matchminer"]


def get_structural_variants(g):
    """
    Performs a string search for the structural variant.

    :param g: Genomic query in
    :return: Genomic query out
    """

    # get the genes.
    hugo = g['TRUE_HUGO_SYMBOL']
    k = hugo.keys()[0]
    genes = hugo[k]

    if not isinstance(genes, list):
        genes = [genes]

    # TODO add synonyms

    # encode as full search criteria.
    sv_clauses = []
    for gene in genes:
        abc = "(.*\W{0}\W.*)|(^{0}\W.*)|(.*\W{0}$)".format(gene)
        sv_clauses.append(re.compile(abc, re.IGNORECASE))

    # add it to filter and remove gene criteria.
    del g['TRUE_HUGO_SYMBOL']
    g['STRUCTURAL_VARIANT_COMMENT'] = {"$in": sv_clauses}

    return g


def clean_query_for_msi(g):
    if 'MMR_STATUS' in g and 'TRUE_HUGO_SYMBOL' in g:
        del g['TRUE_HUGO_SYMBOL']
    return g


def get_cancer_type_match(trial):
    """
    Determines if the trial has criteria to match all solid or all liquid tumors in it.

    :param trial: Entire trial object
    :return: cancer_type_match
    """

    if '_summary' not in trial or 'tumor_types' not in trial['_summary']:
        return 'unknown'

    if '_SOLID_' in trial['_summary']['tumor_types']:
        return 'all_solid'
    elif '_LIQUID_' in trial['_summary']['tumor_types']:
        return 'all_liquid'
    else:
        return 'specific'


def get_coordinating_center(trial):
    """
    Returns the trials' coordinating center

    :param trial: Entire trial object
    """

    if '_summary' not in trial or 'coordinating_center' not in trial['_summary']:
        return 'unknown'
    else:
        return trial['_summary']['coordinating_center']


def oncokb_api_match(db, collection_name):
    """
    Get all trial matched results from OncoKB API
    :param db: mongo database
    :param collection_name: genomic or new_genomic
    :return: matched result object
    matched result object data structure:
    {
        gene: {
            proteinChange: [variants]
        },
        ......
    }
    matched result object example:
    {
      'TP53': {
        'H214L': [
          'Oncogenic Mutations'
        ]
      },
      'BRAF': {
        'V600E': [
          'Oncogenic Mutations',
          'V600',
          'V600E'
        ]
      }
    }
    """

    queries = list()
    annotated_variants = list()
    matched_results = {}
    api_url = 'http://oncokb.org/api/private/utils/match/variant'

    # get genomic info from collection genomic or new_genomic
    genomic_proj = {
        'SAMPLE_ID': 1,
        'TRUE_HUGO_SYMBOL': 1,
        'TRUE_PROTEIN_CHANGE':1
    }
    genomic_collection = db[collection_name]
    genomic_results = list(genomic_collection.find({}, genomic_proj))
    for genomic in genomic_results:
        if 'TRUE_HUGO_SYMBOL' in genomic and genomic['TRUE_HUGO_SYMBOL'] and 'TRUE_PROTEIN_CHANGE' in genomic and genomic['TRUE_PROTEIN_CHANGE']:
            query = {
                "id": genomic['SAMPLE_ID'],
                "hugoSymbol": genomic['TRUE_HUGO_SYMBOL'],
                "alteration": genomic['TRUE_PROTEIN_CHANGE']
            }
            queries.append(query)

    # get genomic node info from collection trial
    steps = list(db.trial.find({'treatment_list.step':{'$exists': 'true', '$ne': []}}, {'_id': 0}));
    for step in steps:
        for arm_match in step['treatment_list']['step']:
            if 'arm' in arm_match and arm_match['arm']:
                for arm in arm_match['arm']:
                    if 'match' in arm and arm['match']:
                        for match in arm['match']:
                            find_genomic_node(match, annotated_variants)
                    if 'dose_level' in arm:
                        for dose in arm['dose_level']:
                            if 'match' in dose and dose['match']:
                                for match in dose['match']:
                                    find_genomic_node(match, annotated_variants)
            if 'match' in arm_match and arm_match['match']:
                for match in arm_match['match']:
                    find_genomic_node(match, annotated_variants)

    body = {
        "oncokbVariants": annotated_variants,
        "queries": queries
    }
    body = json.dumps(body)
    headers = {'Content-type': 'application/json'}
    response = requests.post(api_url, data=body, headers=headers)
    result = response.json()
    for trial_match in result:
        if trial_match['result']:
            protein_change = trial_match['query']['alteration']
            for genomic_alteration in trial_match['result']:
                if genomic_alteration['hugoSymbol'] in matched_results:
                    if protein_change in matched_results[genomic_alteration['hugoSymbol']]:
                        if not (genomic_alteration['alteration'] in matched_results[genomic_alteration['hugoSymbol']][protein_change]):
                            matched_results[genomic_alteration['hugoSymbol']][protein_change].append(genomic_alteration['alteration'])
                    else:
                        matched_results[genomic_alteration['hugoSymbol']][protein_change] = [genomic_alteration['alteration']]
                else:
                    matched_results[genomic_alteration['hugoSymbol']] = {}
                    matched_results[genomic_alteration['hugoSymbol']][protein_change] = [genomic_alteration['alteration']]

    return matched_results


def find_genomic_node(match, node_infos):
    """Find all genomic nodes under 'match' object """

    if 'genomic' in match and match['genomic'] and \
                    'hugo_symbol' in match['genomic'] and match['genomic']['hugo_symbol'] and \
                    'annotated_variant' in match['genomic'] and match['genomic']['annotated_variant']:
        node_infos.append(get_hugo_variant_info(match['genomic']))
    if 'and' in match and match['and']:
        for and_node in match['and']:
            find_genomic_node(and_node, node_infos)
    if 'or' in match and match['or']:
        for or_node in match['or']:
            find_genomic_node(or_node, node_infos)
    return node_infos


def get_hugo_variant_info(genomic_node):
    """Get hugo_symbol and annotated_variant from a genomic node of trials"""

    annotated_variant = {
        "alteration": genomic_node['annotated_variant'],
        "hugoSymbol": genomic_node['hugo_symbol']
    }
    return annotated_variant


def process_cmd(type, uri, file, collection = None, upsert = None, is_json_array = False):
    """
    Generate mongo command line for loading data
    :param type: command line type 'mongorestore' or 'mongoimport'
    :param uri: mongo uri
    :param collection: collection name
    :param file: data file
    :param upsert: store attributes related to upsert
           upsert = {
                 is_upsert: True/False,
                 fields: index used to identify data record
            }
           --upsert: Replace existing documents in the database with matching documents from the import file.
           --upsertFields: Specifies a list of fields for the query portion of the upsert.
    :return: command line string
    """
    cmd = ''
    if 'mlab' in uri:
        user, password, address, dbname = process_mlab_uri(uri)
        cmd = '%s -h %s -d %s -u %s -p %s' % (type, address, dbname, user, password)
        if type == 'mongorestore':
            cmd += file
        elif type == 'mongoimport':
            cmd += ' -c %s --file %s' % (collection, file)
    else:
        cmd = '%s --host localhost:27017 --db matchminer' % type
        if type == 'mongorestore':
            cmd += file
        elif type == 'mongoimport':
            cmd += ' --collection %s --file %s' % (collection, file)

    if not (upsert is None) and type == 'mongoimport':
        if upsert['is_upsert']:
            upsert_fields = ', '.join(str(x) for x in upsert['fields'])
            cmd += ' --upsert --upsertFields %s' % upsert_fields
    if is_json_array:
        cmd += ' --jsonArray'

    return cmd


def process_mlab_uri(uri):
    user_pass = ((uri.split('//', 1)[-1]).split('@', 1)[0]).split(':', 1)
    user = user_pass[0]
    password = user_pass[1]
    address = (uri.split('@', 1)[-1]).split('/', 1)[0]
    dbname = uri.split('/')[-1]
    return user, password, address, dbname


def set_match_method(method):
    global match_method
    match_method = method


def get_match_method():
    return match_method
