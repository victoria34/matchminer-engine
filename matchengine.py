"""Copyright 2016 Dana-Farber Cancer Institute"""

import os
import sys
import json
import time
import yaml
import logging
import argparse
import subprocess
import pandas as pd
import datetime as dt
from pymongo import ASCENDING

from matchengine.engine import MatchEngine
from matchengine.utilities import get_db

MONGO_URI = ""
MONGO_DBNAME = "matchminer"
MATCH_FIELDS = "mrn,sample_id,first_last,protocol_no,nct_id,genomic_alteration,tier,match_type," \
               "trial_accrual_status,match_level,code,internal_id,ord_physician_name,ord_physician_email," \
               "vital_status,oncotree_primary_diagnosis_name,true_hugo_symbol,true_protein_change," \
               "true_variant_classification,variant_category,report_date,chromosome,position," \
               "true_cdna_change,reference_allele,true_transcript_exon,canonical_strand,allele_fraction," \
               "cnv_call,wildtype,_id,oncokb_variant"


class Trial:

    def __init__(self, uri, db):

        self.db = db
        self.load_dict = {
            'yml': self.yaml_to_mongo,
            'bson': self.bson_to_mongo,
            'json': self.json_to_mongo
        }
        user_pass = ((uri.split('//', 1)[-1]).split('@', 1)[0]).split(':', 1)
        self.user = user_pass[0]
        self.password = user_pass[1]
        self.address = (uri.split('@', 1)[-1]).split('/', 1)[0]

    def yaml_to_mongo(self, yml):
        """
        If you specify the path to a directory, all files with extension YML will be added to MongoDB.
        If you specify the path to a specific YML file, it will add that file to MongoDB.

        :param yml: Path to YML file.
        """

        # search directory for ymls
        if os.path.isdir(yml):
            for y in os.listdir(yml):
                ymlpath = os.path.join(yml, y)

                # only add files of extension ".yml"
                if ymlpath.split('.')[-1] != 'yml':
                    continue

                # convert yml to json format
                add_trial(ymlpath, self.db)
        else:
            add_trial(yml, self.db)

    def bson_to_mongo(self, bson):
        """
        If you specify the path to a directory, all files with extension BSON will be added to MongoDB.
        If you specify the path to a specific BSON file, it will add that file to MongoDB.

        :param bson: Path to BSON file.
        """
        cmd = "mongorestore -h %s -d %s -u %s -p %s %s" % (self.address, self.db.name, self.user, self.password, bson)
        subprocess.call(cmd.split(' '))

    def json_to_mongo(self, json):
        """
        If you specify the path to a directory, all files with extension JSON will be added to MongoDB.
        If you specify the path to a specific JSON file, it will add that file to MongoDB.

        :param json: Path to JSON file.
        """
        # --upsert: Replace existing documents in the database with matching documents from the import file.
        # --upsertFields: Specifies a list of fields for the query portion of the upsert.
        cmd = "mongoimport -h %s -d %s -c trial --upsert -u %s -p %s --upsertFields nct_id --file %s" % (self.address, self.db.name, self.user, self.password, json)
        subprocess.call(cmd.split(' '))


class Patient:

    def __init__(self, uri, db):

        self.db = db
        self.load_dict = {
            'csv': self.load_csv,
            'pkl': self.load_pkl,
            'bson': self.load_bson,
            'json': self.load_json
        }
        self.clinical_df = None
        self.genomic_df = None
        user_pass = ((uri.split('//', 1)[-1]).split('@', 1)[0]).split(':', 1)
        self.user = user_pass[0]
        self.password = user_pass[1]
        self.address = (uri.split('@', 1)[-1]).split('/', 1)[0]

    def load_csv(self, clinical, genomic):
        """Load CSV file into a Pandas dataframe"""
        self.clinical_df = pd.read_csv(clinical)
        self.genomic_df = pd.read_csv(genomic, low_memory=False)

    def load_pkl(self, clinical, genomic):
        """Load PKL file into a Pandas dataframe"""
        self.clinical_df = pd.read_pickle(clinical)
        self.genomic_df = pd.read_pickle(genomic)

    def load_bson(self, clinical, genomic):
        """Load bson file into MongoDB"""
        cmd1 = "mongorestore -h %s -d %s -u %s -p %s %s" % (self.address, self.db.name, self.user, self.password, clinical)
        cmd2 = "mongorestore -h %s -d %s -u %s -p %s %s" % (self.address, self.db.name, self.user, self.password, genomic)
        subprocess.call(cmd1.split(' '))
        subprocess.call(cmd2.split(' '))
        return True

    def load_json(self, clinical, genomic):
        """
        If you specify the path to a directory, all files with extension JSON will be added to MongoDB.
        If you specify the path to a specific JSON file, it will add that file to MongoDB.

        :param json: Path to JSON file.

        Note: For the empty fields in genomic data, their values should be null rather than "".
              For the false fields in genomic data, their values should be false rather than "false".
              For the date fields in clinical data, the type of their values should be date object rather than string.
        """

        cmd1 = "mongoimport -h %s -d %s -c clinical -u %s -p %s --file %s --upsert --upsertFields ONCOKB_CLINICAL_ID --stopOnError --jsonArray" % (self.address, self.db.name, self.user, self.password, clinical)
        cmd2 = "mongoimport -h %s -d %s -c genomic -u %s -p %s --file %s --upsert --upsertFields ONCOKB_GENOMIC_ID --stopOnError --jsonArray" % (self.address, self.db.name, self.user, self.password, genomic)
        subprocess.call(cmd1.split(' '))
        subprocess.call(cmd2.split(' '))

        self.db.new_clinical.drop()
        self.db.new_genomic.drop()

        # save patient data into new_genomic and new_collections for matching based on current sent data
        cmd3 = "mongoimport -h %s -d %s -c new_clinical -u %s -p %s --file %s --upsert --upsertFields ONCOKB_CLINICAL_ID --stopOnError --jsonArray" % (self.address, self.db.name, self.user, self.password, clinical)
        cmd4 = "mongoimport -h %s -d %s -c new_genomic -u %s -p %s --file %s --upsert --upsertFields ONCOKB_GENOMIC_ID --stopOnError --jsonArray" % (self.address, self.db.name, self.user, self.password, genomic)
        subprocess.call(cmd3.split(' '))
        subprocess.call(cmd4.split(' '))

        # convert string to date object
        for clinical_collection in [self.db.clinical.find(), self.db.new_clinical.find()]:
            for clinical_item in clinical_collection:
                cols = list()
                keys = clinical_item.keys()
                if 'BIRTH_DATE' in keys:
                    cols.append('BIRTH_DATE')
                if 'REPORT_DATE' in keys:
                    cols.append('REPORT_DATE')
                for col in cols:
                    if type(clinical_item[col]) is not dt.datetime:
                        clinical_item[col] = dt.datetime.strptime(str(clinical_item[col]), '%Y-%m-%d')
                        clinical_item[col] = dt.datetime.strptime(str(clinical_item[col]), '%Y-%m-%d %X')
                        self.db.new_clinical.update({'_id':clinical_item['_id']}, {"$set": {col: clinical_item[col]}}, upsert=False)

        return True

def load(args):
    """
    Sets up MongoDB for matching

    :param args: clinical: Path to csv file containing clinical data. Required fields are:
        - DFCI_MRN (Unique patient identifier)
        - SAMPLE_ID (Unique sample identifier)
        - ONCOTREE_PRIMARY_DIAGNOSIS_NAME (Disease diagnosis)
        - BIRTH_DATE (Date of birth in format 'YYYY-MM-DD 00:00:00.000')

        Suggested additional fields:
        - ORD_PHYSICIAN_NAME
        - ORD_PHYSICIAN_EMAIL
        - REPORT_DATE
        - VITAL_STATUS (alive or deceased)
        - FIRST_LAST (Patient's first and last name)
        - GENDER (Male or Female)

    :param args: genomic: Path to csv file containing genomic data. The following fields are used in matching:
        - SAMPLE_ID (Unique sample identifier)
        - TRUE_HUGO_SYMBOL (Gene name)
        - TRUE_PROTEIN_CHANGE (Specific variant)
        - TRUE_VARIANT_CLASSIFICATION (Variant type)
        - VARIANT_CATEGORY (CNV, MUTATION, or SV)
        - TRUE_TRANSCRIPT_EXON (Exon number <integer>
        - CNV_CALL (Heterozygous deletion, Homozygous deletion, Gain, High Level amplification, or null)
        - WILDTYPE (True or False)

        Suggested additional fields:
        - CHROMOSOME (Chromosome number in format 'chr01')
        - POSITION <integer>
        - TRUE_CDNA_CHANGE
        - REFERENCE_ALLELE
        - CANONICAL_STRAND (- or +)
        - ALLELE_FRACTION <float>
        - TIER <integer>

    :param args: trials: Path to bson trial file.
    """

    db = get_db(args.mongo_uri)
    t = Trial(args.mongo_uri, db)
    p = Patient(args.mongo_uri, db)

    # Add trials to mongo
    if args.trials:
        logging.info('Adding trials to mongo...')
        t.load_dict[args.trial_format](args.trials)

    # Add patient data to mongo
    if args.clinical and args.genomic:
        logging.info('Reading data into pandas...')
        is_bson_or_json = p.load_dict[args.patient_format](args.clinical, args.genomic)

        if not is_bson_or_json:

            # reformatting
            for col in ['BIRTH_DATE', 'REPORT_DATE']:
                try:
                    p.clinical_df[col] = p.clinical_df[col].apply(lambda x: str(dt.datetime.strptime(x, '%Y-%m-%d')))
                except ValueError as exc:
                    if col == 'BIRTH_DATE':
                        print '## WARNING ## Birth dates should be formatted %Y-%m-%d to be properly stored in MongoDB.'
                        print '##         ## Birth dates may be malformed in the database and will therefore not match'
                        print '##         ## trial age restrictions properly.'
                        print '##         ## System error: \n%s' % exc

            p.genomic_df['TRUE_TRANSCRIPT_EXON'] = p.genomic_df['TRUE_TRANSCRIPT_EXON'].apply(
                lambda x: int(x) if x != '' and pd.notnull(x) else x)

            # Add clinical data to mongo
            logging.info('Adding clinical data to mongo...')
            clinical_json = json.loads(p.clinical_df.T.to_json()).values()
            for item in clinical_json:
                for col in ['BIRTH_DATE', 'REPORT_DATE']:
                    if col in item:
                        item[col] = dt.datetime.strptime(str(item[col]), '%Y-%m-%d %X')

            db.clinical.insert(clinical_json)

            # Get clinical ids from mongo
            logging.info('Adding clinical ids to genomic data...')
            clinical_doc = list(db.clinical.find({}, {"_id": 1, "SAMPLE_ID": 1}))
            clinical_dict = dict(zip([i['SAMPLE_ID'] for i in clinical_doc], [i['_id'] for i in clinical_doc]))

            # pd -> json
            if args.trial_format == 'pkl':
                genomic_json = json.loads(p.genomic_df.to_json(orient='records'))
            else:
                genomic_json = json.loads(p.genomic_df.T.to_json()).values()

            # Map clinical ids to genomic data
            for item in genomic_json:
                if item['SAMPLE_ID'] in clinical_dict:
                    item["CLINICAL_ID"] = clinical_dict[item['SAMPLE_ID']]
                else:
                    item["CLINICAL_ID"] = None

            # Add genomic data to mongo
            logging.info('Adding genomic data to mongo...')
            db.genomic.insert(genomic_json)

        # Create index
        logging.info('Creating index...')
        db.genomic.create_index([("TRUE_HUGO_SYMBOL", ASCENDING), ("WILDTYPE", ASCENDING)])

    elif args.clinical and not args.genomic or args.genomic and not args.clinical:
        logging.error('If loading patient information, please provide both clinical and genomic data.')
        sys.exit(1)


def add_trial(yml, db):
    """
    Adds file in YAML format to MongoDB

    :param yml: Path to file
    :param db: MongoDB connection
    """

    with open(yml) as f:
        t = yaml.load(f.read())
        db.trial.insert_one(t)


def export_results(file_format, outpath, uri, db_name):
    """Return csv file containing the match results to the current working directory"""
    user_pass = ((uri.split('//', 1)[-1]).split('@', 1)[0]).split(':', 1)
    user = user_pass[0]
    password = user_pass[1]
    address = (uri.split('@', 1)[-1]).split('/', 1)[0]
    cmd = "mongoexport -h {3} -d {4} -u {5} -p {6} -c trial_match -f {0} --csv -o {2}.{1}".format(MATCH_FIELDS, file_format, outpath, address, db_name, user, password)
    subprocess.call(cmd.split(' '))


def match(args):
    """
    Matches all trials in database to patients

    :param daemon: Boolean flag; when true, runs the matchengine once per 24 hours.
    """

    db = get_db(args.mongo_uri)

    while True:
        me = MatchEngine(db)
        me.find_trial_matches()
        # delete queried data after matching
        db.drop_collection('new_clinical')
        db.drop_collection('new_genomic')

        # exit if it is not set to run as a nightly automated daemon, otherwise sleep for a day
        if not args.daemon:

            # choose output file format
            if args.json_format:
                file_format = 'json'
            elif args.outpath and len(args.outpath.split('.')) > 1:
                file_format = args.outpath.split('.')[-1]
                if file_format not in ['json', 'csv']:
                    file_format = 'csv'
            else:
                file_format = 'csv'

            # choose output path
            if args.outpath:
                outpath = args.outpath.split('.')[0]
            else:
                outpath = './results'

            # export results
            export_results(file_format, outpath, args.mongo_uri, db.name)

            break
        else:
            time.sleep(86400)   # sleep for 24 hours

if __name__ == '__main__':

    param_trials_help = 'Path to your trial data file or a directory containing a file for each trial.' \
                        'Default expected format is YML.'
    param_mongo_uri_help = 'Your MongoDB URI. If you do not supply one it will default to whatever is set to ' \
                           '"MONGO_URI" in your secrets file. ' \
                           'See https://docs.mongodb.com/manual/reference/connection-string/ for more information.'
    param_daemon_help = 'Set to launch the matchengine as a nightly automated process'
    param_clinical_help = 'Path to your clinical file. Default expected format is CSV.'
    param_genomic_help = 'Path to your genomic file. Default expected format is CSV'
    param_json_help = 'Set this flag to export your results in a .json file.'
    param_csv_help = 'Set this flag to export your results in a .csv file. Default.'
    param_outpath_help = 'Destination and name of your results file.'
    param_trial_format_help = 'File format of input trial data. Default is YML.'
    param_patient_format_help = 'File format of input patient data (both clinical and genomic files). Default is CSV.'

    # mode parser.
    main_p = argparse.ArgumentParser()
    subp = main_p.add_subparsers(help='sub-command help')

    # load
    subp_p = subp.add_parser('load', help='Sets up your MongoDB for matching.')
    subp_p.add_argument('-t', dest='trials', help=param_trials_help)
    subp_p.add_argument('-c', dest='clinical', help=param_clinical_help)
    subp_p.add_argument('-g', dest='genomic', help=param_genomic_help)
    subp_p.add_argument('--mongo-uri', dest='mongo_uri', required=False, default=None, help=param_mongo_uri_help)
    subp_p.add_argument('--trial-format',
                        dest='trial_format',
                        default='yml',
                        action='store',
                        choices=['yml', 'json', 'bson'],
                        help=param_trial_format_help)
    subp_p.add_argument('--patient-format',
                        dest='patient_format',
                        default='csv',
                        action='store',
                        choices=['csv', 'pkl', 'bson', 'json'],
                        help=param_patient_format_help)
    subp_p.set_defaults(func=load)

    # match
    subp_p = subp.add_parser('match', help='Matches all trials in database to patients')
    subp_p.add_argument('--mongo-uri', dest='mongo_uri', required=False, default=None, help=param_mongo_uri_help)
    subp_p.add_argument('--daemon', dest="daemon", required=False, action="store_true", help=param_daemon_help)
    subp_p.add_argument('--json', dest="json_format", required=False, action="store_true", help=param_json_help)
    subp_p.add_argument('--csv', dest="csv_format", required=False, action="store_true", help=param_csv_help)
    subp_p.add_argument('-o', dest="outpath", required=False, help=param_outpath_help)
    subp_p.set_defaults(func=match)

    # parse args.
    args = main_p.parse_args()
    args.func(args)
