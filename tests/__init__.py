"""Copyright 2016 Dana-Farber Cancer Institute"""

import os
import json
import yaml
import random
import string
import unittest
import datetime as dt
from bson.objectid import ObjectId

from matchengine.utilities import get_db
from matchengine.engine import MatchEngine

YAML_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/yaml/'))


class TestSetUp(unittest.TestCase):

    def setUp(self):
        """
        Descriptions of test patients

        1: >18, Adrenal Gland, Female, BRAF F346R Mutation
        2: >18, Melanoma, Female, EGFR L858R Mutation
        3: >18, Melanoma, Female, EGFR F346A Mutation
        4: >18, Melanoma, Female, EGFR F346B Mutation
        5: >18, Melanoma, Female, EGFR F000F Mutation
        6: >0.5 && <18, Melanoma, Male, EGFR SV
        7: >0.5 && <18, Glioblastoma, Male, EGFR CNV Hetero del
        8: >0.5 && <18, Glioblastoma, Male, EGFR CNV Gain
        9: >0.5 && <18, Glioblastoma, Male, EGFR CNV Homo del
        10: <0.5, Glioblastoma, Male, EGFR CNV High amp

        Descriptions of test trials
        00-001.yml: dose: EGFR L858R && >=18/_SOLID_
        00-002.yml: arm: EGFR L858R && >=18/_SOLID_
        00-003.yml: step: EGFR L858R && >=18/_SOLID_
        00-004.yml dose: EGFR L858R && >=18/_SOLID_
        00-005.yml 2 doses: EGFR L858R && >=18/_SOLID_
        00-006.yml exon: !13
        """

        self.db = get_db(None)
        for res in ["clinical", "dashboard", "filter", "genomic", "hipaa", "match", "normalize", "oplog"
                    "response", "statistics", "status", "team", "trial", "trial_match", "user"]:
            self.db.drop_collection(res)

        self.me = MatchEngine(self.db)

        self.trials = {}
        self.clinical_id = ObjectId()
        self.mrn = 'TCGA-BH-A1FR'
        self.sample_id = 'TCGA-OR-A5J1'
        self.mrns = [self.mrn] + [self.__random_id() for _ in range(9)]
        self.sample_ids = [self.sample_id] + [self.__random_id() for _ in range(9)]
        self.clinical_ids = [self.clinical_id] + [ObjectId() for _ in range(9)]
        self.static_date = dt.datetime.today().replace(year=2016, month=11, day=3)

        # clinical collection
        self.oncotree_diagnoses = ['Adrenal Gland'] + ['Melanoma'] * 5 + ['Glioblastoma'] * 4
        self.genders = ['Female'] * 5 + ['Male'] * 5
        self.ages = [self.static_date.replace(year=1997)] * 5 + [self.static_date.replace(year=2010)] * 4 + \
                    [self.static_date.replace(month=9)]
        self.clinical = [{
            '_id': clinical_id,
            'ONCOTREE_PRIMARY_DIAGNOSIS_NAME': diagnosis,
            'SAMPLE_ID': sample_id,
            'VITAL_STATUS': 'alive',
            'DFCI_MRN': mrn,
            'REPORT_DATE': self.static_date,
            'BIRTH_DATE': age,
            'GENDER': gender
        } for diagnosis, gender, age, clinical_id, sample_id, mrn in zip(
            self.oncotree_diagnoses, self.genders, self.ages, self.clinical_ids, self.sample_ids, self.mrns)]

        # genomic collection
        self.genes = ['BRAF'] + ['EGFR'] * 9
        self.protein_changes = ['p.F346R', 'p.L858R', 'p.F346A', 'p.F346B', 'p.F000F', None, None, None, None, None]
        self.variant_categories = ['MUTATION'] * 5 + ['SV', 'CNV', 'CNV', 'CNV', 'CNV']
        self.wildtypes = [False] * 10
        self.cnv_calls = [None, None, None, None, None, None,
                          'Heterozygous deletion', 'Gain', 'Homozygous deletion', 'High level amplification']
        self.genomic = [{
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': protein_change,
            'VARIANT_CATEGORY': variant_category,
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'WILDTYPE': wildtype,
            'CLINICAL_ID': _id,
            'CNV_CALL': cnv_call,
            'TRUE_HUGO_SYMBOL': gene,
            'SAMPLE_ID': sample_id,
            'TRUE_TRANSCRIPT_EXON': 19
        } for protein_change, variant_category, wildtype, cnv_call, gene, _id, sample_id in zip(
            self.protein_changes, self.variant_categories, self.wildtypes,
            self.cnv_calls, self.genes, self.clinical_ids, self.sample_ids
        )]

        # test trials
        self.test_trials = ['00-001', '00-002', '00-003']

        # demo match results
        pnos = ['00-001', '00-001', '00-001', '00-002', '00-002', '00-002']
        mlevels = ['arm', 'arm', 'arm', 'dose', 'dose', 'dose']
        iids = ['1', '2', '3', '4', '5', '6']
        galts = ['Alt1', 'Alt2', 'Alt2', 'Alt3', 'Alt3', 'Alt3']
        self.matches = [{
            'mrn': 'SAMPLE1',
            'sample_id': 'SAMPLE1-ID',
            'protocol_no': protocol_no,
            'match_level': match_level,
            'internal_id': internal_id,
            'genomic_alteration': genomic_alteration
        } for protocol_no, match_level, internal_id, genomic_alteration in zip(
            pnos, mlevels, iids, galts
        )]

    def tearDown(self):
        self.db.map.drop()
        self.db.trial_match.drop()

    def add_clinical(self):
        """Add all clinical documents to database needed for unit tests"""
        self.db.clinical.insert_many(self.clinical)

    def add_genomic(self):
        """Add all genomic documents to database needed for unit tests"""
        self.db.genomic.insert_many(self.genomic)

    def add_genomic_for_exon_mutation(self):
        """
        One of the unit tests has to test for the scenario where a trial calls for no mutation in a specific exon
        This method sets up the db for that unit test
        """

        # reset
        self.db.genomic.drop()

        # add genomic collection to db
        variants = ['MUTATION', 'MUTATION', 'CNV']
        clinical_ids = [ObjectId() for _ in range(3)]
        cnv_calls = [None, None, 'Homozygous deletion']
        sample_ids = ['1', '2', '3']                        # 2 and 3 should match; 1 should not match
        exons = [13, 19, 13]
        g = [{
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': 'p.V600E',
            'VARIANT_CATEGORY': variant_category,
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'WILDTYPE': False,
            'CLINICAL_ID': _id,
            'CNV_CALL': cnv_call,
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'SAMPLE_ID': sample_id,
            'TRUE_TRANSCRIPT_EXON': exon
        } for variant_category, _id, cnv_call, sample_id, exon in zip(
            variants, clinical_ids, cnv_calls, sample_ids, exons
        )]

        # add
        self.db.genomic.insert_many(g)

    def add_genomic_for_regex(self):
        """Adds genomic entries to specifically test pymongo $regex operations"""

        # reset
        self.db.genomic.drop()

        # add genomic collection to db
        # p.A0
        muts = ['p.A000Z', 'p.B0_A0B', 'p.B0A', 'p.A0B', 'p.A0fs*6', 'p.A0*', 'p.A0_B12insL']
        g = [{'TRUE_PROTEIN_CHANGE': mut} for mut in muts]
        self.db.genomic.insert_many(g)
        return muts

    def add_trials(self, trials=None):
        """Add all trial documents to database needed for unit tests"""

        if not trials:
            trials=self.test_trials

        # read trials from direcotry
        for yml in os.listdir(YAML_DIR):
            ymlpath = os.path.join(YAML_DIR, yml)

            if yml.split('.')[0] not in trials:
                continue

            # convert yml to json format
            with open(ymlpath) as f:
                t = yaml.load(f.read())
                self.trials[t['protocol_no']] = t

                # add trial to db
                self.db.trial.insert_one(t)

    def add_wildtype(self):
        """Adds a variety of wildtype genomic entries to the database"""

        clinical_ids = [ObjectId() for _ in range(3)]
        wts = [True, True, False]
        sample_ids = ['1', '2', '3']
        g = [{
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': 'p.V600E',
            'VARIANT_CATEGORY': 'MUTATION',
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'WILDTYPE': wt,
            'CLINICAL_ID': _id,
            'CNV_CALL': None,
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'SAMPLE_ID': sample_id,
            'TRUE_TRANSCRIPT_EXON': 13
        } for _id, sample_id, wt in zip(clinical_ids, sample_ids, wts)]

        g.append({
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': 'p.V600E',
            'VARIANT_CATEGORY': 'MUTATION',
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'CLINICAL_ID': ObjectId(),
            'CNV_CALL': None,
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'SAMPLE_ID': '4',
            'TRUE_TRANSCRIPT_EXON': 13
        })

        self.db.genomic.insert_many(g)

    @staticmethod
    def __random_id():
        """Returns a random MRN or SAMPLE_ID formatted "XXXX-XX-XXXX"""
        return 'TCGA-%s-%s' % (
            ''.join(random.choice(string.ascii_uppercase) for _ in range(2)),
            ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
        )

    @staticmethod
    def _debug(t):
        try:
            return json.dumps(t, sort_keys=True, indent=4)
        except TypeError:
            return t
