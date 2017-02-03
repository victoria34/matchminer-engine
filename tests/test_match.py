"""Copyright 2016 Dana-Farber Cancer Institute"""

import os
import json

from tests import TestSetUp

YAML_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/yaml/'))


class TestMatch(TestSetUp):

    def setUp(self):
        super(TestMatch, self).setUp()
        self.add_clinical()
        self.add_genomic()
        self.add_trials()

    def tearDown(self):
        self.db.clinical.drop()
        self.db.genomic.drop()
        self.db.trial.drop()

    def _match(self, match):
        g = self.me.create_match_tree(match)
        results, ginfo = self.me.traverse_match_tree(g)
        return results

    def _find(self, collection, query):
        return set(item['SAMPLE_ID'] for item in list(self.db[collection].find(query, {'SAMPLE_ID': 1})))

    def test_simple_match(self):

        match = {
            'and': [
                {'genomic': {'hugo_symbol': 'EGFR'}},
                {'clinical': {'oncotree_primary_diagnosis': "Non-Small Cell Lung Cancer"}}
            ]
        }

        # turn each trial criterium into a mongo query
        g, _ = self.me.prepare_genomic_criteria(match['and'][0]['genomic'])
        c = self.me.prepare_clinical_criteria(match['and'][1]['clinical'])

        # check for matching patients by "SAMPLE_ID" manually
        genomic_sample_ids = self._find('genomic', g)
        clinical_sample_ids = self._find('clinical', c)
        mongo_results = genomic_sample_ids.intersection(clinical_sample_ids)

        # check that you get the same number of matching patients through the matchengine
        matchengine_results = self._match(match)
        assert sorted(mongo_results) == sorted(matchengine_results)

    def test_dose_level_match(self):

        # get dose level match from yml
        test_inp = self._read_file(os.path.join(YAML_DIR, '00-004.yml'))
        status, data = self.me.validate_yaml_format(test_inp)
        match = data['treatment_list']['step'][0]['arm'][0]['dose_level'][0]['match'][0]

        # turn each trial criterium into a mongo query
        g1, _ = self.me.prepare_genomic_criteria(match['and'][0]['genomic'])
        c1 = self.me.prepare_clinical_criteria(match['and'][1]['or'][0]['clinical'])
        c2 = self.me.prepare_clinical_criteria(match['and'][1]['or'][1]['clinical'])
        c3 = self.me.prepare_clinical_criteria(match['and'][1]['or'][2]['clinical'])
        c4 = self.me.prepare_clinical_criteria(match['and'][1]['or'][3]['clinical'])
        c5 = self.me.prepare_clinical_criteria(match['and'][1]['or'][4]['clinical'])
        c6 = self.me.prepare_clinical_criteria(match['and'][1]['or'][5]['clinical'])

        # check for matching_patients by "SAMPLE_ID" manually
        genomic_ids = self._find('genomic', g1)
        cdb1 = self._find('clinical', c1)
        cdb2 = self._find('clinical', c2)
        cdb3 = self._find('clinical', c3)
        cdb4 = self._find('clinical', c4)
        cdb5 = self._find('clinical', c5)
        cdb6 = self._find('clinical', c6)
        cdbsearch = cdb1.intersection(cdb2).intersection(cdb3).intersection(cdb4).intersection(cdb5).intersection(cdb6)
        dbsearch = set(list(cdbsearch) + list(genomic_ids))

        # check that you get the same number of matching patients through the matchengine
        matchengine = self._match(match)
        assert sorted(dbsearch) == sorted(matchengine)

    def test_system_match(self):
        """
        Loops through all trials in the database, finds the matches at each dose, arm, and step level and
        creates a match tree and finds all matches in the database.
        If the sample id of the given mrn matches that match tree, it is added to a master list tracking which
        trials the mrn matched to.
        """

        # checks that the entire process executes successfully
        self.me.find_trial_matches()

    def test_assess_match(self):

        p = self.mrns[1]
        mrn_map = dict(zip([self.sample_ids[1]], [p]))
        trial_matches = []
        trial = self.db.trial.find_one({'protocol_no': '00-001'})
        trial_segment = trial['treatment_list']['step'][0]['arm'][0]['dose_level'][0]
        match_segment = 'dose'

        # add sample id to trial_matches dictionary
        t = self.me._assess_match(mrn_map, trial_matches, trial, trial_segment, match_segment, 'open')
        mrns = [item['mrn'] for item in t]
        trials = [item['protocol_no'] for item in t if item['mrn'] == p]
        doses = [item['internal_id'] for item in t if item['match_level'] == 'dose' and item['mrn'] == p]
        assert p in mrns, self._debug(mrns)
        assert trials == ['00-001'], self._debug(trials)
        assert doses == ['1'], self._debug(doses)

        # match multiple intratrial doses
        self.add_trials(trials=['00-004'])
        trial = self.db.trial.find_one({'protocol_no': '00-004'})
        trial_segment = trial['treatment_list']['step'][0]['arm'][0]['dose_level'][0]
        t = self.me._assess_match(mrn_map, trial_matches, trial, trial_segment, match_segment, 'open')

        mrns = [item['mrn'] for item in t]
        trials = [item['protocol_no'] for item in t if item['mrn'] == p]
        doses = [item['internal_id'] for item in t if item['match_level'] == 'dose' and item['mrn'] == p]
        assert p in mrns, self._debug(mrns)
        assert trials == ['00-001'], self._debug(trials)
        assert doses == ['1'], self._debug(doses)

        # match multiple intertrial doses
        self.db.trial.drop()
        self.add_trials(trials=['00-005'])
        trial = self.db.trial.find_one({'protocol_no': '00-005'})

        # dose 1
        trial_segment = trial['treatment_list']['step'][0]['arm'][0]['dose_level'][0]
        t = self.me._assess_match(mrn_map, [], trial, trial_segment, match_segment, 'open')

        # dose 2
        trial_segment = trial['treatment_list']['step'][0]['arm'][0]['dose_level'][1]
        t = self.me._assess_match(mrn_map, t, trial, trial_segment, match_segment, 'open')

        mrns = [item['mrn'] for item in t]
        trials = [item['protocol_no'] for item in t if item['mrn'] == p]
        doses = [item['internal_id'] for item in t if item['match_level'] == 'dose' and item['mrn'] == p]
        assert p in mrns, self._debug(mrns)
        assert list(set(trials)) == ['00-005'], self._debug(trials)
        assert doses == ['5', '6'], self._debug(doses)

    @staticmethod
    def _read_file(file):
        fh = open(file, 'r')
        data = ''
        for stream in fh:
            data += stream
        return data
