"""Copyright 2016 Dana-Farber Cancer Institute"""

from tests import TestSetUp


class TestClinical(TestSetUp):

    def setUp(self):
        super(TestClinical, self).setUp()

        # add clinical collection
        self.add_clinical()

    def tearDown(self):
        self.db.clinical.drop()

    def _assert(self, key, val, num):
        c = self.me.prepare_clinical_criteria({key: val})
        num_found = len(list(self.db.clinical.find(c)))
        assert num_found == num, '%s found with query %s' % (num_found, c)

    def test_age_numerical(self):
        self._assert('age_numerical', '>=18', 5)
        self._assert('age_numerical', '<=18', 5)
        self._assert('age_numerical', '<.5', 1)

    def test_oncotree_primary_diagnosis(self):
        self._assert('oncotree_primary_diagnosis', 'Melanoma', 5)
        self._assert('oncotree_primary_diagnosis', 'Glioblastoma', 4)
        self._assert('oncotree_primary_diagnosis', 'Adrenal Gland', 1)

    def test_gender(self):
        self._assert('gender', 'Male', 5)
        self._assert('gender', 'Female', 5)
