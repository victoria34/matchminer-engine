"""Copyright 2016 Dana-Farber Cancer Institute"""

from tests import TestSetUp


class TestGenomic(TestSetUp):

    def setUp(self):
        super(TestGenomic, self).setUp()

        # add genomic collection
        self.add_genomic()

    def tearDown(self):
        self.db.genomic.drop()

    def _assert(self, key, val, num, neg):
        g, check_neg, _ = self.me.prepare_genomic_criteria({key: val})
        num_found = len(list(self.db.genomic.find(g)))
        assert num_found == num, '%s found with query %s' % (num_found, g)
        assert check_neg is neg
        return g

    def test_build_exon(self):

        # basic query
        self._assert('EXON', 19, 10, False)

        # add 3 genomic documents to database
        # 1 is a Mutation in exon 13 and should not match
        # 2 is a Mutation in exon 19 and should match
        # 3 is a CNV in exon 13 and should match
        self.add_genomic_for_exon_mutation()

        # match no mutation in exon 13
        item = {'exon': 13, 'variant_category': '!Mutation'}
        g, neg, _ = self.me.prepare_genomic_criteria(item)
        num_found = list(self.db.genomic.find(g))
        assert neg is True
        assert len(num_found) == 1, len(num_found)
        for found in num_found:
            assert found['SAMPLE_ID'] == '1', found
            assert found['SAMPLE_ID'] != '2' or found['SAMPLE_ID'] != '3', found

    def test_build_hugo_symbol(self):
        self._assert('hugo_symbol', 'EGFR', 9, False)
        self._assert('hugo_symbol', '!BRAF', 1, True)
        self._assert('hugo_symbol', 'BRAF', 1, False)

    def test_build_protein_change(self):
        self._assert('protein_change', 'p.L858R', 1, False)
        self._assert('wildcard_protein_change', 'p.F346', 3, False)
        self._assert('wildcard_protein_change', '!p.F346', 3, True)

        # Add more protein changes to database and test all edge cases
        # regex should not only match the 4th item
        muts = self.add_genomic_for_regex()
        g = self._assert('wildcard_protein_change', 'p.A0', 1, False)

        genomic = list(self.db.genomic.find(g))
        gmuts = [item['TRUE_PROTEIN_CHANGE'] for item in genomic]
        for mut in gmuts:
            assert mut in muts, gmuts
        for item in muts[:3]:
            assert item not in gmuts, '%s\t%s' % (item, gmuts)

    def test_build_variant_category(self):
        self._assert('variant_category', 'Mutation', 5, False)
        self._assert('variant_category', 'Any Variation', 9, False)

    def test_build_variant_classification(self):
        self._assert('variant_classification', 'In_Frame_Del', 10, False)

    def test_build_wildtype(self):

        self.db.genomic.drop()
        self.add_wildtype()

        self._assert('wildtype', 'true', 2, False)
        self._assert('wildtype', 'false', 1, False)

        self._assert('wildcard_protein_change', 'p.', 2, False)
        self._assert('wildcard_protein_change', '!p.V60000000000E', 0, True)

    def test_build_cnv_call(self):
        self._assert('cnv_call', 'High Amplification', 1, False)
        self._assert('cnv_call', 'Homozygous Deletion', 1, False)
        self._assert('cnv_call', 'Heterozygous Deletion', 1, False)
        self._assert('cnv_call', 'Gain', 1, False)

    def test_build_mmr_status(self):

        self.db.genomic.drop()
        self.add_msi()

        self._assert('mmr_status', 'MMR-Proficient', 1, False)
        self._assert('mmr_status', 'MMR-Deficient', 1, False)
        self._assert('ms_status', 'MSI-H', 1, False)
        self._assert('ms_status', 'MSI-L', 1, False)
        self._assert('ms_status', 'MSS', 1, False)
