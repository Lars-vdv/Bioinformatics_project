import unittest
from bio_seq import bio_seq

class TestBioSeq(unittest.TestCase):

    def test_transitions_transversions_ratio(self):
        """Test the static and instance methods for calculating the T/t ratio."""
        seq1_str = "GCAACGCACAACCACCGC"
        seq2_str = "GTTATGCATGACCATAGC"

        # Test the static method with strings
        transitions, transversions, ratio = bio_seq.transitions_transversions_ratio(seq1_str, seq2_str)
        self.assertEqual(transitions, 5)
        self.assertEqual(transversions, 2)
        self.assertAlmostEqual(ratio, 2.5)

        # Test the static method with bio_seq objects
        seq1_obj = bio_seq(seq1_str)
        seq2_obj = bio_seq(seq2_str)
        transitions, transversions, ratio = bio_seq.transitions_transversions_ratio(seq1_obj, seq2_obj)
        self.assertEqual(transitions, 5)
        self.assertEqual(transversions, 2)
        self.assertAlmostEqual(ratio, 2.5)

        # Test the instance method
        transitions, transversions, ratio = seq1_obj.t_t_ratio(seq2_obj)
        self.assertEqual(transitions, 5)
        self.assertEqual(transversions, 2)
        self.assertAlmostEqual(ratio, 2.5)

        # Test case with no transversions
        seq3_str = "AG"
        seq4_str = "GA"
        transitions, transversions, ratio = bio_seq.transitions_transversions_ratio(seq3_str, seq4_str)
        self.assertEqual(transitions, 2)
        self.assertEqual(transversions, 0)
        self.assertEqual(ratio, float('inf'))

if __name__ == '__main__':
    unittest.main()
