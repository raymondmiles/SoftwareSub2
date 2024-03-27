import unittest
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fastaTogenbank import parse_file
from fastaTogenbank import write_genbank
from fastaTogenbank import has_valid_header, is_fasta
from fastaTogenbank import print_fasta_statistics
# These tests can be run as python -m unittest test_script.py -v


class TestParseFile(unittest.TestCase):
    def setUp(self):
        # Create a temporary fasta file for testing
        self.test_fasta = "test.fasta"
        with open(self.test_fasta, 'w') as file:
            file.write(">seq1\nATGCTAGCTAGCTACGATCG\n")
            file.write(">seq2\nATCGATCGATCGATCGATCG\n")

    def tearDown(self):
        # Remove the temporary file after tests
        os.remove(self.test_fasta)

    def test_parse_file_sequence_content(self):
        # Test if the sequences are correctly parsed
        sequences = parse_file(self.test_fasta)
        expected_sequences = ['ATGCTAGCTAGCTACGATCG', 'ATCGATCGATCGATCGATCG']
        for sequence, expected_seq in zip(sequences, expected_sequences):
            self.assertEqual(str(sequence.seq), expected_seq)


class TestWriteGenbank(unittest.TestCase):
    def setUp(self):
        # Create a list of SeqRecord objects to simulate a parsed list
        self.parsed_list = [
            SeqRecord(Seq("ATGCTAGCTAGCTACGATCG"), id="seq1"),
            SeqRecord(Seq("ATCGATCGATCGATCGATCG"), id="seq2"),
            SeqRecord(Seq("ATCGATCGATCGATCGATCGATCG"), id="seq3")
        ]
        for sequence in self.parsed_list:
            sequence.annotations["molecule_type"] = "DNA"

        self.genbank_file = "test_output.gb"

    def tearDown(self):
        # Remove the temporary file after tests
        if os.path.exists(self.genbank_file):
            os.remove(self.genbank_file)

    def test_write_genbank_returns_count(self):
        # Test if the function returns the correct count of SeqRecord objects written
        count = write_genbank(self.parsed_list, self.genbank_file)
        self.assertEqual(count, len(self.parsed_list))

    def test_write_genbank_creates_file(self):
        # Test if the function creates a file at the specified path
        write_genbank(self.parsed_list, self.genbank_file)
        self.assertTrue(os.path.exists(self.genbank_file))
    
    def test_write_genbank_file_content(self):
        # Test if the content of the created file is correct
        write_genbank(self.parsed_list, self.genbank_file)
        with open(self.genbank_file, 'r') as file:
            content = file.read()
            # Check for some expected content in the file
            self.assertIn("LOCUS", content)
            self.assertIn("seq1", content)
            self.assertIn("seq2", content)
            self.assertIn("seq3", content)
    def test_write_genbank_ioerror(self):
        with self.assertRaises(IOError):
            write_genbank(self.parsed_list, "/invalid/path/output.gb")
    def test_write_genbank_unexpectederror(self):
        with self.assertRaises(Exception):
            write_genbank("not a parsed list", "OUTPUTS/output.gb")



    
class TestPrintFastaStatistics(unittest.TestCase):
    def setUp(self):
        self.testing_fasta = "test.fasta"
        with open(self.testing_fasta, 'w') as file:
            file.write(">seq1\nATG\n") #3bp
            file.write(">seq2\nATCGATCGATCG\n") #12bp
            file.write(">seq3\nATCGATCGATCGATCGATCGT\n") #21bp
    def tearDown(self):
        # Remove the temporary file after tests
        os.remove(self.testing_fasta)
    def test_fasta_statistics(self):
        # Test if the sequences are correctly parsed
        statisticsDictionary = print_fasta_statistics(self.testing_fasta)
        self.assertEqual(statisticsDictionary["max_read_length"], 21)
        self.assertEqual(statisticsDictionary["mean_read_length"], 12)
        self.assertEqual(statisticsDictionary["median_read_length"], 12)
        self.assertEqual(statisticsDictionary["min_read_length"], 3)
        self.assertEqual(statisticsDictionary["total_reads"], 3)

class TestTestFastaFormat(unittest.TestCase):
    def setUp(self):
        self.test_fasta = "test.fasta"
        self.not_fasta = "test.not"
        with open(self.test_fasta, 'w') as file:
            file.write(">seq1\nAGGGAGCTAGCTACGATCG\n")
            file.write(">seq2\nATTAATATCGATCGATCGATCG\n")
        with open(self.not_fasta, 'w') as file:
            file.write("seq1\nAGGGAGCTAGCTACGATCG\n")
            file.write("seq2\nATTAATATCGATCGATCGATCG\n")

    def tearDown(self):
        # Remove the temporary file after tests
        os.remove(self.test_fasta)
        os.remove(self.not_fasta)
    def testIsFasta(self):
        self.assertEqual(is_fasta(self.test_fasta), True)
        self.assertEqual(is_fasta(self.not_fasta), False)
        self.assertEqual(has_valid_header(self.test_fasta), True)
        self.assertEqual(has_valid_header(self.not_fasta), False)

# Run the tests
if __name__ == '__main__':
    unittest.main()
