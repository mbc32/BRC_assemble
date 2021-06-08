import unittest
import brc_assemble as brc
dna_translation_table = str.maketrans("ACTG", "TGAC")


class AssembleTestCase(unittest.TestCase):
    def test_read_fastq(self):
        fastq_lines = ['@xxxxx', 'ACTG', '+', '!!!!']
        expected_dict = {'header': '@xxxxx', 'sequence': 'ACTG', 'separator': '+', 'base_quality': '!!!!'}
        observed_dict = brc.process_fastq(fastq_lines)
        print(observed_dict['sequence'])
        self.assertEqual(expected_dict, observed_dict)

    def test_kmers(self):
        sequence = 'AGTTCGAAG'
        expected_kmers = [('AGT', 0), ('GTT', 1), ('TTC', 2), ('TCG', 3), ('CGA', 4), ('GAA', 5), ('AAG', 6)]
        observed_kmers = brc.create_kmers(sequence, 3)
        self.assertEqual(expected_kmers, observed_kmers, msg='The correct kmers was created ')

    def test_Reverse_complement(self):
        expected_dna = 'CGTA'
        observed_dna = brc.reverse_complement('TaCg')
        self.assertEqual(expected_dna, observed_dna, msg='reverse complement DNA is correct')

    def test_kmer_match(self):
        read1 = {'sequence': 'AGTTCGAAG'}
        read2 = {'sequence': 'AAGGATCAA'}
        expected_match = [(6, 0)]
        observed_match = brc.paired_read_kmer_match(read1, read2, 3)
        self.assertEqual(expected_match, observed_match, msg='The correct assemble was created')

    def test_seq_assemble(self):
        expected_seq = 'AGGTTTCCC'
        observed_seq = brc.assemble_sequence({'sequence': 'AGGTTT'}, {'sequence': 'GGGAAA'}, [(3, 0)])
        self.assertEqual(expected_seq, observed_seq, msg='The correct sequence was selected')

    def test_fasta(self):
        sequence = ''.join(['A' for _ in range(50)])
        expected_file = './expected.fasta'
        observed_file = './observed.fasta'
        brc.write_fasta_file(observed_file, sequence, 'test')
        #self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
