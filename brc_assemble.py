import argparse
"""
Module that assembles paired sequence reads. The clone insert length is assumed to be less than the combined read length
so the reads overlap.
"""

dna_translation_table = str.maketrans("ACTG", "TGAC")
fastq_lines_per_record = 4
kmer_size = 8


def process_fastq(fastq_lines):
    """
    Reads FASTQ record which have 4 lines per read: Header(@), sequence, separator(+), base quality.
    :param fastq_lines: list with 4 elements: Header(@), sequence, separator(+), base_quality
    :return: Dictionary with 4 keys: header, sequence, separator, base_quality
    """

    keys = ['header', 'sequence', 'separator', 'base_quality']
    return {key: value for key, value in zip(keys, fastq_lines)}


def reverse_complement(dna):
    """
    create the reverse complement DNA string in uppercase
    :param dna: str
    :return: str
    """
    return dna.upper().translate(dna_translation_table)[::-1]


def create_kmers(sequence, kmer_length):
    """
    Creates kmers of length 'kmer_length' from 'sequence'.
    :param sequence: Any sequence of length > kmer_length
    :param kmer_length: length of the kmer
    :return: list of tuples (kmer, start index)
    """
    if len(sequence) < kmer_length:
        raise ValueError("The Sequence {} must be equal to or longer than kmer_length".format(sequence, kmer_length))

    kmers = list()
    number_of_kmers = len(sequence) - kmer_length + 1

    for i in range(number_of_kmers):
        kmer = sequence[i:i + kmer_length]
        kmers.append((kmer, i))
    return kmers


def paired_read_kmer_match(read_1, read_2, kmer_length):
    """
    Creates kmers from each fastq file and look for matches
    :param read_1: Dictionary with 4 keys: header, sequence, separator, base_quality
    :param read_2: Dictionary with 4 keys: header, sequence, separator, base_quality
    :param kmer_length:
    :return: list of tuples with start position of the kmer match
    """
    kmers1 = create_kmers(read_1['sequence'], kmer_length)
    kmers2 = create_kmers(read_2['sequence'], kmer_length)
    kmer_match = list()
    for kmer1 in kmers1:
        kmer1_seq, kmer1_index = kmer1
        for kmer2 in kmers2:
            kmer2_seq, kmer2_index = kmer2
            count = sum(1 for base1, base2 in zip(kmer1_seq, kmer2_seq) if base1 != base2)

            if count == 0:
                kmer_match.append((kmer1_index, kmer2_index))

    return kmer_match


def assemble_sequence(read_1, read_2, kmer_match):
    """
    Assemble the reads into on sequence using the overlap found by comparing kmers
    :param read_1:
    :param read_2:
    :param kmer_match:
    :return: The longest assembled sequence
    """
    sequences = list()
    for match in kmer_match:
        read_1_index, read_2_index = match
        seq_1 = read_1['sequence'][:read_1_index]
        read_2_rc = reverse_complement(read_2['sequence'])
        seq_2 = read_2_rc[read_2_index:]
        sequences.append(seq_1 + seq_2)
    if sequences:
        return max(sequences, key=len)
    else:
        return None


def write_fasta_file(file_path, sequence, header):
    """
    write fasta file
    :param file_path:
    :param sequence:
    :param header:
    :return:
    """
    with open(file_path, 'a') as file:
        file.write('> ' + header + '\n')
        file.write(sequence + '\n')
        file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1')
    parser.add_argument('--read2')
    parser.add_argument('--output')
    args = parser.parse_args()

    read1_path = args.read1
    read2_path = args.read2
    out_file_path = args.output

    with open(read1_path, 'r') as read1_handle, open(read2_path, 'r') as read2_handle:
        read1_temp_lines = list()
        read2_temp_lines = list()
        for read1_line, read2_line in zip(read1_handle, read2_handle):
            read1_temp_lines.append(read1_line.rstrip())
            read2_temp_lines.append(read2_line.rstrip())
            if len(read1_temp_lines) == fastq_lines_per_record:
                read1_record = process_fastq(read1_temp_lines)
                read2_record = process_fastq(read2_temp_lines)
                read1_temp_lines = []
                read2_temp_lines = []

                kmer_matches = paired_read_kmer_match(read1_record, read2_record, kmer_size)
                assembled_sequence = assemble_sequence(read1_record, read2_record, kmer_matches)
                if assembled_sequence is not None:
                    write_fasta_file(out_file_path, assembled_sequence, read1_record['header'])


