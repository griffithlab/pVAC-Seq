import argparse
import vcf
import re
import sys

def parse_csq_format(vcf_reader):
    info_fields = vcf_reader.infos

    if info_fields['CSQ'] is None:
        sys.exit('Failed to extract format string from info description for tag (CSQ)')
    else:
        csq_header = info_fields['CSQ']
        format_pattern = re.compile('Format: (.*)')
        match = format_pattern.search(csq_header.desc)
        return match.group(1)

def parse_csq_entries(csq_entries, csq_format):
    csq_format_array = csq_format.split('|')

    transcripts = []
    for entry in csq_entries:
        values = entry.split('|')
        transcript = {}
        for key, value in zip(csq_format_array, values):
            transcript[key] = value
        transcripts.append(transcript)

    return transcripts

def position_out_of_bounds(position, sequence):
    return position > len(sequence)-1

#This subroutine is a bit funky but it was designed that way to mirror
#distance_from_end to increase code readability from the caller's perspective
def distance_from_start(position, string):
    return position

def distance_from_end(position, string):
    return len(string) - 1 - position;

def determine_peptide_sequence_length(full_wildtype_sequence_length, peptide_sequence_length, entry):
    actual_peptide_sequence_length = peptide_sequence_length

    #If the wildtype sequence is shorter than the desired peptide sequence
    #length we use the wildtype sequence length instead so that the extraction
    #algorithm below works correctly
    if full_wildtype_sequence_length < actual_peptide_sequence_length:
        actual_peptide_sequence_length = full_wildtype_sequence_length
        print('Wildtype sequence length is shorter than desired peptide sequence length at position (%s, %s). Using wildtype sequence length (%s) instead.' % (entry.CHROM, entry.POS, actual_peptide_sequence_length))

    return actual_peptide_sequence_length

def determine_flanking_sequence_length(full_wildtype_sequence_length, peptide_sequence_length, entry):
    actual_peptide_sequence_length = determine_peptide_sequence_length(full_wildtype_sequence_length, peptide_sequence_length, entry)
    if actual_peptide_sequence_length%2 == 0:
        return (actual_peptide_sequence_length-2) / 2
    else:
        return (actual_peptide_sequence_length-1) / 2

def get_wildtype_subsequence(position, full_wildtype_sequence, wildtype_amino_acid_length, peptide_sequence_length, entry):
    one_flanking_sequence_length = int(determine_flanking_sequence_length(len(full_wildtype_sequence), peptide_sequence_length, entry))
    peptide_sequence_length = 2 * one_flanking_sequence_length + wildtype_amino_acid_length

    # We want to extract a subset from full_wildtype_sequence that is
    # peptide_sequence_length long so that the position ends
    # up in the middle of the extracted sequence.
    # If the position is too far toward the beginning or end of
    # full_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.
    if distance_from_start(position, full_wildtype_sequence) < one_flanking_sequence_length:
        wildtype_subsequence = full_wildtype_sequence[:peptide_sequence_length]
        mutation_position = position
    elif distance_from_end(position, full_wildtype_sequence) < one_flanking_sequence_length:
        start_position = len(full_wildtype_sequence) - peptide_sequence_length
        wildtype_subsequence = full_wildtype_sequence[start_position:]
        mutation_position = peptide_sequence_length - distance_from_end(position, full_wildtype_sequence) - 1
    elif distance_from_start(position, full_wildtype_sequence) >= one_flanking_sequence_length and distance_from_end(position, full_wildtype_sequence) >= one_flanking_sequence_length:
        start_position = position - one_flanking_sequence_length
        end_position   = start_position + peptide_sequence_length
        wildtype_subsequence = full_wildtype_sequence[start_position:end_position]
        mutation_position = one_flanking_sequence_length
    else:
        sys.exit("Something went wrong during the retrieval of the wildtype sequence at position(%s, %s, %s)" % chromsome, start, stop)
    return mutation_position, wildtype_subsequence

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('Generate Variant Sequences', description='')
    parser.add_argument('input_file', type=argparse.FileType('r'), help='input vcf',)
    parser.add_argument('peptide_sequence_length', type=int, help='length of the peptide sequence')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='output FASTA file')

    args = parser.parse_args(args_input)

    peptide_sequence_length = args.peptide_sequence_length

    vcf_reader = vcf.Reader(args.input_file)
    csq_format = parse_csq_format(vcf_reader)

    transcript_count = {}
    pattern = re.compile('([A-Z])(\d+)([A-Z])');
    for entry in vcf_reader:
        transcripts = parse_csq_entries(entry.INFO['CSQ'], csq_format)
        for transcript in transcripts:
            consequence_string = transcript['Consequence']
            if consequence_string is None:
                continue
            consequences = list(set(consequence.lower() for consequence in consequence_string.split('&')))
            full_wildtype_sequence = transcript['WildtypeProtein']

            if 'missense_variant' in consequences:
                wildtype_amino_acid, mutant_amino_acid = transcript['Amino_acids'].split('/')
                if wildtype_amino_acid == '-':
                    position = int(transcript['Protein_position'].split('-', 1)[0])
                    wildtype_amino_acid_length = 0
                else:
                    position = int(transcript['Protein_position']) - 1;
                    wildtype_amino_acid_length = len(wildtype_amino_acid)
            else:
                continue

            if position_out_of_bounds(position, full_wildtype_sequence):
                continue

            if transcript['Feature'] in transcript_count:
                transcript_count[transcript['Feature']] += 1
            else:
                transcript_count[transcript['Feature']] = 1

            mutation_start_position, wildtype_subsequence = get_wildtype_subsequence(position, full_wildtype_sequence, wildtype_amino_acid_length, peptide_sequence_length, entry)
            mutation_end_position = mutation_start_position + wildtype_amino_acid_length
            mutant_subsequence = wildtype_subsequence[:mutation_start_position] + mutant_amino_acid + wildtype_subsequence[mutation_end_position:];
            if mutant_amino_acid is '':
                mutant_amino_acid = '-'
            variant_id = '%s_%s_%s.%s%s%s' % (transcript['SYMBOL'], transcript['Feature'], transcript_count[transcript['Feature']], wildtype_amino_acid, transcript['Protein_position'], mutant_amino_acid)

            for designation, subsequence in zip(['WT', 'MT'], [wildtype_subsequence, mutant_subsequence]):
                args.output_file.writelines('>%s.%s\n' % (designation, variant_id))
                args.output_file.writelines('%s\n' % subsequence)

    args.input_file.close()
    args.output_file.close()

if __name__ == '__main__':
    main()
