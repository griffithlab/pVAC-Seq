import sys
from abc import ABCMeta


class Variant(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        #This should be a parameter that is passed in to determine_fasta_sequences
        self.peptide_sequence_length      = kwargs['peptide_sequence_length']
        #This will go away when the check for epitope length is moved to the IEDB caller
        self.epitope_length               = kwargs['epitope_length']

    def determine_peptide_sequence_length(self, full_wildtype_sequence_length):
        actual_peptide_sequence_length = self.peptide_sequence_length

        #If the wildtype sequence is shorter than the desired peptide sequence
        #length we use the wildtype sequence length instead so that the extraction
        #algorithm below works correctly
        if full_wildtype_sequence_length < actual_peptide_sequence_length:
            actual_peptide_sequence_length = full_wildtype_sequence_length

        return actual_peptide_sequence_length

    def determine_flanking_sequence_length(self, full_wildtype_sequence_length):
        actual_peptide_sequence_length = self.determine_peptide_sequence_length(full_wildtype_sequence_length)
        if actual_peptide_sequence_length%2 == 0:
            return int((actual_peptide_sequence_length-2) / 2)
        else:
            return int((actual_peptide_sequence_length-1) / 2)

    def sequence_is_valid(self, sequence):
        if '*' in sequence:
            return 0

        if 'X' in sequence:
            return 0

        #this check should be done in the IEDB caller
        if len(sequence) < self.epitope_length:
            return 0

        return 1

class SingleTranscriptVariant(Variant):
    def __init__(self, **kwargs):
        Variant.__init__(self, **kwargs)
        self.protein_position             = kwargs['protein_position']
        self.wildtype_amino_acid_sequence = kwargs['wildtype_amino_acid_sequence']

    def position_out_of_bounds(self):
        return self.position > len(self.wildtype_amino_acid_sequence)-1

class InframeVariant(SingleTranscriptVariant):
    def __init__(self, **kwargs):
        SingleTranscriptVariant.__init__(self, **kwargs)
        self.amino_acid_change          = kwargs['amino_acid_change']
        self.wildtype_amino_acid        = self.determine_wildtype_amino_acid()
        self.mutant_amino_acid          = self.determine_mutant_amino_acid()
        self.wildtype_amino_acid_length = self.determine_wildtype_amino_acid_length()
        self.position                   = self.calculate_position()

    #This subroutine is a bit funky but it was designed that way to mirror
    #distance_from_end to increase code readability from the caller's perspective
    def distance_from_start(self, position, string):
        return position

    def distance_from_end(self, position, string):
        return len(string) - 1 - position

    def get_wildtype_subsequence(self):
        one_flanking_sequence_length = self.determine_flanking_sequence_length(len(self.wildtype_amino_acid_sequence))
        peptide_sequence_length = 2 * one_flanking_sequence_length + self.wildtype_amino_acid_length

        # We want to extract a subset from full_wildtype_sequence that is
        # peptide_sequence_length long so that the position ends
        # up in the middle of the extracted sequence.
        # If the position is too far toward the beginning or end of
        # full_wildtype_sequence there aren't enough amino acids on one side
        # to achieve this.
        if self.distance_from_start(self.position, self.wildtype_amino_acid_sequence) < one_flanking_sequence_length:
            wildtype_subsequence = self.wildtype_amino_acid_sequence[:peptide_sequence_length]
            mutation_position = self.position
        elif self.distance_from_end(self.position, self.wildtype_amino_acid_sequence) < one_flanking_sequence_length:
            start_position = len(self.wildtype_amino_acid_sequence) - peptide_sequence_length
            wildtype_subsequence = self.wildtype_amino_acid_sequence[start_position:]
            mutation_position = peptide_sequence_length - self.distance_from_end(self.position, self.wildtype_amino_acid_sequence) - 1
        elif self.distance_from_start(self.position, self.wildtype_amino_acid_sequence) >= one_flanking_sequence_length and self.distance_from_end(self.position, self.wildtype_amino_acid_sequence) >= one_flanking_sequence_length:
            start_position = self.position - one_flanking_sequence_length
            end_position   = start_position + peptide_sequence_length
            wildtype_subsequence = self.wildtype_amino_acid_sequence[start_position:end_position]
            mutation_position = one_flanking_sequence_length
        else:
            sys.exit("ERROR: Something went wrong during the retrieval of the wildtype sequence")

        return mutation_position, wildtype_subsequence

    def determine_fasta_sequences(self):
        if self.position_out_of_bounds():
            return (None, None)

        mutation_start_position, wildtype_subsequence = self.get_wildtype_subsequence()
        mutation_end_position = mutation_start_position + self.wildtype_amino_acid_length
        mutant_subsequence = wildtype_subsequence[:mutation_start_position] + self.mutant_amino_acid + wildtype_subsequence[mutation_end_position:]

        if self.sequence_is_valid(wildtype_subsequence) and self.sequence_is_valid(mutant_subsequence):
            return (wildtype_subsequence, mutant_subsequence)
        else:
            return (None, None)

class InframeDeletionVariant(InframeVariant):
    def determine_wildtype_amino_acid(self):
        return self.amino_acid_change.split('/')[0]

    def determine_mutant_amino_acid(self):
        mutant_amino_acid = self.amino_acid_change.split('/')[1]
        if mutant_amino_acid == '-':
            return ''
        else:
            return mutant_amino_acid

    def determine_wildtype_amino_acid_length(self):
        return len(self.wildtype_amino_acid)

    def calculate_position(self):
        return int(self.protein_position.split('-', 1)[0]) - 1

class InframeInsertionVariant(InframeVariant):
    def determine_wildtype_amino_acid(self):
        return self.amino_acid_change.split('/')[0]

    def determine_mutant_amino_acid(self):
        return self.amino_acid_change.split('/')[1]

    def determine_wildtype_amino_acid_length(self):
        if self.wildtype_amino_acid == '-':
            return 0
        else:
            return len(self.wildtype_amino_acid)

    def calculate_position(self):
        if self.wildtype_amino_acid == '-':
            return int(self.protein_position.split('-', 1)[0])
        else:
            return int(self.protein_position) - 1

class MissenseVariant(InframeVariant):
    def determine_wildtype_amino_acid(self):
        return self.amino_acid_change.split('/')[0]

    def determine_mutant_amino_acid(self):
        return self.amino_acid_change.split('/')[1]

    def determine_wildtype_amino_acid_length(self):
        if self.wildtype_amino_acid == '-':
            return 0
        else:
            return len(self.wildtype_amino_acid)

    def calculate_position(self):
        if self.wildtype_amino_acid == '-':
            return int(self.protein_position.split('-', 1)[0])
        else:
            return int(self.protein_position) - 1

class FrameshiftVariant(SingleTranscriptVariant):
    def __init__(self, **kwargs):
        SingleTranscriptVariant.__init__(self, **kwargs)
        self.downstream_amino_acid_sequence = kwargs['downstream_amino_acid_sequence']
        self.downstream_sequence_length     = kwargs['downstream_sequence_length']
        self.position                       = self.calculate_position()

    def calculate_position(self):
        return int(self.protein_position.split('-', 1)[0]) - 1

    def get_subsequences(self):
        one_flanking_sequence_length = self.determine_flanking_sequence_length(len(self.wildtype_amino_acid_sequence))
        if self.position < one_flanking_sequence_length:
            start_position = 0
        else:
            start_position = self.position - one_flanking_sequence_length
        wildtype_subsequence_stop_position = self.position + one_flanking_sequence_length
        mutation_subsequence_stop_position = self.position
        wildtype_subsequence = self.wildtype_amino_acid_sequence[start_position:wildtype_subsequence_stop_position]
        mutation_start_subsequence = self.wildtype_amino_acid_sequence[start_position:mutation_subsequence_stop_position]
        return wildtype_subsequence, mutation_start_subsequence

    def determine_fasta_sequences(self):
        if self.position_out_of_bounds():
            return (None, None)

        wildtype_subsequence, mutant_subsequence = self.get_subsequences()
        downstream_sequence = self.downstream_amino_acid_sequence

        if self.downstream_sequence_length and len(downstream_sequence) > self.downstream_sequence_length:
            downstream_sequence = downstream_sequence[0:self.downstream_sequence_length]
        mutant_subsequence += downstream_sequence

        if self.sequence_is_valid(wildtype_subsequence) and self.sequence_is_valid(mutant_subsequence):
            return (wildtype_subsequence, mutant_subsequence)
        else:
            return (None, None)

class FusionVariant(Variant):
    def determine_fasta_sequences(self):
        position     = int(self.tsv_entry['fusion_position'])
        sequence     = self.tsv_entry['fusion_amino_acid_sequence']
        one_flanking_sequence_length = self.determine_flanking_sequence_length(len(sequence))
        if position < one_flanking_sequence_length:
            start_position = 0
        else:
            start_position = position - one_flanking_sequence_length

        if variant_type == 'inframe_fusion':
            stop_position = position + one_flanking_sequence_length
            subsequence   = sequence[start_position:stop_position]
        elif variant_type == 'frameshift_fusion':
            subsequence = sequence[start_position:]
            if subsequence.endswith('X'):
                subsequence = subsequence[:-1]
        else:
            return

        if '*' in subsequence:
            return

        if 'X' in subsequence:
            return

        if len(subsequence) < self.epitope_length:
            return

        if subsequence in fasta_sequences:
            fasta_sequences[subsequence].append(self.tsv_entry['index'])
        else:
            fasta_sequences[subsequence] = [self.tsv_entry['index']]

        writer                  = open(self.output_file, 'w')
        key_writer              = open(self.output_key_file, 'w')
        count                   = 1
        for (subsequence, keys) in fasta_sequences.items():
            writer.writelines('>%s\n' % count)
            writer.writelines('%s\n' % subsequence)
            yaml.dump({count: keys}, key_writer, default_flow_style=False)
            count += 1

        writer.close()
        key_writer.close()

