import sys
import os
import unittest
import tempfile
from filecmp import cmp
import py_compile
from pvacseq.lib.variant import *

class VariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'fasta_generator')
        cls.peptide_sequence_length = 21
        cls.epitope_length = 8

    def test_source_compiles(self):
        module = os.path.join(self.executable_dir, 'variant.py')
        self.assertTrue(py_compile.compile(module))

    def test_determine_peptide_sequence_length(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYK",
            amino_acid_change               = "R/H",
        )
        self.assertEqual(variant.determine_peptide_sequence_length(20), 20)
        self.assertEqual(variant.determine_peptide_sequence_length(21), 21)
        self.assertEqual(variant.determine_peptide_sequence_length(22), 21)

    def test_determine_flanking_sequence_length(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYK",
            amino_acid_change               = "R/H",
        )
        self.assertEqual(variant.determine_flanking_sequence_length(22), 10)
        self.assertEqual(variant.determine_flanking_sequence_length(21), 10)
        self.assertEqual(variant.determine_flanking_sequence_length(20), 9)
        self.assertEqual(variant.determine_flanking_sequence_length(19), 9)
        self.assertEqual(variant.determine_flanking_sequence_length(18), 8)

    def test_position_out_of_bounds(self):
        variant1 = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "5",
            wildtype_amino_acid_sequence    = "QELEAALHRD",
            amino_acid_change               = "R/H",
        )
        self.assertFalse(variant1.position_out_of_bounds())
        variant2 = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "10",
            wildtype_amino_acid_sequence    = "QELEAALHRD",
            amino_acid_change               = "R/H",
        )
        self.assertFalse(variant2.position_out_of_bounds())
        variant3 = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "11",
            wildtype_amino_acid_sequence    = "QELEAALHRD",
            amino_acid_change               = "R/H",
        )
        self.assertTrue(variant3.position_out_of_bounds())

    def test_position_out_of_bounds_generates_no_sequences(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYK",
            amino_acid_change               = "R/H",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, None)
        self.assertEqual(mt_sequence, None)

    def test_sequence_is_valid(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYK",
            amino_acid_change               = "R/H",
        )
        self.assertTrue(variant.sequence_is_valid("QELEAALHRD"))
        self.assertFalse(variant.sequence_is_valid("QXLEAALHRD"))
        self.assertFalse(variant.sequence_is_valid("Q*LEAALHRD"))
        self.assertTrue(variant.sequence_is_valid("QELEAALH"))
        self.assertFalse(variant.sequence_is_valid("QELEAAL"))

    def test_sequence_containing_asterisk_generates_no_sequences(self):
        variant = InframeDeletionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "112-114",
            wildtype_amino_acid_sequence    = "MEELDGEPTVTLIPGVNSKKNQMYFDWGPGEMLVCETSFNKKEKSEMVPSCPFIYIIRKDVDVYSQILRKLFNESHGIFLGLQRIDEELTGKSRKSQLVRVSKNYRSVIRACMEEMHQVAIAAKDPANGRQFSSQVSILSAMELIWNLCEILFIEVAPAGPLLLHLLDWVRLHVCEVDSLSADVLGSENPSKHDSFWNLVTILVLQGRLDEARQMLSKEADASPASAGICRIMGDLMRTMPILSPGNTQTLTELELKWQHWHEECERYLQDSTFATSPHLESLLKIMLGDEAALLEQKELLSNWYHFLVTRLLYSNPTVKPIDLHYYAQSSLDLFLGGESSPEPLDNILLAAFEFDIHQVIKECSIALSNWWFVAHLTDLLDHCKLLQSHNLYFGSNMREFLLLEYASGLFAHPSLWQLGVDYFDYCPELGRVSLELHIERIPLNTEQKALKVLRICEQRQMTEQVRSICKILAMKAVRNNRLGSALSWSIRAKDAAFATLVSDRFLRDYCERGCFSDLDLIDNLGPAMMLSDRLTFLGKYREFHRMYGEKRFADAASLLLSLMTSRIAPRSFWMTLLTDALPLLEQKQVIFSAEQTYELMRCLEDLTSRRPVHGESDTEQLQDDDIETTKVEMLRLSLARNLARAIIREGSLEGS",
            amino_acid_change               = "CME/*",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, None)
        self.assertEqual(mt_sequence, None)

    def test_sequence_containing_X_generates_no_sequences(self):
        variant = FrameshiftVariant(
            desired_peptide_sequence_length = 31,
            epitope_length                  = 8,
            protein_position                = "13",
            wildtype_amino_acid_sequence    = "XARLLMQRGRPKSDRLGKIRSLDGVESGVVARASSPSYSRG",
            downstream_amino_acid_sequence  = "SDRPAGEDPESGWCGVGSSGARF",
            downstream_sequence_length      = None,
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, None)
        self.assertEqual(mt_sequence, None)

    def test_resulting_short_fasta_sequence_generates_no_sequences(self):
        variant = FrameshiftVariant(
            desired_peptide_sequence_length = 31,
            epitope_length                  = 8,
            protein_position                = "1",
            wildtype_amino_acid_sequence    = "MANEVQDLLSPRKGGHPPAVKAGGMRISKKQEIGTLERHTKKTGFEKTSAIANVAKIQTLDALNDALEKLNYKFPATVHMAHQKPTPALEKVVPLKRIYIIQQPRKC",
            downstream_amino_acid_sequence  = "IGK",
            downstream_sequence_length      = None,
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, None)
        self.assertEqual(mt_sequence, None)

    def test_distance_from_start_works_as_expected(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 17,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYKVLDVVKAHKDSWPFLEPVDESYAPNYYQIIKAPMDISSMEKKLNGGLYCTKEEFVNDMKTMFRNCRKYNGESSEYTKMSDNLERCFHRAMMKHFPGEDGDTDEEFWIREDEKREKRRSRAGRSGGSHVWTRSRDPEGSSRKQQPMENGGKSLPPTRRAPSSGDDQSSSSTQPPREVGTSNGRGFSHPLHCGGTPSQAPFLNQMRPAVPGTFGPLRGSDPATLYGSSGVPEPHPGEPVQQRQPFTMQPPVGINSLRGPRLGTPEEKQMCGGLTHLSNMGPHPGSLQLGQISGPSQDGSMYAPAQFQPGFIPPRHGGAPARPPDFPESSEIPPSHMYRSYKYLNRVHSAVWNGNHGATNQGPLGPDEKPHLGPGPSHQPRTLGHVMDSRVMRPPVPPNQWTEQSGFLPHGVPSSGYMRPPCKSAGHRLQPPPVPAPSSLFGAPAQALRGVQGGDSMMDSPEMIAMQQLSSRVCPPGVPYHPHQPAHPRLPGPFPQVAHPMSVTVSAPKPALGNPGRAPENSEAQEPENDQAEPLPGLEEKPPGVGTSEGVYLTQLPHPTPPLQTDCTRQSSPQERETVGPELKSSSSESADNCKAMKGKNPWPSDSSYPGPAAQGCVRDLSTVADRGALSENGVIGEASPCGSEGKGLGSSGSEKLLCPRGRTLQETMPCTGQNAATPPSTDPGLTGGTVSQFPPLYMPGLEYPNSAAHYHISPGLQGVGPVMGGKSPASHPQHFPPRGFQSNHPHSGGFPRYRPPQGMRYSYHPPPQPSYHHYQRTPYYACPQSFSDWQRPLHPQGSPSGPPASQPPPPRSLFSDKNAMASLQGCETLNAALTSPTRMDAVAAKVPNDGQNPGPEEEKLDESMERPESPKEFLDLDNHNAATKRQSSLSASEYLYGTPPPLSSGMGFGSSAFPPHSVMLQTGPPYTPQRPASHFQPRAYSSPVAALPPHHPGATQPNGLSQEGPIYRCQEEGLGHFQAVMMEQIGTRSGIRGPFQEMYRPSGMQMHPVQSQASFPKTPTAATSQEEVPPHKPPTLPLDQS",
            amino_acid_change               = "R/H",
        )
        sequence = 'KKLKILGMPFRNIRSILKMVN'
        position = 5
        self.assertEqual(variant.distance_from_start(position, sequence), 5)

    def test_distance_from_end_works_as_expected(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 17,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYKVLDVVKAHKDSWPFLEPVDESYAPNYYQIIKAPMDISSMEKKLNGGLYCTKEEFVNDMKTMFRNCRKYNGESSEYTKMSDNLERCFHRAMMKHFPGEDGDTDEEFWIREDEKREKRRSRAGRSGGSHVWTRSRDPEGSSRKQQPMENGGKSLPPTRRAPSSGDDQSSSSTQPPREVGTSNGRGFSHPLHCGGTPSQAPFLNQMRPAVPGTFGPLRGSDPATLYGSSGVPEPHPGEPVQQRQPFTMQPPVGINSLRGPRLGTPEEKQMCGGLTHLSNMGPHPGSLQLGQISGPSQDGSMYAPAQFQPGFIPPRHGGAPARPPDFPESSEIPPSHMYRSYKYLNRVHSAVWNGNHGATNQGPLGPDEKPHLGPGPSHQPRTLGHVMDSRVMRPPVPPNQWTEQSGFLPHGVPSSGYMRPPCKSAGHRLQPPPVPAPSSLFGAPAQALRGVQGGDSMMDSPEMIAMQQLSSRVCPPGVPYHPHQPAHPRLPGPFPQVAHPMSVTVSAPKPALGNPGRAPENSEAQEPENDQAEPLPGLEEKPPGVGTSEGVYLTQLPHPTPPLQTDCTRQSSPQERETVGPELKSSSSESADNCKAMKGKNPWPSDSSYPGPAAQGCVRDLSTVADRGALSENGVIGEASPCGSEGKGLGSSGSEKLLCPRGRTLQETMPCTGQNAATPPSTDPGLTGGTVSQFPPLYMPGLEYPNSAAHYHISPGLQGVGPVMGGKSPASHPQHFPPRGFQSNHPHSGGFPRYRPPQGMRYSYHPPPQPSYHHYQRTPYYACPQSFSDWQRPLHPQGSPSGPPASQPPPPRSLFSDKNAMASLQGCETLNAALTSPTRMDAVAAKVPNDGQNPGPEEEKLDESMERPESPKEFLDLDNHNAATKRQSSLSASEYLYGTPPPLSSGMGFGSSAFPPHSVMLQTGPPYTPQRPASHFQPRAYSSPVAALPPHHPGATQPNGLSQEGPIYRCQEEGLGHFQAVMMEQIGTRSGIRGPFQEMYRPSGMQMHPVQSQASFPKTPTAATSQEEVPPHKPPTLPLDQS",
            amino_acid_change               = "R/H",
        )
        sequence = 'KKLKILGMPFRNIRSILKMVN'
        position = 5
        self.assertEqual(variant.distance_from_end(position, sequence), 15)

class MissenseVariantTests(unittest.TestCase):
    def test_peptide_sequence_length_17(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 17,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYKVLDVVKAHKDSWPFLEPVDESYAPNYYQIIKAPMDISSMEKKLNGGLYCTKEEFVNDMKTMFRNCRKYNGESSEYTKMSDNLERCFHRAMMKHFPGEDGDTDEEFWIREDEKREKRRSRAGRSGGSHVWTRSRDPEGSSRKQQPMENGGKSLPPTRRAPSSGDDQSSSSTQPPREVGTSNGRGFSHPLHCGGTPSQAPFLNQMRPAVPGTFGPLRGSDPATLYGSSGVPEPHPGEPVQQRQPFTMQPPVGINSLRGPRLGTPEEKQMCGGLTHLSNMGPHPGSLQLGQISGPSQDGSMYAPAQFQPGFIPPRHGGAPARPPDFPESSEIPPSHMYRSYKYLNRVHSAVWNGNHGATNQGPLGPDEKPHLGPGPSHQPRTLGHVMDSRVMRPPVPPNQWTEQSGFLPHGVPSSGYMRPPCKSAGHRLQPPPVPAPSSLFGAPAQALRGVQGGDSMMDSPEMIAMQQLSSRVCPPGVPYHPHQPAHPRLPGPFPQVAHPMSVTVSAPKPALGNPGRAPENSEAQEPENDQAEPLPGLEEKPPGVGTSEGVYLTQLPHPTPPLQTDCTRQSSPQERETVGPELKSSSSESADNCKAMKGKNPWPSDSSYPGPAAQGCVRDLSTVADRGALSENGVIGEASPCGSEGKGLGSSGSEKLLCPRGRTLQETMPCTGQNAATPPSTDPGLTGGTVSQFPPLYMPGLEYPNSAAHYHISPGLQGVGPVMGGKSPASHPQHFPPRGFQSNHPHSGGFPRYRPPQGMRYSYHPPPQPSYHHYQRTPYYACPQSFSDWQRPLHPQGSPSGPPASQPPPPRSLFSDKNAMASLQGCETLNAALTSPTRMDAVAAKVPNDGQNPGPEEEKLDESMERPESPKEFLDLDNHNAATKRQSSLSASEYLYGTPPPLSSGMGFGSSAFPPHSVMLQTGPPYTPQRPASHFQPRAYSSPVAALPPHHPGATQPNGLSQEGPIYRCQEEGLGHFQAVMMEQIGTRSGIRGPFQEMYRPSGMQMHPVQSQASFPKTPTAATSQEEVPPHKPPTLPLDQS",
            amino_acid_change               = "R/H",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "SGGSHVWTRSRDPEGSS")
        self.assertEqual(mt_sequence, "SGGSHVWTHSRDPEGSS")

    def test_peptide_sequence_length_21(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYKVLDVVKAHKDSWPFLEPVDESYAPNYYQIIKAPMDISSMEKKLNGGLYCTKEEFVNDMKTMFRNCRKYNGESSEYTKMSDNLERCFHRAMMKHFPGEDGDTDEEFWIREDEKREKRRSRAGRSGGSHVWTRSRDPEGSSRKQQPMENGGKSLPPTRRAPSSGDDQSSSSTQPPREVGTSNGRGFSHPLHCGGTPSQAPFLNQMRPAVPGTFGPLRGSDPATLYGSSGVPEPHPGEPVQQRQPFTMQPPVGINSLRGPRLGTPEEKQMCGGLTHLSNMGPHPGSLQLGQISGPSQDGSMYAPAQFQPGFIPPRHGGAPARPPDFPESSEIPPSHMYRSYKYLNRVHSAVWNGNHGATNQGPLGPDEKPHLGPGPSHQPRTLGHVMDSRVMRPPVPPNQWTEQSGFLPHGVPSSGYMRPPCKSAGHRLQPPPVPAPSSLFGAPAQALRGVQGGDSMMDSPEMIAMQQLSSRVCPPGVPYHPHQPAHPRLPGPFPQVAHPMSVTVSAPKPALGNPGRAPENSEAQEPENDQAEPLPGLEEKPPGVGTSEGVYLTQLPHPTPPLQTDCTRQSSPQERETVGPELKSSSSESADNCKAMKGKNPWPSDSSYPGPAAQGCVRDLSTVADRGALSENGVIGEASPCGSEGKGLGSSGSEKLLCPRGRTLQETMPCTGQNAATPPSTDPGLTGGTVSQFPPLYMPGLEYPNSAAHYHISPGLQGVGPVMGGKSPASHPQHFPPRGFQSNHPHSGGFPRYRPPQGMRYSYHPPPQPSYHHYQRTPYYACPQSFSDWQRPLHPQGSPSGPPASQPPPPRSLFSDKNAMASLQGCETLNAALTSPTRMDAVAAKVPNDGQNPGPEEEKLDESMERPESPKEFLDLDNHNAATKRQSSLSASEYLYGTPPPLSSGMGFGSSAFPPHSVMLQTGPPYTPQRPASHFQPRAYSSPVAALPPHHPGATQPNGLSQEGPIYRCQEEGLGHFQAVMMEQIGTRSGIRGPFQEMYRPSGMQMHPVQSQASFPKTPTAATSQEEVPPHKPPTLPLDQS",
            amino_acid_change               = "R/H",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "GRSGGSHVWTRSRDPEGSSRK")
        self.assertEqual(mt_sequence, "GRSGGSHVWTHSRDPEGSSRK")

    def test_peptide_sequence_length_31(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 31,
            epitope_length                  = 8,
            protein_position                = "535",
            wildtype_amino_acid_sequence    = "QELEAALHRDDVEFISDLIACLLQGCYQRRDITPQTFHSYLEDIINYRWELEEGKPNPLREASFQDLPLRTRVEILHRLCDYRLDADDVFDLLKGLDADSLRVEPLGEDNSGALYWYFYGTRMYKEDPVQGKSNGELSLSRESEGQKNVSSIPGKTGKRRGRPPKRKKLQEEILLSEKQEENSLASEPQTRHGSQGPGQGTWWLLCQTEEEWRQVTESFRERTSLRERQLYKLLSEDFLPEICNMIAQKGKRPQRTKAELHPRWMSDHLSIKPVKQEETPVLTRIEKQKRKEEEEERQILLAVQKKEQEQMLKEERKRELEEKVKAVEGMCSVRVVWRGACLSTSRPVDRAKRRKLREERAWLLAQGKELPPELSHLDPNSPMREEKKTKDLFELDDDFTAMYKVLDVVKAHKDSWPFLEPVDESYAPNYYQIIKAPMDISSMEKKLNGGLYCTKEEFVNDMKTMFRNCRKYNGESSEYTKMSDNLERCFHRAMMKHFPGEDGDTDEEFWIREDEKREKRRSRAGRSGGSHVWTRSRDPEGSSRKQQPMENGGKSLPPTRRAPSSGDDQSSSSTQPPREVGTSNGRGFSHPLHCGGTPSQAPFLNQMRPAVPGTFGPLRGSDPATLYGSSGVPEPHPGEPVQQRQPFTMQPPVGINSLRGPRLGTPEEKQMCGGLTHLSNMGPHPGSLQLGQISGPSQDGSMYAPAQFQPGFIPPRHGGAPARPPDFPESSEIPPSHMYRSYKYLNRVHSAVWNGNHGATNQGPLGPDEKPHLGPGPSHQPRTLGHVMDSRVMRPPVPPNQWTEQSGFLPHGVPSSGYMRPPCKSAGHRLQPPPVPAPSSLFGAPAQALRGVQGGDSMMDSPEMIAMQQLSSRVCPPGVPYHPHQPAHPRLPGPFPQVAHPMSVTVSAPKPALGNPGRAPENSEAQEPENDQAEPLPGLEEKPPGVGTSEGVYLTQLPHPTPPLQTDCTRQSSPQERETVGPELKSSSSESADNCKAMKGKNPWPSDSSYPGPAAQGCVRDLSTVADRGALSENGVIGEASPCGSEGKGLGSSGSEKLLCPRGRTLQETMPCTGQNAATPPSTDPGLTGGTVSQFPPLYMPGLEYPNSAAHYHISPGLQGVGPVMGGKSPASHPQHFPPRGFQSNHPHSGGFPRYRPPQGMRYSYHPPPQPSYHHYQRTPYYACPQSFSDWQRPLHPQGSPSGPPASQPPPPRSLFSDKNAMASLQGCETLNAALTSPTRMDAVAAKVPNDGQNPGPEEEKLDESMERPESPKEFLDLDNHNAATKRQSSLSASEYLYGTPPPLSSGMGFGSSAFPPHSVMLQTGPPYTPQRPASHFQPRAYSSPVAALPPHHPGATQPNGLSQEGPIYRCQEEGLGHFQAVMMEQIGTRSGIRGPFQEMYRPSGMQMHPVQSQASFPKTPTAATSQEEVPPHKPPTLPLDQS",
            amino_acid_change               = "R/H",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "RRSRAGRSGGSHVWTRSRDPEGSSRKQQPME")
        self.assertEqual(mt_sequence, "RRSRAGRSGGSHVWTHSRDPEGSSRKQQPME")

    def test_mutation_at_relative_end_of_full_sequence(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "132",
            wildtype_amino_acid_sequence    = "MTGELEVKNMDMKPGSTLKITGSIADGTDGFVINLGQGTDKLNLHFNPRFSESTIVCNSLDGSNWGQEQREDHLCFSPGSEVKFTVTFESDKFKVKLPDGHELTFPNRLGHSHLSYLSVRGGFNMSSFKLKE",
            amino_acid_change               = "E/Q",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "SHLSYLSVRGGFNMSSFKLKE")
        self.assertEqual(mt_sequence, "SHLSYLSVRGGFNMSSFKLKQ")

    def test_mutation_at_relative_beginning_of_full_sequence(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "8",
            wildtype_amino_acid_sequence    = "MSRLKRIAGQDLRAGFKAGGRDCGTSVPQGLLKAARKSGQLNLSGRNLSEVPQCVWRINVDIPEEANQNLSFGATERWWEQTDLTKLIISNNKLQSLTDDLRLLPALTVLDIHDNQLTSLPSAIRELENLQKLNVSHNKLKILPEEITNLRNLKCLYLQHNELTCISEGFEQLSNLEDLDLSNNHLTTVPASFSSLSSLVRLNLSSNELKSLPAEINRMKRLKHLDCNSNLLETIPPELAGMESLELLYLRRNKLRFLPEFPSCSLLKELHVGENQIEMLEAEHLKHLNSILVLDLRDNKLKSVPDEIILLRSLERLDLSNNDISSLPYSLGNLHLKFLALEGNPLRTIRREIISKGTQEVLKYLRSKIKDDGPSQSESATETAMTLPSESRVNIHAIITLKILDYSDKQATLIPDEVFDAVKSNIVTSINFSKNQLCEIPKRMVELKEMVSDVDLSFNKLSFISLELCVLQKLTFLDLRNNFLNSLPEEMESLVRLQTINLSFNRFKMLPEVLYRIFTLETILISNNQVGSVDPQKMKMMENLTTLDLQNNDLLQIPPELGNCVNLRTLLLDGNPFRVPRAAILMKGTAAILEYLRDRIPT",
            amino_acid_change               = "A/V",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "MSRLKRIAGQDLRAGFKAGGR")
        self.assertEqual(mt_sequence, "MSRLKRIVGQDLRAGFKAGGR")

    def test_wildtype_sequence_shorter_than_desired_peptide_sequence_length(self):
        variant = MissenseVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "17",
            wildtype_amino_acid_sequence    = "MSANGAVWGRVRSRLRAFP",
            amino_acid_change               = "A/T",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "MSANGAVWGRVRSRLRAFP")
        self.assertEqual(mt_sequence, "MSANGAVWGRVRSRLRTFP")

class InframeInsertionVariantTests(unittest.TestCase):
    def test_amino_acid_replacement(self):
        variant = InframeInsertionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "20",
            wildtype_amino_acid_sequence    = "MLPRVGCPALPLPPPPLLPLLLLLLGASGGGGGARAEVLFRCPPCTPERLAACGPPPVAPPAAVAAVAGGARMPCAELVREPGCGCCSVCARLEGEACGVYTPRCGQGLRCYPHPGSELPLQALVMGEGTCEKRRDAEYGASPEQVADNGDDHSEGGLVENHVDSTMNMLGGGGSAGRKPLKSGMKELAVFREKVTEQHRQMGKGGKHHLGLEEPKKLRPPPARTPCQQELDQVLERISTMRLPDERGPLEHLYSLHIPNCDKHGLYNLKQCKMSLNGQRGECWCVNPNTGKLIQGAPTIRGDPECHLFYNEQQEARGVHTQRMQ",
            amino_acid_change               = "L/LLP",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "LPLPPPPLLPLLLLLLGASGG")
        self.assertEqual(mt_sequence, "LPLPPPPLLPLLPLLLLLGASGG")

    def test_amino_acid_insertion(self):
        variant = InframeInsertionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "287-288",
            wildtype_amino_acid_sequence    = "MSVQNSGWPHQEDSPKPQDPGPPANSDSDSGHLPGEDPEDTHAQGPAVLSLGSLCLDTNQAPNWTGLQTLLQQLPPQDIDERYCLALGEEERAELQLFCARRKQEALGQGVARLVLPKLEGHTCEKCRELLKPGEYGVFAARAGEQRCWHQPCFACQACGQALINLIYFYHDGQLYCGRHHAELLRPRCPACDQLIFSWRCTEAEGQRWHENHFCCQDCAGPLGGGRYALPGGSPCCPSCFENRYSDAGSSWAGALEGQAFLGETGLDRTEGRDQTSVNSATLSRTLLAAAGGSSLQTQRGLPGSSPQQENRPGDKAEAPKGQEQCRLETIRDPKDTPFSTCSSSSDSEPEGFFLGERLPQSWKTPGSLQAEDSNASKTHCTMC",
            amino_acid_change               = "-/L",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "VNSATLSRTLLAAAGGSSLQ")
        self.assertEqual(mt_sequence, "VNSATLSRTLLLAAAGGSSLQ")

class InframeDeletionVariantTests(unittest.TestCase):
    def test_amino_acid_replacement(self):
        variant = InframeDeletionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "495-502",
            wildtype_amino_acid_sequence    = "MTAEDSTAAMSSDSAAGSSAKVPEGVAGAPNEAALLALMERTGYSMVQENGQRKYGGPPPGWEGPHPQRGCEVFVGKIPRDVYEDELVPVFEAVGRIYELRLMMDFDGKNRGYAFVMYCHKHEAKRAVRELNNYEIRPGRLLGVCCSVDNCRLFIGGIPKMKKREEILEEIAKVTEGVLDVIVYASAADKMKNRGFAFVEYESHRAAAMARRKLMPGRIQLWGHQIAVDWAEPEIDVDEDVMETVKILYVRNLMIETTEDTIKKSFGQFNPGCVERVKKIRDYAFVHFTSREDAVHAMNNLNGTELEGSCLEVTLAKPVDKEQYSRYQKAARGGGAAEAAQQPSYVYSCDPYTLAYYGYPYNALIGPNRDYFVKAGSIRGRGRGAAGNRAPGPRGSYLGGYSAGRGIYSRYHEGKGKQQEKGYELVPNLEIPTVNPVAIKPGTVAIPAIGAQYSMFPAAPAPKMIEDGKIHTVEHMISPIAVQPDPASAAAAAAAAAAAAAAVIPTVSTPPPFQGRPITPVYTVAPNVQRIPTAGIYGASYVPFAAPATATIATLQKNAAAAAAMYGGYAGYIPQAFPAAAIQVPIPDVYQTY",
            amino_acid_change               = "AAAAAAAA/A",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "DPASAAAAAAAAAAAAAAVIPTVSTPPP")
        self.assertEqual(mt_sequence, "DPASAAAAAAAVIPTVSTPPP")

    def test_amino_acid_deletion(self):
        variant = InframeDeletionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "750",
            wildtype_amino_acid_sequence    = "MQCLAAALKDETNMSGGGEQADILPANYVVKDRWKVLKKIGGGGFGEIYEAMDLLTRENVALKVESAQQPKQVLKMEVAVLKKLQGKDHVCRFIGCGRNEKFNYVVMQLQGRNLADLRRSQPRGTFTLSTTLRLGKQILESIEAIHSVGFLHRDIKPSNFAMGRLPSTYRKCYMLDFGLARQYTNTTGDVRPPRNVAGFRGTVRYASVNAHKNREMGRHDDLWSLFYMLVEFAVGQLPWRKIKDKEQVGMIKEKYEHRMLLKHMPSEFHLFLDHIASLDYFTKPDYQLIMSVFENSMKERGIAENEAFDWEKAGTDALLSTSTSTPPQQNTRQTAAMFGVVNVTPVPGDLLRENTEDVLQGEHLSDQENAPPILPGRPSEGLGPSPHLVPHPGGPEAEVWEETDVNRNKLRINIGKSPCVEEEQSRGMGVPSSPVRAPPDSPTTPVRSLRYRRVNSPESERLSTADGRVELPERRSRMDLPGSPSRQACSSQPAQMLSVDTGHADRQASGRMDVSASVEQEALSNAFRSVPLAEEEDFDSKEWVIIDKETELKDFPPGAEPSTSGTTDEEPEELRPLPEEGEERRRLGAEPTVRPRGRSMQALAEEDLQHLPPQPLPPQLSQGDGRSETSQPPTPGSPSHSPLHSGPRPRRRESDPTGPQRQVFSVAPPFEVNGLPRAVPLSLPYQDFKRDLSDYRERARLLNRVRRVGFSHMLLTTPQVPLAPVQPQANGKEEEEEEEEDEEEEEEDEEEEEEEEEEEEEEEEEEEEEEEAAAAVALGEVLGPRSGSSSEGSERSTDRSQEGAPSTLLADDQKESRGRASMADGDLEPEEGSKTLVLVSPGDMKKSPVTAELAPDPDLGTLAALTPQHERPQPTGSQLDVSEPGTLSSVLKSEPKPPGPGAGLGAGTVTTGVGGVAVTSSPFTKVERTFVHIAEKTHLNVMSSGGQALRSEEFSAGGELGLELASDGGAVEEGARAPLENGLALSGLNGAEIEGSALSGAPRETPSEMATNSLPNGPALADGPAPVSPLEPSPEKVATISPRRHAMPGSRPRSRIPVLLSEEDTGSEPSGSLSAKERWSKRARPQQDLARLVMEKRQGRLLLRLASGASSSSSEEQRRASETLSGTGSEEDTPASEPAAALPRKSGRAAATRSRIPRPIGLRMPMPVAAQQPASRSHGAAPALDTAITSRLQLQTPPGSATAADLRPKQPPGRGLGPGRAQAGARPPAPRSPRLPASTSAARNASASPRSQSLSRRESPSPSHQARPGVPPPRGVPPARAQPDGTPSPGGSKKGPRGKLQAQRATTKGRAGGAEGRAGAR",
            amino_acid_change               = "E/-",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "EDEEEEEEDEEEEEEEEEEEE")
        self.assertEqual(mt_sequence, "EDEEEEEEDEEEEEEEEEEE")

    def test_range(self):
        variant = InframeDeletionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "771-772",
            wildtype_amino_acid_sequence    = "MGDFAAPAAAANGSSICINSSLNSSLGGAGIGVNNTPNSTPAAPSSNHPAAGGCGGSGGPGGGSAAVPKHSTVVERLRQRIEGCRRHHVNCENRYQQAQVEQLELERRDTVSLYQRTLEQRAKKSGAGTGKQQHPSKPQQDAEAASAEQRNHTLIMLQETVKRKLEGARSPLNGDQQNGACDGNFSPTSKRIRKDISAGMEAINNLPSNMPLPSASPLHQLDLKPSLPLQNSGTHTPGLLEDLSKNGRLPEIKLPVNGCSDLEDSFTILQSKDLKQEPLDDPTCIDTSETSLSNQNKLFSDINLNDQEWQELIDELANTVPEDDIQDLFNEDFEEKKEPEFSQPATETPLSQESASVKSDPSHSPFAHVSMGSPQARPSSSGPPFSTVSTATSLPSVASTPAAPNPASSPANCAVQSPQTPNQAHTPGQAPPRPGNGYLLNPAAVTVAGSASGPVAVPSSDMSPAEQLKQMAAQQQQRAKLMQQKQQQQQQQQQQQQQQQQQQQQQQQQQHSNQTSNWSPLGPPSSPYGAAFTAEKPNSPMMYPQAFNNQNPIVPPMANNLQKTTMNNYLPQNHMNMINQQPNNLGTNSLNKQHNILTYGNTKPLTHFNADLSQRMTPPVANPNKNPLMPYIQQQQQQQQQQQQQQQQQQPPPPQLQAPRAHLSEDQKRLLLMKQKGVMNQPMAYAALPSHGQEQHPVGLPRTTGPMQSSVPPGSGGMVSGASPAGPGFLGSQPQAAIMKQMLIDQRAQLIEQQKQQFLREQRQQQQQQQQQILAEQQLQQSHLPRQHLQPQRNPYPVQQVNQFQGSPQDIAAVRSQAALQSMRTSRLMAQNAGMMGIGPSQNPGTMATAAAQSEMGLAPYSTTPTSQPGMYNMSTGMTQMLQHPNQSGMSITHNQAQGPRQPASGQGVGMVSGFGQSMLVNSAITQQHPQMKGPVGQALPRPQAPPRLQSLMGTVQQGAQSWQQRSLQGMPGRTSGELGPFNNGASYPLQAGQPRLTKQHFPQGLSQSVVDANTGTVRTLNPAAMGRQMMPSLPGQQGTSQARPMVMSGLSQGVPGMPAFSQPPAQQQIPSGSFAPSSQSQAYERNAPQDVSYNYSGDGAGGSFPGLPDGADLVDSIIKGGPGDEWMQELDELFGNP",
            amino_acid_change               = "QQ/-",
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "EQRQQQQQQQQQILAEQQLQQS")
        self.assertEqual(mt_sequence, "EQRQQQQQQQILAEQQLQQS")

class FrameshiftVariantTests(unittest.TestCase):
    def test_feature_truncation(self):
        variant = FrameshiftVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "342",
            wildtype_amino_acid_sequence    = "MMSIKAFTLVSAVERELLMGDKERVNIECVECCGRDLYVGTNDCFVYHFLLEERPVPAGPATFTATKQLQRHLGFKKPVNELRAASALNRLLVLCDNSISLVNMLNLEPVPSGARIKGAATFALNENPVSGDPFCVEVCIISVKRRTIQMFLVYEDRVQIVKEVSTAEQPLAVAVDGHFLCLALTTQYIIHNYSTGVSQDLFPYCSEERPPIVKRIGRQEFLLAGPGGLGMFATVAGISQRAPVHWSENVIGAAVSFPYVIALDDEFITVHSMLDQQQKQTLPFKEGHILQDFEGRVIVATSKGVYILVPLPLEKQIQDLLASRRVEEALVLAKGARRNIPKEKFQVMYRRILQQAGFIQFAQLQFLEAKELFRSGQLDVRELISLYPFLLPTSSSFTRSHPPLHEYADLNQLTQGDQEKMAKCKRFLMSYLNEVRSTEVANGYKEDIDTALLKLYAEADHDSLLDLLVTENFCLLTDSAAWLEKHKKYFALGLLYHYNNQDAAAVQLWVNIVNGDVQDSTRSDLYEYIVDFLTYCLDEELVWAYADWVLQKSEEVGVQVFTKRPLDEQQKNSFNPDDIINCLKKYPKALVKYLEHLVIDKRLQKEEYHTHLAVLYLEEVLLQRASASGKGAEATETQAKLRRLLQKSDLYRVHFLLERLQGAGLPMESAILHGKLGEHEKALHILVHELQDFAAAEDYCLWCSEGRDPPHRQQLFHTLLAIYLHAGPTAHELAVAAVDLLNRHATEFDAAQVLQMLPDTWSVQLLCPFLMGAMRDSIHARRTMQVALGLARSENLIYTYDKMKLKGSSIQLSDKKLCQICQNPFCEPVFVRYPNGGLVHTHCAASRHTNPSSSSPGTRT",
            downstream_amino_acid_sequence  = "RKNFR",
            downstream_sequence_length      = None,
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "LAKGARRNIPKEKFQVMYRR")
        self.assertEqual(mt_sequence, "LAKGARRNIPRKNFR")

    def test_feature_truncation_range(self):
        variant = FrameshiftVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "175-180",
            wildtype_amino_acid_sequence    = "SSAMTSTSLASKLTTLFSTGQAARSGSSSSPISLSTEKETSFLSPTASTSRKTSLFLGPSMARQPNILVHLQTSALTLSPTSTLNMSQEEPPELTSSQTIAEEEGTTAETQTLTFTPSETPTSLLPVSSPTEPTARRKSSPETWASSISVPAKTSLVETTDGTLVTTIKMSSQAAQGNSTWPAPAEETGSSPAGTSPGSPEMSTTLKIMSSKEPSISPEIRSTVRNSPW",
            downstream_amino_acid_sequence  = "VACPSRGDGEQSSRHIPRKPRNVYHSQNHELQGTQHQPRDQVHCEKFSLEDSRNNCSHGDHSGTSHPSVHSPRKWQHQHLSPAHRNHITNQVTNRKYVGYRKGLPLPIPT",
            downstream_sequence_length      = 100,
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "VTTIKMSSQAAQGNSTWPAP")
        self.assertEqual(mt_sequence, "VTTIKMSSQAVACPSRGDGEQSSRHIPRKPRNVYHSQNHELQGTQHQPRDQVHCEKFSLEDSRNNCSHGDHSGTSHPSVHSPRKWQHQHLSPAHRNHITNQVTNRKYVGY")

    def test_feature_elongation(self):
        variant = FrameshiftVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "142",
            wildtype_amino_acid_sequence    = "MLYSPGPSLPESAESLDGSQEDKPRGSCAEPTFTDTGMVAHINNSRLKAKGVGQHDNAQNFGNQSFEELRAACLRKGELFEDPLFPAEPSSLGFKDLGPNSKNVQNISWQRPKDIINNPLFIMDGISPTDICQGILGDCWLLAAIGSLTTCPKLLYRVVPRGQSFKKNYAGIFHFQIWQFGQWVNVVVDDRLPTKNDKLVFVHSTERSEFWSALLEKAYAKLSGSYEALSGGSTMEGLEDFTGGVAQSFQLQRPPQNLLRLLRKAVERSSLMGCSIEVTSDSELESMTDKMLVRGHAYSVTGLQDVHYRGKMETLIRVRNPWGRIEWNGAWSDSAREWEEVASDIQMQLLHKTEDGEFWMSYQDFLNNFTLLEICNLTPDTLSGDYKSYWHTTFYEGSWRRGSSAGGCRNHPGTFWTNPQFKISLPEGDDPEDDAEGNVVVCTCLVALMQKNWRHARQQGAQLQTIGFVLYAVPKEFQNIQDVHLKKEFFTKYQDHGFSEIFTNSREVSSQLRLPPGEYIIIPSTFEPHRDADFLLRVFTEKHSESWELDEVNYAEQLQEEKVSEDDMDQDFLHLFKIVAGEGKEIGVYELQRLLNRMAIKFKSFKTKGFGLDACRCMINLMDKDGSGKLGLLEFKILWKKLKKWMDIFRECDQDHSGTLNSYEMRLVIEKAGIKLNNKVMQVLVARYADDDLIIDFDSFISCFLRLKTMFTFFLTMDPKNTGHICLSLEQWLQMTMWG",
            downstream_amino_acid_sequence  = "LAAGCHRLPYHLPQTAIPRGAQRTELQEKLCWHLPFSDLAVWTVGERGGR",
            downstream_sequence_length      = None,
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "CQGILGDCWLLAAIGSLTTC")
        self.assertEqual(mt_sequence, "CQGILGDCWLLAAGCHRLPYHLPQTAIPRGAQRTELQEKLCWHLPFSDLAVWTVGERGGR")

    def test_range(self):
        variant = FrameshiftVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "162-168",
            wildtype_amino_acid_sequence    = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSRPRQVDHLRSGVQDSLANIAKSHLY",
            downstream_amino_acid_sequence  = "T",
            downstream_sequence_length      = None,
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "PPGTRVRAMAIYKQSQHMTE")
        self.assertEqual(mt_sequence, "PPGTRVRAMAT")

    def test_downstream_sequence_length_limit(self):
        variant = FrameshiftVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            protein_position                = "142",
            wildtype_amino_acid_sequence    = "MLYSPGPSLPESAESLDGSQEDKPRGSCAEPTFTDTGMVAHINNSRLKAKGVGQHDNAQNFGNQSFEELRAACLRKGELFEDPLFPAEPSSLGFKDLGPNSKNVQNISWQRPKDIINNPLFIMDGISPTDICQGILGDCWLLAAIGSLTTCPKLLYRVVPRGQSFKKNYAGIFHFQIWQFGQWVNVVVDDRLPTKNDKLVFVHSTERSEFWSALLEKAYAKLSGSYEALSGGSTMEGLEDFTGGVAQSFQLQRPPQNLLRLLRKAVERSSLMGCSIEVTSDSELESMTDKMLVRGHAYSVTGLQDVHYRGKMETLIRVRNPWGRIEWNGAWSDSAREWEEVASDIQMQLLHKTEDGEFWMSYQDFLNNFTLLEICNLTPDTLSGDYKSYWHTTFYEGSWRRGSSAGGCRNHPGTFWTNPQFKISLPEGDDPEDDAEGNVVVCTCLVALMQKNWRHARQQGAQLQTIGFVLYAVPKEFQNIQDVHLKKEFFTKYQDHGFSEIFTNSREVSSQLRLPPGEYIIIPSTFEPHRDADFLLRVFTEKHSESWELDEVNYAEQLQEEKVSEDDMDQDFLHLFKIVAGEGKEIGVYELQRLLNRMAIKFKSFKTKGFGLDACRCMINLMDKDGSGKLGLLEFKILWKKLKKWMDIFRECDQDHSGTLNSYEMRLVIEKAGIKLNNKVMQVLVARYADDDLIIDFDSFISCFLRLKTMFTFFLTMDPKNTGHICLSLEQWLQMTMWG",
            downstream_amino_acid_sequence  = "LAAGCHRLPYHLPQTAIPRGAQRTELQEKLCWHLPFSDLAVWTVGERGGR",
            downstream_sequence_length      = 20,
        )
        wt_sequence, mt_sequence = variant.determine_fasta_sequences()
        self.assertEqual(wt_sequence, "CQGILGDCWLLAAIGSLTTC")
        self.assertEqual(mt_sequence, "CQGILGDCWLLAAGCHRLPYHLPQTAIPRG")

class InframeFusionVariantTests(unittest.TestCase):
    def test(self):
        variant = InframeFusionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            fusion_position                 = "22",
            fusion_amino_acid_sequence      = "QALDENMDLLEGITGFEDSVRKSSIPKNVFLALHEKLYIMLKGKMGTVNLHQFTGQLTEELHEQLENLGTHGTMDLNNLV",
        )
        sequence = variant.determine_fasta_sequences()
        self.assertEqual(sequence, "ITGFEDSVRKSSIPKNVFLA")

class FrameshiftFusionVariantTests(unittest.TestCase):
    def test(self):
        variant = FrameshiftFusionVariant(
            desired_peptide_sequence_length = 21,
            epitope_length                  = 8,
            fusion_position                 = "37",
            fusion_amino_acid_sequence      = "LLQAATNYHNGHTGQLSAITVFLLFGGSLARIFTSIQKTPLHAQLHPHELVCFFHPENPGCTGEGRRLLQLLLQEAX",
        )
        sequence = variant.determine_fasta_sequences()
        self.assertEqual(sequence, "SLARIFTSIQKTPLHAQLHPHELVCFFHPENPGCTGEGRRLLQLLLQEA")

#Test for fusion position at relative beginning or end of fusion sequence

if __name__ == '__main__':
    unittest.main()
