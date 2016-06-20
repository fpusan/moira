#NOTE THAT THIS IS MEANT TO BE IMPORTED AS A MODULE FROM THE PARENT DIRECTORY

import unittest

import moira

bernoulliOK = False
try:
    import bernoulli
    bernoulliOK = True
except ImportError:
    print 'bernoulli module not found.'
    
nwOK = False
try:
    import nw_align
    nwOK = True
except ImportError:
    print 'nw_align module not found'


import gzip
import bz2
import os

try: 
    os.makedirs('test/test_output')
except OSError:
    if not os.path.isdir('test/test_output'):
        raise


class Arguments():
    def __init__(self, **kwargs):     
        for key, value in kwargs.items():
              setattr(self, key, value)

class TestErrorCalculationAlgorithms(unittest.TestCase):
    def testPbPython(self):
        self.assertEqual(moira.calculate_errors_PB(testSeq1, testQual1, 0.005), (6.446879136706666, 0))
    @unittest.skipIf(not bernoulliOK, 'Module not present')
    def testPbC(self):
        self.assertEqual(bernoulli.calculate_errors_PB(testSeq1, testQual1, 0.005), (6.446879136706666, 0))
    def testPoisson(self):
        self.assertEqual(moira.calculate_errors_poisson(testSeq1, testQual1, 0.005), (6.932519986616133, 0))
        

class TestContigConstructor(unittest.TestCase):
    def testReverseComplement(self):
        self.assertEqual(moira.reverse_complement(testSeq2, testQual2), testRC2)
    def testNwPython(self):
        self.assertEqual(moira.nw_align(testSeq1, moira.reverse_complement(testSeq2), args.match, args.mismatch, args.gap), test_aligned)
    @unittest.skipIf(not nwOK, 'Module not present')
    def testNwC(self):
        self.assertEqual(nw_align.nw_align(testSeq1, moira.reverse_complement(testSeq2), args.match, args.mismatch, args.gap), test_aligned)
    def testContig(self):
        testSeq2RC, testQual2RC = moira.reverse_complement(testSeq2, testQual2)
        aligned1, aligned2 = moira.nw_align(testSeq1, testSeq2RC, args.match, args.mismatch, args.gap)[:2]
        self.assertEqual(moira.make_contig(aligned1, testQual1, aligned2, testQual2RC, args.insert, args.deltaq, args.consensus_qscore, args.qscore_cap, args.trim_overlap), test_contig)

        
class TestProcessing(unittest.TestCase):
    def testProcessForwardSeq(self):
        args.truncate = 200
        args.paired = False
        self.assertEqual(moira.process_data('foo', testSeq1, testQual1, None, None, args), test_ForwardProcess)
    def testProcessPairedSeq(self):
        args.truncate = 200
        args.paired = True
        self.assertEqual(moira.process_data('foo', testSeq1, testQual1, testSeq2, testQual2, args), test_PairedProcess)


class TestFullPipeline(unittest.TestCase):
    def testProcessForwardDataset(self):
        args.truncate = None
        args.paired = False
        args.forward_fastq = 'test/test1.fastq'
        args.reverse_fastq = None
        args.output_prefix = 'test/test_output/forward'
        args.output_compression = 'none'
        moira.main(args)
        self.assertEqual(open('test/test_output/forward.qc.good.fasta').read(), open('test/test_results/forward.qc.good.fasta').read())
        self.assertEqual(open('test/test_output/forward.qc.good.qual').read(), open('test/test_results/forward.qc.good.qual').read())
        self.assertEqual(open('test/test_output/forward.qc.good.names').read(), open('test/test_results/forward.qc.good.names').read())
        self.assertEqual(open('test/test_output/forward.qc.bad.fasta').read(), open('test/test_results/forward.qc.bad.fasta').read())
        self.assertEqual(open('test/test_output/forward.qc.bad.qual').read(), open('test/test_results/forward.qc.bad.qual').read())
        self.assertEqual(open('test/test_output/forward.qc.bad.names').read(), open('test/test_results/forward.qc.bad.names').read())
    def testProcessPairedDataset(self):
        args.truncate = None
        args.paired = True
        args.forward_fastq = 'test/test1.fastq'
        args.reverse_fastq = 'test/test2.fastq'
        args.output_prefix = 'test/test_output/paired'
        args.output_compression = 'none'
        moira.main(args)
        self.assertEqual(open('test/test_output/paired.qc.good.fasta').read(), open('test/test_results/paired.qc.good.fasta').read())
        self.assertEqual(open('test/test_output/paired.qc.good.qual').read(), open('test/test_results/paired.qc.good.qual').read())
        self.assertEqual(open('test/test_output/paired.qc.good.names').read(), open('test/test_results/paired.qc.good.names').read())
        self.assertEqual(open('test/test_output/paired.qc.bad.fasta').read(), open('test/test_results/paired.qc.bad.fasta').read())
        self.assertEqual(open('test/test_output/paired.qc.bad.qual').read(), open('test/test_results/paired.qc.bad.qual').read())
        self.assertEqual(open('test/test_output/paired.qc.bad.names').read(), open('test/test_results/paired.qc.bad.names').read())
    def testCompression(self):
        args.truncate = None
        args.paired = True
        args.forward_fastq = 'test/test1.fastq.gz'
        args.reverse_fastq = 'test/test2.fastq.bz2'
        args.output_prefix = 'test/test_output/paired'
        args.output_compression = 'gz'
        moira.main(args)
        self.assertEqual(gzip.GzipFile('test/test_output/paired.qc.good.fasta.gz').read(), open('test/test_results/paired.qc.good.fasta').read())        
        args.output_compression = 'bz2'
        moira.main(args)
        self.assertEqual(bz2.BZ2File('test/test_output/paired.qc.good.fasta.bz2').read(), open('test/test_results/paired.qc.good.fasta').read())




testSeq1 = 'CCTACGGGTGGCAGCAGTAGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGCGCGATGAAGGCCTTCGGGTCGTAAAGCGCTTTTTGAGGAGATGAGGAAGGACAGTATCCTCAGAATAAGGATCGGCTAACTACGTGCCAGCAGCCGCGGTAACACGTAGGATCCGAGCGTTATCCGAATTTACTGGGCGTAAAGCGCGTGCAGGCGGTTTGGTAAGTTGGATGTGAAAGCTCCTGGCTCAACTGGGAGAGGCCGTTCAAAACTACCAGACTCGAGGGTGGTTGG'
testSeq2 = 'GANTACAGGGGTATCTAATCCCGTTCGCTCCCCTAGCTTTCGCGTCTTAGCGTCAGGAATGGTCCAGGAGGCCGCCTTCGCCTCTGGTGTTCCTCCCGATATCTACGCATTTCACTACTACACCGGGAATTCCACCTCCCTCTACCACCCTCGAGTCTGGTAGTTTTGAACGGCCTCTCCCAGTTGAGCCAGGAGCTTTCACATCCAACTTACCAAACCGCCTGCACGCGCTTTACGCCCAGTAAATTCGGCTAACGATCTGATCCTACGTGTGCCGTCGGCTGCGGGCACATAGTTAGCG'
testQual1 = '68BCCGGGGBADCBBCFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFEGGGGGGGGGGGGGGGFGGFGGGGGGGFGGGGGGFGFEGFGGGGGGGGGGGGGGGGGGGGGGGCGGGGGGGGGGGGGGGGGD7>CEFGGFCG@@FEEF9D8EGGGEF<CGGCEEG:EEF?9FEGGGGGGGCFFCEEEEFF79E>EG2CFD38EGGGG>)<AC>F73)<CC>:D<EFC<EF8797:?F,<CGF6?5*;.77);8672:CGG6CF+(*0<(5:>>6:B(((('
testQual2 = 'CC#8AFGGGGGGGGGGGGGGEGGGGGGGGGGGGGGFGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGFEFFGGEFCDFDGGFEGD>FFFFGGGFGGGEFGCFFFFGCEAFGGFFFFGFGFFGFEEDGG6=@><FGG708;2EEC;FDEABDGEEF4*,8<>+**/*,1*;?AF+09BE<F(8=6@CEGG=3A:>0*=;;AD;;690(38AB@596148>+24:<(((436+****()65444:><273-))(-(4(7531((((,,()..6>)6)('
testQual1 = [ord(qual) - 33 for qual in testQual1]
testQual2 = [ord(qual) - 33 for qual in testQual2]
testRC2 = ('CGCTAACTATGTGCCCGCAGCCGACGGCACACGTAGGATCAGATCGTTAGCCGAATTTACTGGGCGTAAAGCGCGTGCAGGCGGTTTGGTAAGTTGGATGTGAAAGCTCCTGGCTCAACTGGGAGAGGCCGTTCAAAACTACCAGACTCGAGGGTGGTAGAGGGAGGTGGAATTCCCGGTGTAGTAGTGAAATGCGTAGATATCGGGAGGAACACCAGAGGCGAAGGCGGCCTCCTGGACCATTCCTGACGCTAAGACGCGAAAGCTAGGGGAGCGAACGGGATTAGATACCCCTGTANTC', [7, 8, 21, 8, 29, 21, 13, 13, 8, 7, 11, 11, 7, 7, 7, 7, 16, 18, 20, 22, 7, 19, 7, 12, 7, 8, 8, 12, 18, 22, 17, 27, 29, 25, 19, 19, 19, 20, 21, 8, 7, 9, 9, 9, 9, 10, 21, 18, 19, 7, 7, 7, 27, 25, 19, 17, 10, 29, 23, 19, 16, 21, 24, 20, 31, 33, 32, 23, 18, 7, 15, 24, 21, 26, 26, 35, 32, 26, 26, 28, 9, 15, 29, 25, 32, 18, 28, 38, 38, 36, 34, 31, 21, 28, 23, 7, 37, 27, 36, 33, 24, 15, 10, 37, 32, 30, 26, 9, 16, 11, 9, 14, 9, 9, 10, 29, 27, 23, 11, 9, 19, 37, 36, 36, 38, 35, 33, 32, 36, 35, 37, 26, 34, 36, 36, 17, 26, 23, 15, 22, 38, 38, 37, 27, 29, 31, 28, 21, 38, 38, 35, 36, 36, 37, 38, 37, 37, 38, 37, 38, 37, 37, 37, 37, 38, 38, 37, 32, 36, 34, 38, 37, 37, 37, 37, 34, 38, 37, 36, 38, 38, 38, 37, 38, 38, 38, 37, 37, 37, 37, 29, 35, 38, 36, 37, 38, 38, 35, 37, 35, 34, 37, 36, 38, 38, 37, 37, 36, 37, 38, 38, 38, 38, 35, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 36, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 32, 23, 2, 34, 34])
test_aligned = ('CCTACGGGTGGCAGCAGTAGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGCGCGATGAAGGCCTTCGGGTCGTAAAGCGCTTTTTGAGGAGATGAGGAAGGACAGTATCCTCAGAATAAGGATCGGCTAACTACGTGCCAGCAGCCG-CGGTAACACGTAGGATCCGAGCGTTATCCGAATTTACTGGGCGTAAAGCGCGTGCAGGCGGTTTGGTAAGTTGGATGTGAAAGCTCCTGGCTCAACTGGGAGAGGCCGTTCAAAACTACCAGACTCGAGGGTGGTTG-G-------------------------------------------------------------------------------------------------------------------------------------------', '--------------------------------------------------------------------------------------------------------------------------------------------CGCTAACTATGTGCCCGCAGCCGACGG-CACACGTAGGATCAGATCGTTAGCCGAATTTACTGGGCGTAAAGCGCGTGCAGGCGGTTTGGTAAGTTGGATGTGAAAGCTCCTGGCTCAACTGGGAGAGGCCGTTCAAAACTACCAGACTCGAGGGTGGTAGAGGGAGGTGGAATTCCCGGTGTAGTAGTGAAATGCGTAGATATCGGGAGGAACACCAGAGGCGAAGGCGGCCTCCTGGACCATTCCTGACGCTAAGACGCGAAAGCTAGGGGAGCGAACGGGATTAGATACCCCTGTANTC', 13431)
test_contig = ('CCTACGGGTGGCAGCAGTAGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGCGCGATGAAGGCCTTCGGGTCGTAAAGCGCTTTTTGAGGAGATGAGGAAGGACAGTATCCTCAGAATAAGGATCGGCTAACTACGTGCCAGCAGCCGCGGTAACACGTAGGATCCGAGCGTTATCCGAATTTACTGGGCGTAAAGCGCGTGCAGGCGGTTTGGTAAGTTGGATGTGAAAGCTCCTGGCTCAACTGGGAGAGGCCGTTCAAAACTACCAGACTCGAGGGTGGTAGAGGGAGGTGGAATTCCCGGTGTAGTAGTGAAATGCGTAGATATCGGGAGGAACACCAGAGGCGAAGGCGGCCTCCTGGACCATTCCTGACGCTAAGACGCGAAAGCTAGGGGAGCGAACGGGATTAGATACCCCTGTANTC', [21, 23, 33, 34, 34, 38, 38, 38, 38, 33, 32, 35, 34, 33, 33, 34, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 36, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 37, 38, 37, 36, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 34, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 35, 22, 29, 34, 36, 37, 38, 38, 37, 34, 38, 31, 31, 37, 36, 36, 37, 24, 35, 27, 36, 38, 38, 38, 36, 37, 27, 34, 38, 38, 34, 36, 36, 38, 25, 36, 36, 37, 30, 24, 37, 36, 38, 38, 38, 38, 38, 38, 38, 34, 37, 37, 34, 36, 36, 36, 36, 37, 37, 24, 24, 36, 29, 36, 38, 26, 34, 37, 35, 18, 29, 36, 38, 38, 38, 38, 38, 36, 34, 32, 34, 29, 37, 22, 37, 27, 36, 34, 34, 29, 25, 37, 32, 36, 37, 34, 27, 36, 37, 23, 22, 24, 22, 29, 30, 37, 11, 27, 34, 38, 37, 36, 38, 35, 33, 32, 36, 35, 37, 26, 34, 36, 36, 22, 26, 25, 34, 38, 38, 38, 37, 37, 29, 31, 28, 21, 38, 38, 35, 36, 36, 37, 38, 37, 37, 38, 37, 38, 37, 37, 37, 37, 38, 38, 37, 32, 36, 34, 38, 37, 37, 37, 37, 34, 38, 37, 36, 38, 38, 38, 37, 38, 38, 38, 37, 37, 37, 37, 29, 35, 38, 36, 37, 38, 38, 35, 37, 35, 34, 37, 36, 38, 38, 37, 37, 36, 37, 38, 38, 38, 38, 35, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 36, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 32, 23, 2, 34, 34], 162, 3, 8)
test_ForwardProcess = ('foo', 'CCTACGGGTGGCAGCAGTAGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGCGCGATGAAGGCCTTCGGGTCGTAAAGCGCTTTTTGAGGAGATGAGGAAGGACAGTATCCTCAGAATAAGGATCGGCTAACTACGTGCCAGCAGCCGCGGTAACACGTAGGATCCGAGCGTTATCCGAATTTAC', [21, 23, 33, 34, 34, 38, 38, 38, 38, 33, 32, 35, 34, 33, 33, 34, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 36, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 37, 38, 37, 36, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 34, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 35, 22, 29, 34, 36, 37, 38, 38, 37, 34, 38, 31, 31, 37, 36, 36, 37, 24, 35, 23, 36, 38, 38, 38, 36, 37, 27, 34, 38, 38, 34, 36, 36, 38, 25, 36, 36, 37, 30, 24, 37, 36, 38, 38, 38, 38, 38, 38], 0.9685179556745876, 0, 0, 0)
test_PairedProcess = ('foo', 'CCTACGGGTGGCAGCAGTAGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGCGCGATGAAGGCCTTCGGGTCGTAAAGCGCTTTTTGAGGAGATGAGGAAGGACAGTATCCTCAGAATAAGGATCGGCTAACTACGTGCCAGCAGCCGCGGTAACACGTAGGATCCGAGCGTTATCCGAATTTAC', [21, 23, 33, 34, 34, 38, 38, 38, 38, 33, 32, 35, 34, 33, 33, 34, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 36, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 37, 38, 37, 36, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 34, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 35, 22, 29, 34, 36, 37, 38, 38, 37, 34, 38, 31, 31, 37, 36, 36, 37, 24, 35, 27, 36, 38, 38, 38, 36, 37, 27, 34, 38, 38, 34, 36, 36, 38, 25, 36, 36, 37, 30, 24, 37, 36, 38, 38, 38, 38, 38, 38], 0.9643903629780557, 162, 3, 8)

args = Arguments(alpha = 0.005, match = 1, gap = -2, mismatch = -1, insert = 20, deltaq = 6, consensus_qscore = 'best',
                 paired = True, truncate = 200, only_contig = False, error_calc = 'poisson_binomial', ambigs = 'treat_as_errors',
                 round = False, silent = True, nowarnings = False, doc = False, uncert = 0.01, maxerrors = None, processors = 1,
                 forward_fasta = None, forward_quals = None, reverse_fasta = None, reverse_quals = None,
                 forward_fastq = None, reverse_fastq = None, output_format = 'fasta', collapse = True, pipeline = 'mothur', fastq_offset = 33,
                 relabel = None, output_compression = 'none', qscore_cap = 40, min_overlap = None, trim_overlap = False)


if __name__ == '__main__':
    unittest.main()
