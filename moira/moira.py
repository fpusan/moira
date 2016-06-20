#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
--------------------------------------------------------------------------------------------------------------------------------------------

Quality-filter raw sequence reads using the Poisson binomial filtering algorithm.

The moira.py script and the Poisson binomial filtering algorithm are now described in the following paper:

    Puente-Sánchez, F., Aguirre, J., & Parro, V. (2016).
    A novel conceptual approach to read-filtering in high-throughput amplicon sequencing studies.
    Nucleic acids research, 44 (4): e40.

which can be accessed online at: http://nar.oxfordjournals.org/content/44/4/e40


REQUIREMENTS:

- Expects that input sequences (single or paired) and qualities are in the same order.
- Expects that sequences and qualities are stored only in one line (i.e. >header\\nsequence\\n>header2\\nsequence2).
- Requires numpy.
- OPTIONAL: Requires the bernoulli library, which includes the C implementation of the Poisson binomial filtering algorithm.
  If not available, the script will automaticaly switch to the pure
  python implementation.
- OPTIONAL: Requires the nw_align library, which includes the cython implementation of the Needleman-Wunsch alignment algorithm.
  If not present, the script will automatically switch to the pure python implementation.


USAGE:

    - Make contigs from paired reads (fastq+qual) without quality-filtering:

        moira.py --forward_fasta=<FILE> --forward_qual=<FILE> --reverse_fasta=<FILE> --reverse_qual=<FILE> --paired --only_contig

    - Make contigs from paired reads (fastq) and perform quality-filtering, output results in fastq format:

        moira.py --forward_fastq=<FILE> --reverse_fastq=<FILE> --paired --consensus_qscore posterior --output_format fastq

    - Quality-filter already assembled contigs or single reads:

        moira.py --forward_fasta=<FILE> --forward_qual=<FILE>

    - Input files can be compressed with gzip or bzip2.


OUTPUT:

    - If quality control is being performed, files will be generated with both the sequences that passed the QC and the ones 
      that didn't. A brief report will be included on the headers of the contigs that didn't pass the QC.

        <PREFIX>.qc.good.fasta
        <PREFIX>.qc.good.qual
        <PREFIX>.qc.bad.fasta
        <PREFIX>.qc.bad.qual

    - Else, only two files will be generated.

        <PREFIX>.contigs.fasta
        <PREFIX>.contigs.qual

    - The default prefix is the forward input name without the extension. A custom prefix can be specified with the --output_prefix option.
    - If --output_format is set to "fastq", fastq files will be generated instead of fasta + qual files.
    - If identical sequences are being collapsed, mothur-formatted name files (or UPARSE formatted sequence headers) will also be generated.
    - moira.py will replace ':' for '_' in sequence names for compatibility with the mothur pipeline.


PARAMETERS:

        - Needleman-Wunsch aligner parameters:

                --match (default 1): match score
                --gap (default -2): gap penalty
                --mismatch (default -1): mismatch penalty

        - Contig constructor parameters:
                
                --insert (default 20): quality above which a base will be used for filling a complementary gap or ambiguity.
                --deltaq (default 6): minimum quality difference allowed between two mismatched bases for not including an N 
                                      in the consensus sequence.
                --consensus_qscore (default 'best')
                        best: use the best quality on each position of the alignment as the consensus quality score (Unless an 
                              ambiguity is introduced in that position by the contig constructor. In that case, the reported quality
                              score will be always 2).
                        sum: in matching bases, consensus quality score will be the sum of the qualities of both reads in that
                             position of the alignment.
                        posterior: use Edgar & Flyvbjerg's (2015) method for calculating consensus quality scores. The insert and deltaq
                             parameters will be ignored. Ambiguities will be introduced in gaps, or if two mismatched bases have exactly
                             the same quality score. In that case, the reported quality score will be always 2.
                --qscore_cap (default 40): Maximum consensus quality score to report. Higher consensus quality scores will be trimmed to
                             the value of --qscore_cap. Setting it to 0 will remove the cap.
                --trim_overlap (default False): trim the contig to the overlapping region.

        - Quality-filtering parameters:

                --collapse (default True): if True, identical sequences will be collapsed before quality control, and the one with
                        the best quality will be used as a representative of the whole group.

                --error_calc (default 'poisson_binomial'): algorithm used for error calculation.
                        poisson_binomial: calculate the Poisson binomial distribution (sum of bernoulli random variables).
                        poisson_binomial_py: calculate the Poisson binomial distribution (sum of bernoulli random variables), using the Python implementation.
                        poisson: approximating sum of bernoulli random variables to a poisson distribution.
                        bootstrap: numerical generation of an error distribution (deprecated).

                --ambigs (default treat_as_error): handling of ambiguous positions during quality checking.
                        treat_as_error: will consider than ambiguities always result in a misread base.
                        disallow: will discard sequences with ambiguities.
                        ignore: will ignore ambiguities.

                --truncate (default None): truncate sequences to a fixed length before quality control. Discard sequences smaller than that length.

                --min_overlap (default None): discard contigs with less than the specified overlap length. Will be ignored if using single-end reads.

                --round: Round down the predicted errors to the nearest integer prior to filtering.

                --uncert (default 0.01): maximum divergence of the observed sequence from the original one due to sequencing errors.

                --maxerrors (no default value): maximum errors allowed in the sequence.
                        Will override --uncert if specified as a parameter.

                --alpha (default 0.005): probability of underestimating the actual errors of a sequence.

                --bootstrap (default 100): number of replicates per position used for error calculation by the bootstrap method.
        
        - Other:

                --paired: input files are paired end files and will be assembled into contigs.
                --only_contigs: assemble contigs, don\'t do quality control.
                --relabel (default False): if a prefix string is introduced, sequential labels will be generated for the sequences,
                                           with the format <prefix>N, where N=1,2,3 etc.
                --output_format (default fasta):
                        fasta: output files in fasta + qual format.
                        fastq: output files in fastq format.
                --output_prefix: Prefix for the output file names.
                --output_compression (default none): Compression applied to output files. Can be 'none', 'bz2', or 'gzip'.
                --pipeline (default mothur):
                        mothur: output for collapsed sequences will be in mothur\'s fasta + names format.
                        USEARCH: output for collapsed sequences will be in a single fasta file, with abundance information stored in the sequence header.
                --fastq_offset (default 33): ASCII/qscore encoding.
                --processors (default 1): number of processes to use.
                --silent: Do not print welcome, progress and goodbye messages. Warnings will still be printed.

COMMENTS:

        - Alignment parameters are set to replicate mothur's default implementation of the Needleman-Wunsch algorithm.

        - The 'insert' and 'deltaq' parameters from mothur make.contig are also reproduced. They are set at their default values.
          More details can be found at www.mothur.org/wiki/Make.contigs. They will be applied only if the --consensus_qscore parameter is
          set to 'best' or 'sum'.

        - We now provide the option to calculate posterior error probabilities for reporting consensus quality scores, as described by Edgar &
          Flyvbjerg (2015). This can be achieved by adding the '--consensus_qscore posterior' or '-q posterior' flags, as shown in the example above.
          While this method is theoretically sound, it only improved the results obtained by using our default '-q best' option for some of our mock
          community samples.

        - Approximating the sum of Bernoulli random variables to a poisson distribution is quicker than calculating 
          their exact sum (Poisson binomial distribution). It proves specially useful for long reads (>500 nt).
          That said, the Poisson binomial filtering algorithm is also implemented in C and even the python implementation is quick enough
          for processing large datasets. The bootstrap method (--error_calc bootstrap) is a numerical algorithm for performing the sum of bernoulli
          random variables. It is only included for testing purposes.

        - Quality-filtering will discard the contigs expected to have more than 'alpha' chances of diverging from the original 
          sequence more than the value specified by the 'uncert' param. That means that, during distance calculation between two
          given sequences, the observed distance will be at most 'dist + 2*uncert', where 'dist' is the original distance between
          those sequences without sequencing errors. Thus, a good rule of thumb would be considering the effective OTU clustering 
          distance to be actually 'OTUdist - 2*uncert', where OTUdist is the distance used for clustering the observed sequences.


Distributed under the Modified BSD license.

--------------------------------------------------------------------------------------------------------------------------------------------


"""


__author__ = 'Fernando Puente-Sánchez'
__email__ = 'fpusan@gmail.com'
__version__ = '1.3.0'
__date__ = '20-Jun-2016'
__license__ = 'BSD-3'
__copyright__ = 'Copyright 2013-2016 Fernando Puente-Sánchez'

BSD3_LICENSE = """

    Copyright (c) 2015, Fernando Puente Sánchez
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither the name of moira, Poisson binomial filtering and its Poisson
      approximation, nor the names of its contributors may be used to endorse or
      promote products derived from this software without specific prior
      written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

CITATION = """
Puente-Sánchez F, Aguirre J, Parro V (2015),
A novel conceptual approach to read-filtering in high-throughput amplicon sequencing studies.
Nucleic acids research, gkv1113.
"""


#Python stdlib imports.
import sys
import os
import time
import math
import argparse
from multiprocessing import Pool
import traceback
import resource
import io
import gzip
import bz2
from itertools import izip

#Non stdlib imports.
numpy = True
try:
    from numpy.random import random
    from numpy import percentile
    from numpy import isnan
except:
    numpy = False

Cbernoulli = True
try:
    import bernoulli
except:
    Cbernoulli = False

Cy_nw_align = True 
try:
    import nw_align as nw
except:
    Cy_nw_align = False
    

################################################################################################################
#Main program flow.
################################################################################################################

def main(args):
    """Get things done."""

    #Be nice and say hello.
    if not args.silent:
        print
        print '-------------------------------------------------------------------------------'
        print
        print 'moira.py v%s'%__version__
        print 'Last updated: %s'%__date__
        print
        print 'By Fernando Puente-Sánchez'
        print '- Department of Molecular Evolution'
        print '  Centro de Astrobiología (CSIC-INTA), Instituto Nacional de Técnica Aeroespacial'
        print '- Systems Biology Program'
        print '  Centro Nacional de Biotecnología (CSIC)'
        print 'fpusan@gmail.com'
        print 'http://github.com/fpusan/moira'
        print 
        print 'When using, please cite:'
        print CITATION
        print 'Distributed under the Revised BSD License.'
        print
        print '-------------------------------------------------------------------------------'
        print

    #Parse arguments for command line and check everything's ok.
    if not check_arguments(args):
        return 1

    ###############
    #Open and create the necessary files.
    open_write, suffix = {'none' : (open, ''), 'gz': (gzip.GzipFile, '.gz'), 'bz2': (bz2.BZ2File, '.bz2')}[args.output_compression]
    if args.output_prefix:
        output_name = args.output_prefix
    else:
        if args.forward_fastq:
            output_name = '.'.join(args.forward_fastq.split('.')[:-1])
        else:
            output_name = '.'.join(args.forward_fasta.split('.')[:-1])

    try:
        if args.forward_fastq:
            forward_fastq_data = open_input(args.forward_fastq)
        else:
            forward_fasta_data = open_input(args.forward_fasta)
            forward_qual_data = open_input(args.forward_qual)


        if args.paired:
            if args.reverse_fastq:
                reverse_fastq_data = open_input(args.reverse_fastq)
            else:
                reverse_fasta_data = open_input(args.reverse_fasta)
                reverse_qual_data = open_input(args.reverse_qual)
        else:
            reverse_fastq_data, reverse_fasta_data, reverse_qual_data = None, None, None
    

        if args.only_contig:
            if args.output_format == 'fastq':
                contig_output = open_write('%s.contigs.fastq%s'%(output_name, suffix), 'w')
                qual_output = None
            else:
                contig_output = open_write('%s.contigs.fasta%s'%(output_name, suffix), 'w')
                qual_output = open_write('%s.contigs.qual%s'%(output_name, suffix), 'w')
            if args.collapse and args.pipeline == 'mothur':
                names_output = open_write('%s.contigs.names%s'%(output_name, suffix), 'w')
            else:
                names_output = None
            if not args.truncate and not args.min_overlap: 
                bad_contig_output, bad_qual_output, bad_names_output = None, None, None
            else: #Not doing qc par se, but filtering due to length and/or minoverlap.
                if args.output_format == 'fastq':
                    bad_contig_output = open_write('%s.bad.contigs.fastq%s'%(output_name, suffix), 'w')
                    bad_qual_output = None
                else:
                    bad_contig_output = open_write('%s.bad.contigs.fasta%s'%(output_name, suffix), 'w')
                    bad_qual_output = open_write('%s.bad.contigs.qual%s'%(output_name, suffix), 'w')
                if args.collapse and args.pipeline == 'mothur':
                    bad_names_output = open_write('%s.bad.contigs.names%s'%(output_name, suffix), 'w')
                else:
                    bad_names_output = None
                
        else:
            if args.output_format == 'fastq':
                contig_output = open_write('%s.qc.good.fastq%s'%(output_name, suffix), 'w')
                qual_output = None
                bad_contig_output = open_write('%s.qc.bad.fastq%s'%(output_name, suffix), 'w')
                bad_qual_output = None
            else:
                contig_output = open_write('%s.qc.good.fasta%s'%(output_name, suffix), 'w')
                qual_output = open_write('%s.qc.good.qual%s'%(output_name, suffix), 'w')
                bad_contig_output = open_write('%s.qc.bad.fasta%s'%(output_name, suffix), 'w')
                bad_qual_output = open_write('%s.qc.bad.qual%s'%(output_name, suffix), 'w')
            if args.collapse and args.pipeline == 'mothur':
                names_output = open_write('%s.qc.good.names%s'%(output_name, suffix), 'w')
                bad_names_output = open_write('%s.qc.bad.names%s'%(output_name, suffix), 'w')
            else:
                names_output, bad_names_output = None, None
                

        if args.paired:
            contig_report_output = open_write('%s.contigs.report%s'%(output_name, suffix), 'w')
            contig_report_output.write('header\tn_seqs\toverlap_length\tgaps\tmismatches\n')
        else:
            contig_report_output = None

        
    except IOError, e:
        print e
        print
        return 1
    ###############

    #Now that files are ready, start with the actual work..
    try:
        #Get the number of seqs:
        if args.forward_fastq:
            data = forward_fastq_data
        else:
            data = forward_fasta_data
        input_seqs = 0.0
        for line in data:
            if args.forward_fastq:
                input_seqs += 0.25
            elif line.startswith('>'):
                input_seqs += 1
            else:
                pass
        data.seek(0)
        ###############
    
        #Prepare for main loop:
        if args.processors > 1:
            pool = Pool(args.processors)
        if args.forward_fastq:
            parse_seqs = parse_fastq(forward_fastq_data, reverse_fastq_data, args.fastq_offset)
        else:
            parse_seqs = parse_fasta_and_qual(forward_fasta_data, forward_qual_data,
                                              reverse_fasta_data, reverse_qual_data)
        processed_seqs = 0
        discarded_errors = 0.0
        discarded_minlength = 0.0
        discarded_minoverlap = 0.0
        if args.collapse:
            uniques = {} #{sequence{rep_header: header, rep_errors: expected_errors, rep_quals: quals, names_info: [headers]}
        else:
            pass
        time_start = time.time()

        #Go for it:
        while 1:
    
            #Timer.
            if processed_seqs % 2*args.processors  == 0 and processed_seqs != 0:
                memused = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
                elapsed = time.time() - time_start            
                percent = (processed_seqs / input_seqs) * 100
                estimated_left = (elapsed / (processed_seqs / input_seqs)- elapsed) / 60
                msg='Using %.2f MB. %d sequences processed (%.2f%%) in %.1f seconds, %.1f minutes to finish.\r'%(memused, processed_seqs, 
                                                                                                     percent, elapsed, estimated_left)
                if not args.silent:
                    print msg,
            
            #Read data and launch asynchronous processes:
            results = []
            if args.processors > 1:
                for i in range(args.processors):
                    try:
                        header, forward_sequence, forward_quals, reverse_sequence, reverse_quals = parse_seqs.next()
                        results.append(pool.apply_async(process_data,
                                        (header, forward_sequence, forward_quals, reverse_sequence, reverse_quals, args)))
                    except StopIteration:
                        break
            else: #Avoid using the pool (Cprofile does not like the multiprocessing module).
                try:
                    header, forward_sequence, forward_quals, reverse_sequence, reverse_quals = parse_seqs.next()
                    results.append(process_data(header, forward_sequence, forward_quals, reverse_sequence, reverse_quals, args))
                except StopIteration:
                    pass
            
            #Exit loop if we've finished:
            if not results:
                if not args.silent:
                    print '%d sequences processed in %.1f seconds.'%(processed_seqs, (time.time() - time_start)) + ' '*80
                break
        
            #Retrieve results from asynchronous processes:
            if args.processors > 1:
                results = [result.get() for result in results]
            for header, contig, contig_quals, expected_errors, overlap_length, gaps, mismatches in results:
                if isnan(expected_errors): #Should not happen anymore, but we're still testing just in case.
                    raise ReturnedNaNError(header)

                if args.collapse:
                    #Collapse unique sequences and choose the representative with the greater quality.
                    if contig not in uniques:
                        uniques[contig] = {'rep_header': header, 'rep_errors': expected_errors,
                                           'rep_quals': contig_quals, 'names_info': [header],
                                           'overlap_length': overlap_length, 'gaps': gaps, 'mismatches': mismatches}
                    else:
                        if expected_errors < uniques[contig]['rep_errors']:
                            uniques[contig]['rep_header'] = header
                            uniques[contig]['rep_errors'] = expected_errors
                            uniques[contig]['rep_quals'] = contig_quals
                            uniques[contig]['names_info'].insert(0, header)
                            uniques[contig]['overlap_length'] = overlap_length
                            uniques[contig]['gaps'] = gaps
                            uniques[contig]['mismatches'] = mismatches
                        else:
                            uniques[contig]['names_info'].append(header)                         
                else:
                    #Directly check if the sequence has passed the filter and write it.
                    filtering_results = write_results(processed_seqs, header, contig, contig_quals, expected_errors, None,
                                                      overlap_length, gaps, mismatches,
                                                      args,
                                                      contig_output, qual_output, names_output,
                                                      bad_contig_output, bad_qual_output, bad_names_output, contig_report_output)
                    discarded_errors += filtering_results[0]
                    discarded_minlength += filtering_results[1]
                    discarded_minoverlap += filtering_results[3]
                                   
                processed_seqs += 1                
        ###############

        #If we are collapsing unique sequences: after we've finished reading the file, sort by abundance, filter and write the collapsed sequences.
        if args.collapse:
            sorted_uniques = sorted(uniques, key = lambda sequence: len(uniques[sequence]['names_info']), reverse = True)
            for index, sequence in enumerate(sorted_uniques, start=1):
                values = uniques[sequence]
                filtering_results = write_results(index, values['rep_header'], sequence,
                                                  values['rep_quals'], values['rep_errors'],
                                                  values['names_info'],
                                                  values['overlap_length'], values['gaps'], values['mismatches'],
                                                  args,
                                                  contig_output, qual_output, names_output,
                                                  bad_contig_output, bad_qual_output, bad_names_output, contig_report_output)
                discarded_errors += filtering_results[0]
                discarded_minlength += filtering_results[1]
                discarded_minoverlap += filtering_results[2]
        ###############
                                   

        #Be nice and say goodbye.
        if not args.silent:
            remaining_seqs = processed_seqs - discarded_errors - discarded_minlength - discarded_minoverlap
            print '- Kept %d (%.2f%%) of the original sequences.'%(remaining_seqs, ((remaining_seqs / processed_seqs) * 100))
            if args.truncate:
                print '- %d (%.2f%%) of the original sequences were discarded due to length < %s.'%(discarded_minlength, (discarded_minlength / processed_seqs) * 100,
                                                                                                    args.truncate)
            if args.paired and args.min_overlap:
                print '- %d (%.2f%%) of the original sequences were discarded due to paired-end reads having an overlap length < %s.'%(discarded_minoverlap,
                                                                                                        (discarded_minoverlap / processed_seqs) * 100, args.min_overlap)

            print '- %d (%.2f%%) of the original sequences were discarded due to low quality.\n'%(discarded_errors, (discarded_errors / processed_seqs) * 100)
            print 'The following output files were generated:'
            if args.only_contig:
                if args.output_format == 'fastq':
                    print '%s.contigs.fastq%s'%(output_name, suffix)
                    if args.truncate or args.min_overlap:
                        print '%s.bad.contigs.fastq%s'%(output_name, suffix)
                else:
                    print '%s.contigs.fasta%s'%(output_name, suffix)
                    print '%s.contigs.qual%s'%(output_name, suffix)     
                    if args.truncate or args.min_overlap:
                        print '%s.bad.contigs.fasta%s'%(output_name, suffix)
                        print '%s.bad.contigs.qual%s'%(output_name, suffix)
                if args.collapse and args.pipeline == 'mothur':
                    print '%s.contigs.names%s'%(output_name, suffix)
                    if args.truncate or args.min_overlap:
                        print '%s.bad.contigs.names%s'%(output_name, suffix)
            else:
                if args.output_format == 'fastq':
                    print '%s.qc.good.fastq%s'%(output_name, suffix)
                else:
                    print '%s.qc.good.fasta%s'%(output_name, suffix)
                    print '%s.qc.good.qual%s'%(output_name, suffix)
                if args.collapse and args.pipeline == 'mothur':
                    print '%s.qc.good.names%s'%(output_name, suffix)
                if args.output_format == 'fastq':
                    print '%s.qc.bad.fastq%s'%(output_name, suffix)
                else:   
                    print '%s.qc.bad.fasta%s'%(output_name, suffix)
                    print '%s.qc.bad.qual%s'%(output_name, suffix)
                if args.collapse and args.pipeline == 'mothur':
                    print '%s.qc.bad.names%s'%(output_name, suffix)
    
            if args.paired:
                print '%s.contigs.report%s\n'%(output_name, suffix)
            else:
                print

        ###############

    #Tidy your room after you play (even if you kinda broke your toys).
    finally:
        try:
            if args.only_contig:
                contig_output.close() 
                qual_output.close()
                if args.collapse:
                    names_output.close()
            else:
                contig_output.close()
                qual_output.close()
                bad_contig_output.close()
                bad_qual_output.close()
                if args.collapse and args.pipeline == 'mothur':
                    names_output.close()
                    bad_names_output.close()
        except:
            pass

    ###############

        
def parse_arguments():
    """Parse the command line arguments and return an object containing them."""
    
    def str2bool(value):
        return value.lower() in ("yes", "true", "t", "1")

    parser = argparse.ArgumentParser(description = 'Perform quality filtering on a set of sequences.')
    
    general = parser.add_argument_group('General options')
    general.add_argument('-ff', '--forward_fasta', type = str,
                        help = 'Forward fasta file (can be gzip or bzip2 compressed).')
    general.add_argument('-fq', '--forward_qual', type = str,
                        help = 'Forward qual file (can be gzip or bzip2 compressed).')
    general.add_argument('-rf', '--reverse_fasta', type = str,
                        help = 'Reverse fasta file (can be gzip or bzip2 compressed).')
    general.add_argument('-rq', '--reverse_qual', type = str,
                        help = 'Reverse qual file (can be gzip or bzip2 compressed).')
    general.add_argument('-ffq', '--forward_fastq', type= str,
                        help = 'Forward fastq file (can be gzip or bzip2 compressed).')
    general.add_argument('-rfq', '--reverse_fastq', type = str,
                        help = 'Forward fastq file (can be gzip or bzip2 compressed).')
    general.add_argument('-l', '--relabel', type = str,
                        help = 'Generate sequential labels for the ordered sequences, with the specified string at the beginning.') 
    general.add_argument('-o', '--output_format', type = str, default = 'fasta',
                        choices = ('fasta' , 'fastq'),
                        help = 'Output format: fasta with qual (and mothur name file if --collapse) or fastq.')
    general.add_argument('-pi', '--pipeline', type = str, default = 'mothur',
                        choices = ('mothur', 'USEARCH'),
                        help = 'Make the output format compatible with the indicated analysis pipeline.')
    general.add_argument('-op', '--output_prefix', type = str,
                        help = "Prefix for the output files")
    general.add_argument('-oc', '--output_compression', type = str, default = 'none',
                        choices = ('none', 'gz', 'bz2'),
                        help = "Prefix for the output files")
    general.add_argument('-p', '--processors', type = int, default = 1,
                        help = 'Number of processors to be used.')
    general.add_argument('--paired', action = 'store_true',
                        help = 'Assemble paired-end reads and perform quality control on the resulting contig.')
    general.add_argument('-fo', '--fastq_offset', type = int, default = 33)
    general.add_argument('--only_contig', action = 'store_true',
                        help = 'Assemble contigs but don\'t perform quality control.')
    general.add_argument('--silent', action = 'store_true',
                        help = 'Do not print welcome, progress and goodbye messages. Warnings will still be printed.')
    general.add_argument('--nowarnings', action = 'store_true',
                        help = 'Do not print warning messages.')
    general.add_argument('--doc', action = 'store_true',
                        help = 'Print full documentation.')
 
    constructor = parser.add_argument_group('Contig construction options')
    constructor.add_argument('-m', '--match', type = int, default = 1,
                        help = 'Needleman-Wunsch aligner match score.')
    constructor.add_argument('-x', '--mismatch', type = int, default = -1,
                        help = 'Needleman-Wunsch aligner mismatch penalty.')
    constructor.add_argument('-g', '--gap', type = int, default = -2,
                        help = 'Needleman-Wunsch aligner gap penalty.')
    constructor.add_argument('--trim_overlap', action = 'store_true',
                        help = 'Trim the contig to the overlapping region.')
    constructor.add_argument('-i', '--insert', type = int, default = 20,
                        help = 'Contig constructor insert threshold.')
    constructor.add_argument('-d', '--deltaq', type = int, default = 6,
                        help = 'Contig constructor mismatch correction deltaq threshold.')
    constructor.add_argument('-q', '--consensus_qscore', type = str, default = 'best',
                        choices = ('best', 'sum', 'posterior'),
                        help = 'Contig constructor consensus qscore: always report best qscore, sum qscores for matching bases, report posterior qscores.')
    constructor.add_argument('-z', '--qscore_cap', type = int, default = 40,
                        help = 'Maximum consensus quality score reported by the contig constructor. Use 0 for removing the qscore cap.')
    
    filtering = parser.add_argument_group('Sequence filtering options')
    filtering.add_argument('-c', '--collapse', type = str2bool, default = 'True',
                        help = 'Collapse identical sequences before quality control.')
    filtering.add_argument('-t', '--truncate', type = int,
                        help = 'Truncate sequences to a fixed length before quality control. Discard smaller sequences.')
    filtering.add_argument('-mo', '--min_overlap', type = int,
                        help = 'Discard contigs with less than the specified overlap length.')
    filtering.add_argument('-e', '--error_calc', type = str, default = 'poisson_binomial',
                        choices = ('poisson_binomial', 'poisson_binomial_py', 'poisson', 'bootstrap'),
                        help = 'Error calculation method.')
    filtering.add_argument('-n', '--ambigs', type = str, default = 'treat_as_errors',
                        choices = ('disallow', 'ignore', 'treat_as_errors'),
                        help = 'Treatment of ambiguities: remove sequences with ambigs, ignore ambigs, treat ambigs as errors.')
    filtering.add_argument('-r', '--round', action = 'store_true',
                        help = 'Round down the predicted number of errors to their nearest integer prior to filtering.')
    err_uncert = filtering.add_mutually_exclusive_group()
    err_uncert.add_argument('-u', '--uncert', type = float, default = 0.01,
                        help = 'Maximum allowed uncertainty (errors / sequence length).')
    err_uncert.add_argument('-me', '--maxerrors', type = float,
                        help = 'Maximum allowed errors per sequence.')
    filtering.add_argument('-a', '--alpha', type = float, default = 0.005,
                        help = 'Alpha cutoff value for the error distributions.')
    filtering.add_argument('-b', '--bootstrap', type = int, default = 100,
                        help = 'Number of replicates to use with the bootstrap method')
    
    args = parser.parse_args()
    
    return args


def check_arguments(args):
    """Check that all the arguments are ok and print warning messages if needed."""
    print
    ok = True
    #Check for arguments asking for usage or other information.
    if args.doc:
        print __doc__
        ok = False
    if not ok:
        return False
    #Check for missing libraries that make the program unable to run.
    elif not numpy:
        print '\nFailed to import numpy.'
        print 'For more info type moira.py -h or moira.py --doc .\n'
        return False
    else:
        #Check for missing arguments.
        if args.only_contig:
            args.paired = True
        else:
            if not args.forward_fastq and (not args.forward_fasta or not args.forward_qual):
                if not args.nowarnings:
                    print '- You must at least provide one fastq file, or a fasta and quality files.'
                ok = False
            if args.paired:
                if not args.reverse_fastq and (not args.reverse_fasta or not args.reverse_qual):
                    if not args.nowarnings:
                        print '- You must provide one reverse fastq file, or reverse fasta and quality files.'
                    ok = False
        #Check for arguments with wrong values (type checking was performed by argparse).
        if args.match < 0:
            if not args.nowarnings:
                print '- Needleman-Wunsch match score must be a non-negative integer.'
            ok = False
        if args.mismatch > 0:
            if not args.nowarnings:
                print '- Needleman-Wunsch mismatch penalty must be a non-positive integer.'
            ok = False
        if args.gap > 0:
            if not args.nowarnings:
                print '- Needleman-Wunsch gap penalty must be a non-positive integer.'
            ok = False
        if args.insert < 1:
            if not args.nowarnings:
                print '- The contig constructor insert parameter must be a positive integer.'
            ok = False
        if args.deltaq < 1:
            if not args.nowarnings:
                print '- The contig constructor deltaq parameter must be a positive integer.'
            ok = False
        if not 0 < args.uncert <= 1:
            if not args.nowarnings:
                print '- The uncert parameter must be between 0 (not included) and 1.'
            ok = False
        if args.maxerrors != None and args.maxerrors <= 0:
            if not args.nowarnings:
                print '- The maxerrors parameter must be greater than 0.'
            ok = False
        if not 0 < args.alpha < 1:
            if not args.nowarnings:
                print '- The alpha parameter must be between 0 (not included) and 1.'
            ok = False
        if args.truncate and args.truncate <= 0:
            if not args.nowarnings:
                print '- The truncate parameter must be greater than 0.'
            ok = False
        if args.min_overlap != None and args.min_overlap <= 0:
            if not args.nowarnings:
                print '- The min_overlap parameter must be greater than 0.'
            ok = False
        if ok:
            #Check for arguments with wrong values that can be harmlessly changed back to their defaults.
            if args.processors < 1:
                if not args.nowarnings:
                    print '- Processors must be a non-zero positive integer. The default value of 1 will be used.'
                args.processors = 1
            #Check for missing libraries that can be substituted by their python implementation.
            if args.paired and not Cy_nw_align:
                if not args.nowarnings:
                    print '\nCython implementation of Needlemann-Wunsch aligner is not present.'
                    print 'Will use pure python implementation instead.\n'
            if args.error_calc == 'poisson_binomial' and not Cbernoulli:
                if not args.nowarnings:
                    print '\nC implementation of Poisson binomial filtering algorithm is not present.'
                    print 'Will use pure python implementation instead.\n'
            #Print other useful info:
            if (args.reverse_fasta or args.reverse_fastq) and not args.paired:
                if not args.nowarnings:
                    print 'You provided a reverse sequence file, but not the --paired flag. Note that only the forward file will be processed.'
                    print
            if (args.min_overlap) and not args.paired:
                if not args.nowarnings:
                    print 'You specified a value for --min_overlap, but not the --paired flag. Note that contigs will not be assembled.'
                    print
            if args.error_calc == 'bootstrap':
                if not args.nowarnings:
                    print 'The bootstrap method is only included for testing and nostalgia. Mainly the second, at this point.'
                    print 'If your purpose falls outside of these two categories, please consider switching to "-e poisson_binomial" or "-e poisson".'
                    print
            return ok
        else:
            if not args.nowarnings:
                print '\nFor more info type moira.py -h or moira.py --doc.\n'
            return False


def process_data(header, forward_sequence, forward_quals, reverse_sequence, reverse_quals, args):
    """Assemble(if needed), truncate(if needed) and quality check the sequences."""

    #Try/except with traceback so exceptions in worker processes get properly reported.
    try:
        if args.paired:
            #Align and assemble:
            assert reverse_sequence and reverse_quals
            reverse_sequence, reverse_quals = reverse_complement(reverse_sequence, reverse_quals)
            if Cy_nw_align:
                forward_aligned, reverse_aligned, score = nw.nw_align(forward_sequence, reverse_sequence,
                                                                 args.match, args.mismatch, args.gap)
            else:
                forward_aligned, reverse_aligned, score = nw_align(forward_sequence, reverse_sequence,
                                                                 args.match, args.mismatch, args.gap)
            contig, contig_quals, overlap_length, gaps, mismatches = make_contig(forward_aligned, forward_quals, reverse_aligned, reverse_quals,
                                                                                 args.insert, args.deltaq,
                                                                                 args.consensus_qscore, args.qscore_cap, args.trim_overlap)
            ### 
        else:
            overlap_length, gaps, mismatches = 0, 0, 0
            contig, contig_quals = forward_sequence, forward_quals
        if args.truncate:
            contig, contig_quals = contig[:args.truncate], contig_quals[:args.truncate]
 
        if args.only_contig:
            expected_errors = 0

        else:
            #Having qvalues of 0 (happened in artificial datasets) will lead to divisions by zero in the PB error calculation.
            contig_quals = [qual if qual > 0 else 1 for qual in contig_quals]
            if args.error_calc in ('poisson_binomial', 'poisson_binomial_py'):
                if args.error_calc == 'poisson_binomial' and Cbernoulli:
                    expected_errors, Ns = bernoulli.calculate_errors_PB(contig, contig_quals, args.alpha)
                    if isnan(expected_errors): #Happens randomly when using the C implementation with very short sequences.
                        expected_errors, Ns = calculate_errors_PB(contig, contig_quals, args.alpha) #So we fall back to the python version.
                else:
                    expected_errors, Ns = calculate_errors_PB(contig, contig_quals, args.alpha)
            elif args.error_calc == 'poisson':
                expected_errors, Ns = calculate_errors_poisson(contig, contig_quals, args.alpha)
            else:#qmode == 'bootstrap':
                expected_errors, Ns = calculate_errors_bootstrap(contig, contig_quals, args.alpha, args.bootstrap)
        
            if args.ambigs == 'treat_as_errors':
                expected_errors = expected_errors + Ns

        if args.round:
            expected_errors = math.floor(expected_errors)

        return header, contig, contig_quals, expected_errors, overlap_length, gaps, mismatches

    except KeyboardInterrupt: #Let KeyboardInterrupt be raised only in the main function, b/c it hangs the program when catched inside a multiprocessing pool.
        pass
    except:
        traceback.print_exc()
        raise
       
 
def write_results(index, header, sequence, quals, expected_errors, names_info,
                  overlap_length, gaps, mismatches, args,
                  contig_output, qual_output, names_output,
                  bad_contig_output, bad_qual_output, bad_names_output,
                  contig_report_output):

    """
    Decide if we're keeping or discarding the sequence, and write it to the corresponding output files.
    Will return a tuple containing the number of sequences that got discarded due to errors and due to small length.
    (We can discard more than one sequence in a single call, if --collapse = True)
    """

    if args.relabel:
        header = '%s%d'%(args.relabel, index)


    if args.pipeline == 'USEARCH':
        if names_info:
            size = len(names_info)
        else:
            size = 1
        header = header + ';ee=%.2f;size=%d;'%(expected_errors, size)


    if args.paired:
        contig_report_output.write('%s\t%s\t%s\t%s\t%s\n'%(header, len(names_info), overlap_length, gaps, mismatches))


    if args.truncate and len(sequence) < args.truncate:
        if args.output_format == 'fastq':
            bad_contig_output.write('@%s\tlength below %s\n%s\n+\n%s\n'%(header, args.truncate, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
        else:
            bad_contig_output.write('>%s\tlength below %s\n%s\n'%(header, args.truncate, sequence))
            bad_qual_output.write('>%s\tlength below %s\n%s\n'%(header, args.truncate, ' '.join(map(str, quals))))
        if args.collapse:
            if args.pipeline == 'mothur':
                bad_names_output.write('%s\t%s\n'%(header.lstrip('>'), ','.join(names_info)))
            return 0, len(names_info), 0
        else:
            return 0, 1, 0


    elif args.min_overlap and overlap_length < args.min_overlap:
        if args.output_format == 'fastq':
            bad_contig_output.write('@%s\toverlap length below %s\n%s\n+\n%s\n'%(header, args.truncate, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
        else:
            bad_contig_output.write('>%s\toverlap length below %s\n%s\n'%(header, args.min_overlap, sequence))
            bad_qual_output.write('>%s\toverlap length below %s\n%s\n'%(header, args.min_overlap, ' '.join(map(str, quals))))
        if args.collapse:
            if args.pipeline == 'mothur':
                bad_names_output.write('%s\t%s\n'%(header.lstrip('>'), ','.join(names_info)))
            return 0, 0, len(names_info)
        else:
            return 0, 0, 1

    
    elif args.only_contig:
        if args.output_format == 'fastq':
            contig_output.write('@%s\n%s\n+\n%s\n'%(header, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
        else:
            contig_output.write('>%s\n%s\n'%(header, sequence))
            qual_output.write('>%s\n%s\n'%(header, ' '.join(map(str, quals))))
        if args.collapse and args.pipeline == 'mothur':
            names_output.write('%s\t%s\n'%(header.lstrip('>'), ','.join(names_info)))
        return 0, 0, 0

        
    elif 'N' in sequence and args.ambigs == 'disallow':
        if args.output_format == 'fastq':
            bad_contig_output.write('@%s\tcontains ambiguities\n%s\n+\n%s\n'%(header, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
        else:
            bad_contig_output.write('>%s\tcontains ambiguities\n%s\n'%(header, sequence))
            bad_qual_output.write('>%s\tcontains ambiguities\n%s\n'%(header, ' '.join(map(str, quals))))
        if args.collapse:
            if args.pipeline == 'mothur':
                bad_names_output.write('%s\t%s\n'%(header, ','.join(names_info)))
            return len(names_info), 0, 0
        else:
            return 1, 0, 0


    elif args.maxerrors: #maxerrors mode.
        if expected_errors <= args.maxerrors:
            if args.output_format == 'fastq':
                contig_output.write('@%s\n%s\n+\n%s\n'%(header, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
            else:
                contig_output.write('>%s\n%s\n'%(header, sequence))
                qual_output.write('>%s\n%s\n'%(header, ' '.join(map(str, quals))))
            if args.collapse and args.pipeline == 'mothur':
                names_output.write('%s\t%s\n'%(header, ','.join(names_info)))
            return 0, 0, 0
        else:
            if args.output_format == 'fastq':
                bad_contig_output.write('@%s\terrors > %.2f\n%s\n+\n%s\n'%(header, args.maxerrors, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
            else: 
                bad_contig_output.write('>%s\terrors > %.2f\n%s\n'%(header, args.maxerrors, sequence))
                bad_qual_output.write('>%s\terrors > %.2f\n%s\n'%(header, args.maxerrors, ' '.join(map(str, quals))))
            if args.collapse:
                if args.pipeline == 'mothur':
                    bad_names_output.write('%s\t%s\n'%(header.lstrip('>'), ','.join(names_info)))
                return len(names_info), 0, 0
            else:
                return 1, 0, 0


    else: #uncert mode.
        if expected_errors <= len(sequence) * args.uncert:
            if args.output_format == 'fastq':
                contig_output.write('@%s\n%s\n+\n%s\n'%(header, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
            else:
                contig_output.write('>%s\n%s\n'%(header, sequence))
                qual_output.write('>%s\n%s\n'%(header, ' '.join(map(str, quals))))
            if args.collapse and args.pipeline == 'mothur':
                names_output.write('%s\t%s\n'%(header, ','.join(names_info)))
            return 0, 0, 0
        else:
            if args.output_format == 'fastq':
                bad_contig_output.write('@%s\tuncert > %.3f\n%s\n+\n%s\n'%(header, args.uncert, sequence, ''.join([chr(qual + args.fastq_offset) for qual in quals])))
            else:
                bad_contig_output.write('>%s\tuncert > %.3f\n%s\n'%(header, args.uncert, sequence))
                bad_qual_output.write('>%s\tuncert > %.3f\n%s\n'%(header, args.uncert, ' '.join(map(str, quals))))
            if args.collapse:
                if args.pipeline == 'mothur':
                    bad_names_output.write('%s\t%s\n'%(header, ','.join(names_info)))
                return len(names_info), 0, 0
            else:
                return 1, 0, 0


class ReturnedNaNError(Exception):
    """This error should not happen anymore, but just in case we're still checking for it"""
    def __init__(self, header):
        self.header = header
        self.__str__()
    def __str__(self):
        return 'Error calculation returned NaN for sequence %s. If using a C implementation, try switching to the python one instead.'%self.header
#
#
################################################################################################################



################################################################################################################
#moira utils library
#
#The functions below are meant to be part of an independent library for sequence processing.
#They are included here so this script can be used as a stand-alone.
#
################################################################################################################

class UnpairedFilesError(Exception):
    def __init__(self, lffilename, lqfilename, rffilename = None, rqfilename = None):
        self.lffilename = lffilename
        self.lqfilename = lqfilename
        self.rffilename = rffilename
        self.rqfilename = rqfilename
    def __str__(self):
        if not self.rffilename and not self.rqfilename:
            return 'You must provide at least a forward fasta and quality file'
        else:
            return 'If reading from paired-end files, you must provide fasta and quality files for both the forward and the reverse reads'


class NameMismatchError(Exception):
    def __init__(self, lfheader, lqheader, rfheader = None, rqheader = None):
        self.lfheader = lfheader
        self.lqheader = lqheader
        self.rfheader = rfheader
        self.rqheader = rqheader
        self.__str__()
    def __str__(self):
        if not self.rfheader:
            return 'Fasta header does not match Qfile header. Offending headers were: FASTA: %s   QFILE: %s'%(repr(self.lfheader), repr(self.lqheader))
        elif not self.lqheader:
            return 'Header mismatch. Offending headers were: forward_fastq: %s,   reverse_fastq %s'%(repr(self.lfheader), repr(self.rfheader))
        else:
            return 'Header mismatch. Offending headers were: forward_fasta: %s   forward_qual %s   reverse_fasta %s   reverse_qual %s'%(
                    repr(self.lfheader), repr(self.lqheader), repr(self.rfheader), repr(self.rqheader))


class LengthMismatchError(Exception):
    def __init__(self, fheader = None, ffilename = None, qfilename = None):
        self.fheader = fheader
        self.ffilename = ffilename
        self.qfilename = qfilename
        self.__str__()
    def __str__(self):
        if None not in (self.fheader, self.ffilename, self.qfilename):
            return 'Error reading sequence %s in files %s and %s. Sequence length and quality length do not match'%(repr(self.fheader), repr(self.ffilename),
                                                                                                                     repr(self.qfilename))
        elif not self.qfilename:
            return 'Error reading sequence %s in file %s. Sequence length and quality length do not match'%(repr(self.fheader), repr(self.ffilename))
        else:
            return 'Sequence and qualities are of different lengths.'


class EmptySeqError(Exception):
    def __init__(self, fheader, filename):
        self.fheader = fheader
        self.filename = filename
        self.__str__()
    def __str__(self):
        return 'Error reading file %s. Sequence %s was empty'%(repr(self.filename), repr(self.fheader))


class EmptyQualError(Exception):
    def __init__(self, qheader, filename):
        self.qheader = qheader
        self.filename = filename
        self.__str__()
    def __str__(self):
        return 'Error reading file %s. Quality %s was empty'%(repr(self.filename), repr(self.qheader))

    
def open_input(filename):
    """
    Make a best-effort guess as to how to open the input file.
    Deals with gzip and bzip2 compressed files.

    Function modified from screed's open_reader
    """
    file_signatures = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
    }
    try:
        bufferedfile = io.open(file=filename, mode='rb', buffering=8192)
    except IOError as e:
        print e
        print
        sys.exit(1)
    num_bytes_to_peek = max(len(x) for x in file_signatures)
    file_start = bufferedfile.peek(num_bytes_to_peek)
    compression = None
    for signature, ftype in file_signatures.items():
        if file_start.startswith(signature):
            compression = ftype
            break

    if compression is 'bz2':
        return bz2.BZ2File(filename=filename)
    elif compression is 'gz':
        if not bufferedfile.seekable():
            raise IOError('Unable to stream gzipped data.')
        return gzip.GzipFile(filename=filename)
    else:
        return bufferedfile


def parse_fasta_and_qual(forward_fasta_data, forward_qual_data, reverse_fasta_data = None, reverse_qual_data = None):

    """
    Parse fasta and qual files.
    If paired files are provided, will check if the headers match.
    Expects sequences and qualities to be stored in a single line.
    """

    #If reverse_fasta_data and reverse_qual_data are both None it will work in SingleEnd mode.
    if not forward_fasta_data or not forward_qual_data:
        raise UnpairedFilesError(forward_fasta_data, forward_qual_data, reverse_fasta_data, reverse_qual_data)
    if not reverse_fasta_data and reverse_qual_data or reverse_fasta_data and not reverse_qual_data:
        raise UnpairedFilesError(forward_fasta_data, forward_qual_data, reverse_fasta_data, reverse_qual_data)
 
    while 1:

        forward_fasta_header = forward_fasta_data.readline()
        forward_qual_header = forward_qual_data.readline()
        if reverse_fasta_data: reverse_fasta_header = reverse_fasta_data.readline()
        if reverse_qual_data: reverse_qual_header = reverse_qual_data.readline()

        if reverse_qual_data and reverse_fasta_data:
            if not forward_fasta_header and not reverse_fasta_header and not forward_qual_header and not reverse_qual_header:
                break
        else:
            if not forward_fasta_header and not forward_qual_header:
                break

        forward_fasta_header = forward_fasta_header.strip().replace('\t', ' ').split(' ')[0].lstrip('>').replace(':', '_')
        forward_sequence = forward_fasta_data.readline().strip()
        forward_qual_header = forward_qual_header.strip().replace('\t', ' ').split(' ')[0].lstrip('>').replace(':', '_')
        forward_quals = map(int, forward_qual_data.readline().strip().replace('\t', ' ').split(' '))
        if reverse_fasta_data: reverse_fasta_header = reverse_fasta_header.replace('\t', ' ').strip().split(' ')[0].lstrip('>').replace(':', '_')
        if reverse_fasta_data: reverse_sequence = reverse_fasta_data.readline().strip()
        if reverse_qual_data: reverse_qual_header = reverse_qual_header.strip().replace('\t', ' ').split('\t')[0].split(' ')[0].lstrip('>').replace(':', '_')
        if reverse_qual_data: reverse_quals = map(int, reverse_qual_data.readline().strip().replace('\t', ' ').split(' '))
   
        if reverse_fasta_data and reverse_qual_data and len(set([forward_fasta_header, reverse_fasta_header, forward_qual_header, reverse_qual_header])) != 1:
            raise NameMismatchError(forward_fasta_header, forward_qual_header, reverse_fasta_header, reverse_qual_header)
        if not reverse_fasta_data and not reverse_qual_data and len(set([forward_fasta_header, forward_qual_header])) != 1:
            raise NameMismatchError(forward_fasta_header, forward_qual_header)
        if not forward_sequence:
            raise EmptySeqError(forward_fasta_header, forward_fasta_data.name)
        if not forward_quals:
            raise EmptyQualError(forward_qual_header, forward_qual_data.name)
        if reverse_fasta_data and not reverse_sequence:
            raise EmptySeqError(reverse_fasta_header, reverse_fasta_data.name)
        if reverse_qual_data and not reverse_quals:
            raise EmptyQualError(reverse_qual_header, reverse_qual_data.name)
        if len(forward_sequence) != len(forward_quals):
            raise LengthMismatchError(forward_fasta_header, forward_fasta_data.name, forward_qual_data.name)
        if reverse_fasta_data and reverse_qual_data and len(reverse_sequence) != len(reverse_quals):
            raise LengthMismatchError(reverse_fasta_header, forward_fasta_data.name, forward_qual_data.name) 
        if reverse_fasta_data and reverse_qual_data:
            yield forward_fasta_header, forward_sequence, forward_quals, reverse_sequence, reverse_quals
        else:
            yield forward_fasta_header, forward_sequence, forward_quals, None, None


def parse_fastq(forward_fastq_data, reverse_fastq_data = None, fastq_offset = 33):
    """
    Parse fastq files.
    If paired files are provided, will check if the headers match.
    """

    if reverse_fastq_data:
        data = izip(forward_fastq_data, reverse_fastq_data)
    else:
        data = forward_fastq_data

    forward_buffer = []
    reverse_buffer = []
    
    for line in data:
        
        if reverse_fastq_data:
            forward_buffer.append(line[0].strip())
            reverse_buffer.append(line[1].strip())
        else:
            forward_buffer.append(line.strip())
            
        if len(forward_buffer) == 4:
            forward_header = forward_buffer[0].replace('\t', ' ').split(' ')[0].lstrip('@').replace(':', '_')
            forward_sequence = forward_buffer[1]
            forward_quals = [ord(x) - fastq_offset for x in forward_buffer[3]]
            if not forward_sequence:
                raise EmptySeqError(forward_header, forward_fastq_data.name)
            if not forward_quals:
                raise EmptyQualError(forward_header, forward_fastq_data.name)              
            if len(forward_sequence) != len(forward_quals):
                raise LengthMismatchError(forward_header, forward_fastq_data.name)
            forward_buffer = []
            
            if reverse_fastq_data:
                reverse_header = reverse_buffer[0].replace('\t', ' ').split(' ')[0].lstrip('@').replace(':', '_')
                reverse_sequence = reverse_buffer[1]
                reverse_quals = [ord(x) - fastq_offset for x in reverse_buffer[3]]
                if not reverse_sequence:
                    raise EmptySeqError(forward_header, forward_fastq_data.name)
                if not reverse_quals:
                    raise EmptyQualError(forward_header, forward_fastq_data.name)              
                if len(reverse_sequence) != len(reverse_quals):
                    raise LengthMismatchError(reverse_header, reverse_fastq_data.name)

                if forward_header != reverse_header:
                    raise NameMismatchError(forward_fastq_header, None, reverse_fastq_header, None)
                reverse_buffer = []

                yield forward_header, forward_sequence, forward_quals, reverse_sequence, reverse_quals

            else:
                yield forward_header, forward_sequence, forward_quals, None, None


def reverse_complement(sequence, quals = None):
    """Returns the reverse complement of a sequence."""
    
    complementary_matrix = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 
                            'W':'W', 'S':'S', 'R':'Y', 'Y':'R', 'M':'K', 
                            'K':'M', 'B':'V', 'V':'B', 'D':'H', 'H':'D',
                            '-':'-', '.':'.'}
    ###Check parameters:
    sequence = str(sequence)
    if quals: 
        quals = list(quals)
        map(int, [x for x in quals if x not in ('.', '-')]) #Should work even for aligned quality scores.
    if quals and len(sequence.replace('-', '').replace('.', '')) != len([x for x in quals if x not in ('.', '-')]):
        raise LengthMismatchError
    ###

    reverse = sequence[::-1] #Reverse string
    reverse_complement = []
    for base in reverse:
        try:
            reverse_complement.append(complementary_matrix[base])
        except KeyError:
            raise ValueError('"%s" is not a recognizable IUPAC-coded base.'%base)

    if quals:
        quals.reverse()
        return ''.join(reverse_complement), quals
    else:
        return ''.join(reverse_complement)


def nw_align(seq_1, seq_2, match, mismatch, gap, refine_overlap = True, verbose = False):
    """
    Globally align two sequences using the Needleman-Wunsch algorithm.
    Note that there is no gap_extend penalty in NW. Mothur charges the gap_open penalty for every gapped position.
    If refine_overlap=True, mothur's solution for aligning sequences with big non-overlapping regions will be applied.
    - First row and column will be initialized to 0 instead of applying the gap penalty.
    - Last row or column will be modified with by the nw_overlap function to insert the correct number of gaps at 3'.
    """

    #Populate matrices.
    seq_1 =' %s'%seq_1 #We need an extra space at the beginning for building the matrix.
    seq_2 =' %s'%seq_2
    score_matrix = []
    pointer_matrix = [] #Will contain tuples with the displacement in rows and columns (i.e. (0, -1) for going left).
    pointers = {(0, 0) : 'x', (-1, -1) : 'd', (0, -1) : 'l', (-1, 0) : 'u'}

    #Start with first row and column.
    for i in range(len(seq_1)):
        if refine_overlap:
            score_matrix.append([0])
        else:
            score_matrix.append([i * gap])
        pointer_matrix.append([(-1, 0)])
    for j in range(1, len(seq_2)):
        if refine_overlap:
            score_matrix[0].append(0)
        else:
            score_matrix[0].append(j * gap)
        pointer_matrix[0].append((0, -1))

    #Then the rest.
    for i in range(1, len(seq_1)):
        for j in range(1, len(seq_2)):
            if seq_1[i] == seq_2[j]:
                diag_score = score_matrix[i-1][j-1] + match
            else:
                diag_score = score_matrix[i-1][j-1] + mismatch

            up_score = score_matrix[i-1][j] + gap
            left_score = score_matrix[i][j-1] + gap

            if diag_score >= up_score:
                if diag_score >= left_score:
                    score_matrix[i].append(diag_score)
                    pointer_matrix[i].append((-1, -1))
                else:
                    score_matrix[i].append(left_score)
                    pointer_matrix[i].append((0, -1))
            else:
                if up_score >= left_score:
                    score_matrix[i].append(up_score)
                    pointer_matrix[i].append((-1, 0))
                else:
                    score_matrix[i].append(left_score)
                    pointer_matrix[i].append((0, -1))

    if refine_overlap:
        nw_overlap(score_matrix, pointer_matrix)

    if verbose:
        print 'Needleman-Wunsch matrix.'
        if refine_overlap:
            print '(Mothur overlap corrections were applied)'
        print '\t' + '\t\t'.join(seq_2)
        for i in range(len(score_matrix)):
            print seq_1[i], '\t'.join([str((score_matrix[i][j], pointers[pointer_matrix[i][j]])) for j in range(len(score_matrix[0]))])
        print
 
    #Traceback:
    score = 0
    seq_1_aligned = []
    seq_2_aligned = []
    i, j = len(seq_1) - 1, len(seq_2) - 1

    while i > 0 or j > 0:

        score += score_matrix[i][j]
        pointer = pointer_matrix[i][j]
        if pointer[0] == -1:
            seq_1_aligned.append(seq_1[i])
        else:
            seq_1_aligned.append('-')
        if pointer[1] == -1:
            seq_2_aligned.append(seq_2[j])
        else:
            seq_2_aligned.append('-')

        i += pointer[0] #Note that these values are 0 or negative.
        j += pointer[1]

    seq_1_aligned = ''.join(reversed(seq_1_aligned)) #Since we filled it backwards.
    seq_2_aligned = ''.join(reversed(seq_2_aligned))

    return seq_1_aligned, seq_2_aligned, score


def nw_overlap(score_matrix, pointer_matrix):

    """
    Clean up the 3' part of the alignment as mothur does, in order to avoid having a lot of scattered
    bases in the 3' region.
    
    Their solution was to look for the cells with the highest scores in the last row and column, and force gaps
    by making previous cells in those line/column point towards them.

    According to mothur's documentation, issues in the 5' region are solved by themselves during traceback.

    Note that they also initialize the first row and column of the NW matrix to 0 instead of applying gap penalties.
    """
 
    best_in_last_column_score = None
    best_in_last_column_index = None
    for i in range(len(score_matrix)):
        cell_score = score_matrix[i][-1]
        if cell_score >= best_in_last_column_score:
            best_in_last_column_score = cell_score
            best_in_last_column_index = i

    best_in_last_row_score = None
    best_in_last_row_index = None
    for j in range(len(score_matrix[0])):
        cell_score = score_matrix[-1][j]
        if cell_score >= best_in_last_row_score:
            best_in_last_row_score = cell_score
            best_in_last_row_index = j
    
    #Decide which sequence will get the gaps:
    if (best_in_last_column_index, best_in_last_row_index) == ((len(score_matrix) - 1), (len(score_matrix[0]) - 1)):
        pass
    
    elif best_in_last_column_score > best_in_last_row_score:
        for i in range(len(score_matrix) -1, best_in_last_column_index, -1):
            pointer_matrix[i][-1] = -1, 0 
    else:
        for j in range(len(score_matrix[0]) - 1, best_in_last_row_index, -1):
            pointer_matrix[-1][j] = 0, -1


def make_contig(forward_aligned, forward_quals, reverse_aligned, reverse_quals, insert, deltaq, consensus_qscore, qscore_cap, trim_overlap):

    """
    Build a contig from two aligned reads and their corresponding quality scores.

    Expects qualities to be lists of integers.
    
    Internal gaps will be filled with their complementary base if its quality score is equal or greater than the insert parameter.
    
    In case of mismatch, the base with the highest score will be kept, unless the quality score difference between both bases
    is lesser than the deltaq parameter. In that case, an N will be inserted in that position.
    
    If the sum_qscores parameter is equal to True, the sum of both quality scores will be returned for matching bases. If not, only
    the greatest quality score will be returned.
    """

    def qual2prob(qual):
        return 10 ** (qual / (-10.0))

    def prob2qual(prob):
        return int(math.floor(-10 * math.log10(prob)))

    ###Check parameters:
    forward_aligned = str(forward_aligned)
    forward_quals = map(int, list(forward_quals))
    reverse_aligned = str(reverse_aligned)
    reverse_quals = map(int, list(reverse_quals))
    insert = int(insert)
    deltaq = int(deltaq)
    if consensus_qscore not in ('best', 'sum', 'posterior'):
        raise ValueError('consensus_qscore must be "best", "sum" or "posterior".')
    if len(forward_aligned.replace('-', '')) != len(forward_quals):
        raise LengthMismatchError
    if len(reverse_aligned.replace('-', '')) != len(reverse_quals):
        raise LengthMismatchError
    if insert <= 0:
        raise ValueError('insert must be a positive integer')
    if deltaq <= 0:
        raise ValueError('deltaq must be a positive integer')
    if qscore_cap < 0:
        raise ValueError('qscore_cap must be a non-negative integer')
    ###

    #Fit the qualities into the alignment.
    forward_quals_aligned = []
    unaligned_position = 0
    for base in forward_aligned:
        if base == '-':
            forward_quals_aligned.append(base)
        else:
            forward_quals_aligned.append(forward_quals[unaligned_position])
            unaligned_position += 1
 
    reverse_quals_aligned = []
    unaligned_position = 0
    for base in reverse_aligned:
        if base == '-':
            reverse_quals_aligned.append(base)
        else:
            reverse_quals_aligned.append(reverse_quals[unaligned_position])
            unaligned_position += 1


    #Determine the start and the end of both sequences.
    for index, base in enumerate(forward_aligned):
        if base != '-':
            forward_start = index
            break
    for index, base in enumerate(reverse_aligned):
        if base != '-':
            reverse_start = index
            break

    for index, base in reversed(list(enumerate(forward_aligned))):
        if base != '-':
            forward_end = index
            break
    for index, base in reversed(list(enumerate(reverse_aligned))):
        if base != '-':
            reverse_end = index
            break

    #Determine the overlapping region. Will trust the alignment and consider overlap starts at the first non-gap character.
    if forward_start < reverse_start:
        overlap_start = reverse_start
        overlap_end = forward_end
        seqs_reversed = False
    else:
        overlap_start = forward_start
        overlap_end = reverse_end
        seqs_reversed = True

    overlap_length = overlap_end - overlap_start
    
    contig = []
    contig_quals = []
    mismatches = 0
    gaps = 0

    for position in range(len(forward_aligned)):
        if position < overlap_start:
            if not trim_overlap:
                if seqs_reversed:
                    contig.append(reverse_aligned[position])
                    contig_quals.append(reverse_quals_aligned[position])
                else:
                    contig.append(forward_aligned[position])
                    contig_quals.append(forward_quals_aligned[position])
        elif position > overlap_end:
            if not trim_overlap:
                if seqs_reversed:
                    contig.append(forward_aligned[position])
                    contig_quals.append(forward_quals_aligned[position])
                else:
                    contig.append(reverse_aligned[position])
                    contig_quals.append(reverse_quals_aligned[position])
        else:
            if forward_aligned[position] == '-':
                gaps += 1
                if consensus_qscore == 'posterior':
                    contig.append('N')
                    contig_quals.append(2)
                elif int(reverse_quals_aligned[position]) > insert:
                    contig.append(reverse_aligned[position])
                    contig_quals.append(reverse_quals_aligned[position])
                else:
                    pass

            elif reverse_aligned[position] == '-':
                gaps += 1
                if consensus_qscore == 'posterior':
                    contig.append('N')
                    contig_quals.append(2)
                elif int(forward_quals_aligned[position]) > insert:
                    contig.append(forward_aligned[position])
                    contig_quals.append(forward_quals_aligned[position])
                else:
                    pass

            elif forward_aligned[position] == reverse_aligned[position]:
                contig.append(forward_aligned[position])
                if consensus_qscore == 'sum':
                    contig_quals.append(forward_quals_aligned[position] + reverse_quals_aligned[position])
                elif consensus_qscore == 'posterior':
                    p1, p2 = qual2prob(forward_quals_aligned[position]), qual2prob(reverse_quals_aligned[position])
                    post_prob = (p1 * p2 / 3) / (1 - p1 - p2 + (4 * p1 * p2 / 3))
                    contig_quals.append(prob2qual(post_prob))
                elif forward_quals_aligned[position] >= reverse_quals_aligned[position]:
                    contig_quals.append(forward_quals_aligned[position])
                else:
                    contig_quals.append(reverse_quals_aligned[position])

            else:
                mismatches += 1
                if consensus_qscore in ('best', 'sum'):
                    if abs(forward_quals_aligned[position] - reverse_quals_aligned[position]) < deltaq:
                        contig.append('N')
                        contig_quals.append(2)
                    else:
                        if forward_quals_aligned[position] >= reverse_quals_aligned[position]:
                            contig.append(forward_aligned[position])
                            contig_quals.append(forward_quals_aligned[position])
                        else:
                            contig.append(reverse_aligned[position])
                            contig_quals.append(reverse_quals_aligned[position])
                else: #consensus_qscore == 'posterior'
                    if forward_quals_aligned[position] == reverse_quals_aligned[position]:
                        contig.append('N')
                        contig_quals.append(2)
                    else:
                        if forward_quals_aligned[position] > reverse_quals_aligned[position]:
                            p1, p2 = qual2prob(forward_quals_aligned[position]), qual2prob(reverse_quals_aligned[position])
                            contig.append(forward_aligned[position])
                        else:
                            p2, p1 = qual2prob(forward_quals_aligned[position]), qual2prob(reverse_quals_aligned[position])
                            contig.append(reverse_aligned[position])
                        post_prob = p1 * (1 - p2 / 3) / (p1 + p2 - (4 * p1 * p2 / 3))
                        contig_quals.append(prob2qual(post_prob))
                        
    if qscore_cap:
        contig_quals = [qual if qual < qscore_cap else qscore_cap for qual in contig_quals]
    
    return ''.join(contig), contig_quals, overlap_length, gaps, mismatches


def calculate_errors_PB(sequence, quals, alpha):
    """
    Calculate the errors on a sequence using the Poisson binomial method (sum of bernoulli random variables).
    Expects quals to be a list of integers containing Phred quality scores.
    """

    def prob_j_errors(p, j, n): #Where p is the error probability, j is the number of errors and n the number of observations.
        if j > n:
            return 0 #The formula below would end returning 0.
        else:
            per_position_accum_probs = [(1 - p) ** n]   #For j == 0:
            for j in range (1, j+1):                    #For j >= 1:
                per_position_accum_probs.append(((n - j + 1) / float(j)) * (p / (1 - p)) * per_position_accum_probs[j - 1])
        return per_position_accum_probs[-1]

    def sum_of_binomials(j, k): #Where j is the number of errors and k is the position in the sequence.
        probability = 0
        for i in range(j+1):
            probability += prob_j_errors(error_probs[k], i, n) * per_position_accum_probs[j-i][k-1]            
            #Where error_probs[k] is the error probability of the k-th position.
            #Where per_position_accum_probs[j-i][k-1] is the probability that all the bases
            #from position 0 to k-1 had a total of j-i errors.
        return probability                                         


    ###Check parameters:
    sequence = str(sequence)
    quals = map(int, list(quals))
    alpha = float(alpha)

    if len(sequence) != len(quals):
        raise LengthMismatchError()
    if alpha <= 0 or alpha > 1:
        raise ValueError('Alpha must be between 0 (not included) and 1.')
    ###


    n = 1 #Bernouilli distribution. 
    Ns = 0
    error_probs = []

    for base, qscore in zip(sequence, quals):
        if qscore < 0:
            raise ValueError('Qualities must have positive values.')
        if base == 'N':
            Ns += 1
        else:
            error_probs.append(10**(qscore / -10.0))

    expected_errors = 0
    accumulated_probs = [0] #[0] so the first while loop is executed
    per_position_accum_probs = []
    
    if error_probs:
        while 1:
            per_position_accum_probs.append([])
            for k in range(len(error_probs)):
                if k == 0:
                    per_position_accum_probs[expected_errors].append(prob_j_errors(error_probs[k], expected_errors, n))
                else:
                    per_position_accum_probs[expected_errors].append(sum_of_binomials(expected_errors, k))
            probability = per_position_accum_probs[-1][-1]
            accumulated_probs.append(accumulated_probs[-1] + probability)
            if accumulated_probs[-1] > (1 - alpha):
                break
            else:
                expected_errors += 1
                
        expected_errors = interpolate(expected_errors - 1, accumulated_probs[-2], expected_errors, accumulated_probs[-1], alpha)

    else:
        expected_errors = 0

    return expected_errors, Ns


def calculate_errors_poisson(sequence, quals, alpha):
    """
    Calculate the errors in a sequence approximating the sum of bernouilli random variables to a poisson distribution.
    Expects quals to be a list of integers containing Phred quality scores.
    """
    
    ###Check parameters:
    sequence = str(sequence)
    quals = map(int, list(quals))
    alpha = float(alpha)

    if len(sequence) != len(quals):
        raise LengthMismatchError()
    if alpha <= 0 or alpha > 1:
        raise ValueError('Alpha must be between 0 (not included) and 1.')
    ###
   
    Lambda = 0
    Ns = 0

    for base, qscore in zip(sequence, quals):
        if qscore < 0:
            raise ValueError('Qualities must have positive values.')
        if base == 'N':
            Ns += 1
        else:
            Lambda += 10**(qscore / -10.0)

    accumulated_probs = [0]
    expected_errors = 0

    while 1:
        probability = (math.exp(-Lambda) * (Lambda ** expected_errors)) / (math.factorial(expected_errors))
        accumulated_probs.append(accumulated_probs[-1] + probability)
        
        if accumulated_probs[-1] > (1 - alpha):
            break
        else:
            expected_errors += 1

    expected_errors = interpolate(expected_errors - 1, accumulated_probs[-2], expected_errors, accumulated_probs[-1], alpha)
    
    return expected_errors, Ns
        

def calculate_errors_bootstrap(sequence, quals, alpha, bootstrap):
    """
    Calculate the errors in a sequence using the bootstrap method.
    Expects quals to be a list of integers containing Phred quality scores.
    """

    ###Check parameters:
    sequence = str(sequence)
    quals = map(int, list(quals))
    alpha = float(alpha)
    bootstrap = int(bootstrap)

    if len(sequence) != len(quals):
        raise LengthMismatchError()
    if bootstrap < 1:
        raise ValueError('Bootstrap must be a positive integer.')
    if alpha <= 0 or alpha > 1:
        raise ValueError('Alpha must be between 0 (not included) and 1.')
    ###

    results = []
    
    for i in range(bootstrap):
        errors = 0
        Ns = 0
        for base, qscore in zip(sequence, quals):
            if qscore < 0:
                raise ValueError('Qualities must have positive values.')
            if base == 'N':
                Ns += 1
            else:
                if random() <= 10 ** (qscore / (-10.0)): #Error probability for a given Sanger quality value.
                    errors += 1
        
        results.append(errors)

    expected_errors = percentile(results, (1 - alpha) * 100)
    
    return expected_errors, Ns


def interpolate(errors1, prob1, errors2, prob2, alpha):
    """
    Perform a linear interpolation in the errors distribution to return the number of errors that has an accumulated
    probability of 1 - alpha.
    """

    result = errors1 + ((errors2 - errors1) * ((1 - alpha) - prob1) / (prob2 - prob1))
    if result < 0:
        #Happens only for very-short high qual sequences in which the probability of having 0 errors is higher than 1 - alpha.
        result = 0
    return result


#
#
################################################################################################################
    
if __name__ == '__main__':
    sys.exit(main(parse_arguments()))
