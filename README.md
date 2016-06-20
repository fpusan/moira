moira
=====

Quality-filter raw sequence reads using the Poisson binomial filtering algorithm.
The moira.py script and the Poisson binomial filtering algorithm are now described in the following paper:

    Puente-Sánchez, F., Aguirre, J., & Parro, V. (2016).
    A novel conceptual approach to read-filtering in high-throughput amplicon sequencing studies.
    Nucleic acids research, 44 (4): e40.

which can be accessed online at: http://nar.oxfordjournals.org/content/44/4/e40


INSTALLATION INSTRUCTIONS:

- moira is available as a pip-installable python package. To install it just run *pip install moira* in a terminal. Alternatively, download and unzip the distribution package and run *python setup.py install test*.

- The moira.py contains the python implementation of the Poisson binomial algorithm. It will perform as a standalone script as described here.

- bernoullimodule.c contains the C implementation of the Poisson binomial filtering algorithm. It's written to work as a python extension module, and will speed up moira.py if compiled as a shared library. An example command line for compiling it would be:

        gcc -fpic -shared -I /usr/include/python2.7/ -o bernoulli.so bernoullimodule.c

- nw_align.pyx contains a cython implementation of the Needleman-Wunsch aligner, in order to speed up contig construction from paired reads. An example command line for compiling it would be:

        cython nw_align.pyx
        gcc -fpic -shared -I /usr/include/python2.7/ -o nw_align.so nw_align.c

- Manual compilation of C extensions is only needed if not using pip or the setup.py install script.

- Testing can be performed by downloading the distribution package (or cloning the github repository) and running *python setup.py test*. A correct installation should pass all tests. In case that the C modules have not been compiled some of the tests will be skipped. Note that the moira.py script will nevertheless be able to run, but it will do so at slower speeds. If the program is not working in your platform (or everything works, but you have questions anyway) please do not hesitate to contact me at fpusan@gmail.com.


REQUIREMENTS:

- Expects that input sequences (single or paired) and qualities are in the same order.
- Expects that sequences and qualities are stored only in one line (i.e. >header\\nsequence\\n>header2\\nsequence2).
- Requires numpy.


USAGE:

  - Make contigs from paired reads (fasta + qual) without quality-filtering:

        moira.py --forward_fasta=<FILE> --forward_qual=<FILE> --reverse_fasta=<FILE> --reverse_qual=<FILE> --paired --only_contig

  - Make contigs from paired reads (fastq) and perform quality-filtering, output results in fastq format:

        moira.py --forward_fastq=<FILE> --reverse_fastq=<FILE> --paired --consensus_qscore posterior --output_format fastq

  - Quality-filter already assembled contigs or single reads:

        moira.py --forward_fasta=<FILE> --forward_qual=<FILE>

  - Input files can be compressed with gzip or bzip2. This will be automatically detected by the script.


OUTPUT:

  - If quality control is being performed, files will be generated with both the sequences that passed the QC and the ones that didn't. A brief report will be included on the headers of the contigs that didn't pass the QC.

        <PREFIX>.qc.good.fasta
        <PREFIX>.qc.good.qual
        <PREFIX>.qc.bad.fasta
        <PREFIX>.qc.bad.qual

  - Else, only two files will be generated.

        <PREFIX>.contigs.fasta
        <PREFIX>.contigs.qual
   
  - The default prefix is the forward input name without the extension. A custom prefix can be specified with the --output_prefix option. 
  - If --output_format is set to "fastq", fastq files will be generated instead of fasta + qual files.
  - Output can also be compressed in the gzip or bz2 formats via the --output_compression parameter.
  - If identical sequences are being collapsed, mothur-formatted name files (or USEARCH formatted sequence headers) will also be generated.
  - moira.py will replace ':' for '_' in sequence names for compatibility with the mothur pipeline.


PARAMETERS:

  - Needleman-Wunsch aligner parameters:
    - --match (default 1): match score
    - --gap (default -2): gap penalty
    - --mismatch (default -1): mismatch penalty

  - Contig constructor parameters:
 
   - --insert (default 20): quality above which a base will be used for filling a complementary gap or ambiguity.
   - --deltaq (default 6): minimum quality difference allowed between two mismatched bases for not including an N in the consensus sequence.
   - --consensus_qscore (default 'best')
     - best: use the best quality on each position of the alignment as the consensus quality score (Unless an ambiguity is introduced in that position by the contig constructor. In that case, the reported quality score will be always 2.
     - sum: in matching bases, consensus quality score will be the sum of the qualities of both reads in that position of the alignment.
     - posterior: use Edgar & Flyvbjerg's (2015) method for calculating consensus quality scores. The insert and deltaq parameters will be ignored. Ambiguities will be introduced in gaps, or if two mismatched bases have exactly the same quality score. In that case, the reported quality score will be always 2.
   - --qscore_cap (default 40): Maximum consensus quality score to report. Higher consensus quality scores will be trimmed to the value of --qscore_cap. Setting it to 0 will remove the cap.
   - --trim_overlap (default False): trim the contig to the overlapping region.

  - Quality-filtering parameters:

    - --collapse (default True): if True, identical sequences will be collapsed before quality control, and the one with the best quality will be used as a representative of the whole group.
    - --error_calc (default 'poisson_binomial'): algorithm used for error calculation.
      - poisson_binomial: calculate the Poisson binomial distribution (sum of bernoulli random variables).
      - poisson: approximating sum of bernoulli random variables to a poisson distribution.
      - bootstrap: numerical generation of an error distribution (deprecated).
    - --ambigs (default treat_as_error): handling of ambiguous positions during quality checking.
      - treat_as_error: will consider than ambiguities always result in a misread base.
      - disallow: will discard sequences with ambiguities.
      - ignore: will ignore ambiguities.
    - --truncate (default None): truncate sequences to a fixed length before quality control. Discard sequences smaller than that length.
    - --min_overlap (default None): discard contigs with less than the specified overlap length. Will be ignored if using single-end reads.
    - --round: round down the predicted errors to the nearest integer prior to filtering.
    - --uncert (default 0.01): maximum divergence of the observed sequence from the original one due to sequencing errors.
    - --maxerrors (no default value): Maximum errors allowed in the sequence. Will override --uncert if specified as a parameter.
    - --alpha (default 0.005): probability of underestimating the actual errors of a sequence.
    - --bootstrap (default 100): number of replicates per position used for error calculation by the bootstrap method.
        
  - Other:

    - --paired: input files are paired end files and will be assembled into contigs.
    - --only_contigs: assemble contigs, don\'t do quality control.
    - --relabel (default False): if a prefix string is introduced, sequential labels will be generated for the sequences,with the format <prefix>N, where N=1,2,3...
    - --output_format (default fasta):
      - fasta: output files in fasta + qual format.
      - fastq: output files in fastq format.
    - --output_prefix: Prefix for the output file names.
    - --output_compression (default none): Compression applied to output files. Can be 'none', 'bz2', or 'gzip'.
    - --pipeline (default mothur):
     - mothur: output for collapsed sequences will be in mothur\'s fasta + names format.
     - USEARCH: output for collapsed sequences will be in a single fasta file, with abundance information stored in the sequence header.
    - --fastq_offset (default 33): ASCII/qscore encoding.
    - --processors (default 1): number of processes to use.
    - --silent: Do not print welcome, progress and goodbye messages. Warnings will still be printed.


COMMENTS:

   - Alignment parameters are set to replicate mothur's default implementation of the Needleman-Wunsch algorithm.

   - The 'insert' and 'deltaq' parameters from mothur make.contig are also reproduced. They are set at their default values. More details can be found at www.mothur.org/wiki/Make.contigs. They will be applied only if the --consensus_qscore parameter is set to 'best' or 'sum'.

   - We now provide the option to calculate posterior error probabilities for reporting consensus quality scores, as described by Edgar & Flyvbjerg (2015). This can be achieved by adding the '--consensus_qscore posterior' or '-q posterior' flags, as shown in the example above. While this method is theoretically sound, it only improved the results obtained by using our default '-q best' option for some of our mock community samples.

   - Approximating the sum of bernoulli random variables to a poisson distribution is quicker than calculating their exact sum (Poisson binomial distribution). It proves specially useful for long reads (>500 nt). That said, the Poisson binomial filtering algorithm is also implemented in C and even the python implementation is quick enough or processing large datasets. The bootstrap method (--error_calc bootstrap) is a numerical algorithm for performing the sum of bernoulli random variables. It is only included for testing purposes.

   - Quality-filtering will discard the contigs expected to have more than 'alpha' chances of diverging from the original sequence more than the value specified by the 'uncert' param. That means that, during distance calculation between two given sequences, the observed distance will be at most 'dist + 2\*uncert', where 'dist' is the original distance between those sequences without sequencing errors. Thus, a good rule of thumb would be considering the effective OTU clustering distance to be actually 'OTUdist - 2\*uncert', where OTUdist is the distance used for clustering the observed sequences.


Distributed under the Modified BSD license.
