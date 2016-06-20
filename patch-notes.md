- v1.3.0
  - Fixed a bug in which the '>' character was not added to the beginning of sequence headers when using the --only_contig option.
  - Fixed a bug in the bernoulli library that caused error calculation for very short sequences to randomly return NaNs. This bug was not passing silently, so the results obtained with the previous versions of moira can still be trusted.
  - Added the option to keep only the overlapping region when assembling paired-end reads. This can be achieved with the --trim_overlap flag.
  - Added the option to discard contigs with a small overlapping region, via the --min_overlap parameter.

- v1.2.1
  - Fixed a bug that was preventing the --maxerrors option to work (contributed by Christopher Thornton, University of Utah).

- v1.2.0
  - Added support for posterior consensus qscore calculation during contig construction, as described by Edgar & Flyvbjerg (2015). This can be achieved by using the "--paired" and "-q posterior" flags.
  - Added the --qscore_cap option for limiting the maximum consensus qscore reported during contig construction. It's default value is 40, which should result in backwards compatibility when using the old consensus qscore calculation options. It can be nonetheless disabled by setting it to 0.
  - The documentation and examples have been updated to reflect these changes.
  - Optimized the cython code for the Needleman-Wunsch aligner. This should result in a noticeable speed boost when running moira.py on paired-end reads.
  - The script will now print warning messages when using deprecated options.
  - Fixed an error that raised when using the -e bootstrap flag.

- v1.1.0
  - Corrected license information in the moira.py welcome message.
  - Fixed a bug in which the script was not falling back to the pure python implementation of the Needleman-Wunsch algorithm if the cython version was not present.
  - moira.py will now warn the user if reverse files are provided without adding the --paired flag.
  - Fixed a bug in which an additional '>' was inserted at the beginning of sequence headers when using the --relabel option (contributed by Jing Wang, Shanghai Jiao Tong University).
  - Added support for gzip and bzip2 compressed input files. The script automatically detects the input file type (contributed by Christopher Thornton, University of Utah).
  - Added the option for choosing a custom output prefix, via the -op, --output_prefix flags (contributed by Christopher Thornton, University of Utah).
  - Added the option of producing gzip or bz2 compressed output files (--output_compression).
  - Added an option for running without printing welcome, progress or goodbye messages, via the --silent flag. Warning messages, if any, will still be printed.
  - Testing can now be directly performed via the setup.py script (python setup.py test).

- v1.0.2b.
  - Corrected several typos in documentation strings.

- v1.0.2.
  - Updated citation info, ReadMe and moira.py usage information.
  - Updated info in setup.py
  - Added the script README_TEST.sh to facilitate checking that everything is working fine.
  - Corrected an error in the bernoulli C extension module (which performs high-speed Poisson binomial filtering) that prevented it from loading.
  - The bernoulli C extension module now internally checks that the values passed to it are correct, and returns exceptions accordingly.
	  - In previous versions it performed basic type checking during initial argument tuple parsing, but value checking had to be performed in python before calling the calculate_errors_PB C function.
	  - Now having an alpha not between 0 and 1 or a having sequences and qualities of different lenghts will raise a value error.
	  - Similarly, having non-integer elements in the qualities list will raise a type error.
	  - Having qvalues of 0 (happened in some artificial datasets) will lead to divisions by zero in the PB error calculation. The C extension now substitutes them for Q=1 in order to avoid this.
	  - This should facilitate the use of the bernoulli C extension module by other python scripts.

- v1.0.1.
  - The license was updated.
  - Fixed an error that was causing the calculate_errors_PB function from the moira.py script to raise an exception when running on a sequence with only Ns.
  - Fixed an error that was causing the calculate_errors_PB function from the bernoulli library to segfault when running on a sequence with only Ns.

- v1.0.0.
  - Initial release.
