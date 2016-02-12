- v1.0.3
  - Corrected license information in the moira.py welcome message.
  - moira.py will now warn the user if reverse files are provided without adding the --paired flag.

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
