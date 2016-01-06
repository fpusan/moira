#Test with:

moira.py -ffq test1.fastq -rfq test2.fastq --paired
python test_moira.py
cat README_TEST.sh

################################################################################
#
#The expected output files are stored in the test_results folder.
#The output message should be:
#
# 1000 sequences processed in 8.4 seconds.                                                                                
#   - Kept 924 (92.40%) of the original sequences.
#   - 76 (7.60%) of the original sequences were discarded due to low quality.
#
################################################################################



