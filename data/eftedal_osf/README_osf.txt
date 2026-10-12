
### Description of files ###
In our main text, we present 82 separate correlations between relatives on the averaged score on 
the three tests taken in the 8th grade. 
In this repository, there are files containing the corresponding correlations for several more phenotypes. 
There are eight different national tests each student takes as part of going through the Norwegian 
school system: there are three subjects – Norwegian reading comprehension, English reading 
comprehension, and Mathematics – and there are separate tests for the 5th, 8th, and 9th grade. 
As there is no English reading comprehension test in the 9th grade, this adds up to eight. We offer 
files containing correlations for the 83 categories of relatives for each single test. Additionally, we 
have files with averaged scores on the three tests in the 5th grade, the three tests in the 8th grade, and 
the two tests in the 9th grade (with separate files for correlations with and without exclusion of 
students who did not take all three tests, which can be found in the "remove missing" folder). We also 
have files with just the average of the math and norwegian tests in the 5th and 8th grade, to make correlations 
more directly comparable across the 5th, 8th, and 9th grades (since there is no english test in the 9th grade). 
To allow researchers some freedom to choose how to deal with effects of country of origin of adoptees, we 
offer files with correlations between adoptive relatives in their uncorrected form, which can then be combined 
with the other files if desired (files with these correlations are in the "ado_uncorr"-folder).
“analysis_ beyond additive” is an R-script for performing our analyses.
Files with single subject correlations are named “math”, “read”, or “eng”, depending on whether 
they are for math, Norwegian reading comprehension or English reading comprehension, respectively. 
Whether tests are for 5th, 8th, or 9th grade is indicated by a number following the name. “eng8” is 
then correlations for the english tests in the 8th grade, for example. Files containing correlations for 
averaged scores for several subjects, are named either “MRE” if they are for all three subjects, or “MR” if 
they are for math and reading only. Files with the correlations for adoptive relatives that are uncorrected 
for country of origin follow the same naming scheme but with “_ado_uncorr” added at the end.

The files have the following columns:
- corr: The estimated correlation
- SE: The standard error of the estimate of the correlation
- rel: The category of relative the correlation is for. The categories are:
	- mz: monozygotic twins 
	- dz: dizygotic twins
	- uz: twins of unknown zygosity	
	- sibs: siblings who are not twins
	- half sibs: half siblings
	- 1st cousins t/uz: 1st cousins related through same-sex twins of unknown zygosity
	- 1st cousins: 1st cousins not related through same-sex twins
	- half cousins: Half 1st cousins
	- 2nd cousins t/uz: 2nd cousins related through same-sex twins of unknown zygosity
	- 2nd cousins: 2nd cousins not related through same-sex twins		
	- step sibs: Pairs not related by blood who have one or more half siblings in common
	- cocousins: Pairs not related by blood who have one or more cousins in common
	- adoptive sibs: Siblings where at least one is adopted from outside Norway
	- adoptive 1st cousins: 1st cousins where at least one is adopted from outside Norway
	- adoptive 2nd cousins: 2nd cousins where at least one is adopted from outside Norway
	- adoptive cocousins: cocousins where at least one is adopted from outside Norway
- type: Takes the value maternal for relatives related through their mothers, paternal for relatives related through their fathers, and other for pairs of relatives who cannot be categorized as maternal or paternal
- pair: Valued ssf for same-sex female, ssm for same-sex male, and os for opposite sex pairs
- half: Indicator for half relatives, valued 1 for half siblings and half 1st cousins, and 0 otherwise
- npairs: Number of pairs of relatives in the category for that row
- weight: Inverse of the squared standard error
- ns: For biological relatives, this is their degree of relatedness (for ss twins and relatives through ss twins, where relatedness depends on an unknown zygosity, it takes the average of the two alternatives). For relatives-in-law, the value represents the combined distance to their connecting pair of mates.
- ril: Indicator for relatives-in-law, valued 1 for step siblings and co-cousins, and 0 otherwise
- ado: Indicator for adoptive relatives, valued 1 for adoptive siblings, adoptive 1st cousins, adoptive 2nd cousins, and adoptive co-cousins, and 0 otherwise.
- log_corr: The natural logarithm of the values in the corr column
