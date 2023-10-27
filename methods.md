### Comparing our methods to Jorrin and Imperial 2015 paper:

### Methods we have done so far: Filter reads with Trimmomatic (Kbase), align reads to Rlv reference genome, transform output with Sam Tools, detect SNPs (They used varscan, we used PoPoolation)

### Methods they did, which we should do next: Filter SNPs bassed on requirements (FST > 0.1, p < 0.05), calculate euclidean distances and multidimensional scaling based on output

See below the line for details, but: 
### We must also calculate coverage for filtering and QC steps, check correlatedness between variables measured, and plot them to determine appropriate model to use for the data overall
- this paper's definition of coverage: each 'genome' from 100 'isolates' was sequenced to [2 - 3x](https://hyp.is/VTN0oHO9Ee6ol3_8rdRlrw/apsjournals.apsnet.org/doi/pdfdirect/10.1094/MPMI-09-14-0296-FI?hmac=1698296604-O35rEBBIyldRBPjgthHoXkAf28CuoP796%2FYc42kqZZI%3D) coverage. [This](https://hyp.is/VTN0oHO9Ee6ol3_8rdRlrw/apsjournals.apsnet.org/doi/pdfdirect/10.1094/MPMI-09-14-0296-FI?hmac=1698296604-O35rEBBIyldRBPjgthHoXkAf28CuoP796%2FYc42kqZZI%3D) wording is unclear to me, but I *think* they mean each host plant got a sample cultured, and they got pooled data from the culture. Repeat 100 times. SO, that would be analagous to us taking pooled data from a host plant(s) for a treatment, and repeating until we do all mutations (24 times)?
	- in that case, coverage is going to be calculated **per pool**. We do not have control over this figure at this point, but we must discover it
 	- they also calculated another coverage figure- how much [average coverage was there for sequences ** present in the reference genome? **](https://hyp.is/VTN0oHO9Ee6ol3_8rdRlrw/apsjournals.apsnet.org/doi/pdfdirect/10.1094/MPMI-09-14-0296-FI?hmac=1698296604-O35rEBBIyldRBPjgthHoXkAf28CuoP796%2FYc42kqZZI%3D)
  	- they looked at differences on a plasmid by plasmid, and regional basis. Differences in coverage of such regions was something they used [to infer *proportion of actual plasmids* present](](https://hyp.is/sEuT1nO9Ee6PV19rJnJbYw/apsjournals.apsnet.org/doi/pdfdirect/10.1094/MPMI-09-14-0296-FI?hmac=1698296604-O35rEBBIyldRBPjgthHoXkAf28CuoP796%2FYc42kqZZI%3D)), it seems? As in, lower coverage = less plasmids were in the bacterial cells? (check necessary: how variable are plasmid numbers in bacteria of this type?) contrast with what they say [about regional/genic differences](https://hyp.is/9WVbZHO9Ee6xYmfHC6gAOA/apsjournals.apsnet.org/doi/pdfdirect/10.1094/MPMI-09-14-0296-FI?hmac=1698296604-O35rEBBIyldRBPjgthHoXkAf28CuoP796%2FYc42kqZZI%3D)
  	- this paper does what many others do: make a histogram of the coverage [illumina reference](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage.html), both for the whole genome and region by region. (we probably want to hone in on regions-- do we expect to see our differences manifest there?

### Their framing of the biology of the system (population genotype differences reflecting selection by host plant) is part of what we're doing, but not all of it. We might want to think about ways to justify in the analysis/plotting how to distinguish between when things might be top-down (selection from plant) vs. bottom-up (variables influence which bacteria are recruited). If not, we can speculate and compare/contrast with this paper's methodology. This will probably involve looking closely at regional differences. (could we perhaps even speculate on functional differences, if they are linked to functional genes???????????)

-----------------------------------
basic things to keep in mind:

### Coverage
the coverage of your reads affects the allele frequencies you can detect, and CHANGES in allele frequencies you can detect, in your dataset. 
SO, we must estimate coverage to be more informed on how we're going to filter our data, how we can do our statistics, etc. 
- The postdoc I talked to menioned checking depth per individual. That is a method for organismal data that surely can't be the marker for microbial data. the 2015 paper must address that(??), or otherwise just include a method that makes sense for microbes
- coverage = reads x read length, (and if paired, x2) and can maybe scale that against genome size

### FST
- I believe there is some kind of fancy FST that can be done for pooled sequencing data

### Modeling possibilities
Book I am referencing: Statistics for terrified biologists

MIXED LINEAR MODEL 
- linear models assume independent predictors. (I have a note: "ways to get around it: look up "feature selection"")
- correlations between variables guide from Christopher Fiscus: https://benchling.com/s/etr-Yl2XOKr7iGa5lMIPylqN/edit
		- If I remember correctly, the idea is to end up with the matrix at the end with the normal-looking distributions things. 
	 	  the idea is that you plot them against each other and display spearman's correlation coefficient to see if they count as 'independent' ( a requirement)
		  and the distributions show you if the variable (against itself) is roughly normally distributed or not (another requirement for linear models
		- suggestions: gally in ggplot for matrix of correlation values (again, like figure at end) You can also transform the variables if they are 
		  NOT normally distributed, with log or sqrt(log), but this risks losing variance in the data. 
		- corrplot. It requires a symmetrical matrix, can use geom.tile to do it in ggplot

RANDOM FOREST
- a ML approach
- can be used if there are a TON of variables, or if what you're trying to model is nonlinear! (cubic, quadratic, etc). 





