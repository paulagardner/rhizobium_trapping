Comparing our methods to Jorrin and Imperial 2015 paper:

Methods we have done so far: Filter reads with Trimmomatic (Kbase), align reads to Rlv reference genome, transform output with Sam Tools, detect SNPs (They used varscan, we used PoPoolation)

Methods they did, which we should do next: Filter SNPs bassed on requirements (FST > 0.1, p < 0.05), calculate euclidean distances and multidimensional scaling based on output

See below the line for details, but: 
We must also calculate coverage for filtering and QC steps, check correlatedness between variables measured, and plot them to determine appropriate model to use for the data overall

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





