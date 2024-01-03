
# WEEK 11 ####

# Homework 9:

# You are looking for an E. coli gene whose mRNA expression level (1) significantly changes 
# upon each of heat stress and chemical X treatments and (2) are significantly different from the 
# additive effects of the two treatments when the two treatments are simultaneously given. Based 
# on literature, you selected the gene AX3 as a candidate. You measured the gene AX3 expression 
# levels upon heat stress, chemical X, and both treatments and no treatment control in three 
# independent experiments. The data is in “AX3.csv” (which is hyperlinked here). AX3.csv has 
# three columns, treatment (none, chemX, heat, and both), experiment (exp1, 2, and 3), and exp.v 
# (expression level values).  Do you conclude that AX3 has the expression characteristics you are 
# looking for? Don’t forget to run diagnostics of the model.
# Hints: 
# 1. What are the right factors for a model when the above expression characteristics are 
# asked? Just the “treatment” factor does not work.
# 2. Do you need to transform the expression level values?

exp.dat = read.table('AX3.csv',sep=',',header=T)
summary(exp.dat)
head(exp.dat)
boxplot(exp.v ~ treatment,data=exp.dat)
boxplot(log(exp.v) ~ treatment,data=exp.dat)

chem.x = factor(exp.dat$treatment=='chemX' |
                  exp.dat$treatment =='both')

heat.t = factor(exp.dat$treatment == 'heat' |
                  exp.dat$treatment == 'both')

chem.x = factor(chem.x,labels = c('No','Yes'))
heat.t = factor(heat.t,labels = c('No','Yes'))

exp.dat = data.frame(exp.dat,chem.x,heat.t)
head(exp.dat)

lm.ax3 = lm(log(exp.v)~chem.x+heat.t+chem.x:heat.t+experiment,data=exp.dat)
summary(lm.ax3)

opar=par(mfrow=c(2,2))
plot(lm.ax3,which=1:4)
par(opar)

p.val = summary(lm.ax3)$coef[c(2,3,6),4]
p.adjust(p.val,method='holm')






# WEEK 9 ####

# Homework 8:

# (1)   Plot the number of candies versus the age (2 points). Add the linear regression line
# (2 points), the 99% confidence interval of the mean estimates, and the 99% prediction interval
# of the data (2 points) in different colors. Also include legends explaining the color codes
# (1 point). Hint: you should set the ylim values in the first plot to make space for the legends.

# (2)   Do you conclude that the number of candies the children collected is affected by their 
# ages? Why or why not? (3 points)

# Load data
candy = read.table('candies.csv',header=T,sep=',')
dim(candy)
head(candy)

# Plot number of candies versus age
plot(candies ~ age, data=candy,
     xlab="Age",ylab="Number of Candies",
     ylim=c(90,210))

# Add linear regression line
lm.candy = lm(candies ~ age,data=candy)

# Ages 3-15
x5 = seq(3,15,0.01)

# Predicting points for every number between 3 and 15
new.d = data.frame(age=x5)

# Add 99% confidence interval of mean estimates
?predict
pr.c = predict(lm.candy,newdata=new.d,int='c',level=0.99)

# Add 99% prediction interval of the data
pr.p = predict(lm.candy,newdata = new.d,int='p',level=0.99)
pr.mat = cbind(pr.c,pr.p[,2:3])

?legend()
?matplot
matplot(new.d,pr.mat,type='l',add=T,
        col=c('red','blue','blue','green','green'),
        lty=c(1,2,2,3,3))

# Add legend
legend('topright',legend=c('Regression line',
                           '99% conf int of model',
                           '99% pred int of data'),
       col=c('red','blue','green'), lty=c(1,2,3),
       cex=0.7)

summary(lm.candy)
# Since the p value for the age coefficient is large (0.948), the age value does not
# significantly affect the number of candies.

cor.test(candy$age,candy$candies)
# get same p value!

# WEEK 8 ####

# Produced: 334 red, 278 yellow, and 289 orange peppers 
# Sold: 28 red, 14 yellow, and 19 orange peppers

peppers <- matrix(c(334,278,289,28,14,19),ncol=2)

fisher.test(peppers)

# Homework 7

# You are analyzing a barley population consisting of 150 recombinant inbred lines (RILs; L001 to L150)
# derived from two parent inbred lines, Velvet (V) and Wisconsin No. 38 (W). Note that in RILs, 
# practically all the loci are homozygous. You are trying to find a DNA marker that can be used to breed
# a gene that confers resistance to a fungal disease. V is resistant to the disease, and W is not. You
# tested two DNA markers m363 and m612, both of which can distinguish two parent alleles. The table,
# “barleyRIL.txt”,   Download The table, “barleyRIL.txt”, contains the genotypes at each marker (either
# V or W alleles) and the disease phenotype (“S” for susceptible and “R” for resistant) for the RILs
# and the parents. Do you conclude that either (or both or neither) of the markers can be used for 
# breeding the disease resistance? Make it clear how you got to the conclusion.

barley <- read.delim('barleyRIL.txt',header=T)
head(barley)

count.t.m363=table(barley$m363,barley$phenotype)
count.t.m612=table(barley$m612,barley$phenotype)





fisher.test(count.t.m363,alternative='greater')

fisher.test(count.t.m612,alternative='greater')

# Computer Lab Problems 8(1)

# 1. (i)
seed <- read.delim('seed.txt')
head(seed)
lm.seed = lm(distance ~ diameter,data=seed)
#    (ii) intercept = 5.706 (when diameter = 0, distance = 5.706)
#         diameter coefficient = 0.8321

#    (iii) 
lm.df = data.frame(fitted = lm.seed$fitted,
                   resid = lm.seed$resid)

lm.df$c.res = apply(as.matrix(lm.df),1,sum)
lm.df$r.res = seed$distance
identical(round(lm.df$c.res,digits = 5),
          round(lm.df$r.res,digits = 5))

# 2.
library(ISwR)
?thuesen

lm.sv = lm(short.velocity ~ blood.glucose,data=thuesen)
lm.sv$coef
plot(lm.sv$fitted,na.omit(thuesen$short.velocity),
     xlab='Fitted value for short.velocity',
     ylab='Actual short.velocity value')
abline(0,1,col='red')
# How does this line show y = x?


# WEEK 7 ####

# Homework 6

# You are an apple wholesale buyer. Larger “Fireside” apples (larger by weight per apple)
# command higher retail prices per kg. You got offers of the same price per kg for Fireside
# apples from Farms A and B. It is known that the weights of apples from each Farm have a
# normal distribution with the standard deviation of 15 grams. You weighed each of 10 randomly
# chosen apples from each farm and ran a t-test with a 99% confidence level to compare the mean
# apple weights from the two farms. Unknown to you, the average weight of apples from Farm A is
# 6 grams lower than that from Farm B.

# (1)   What is the chance that you can make the right call that the mean weight of apples
#       from Farm B is higher than that from Farm A? (5 points)
?power.t.test
power.t.test(n=10,delta=6,sd=15,sig.level=.01)


# (2)   If you want to have the chance to make the right call higher than 70%, how many apples
#       from each farm at minimum do you need to weigh? (5 points)
power.t.test(delta=6,sd=15,sig.level=0.01,power=.7,type='two.sample')


# Computer lab problems 7(1)
# 1
samp.mat = c()
for (i in 1:1e4) {
  pop1.mean = runif(1,-1,3)
  samp1 = rnorm(5,mean=pop1.mean,sd=1.5)
  samp.mat = rbind(samp.mat,samp1)
}

p.vals = apply(samp.mat,1,function(x) t.test(x,mu=0)$p.value)

sum(p.adjust(p.vals,method='bonferroni') < 0.05)/1e4

sum(p.adjust(p.vals,method='holm') < 0.05)/1e4

sum(p.adjust(p.vals,method='BH') < 0.05)/1e4

hist(p.vals,freq=F)

# 2
binom.test(62,n=62+172,p=3/(3+13))

# Computer lab problems 7(2)
# 1
contingency.mat = matrix(c(43,100-43,92,200-92),ncol=2)
fisher.test(contingency.mat)

# 2
contingency.mat = matrix(c(21,86-21,231-21,32344-86-231+21),ncol=2)
fisher.test(contingency.mat)



# WEEK 6 ####

# HW 5: You got a developmental mutant of zebrafish. Many mutant individuals have shorter
# tails than wild type, but some mutant individuals have long tails. You measured the tail
# length (in mm) of multiple individuals from the mutant and wild type. The data is in the
# list “tail.length” in “fish.tail.RData” file. Test whether the mutant has a significantly
# shorter mean tail length than wild type. Note: (1) the distribution of the tail length of
# the mutant is not symmetric; (2) Since you decided to ask whether it is shorter after
# measuring some of them, it should be a 2-sided test. However, a 2-sided test could be
# complicated (by a bootstrap test), so it is OK if you declare to do a 1-sided test and
# do a 1-sided test (merely for the sake of the assignment).


# [3272] Test it using the best and justifiable test method among the non-resampling-based
# methods we learned.

# [5272] Test it using a bootstrap test.

summary(tail.length)
boxplot(tail.length)

# Bootstrap test:
ft.mut = tail.length$mutant
ft.wt = tail.length$ wild.type

boot.mean.diff = c()
for(i in 1:1e4) {
  ft.mb = sample(ft.mut,replace=T)
  ft.wb = sample(ft.wt,replace=T)
  mean.diff = mean(ft.wb) - mean(ft.mb)
  boot.mean.diff = c(boot.mean.diff,mean.diff)
}
plot(density(boot.mean.diff),
     xlab='Mean tail length difference (wt - mut')
abline(v=0,col='red')

# 1 sided test
sum(boot.mean.diff < 0) / 1e4

# For an alternative, log transformation:
boxplot(log(ft.mut),log(ft.wt))

shapiro.test(ft.mut)
shapiro.test(ft.wt)

shapiro.test(log(ft.mut))
shapiro.test(log(ft.wt))

t.test(log(ft.mut),log(ft.wt))

# WEEK 4 ####

# Make a plot that demonstrates that a Poisson distribution is a good approximation
# of a binomial distribution when n >> 0 and p << 1. To do this, choose an appropriate
# set of n and p values, like slide 56 of Ppt week3(1). (But note that both distributions
# in this assignment are discrete distributions.) Make a single file that contains: all
# the 'code' you used for this (it is fine to copy all the related from the console);
# the plot you made as a pasted picture (5 points); including a legend + plot title + 
# labeled axes that were inserted into your graphic using R (2 points); and an 
# accompanying text explanation of why this plot demonstrates that a Poisson distribution
# is a good approximation of a binomial distribution when n >> 0 and p << 1. Make the 
# explanation understandable to people who only see your plot and explanation - assume that
# those people do not see your R scripts (3 points).

x = 0:15
plot(x,dbinom(x,size=4000,p=0.001),type='h',lwd=3,
     main="Comparison of binomial and poisson distributions",
     xlab="X",ylab="Density")
points(x+0.2,dpois(x,lambda=4),type="h",lwd=3,col="red")
legend(6,0.2,legend=c("Binomial,n=4000,p=0.001","Poisson,lambda=4"),
       col=c("black","red"),lwd=3)

# Computer lab problems 4(1)
# 1. What is the t-value for a 5% upper tail when the degree of freedom is 4?
qt(0.95,4)

# 2. You were interested in whether the yield of a new potato variety is 
# different from the industry standard. The harvested potato weight per plant
# for 6 plants were 1.77, 1.20, 1.51, 1.42, 1.70, and 1.39 kg. You assume that
# the potato weight per plant normally distributes based on your previous experience.
# Do you conclude that the mean yield of this new variety is significantly higher than
# the industry standard of 1.32 kg per plant?

# i.	 Calculate the t-value
samp.vec = c(1.77, 1.20, 1.51, 1.42, 1.70, 1.39)
t.val = ( mean(samp.vec) - 1.32 ) /
  (sd(samp.vec)/sqrt(length(samp.vec)))
t.val

# ii.  What is the null hypothesis? Is this a one-sided or two-sided test case?
#      The mean of the new variety is equal to the industry standard.

# iii. Use (i) and (ii) to calculate the p-value that the null hypothesis is correct. 
#      Don’t use the t.test function.
2*( 1-pt(t.val, df=length(samp.vec) - 1) )

# iv.	 Do (iii) using the t.test function. Was the calculation in (iii) correct?
t.test(samp.vec,mu=1.32)

# v.   What’s your conclusion? Don’t forget to specify the confidence level (or the rejection rate).
# vi.	 What the standard weight per plant value should be to have the p-value of 0.01 with these
#      observed weight per plant values? Calculate them for both cases of higher and lower than 
#      the standard. Hint. What are the t-values should be to have the p-value of 0.01? 
#      Then think how the t-value is calculated. You can write a simple equation by having
#      the standard weight per plant value as x.
upper = qt(0.995, df=length(samp.vec) - 1)
upper *  sd(samp.vec)/sqrt(length(samp.vec)) + mean(samp.vec)
-upper *  sd(samp.vec)/sqrt(length(samp.vec)) + mean(samp.vec)
t.test(samp.vec, mu=upper *  sd(samp.vec)/sqrt(length(samp.vec)) + mean(samp.vec))
t.test(samp.vec, mu=-upper *  sd(samp.vec)/sqrt(length(samp.vec)) + mean(samp.vec))

# 3i. Draw a histogram with the react vector object in the ISwR package. Does it look like a sample
#     from a unimodal and symmetric population?
library(ISwR)
hist(react)
# ii. Using a t-test, test whether the mean of the population is significantly different from 0.
#     Is it significantly different from -1?
#     a. Do it using a t-distribution function.

t.val0 =  mean(react)  /
  (sd(react)/sqrt(length(react)))
t.val0

t.val_1 =  ( mean(react) + 1)  /
  (sd(react)/sqrt(length(react)))
t.val_1

2*(1-pt(abs(t.val0), df=length(react)-1))
#     b. Do it using the t.test function.
t.test(react,mu=0)
t.test(react,mu=-1)

# Computer lab problems week 4(2)

# 1.	Are the mean heights of the boys and girls significantly different? If so, 
#     which mean height is higher?
bg.h = read.table("height30.30.txt",header=T)
mean.diff = mean(bg.h$boy) - mean(bg.h$girl)
sedm = sqrt( sd(bg.h$boy)^2/length(bg.h$boy) + sd(bg.h$girl)^2/length(bg.h$girl) )
deg.f = length(bg.h$boy) + length(bg.h$girl) - 2
2*( 1-pt(abs(mean.diff)/sedm, deg.f) )
sign(mean.diff)

t.test(bg.h$boy,bg.h$girl)

# 2.	The time records for two 100-meter dash runners (U.B. and T.G.) in the last
#     5 competitions were: U.B., 9.66, 9.71, 9.71, 9.68, and 9.71 seconds; T.G., 
#     9.71, 9.76, 9.68, 9.77, and 9.70. Assume that their (true) times distribute
#     normally with the same SD. Do you conclude one is significantly faster than the other?

ub = c(9.66, 9.71, 9.71, 9.68, 9.71)
tg = c(9.71, 9.76, 9.68, 9.77, 9.70)
mean.diff = mean(ub) - mean(tg)
sedm = sqrt(var(ub)/length(ub) + var(tg)/length(tg))
deg.f = length(ub) + length(tg) - 2
2*(1-pt(abs(mean.diff/sedm), df=deg.f))
t.test(ub,tg)

# 3.	You want to compare the weights of eggs of a single sea bird species from two different
#     sites (sites A and B). You measured the weights (in grams) of 40 randomly chosen eggs in
#     each site. The data are in “sea.bird.egg.csv” in the Canvas site. Assume the weights for
#     each site distribute normally with the same SD.
egg.dat = read.table('sea.bird.egg.csv',header=T, sep=',')
head(egg.dat)
#     What is the sample mean difference between the sites?
#     Are the egg weights from two sites significantly different? If so, eggs from which site are heavier?
mean.diff = mean(egg.dat$site.A) - mean(egg.dat$site.B)  
sedm = sqrt(var(egg.dat$site.A)/length(egg.dat$site.A) +
              var(egg.dat$site.B)/length(egg.dat$site.B))
deg.f = length(egg.dat$site.A) + length(egg.dat$site.B) -2
2*pt(mean.diff/sedm, df=deg.f)
# 4.	[5272] Let’s simulate a 2-sample t-test case for no difference in the mean. Two normal
#     populations have the same mean value of 24 and the same SD value of 4.3. Make 1e4 samples
#     of size 7 for each, and run 1e4 tests between samples from the populations. How many 
#     percent of the tests resulted in positive with a 95% confidence level? “Positive” means
#     that the test result shows significantly difference.
samp.mat = matrix(rnorm(1e4*7*2, mean = 24, sd=4.3), ncol=7*2)
deg.f = 7*2 -2
p.vals = c()
for (samp.row in 1:nrow(samp.mat)) {
  vals = samp.mat[samp.row,]
  vals.A = vals[1:7]
  vals.B = vals[8:14]
  mean.diff = mean(vals.A) - mean(vals.B)
  sedm = sqrt(var(vals.A)/7 + var(vals.B)/7)
  p.v = 2*(1-pt(abs(mean.diff/sedm), df=deg.f))
  p.vals = c(p.vals, p.v)
}
sum(p.vals < 0.05) / length(p.vals)
# WEEK 3 ####

# HW 3: You have a F2 population of a plant species (which is a diploid) resulted
#       from a cross between a plant with red flowers and a plant with white flowers.
#       You hypothesized that: the flower color is controlled by a single gene; each
#       of the parents is homozygous for the gene; the red allele is dominant over the
#       white allele. If your hypothesis is correct: when you scored 20 F2 plants, what
#       is the probability to have the number of plants with white flowers is 2 or fewer?

# The first part is biology. You need to know that with this hypothesis the expected average
# ratio between red-flower and white-flower F2 plants are 3:1 (i.e., 0.75:0.25).

# Then we are dealing with a binomial distribution with n=20, p=0.25.The probability to get
# no white-flower plant is:

p0 = choose(20,0)*(0.25^0)*(0.75^20)

# The probability to get one white-flower plant is:
p1 = choose(20,1)*(0.25^1)*(0.75^19)

# The probability to get two white-flower plants is:
p2 = choose(20,2)*(0.25^2)*(0.75^18)

# The probability we need is p0+p1+p2:
p0+p1+p2

sum(choose(20, 0:2) * .25^(0:2) * .75^(20-0:2) )

# We know how to do it in a much simpler way, right?
pbinom(2,20,prob=0.25)
?pbinom
# Using simulation, this assignment can be done like this.
white.count=c()
for (i in 1:10000) {
  wh.numb = sample(c("R","W"),20,replace=T,prob=c(0.75,0.25))
  white.count = c(white.count, sum(wh.numb=="W"))
  }
sum(white.count<=2)/10000[1]


# Computer problems:
# 1. i. Plot the probability density curve for the normal distribution
#       with mean = 7 and sd = 2.5 in the range of [0,14]
x = seq(0,14,0.1)
plot(x,dnorm(x,mean=7,sd=2.5))

#    ii. What is the probability to get the random variable values smaller
#        than 0 or larger than 4 with this normal distribution?
pnorm(0,mean=7,sd=2.5) + (1 - pnorm(4,mean=7,sd=2.5))

#    iii. Get an approximate answer to (ii) using a simulation based on
#         10,000 randomly sampled numbers.
rand.norm = rnorm(10000,mean=7,sd=2.5)
sum(rand.norm < 0 | rand.norm > 4) / length(rand.norm)

# 2. i. Calculate the probability of having 12 or more heads in 15 coin
#       tosses when the probability to have a head is 0.45.
1 - pbinom(11,size=15,prob=0.45)

#    ii. Get an approximate answer to (i) using 10,000 random numbers.
rand.binom=rbinom(10000,size=15,prob=0.45)
sum(rand.binom>=12)/length(rand.binom)

# 3.	This is how generating random numbers for a particular distribution works.
#     i.	 Make the vector object rand.numb that have 5000 random numbers sampled
#          from a uniform distribution from 0 to 1. This is to randomly sample the
#          probability values (the cumulative distribution values). Thus, they are
#          from 0 to 1.
rand.numb  = runif(5000,min=0,max=1)

#     ii.	 For a normal distribution: apply the quantile function qnorm (the inverse
#          of the cumulative distribution function) to rand.numb. Let’s use mean = 10,
#          sd = 3. Name the resulting vector rand.norm.
rand.norm = qnorm(rand.numb,mean=10,sd=3)

#     iii. Plot a density histogram of rand.norm. 
hist(rand.norm,freq=F)

#     iv.  Calculate the mean and sd of rand.norm. Are they as expected?
mean(rand.norm)
sd(rand.norm)

#     v.   For a binomial distribution: apply the quantile function rbinom to rand.numb.
#          Let’s use size = 12, prob = 0.3. Name the resulting vector rand.binom.
rand.binom=rbinom(rand.norm,size=12,prob=0.3)

#     vi.	 Make a density histogram of rand.binom. You should use the method with the 
#          table function and the plot function. Make sure to make the values the probabilities. 
rb.tab=table(rand.binom)
rb.tab=rb.tab/sum(rb.tab)
rb.tab

plot(as.integer(names(rb.tab)),as.numeric(rb.tab),
     type='h',
     xlab='x',ylab='Point probability')

#     vii. Try to include the point probability of the binomial population from which 
#          these random numbers were sampled from. Use the points function with type='h'
#          and col='red' arguments, and for the x-axis value give an offset of 0.2 to avoid
#          the lines from (vi) overlapping. Are the rand.binom density and the population
#          density well corresponding?
x=0:20
points(x+0.2,dbinom(x,size=12,prob=0.3),
       type='h',col='red')

#    viii. Calculate the mean and sd of rand.binom. Are these close to the expected values? What are the expected values if they are truly random numbers sampled from the binomial population? Look them up in Wikipedia (or something else).
mean(rand.binom); sd(rand.binom)

# 4.   Let’s simulate the central limit theorem with a Poisson distribution. Yes, the central
#      limit theorem works with a discrete distribution as well.
# i.   Make the vector object poisson.sampling that have 1e5 random numbers sampled from
#      a Poisson distribution with the mean of 3.
poisson.sampling = rpois(1e5,lambda=3)

# ii. Make the vector poisson.sampling into a matrix of 20000 rows, and call the matrix
#     poisson.sampling.mat. These are random numbers, so it doesn’t matter if you fill the
#     matrix by row or by column. Now we can consider that this is a matrix of 20000 samples
#     of size 5. Check the dimension of the matrix.
poisson.sampling.mat=matrix(poisson.sampling,nrow=20000)
dim(poisson.sampling.mat)

# iii	We haven’t learned the apply function yet, but you can use it in the following way.
mean.vals = apply(poisson.sampling.mat, 1, mean)

#     The apply function takes the values from each row of the matrix poisson.sampling.mat
#     and calculates and returns their mean value. It does this from row 1 to row 20000 one
#     by one. So, the mean.vals is 20000 sample means when the size of each sample was 5. 
#     In the apply function, the second argument “1” means that the apply function applies the
#     specified function to each row (it would do each column if you put “2” for the second 
#     argument); the third argument is the function to apply – in this case, the mean function.

# iv. Plot a density histogram of mean.vals. Make the column width of the histogram 0.2. 
#     Make sure that you plot the density.
range(mean.vals)
hist(mean.vals,freq=F,breaks=seq(0,7,0.2))

# v.	Superimpose the density of a normal distribution with mean = 3 and sd = sqrt(3/5)
#     in red. To superimpose a curve on an existing plot, use the lines function. You 
#     should look up the lines function how to give x- and y-axis values in arguments. 
#     For the x-axis values, use seq(1, 5, 0.01). Add an argument “col = 'red'” to make the curve red.
x1=seq(0,7,0.01)
lines(x1,dnorm(x1,mean=3,sd=sqrt(3/5)),col='red')

# vi.	Superimpose the Poisson population in blue in the plot as well. Note that the Poisson 
#     distribution is a discrete distribution. You should use the points function with type='h' ,
#     col='blue' arguments. You may want to make the line thicker. What about including lwd=3 argument?
x2=0:7  
points(x2,dpois(x2,lambda=3),
       type='h',col='blue',lwd=3)

# Part 3(2)
#1.   Let’s look at the estimated population variance when the population is a normal population with
#     mean = 10, sd = 3.

#i.   Make a matrix of 1e4 samples of size 3.
samp.size=3
true.pop.var=3^2
samp.mat = matrix(rnorm(1e4*samp.size,mean=10,sd=sqrt(true.pop.var)),ncol=samp.size)

#ii.  Calculate the estimated population variance from each sample, and then calculate
#     the mean of the estimated population variances across the samples.
pop.var=apply(samp.mat,1,var)
mean(pop.var)

#iii. What is the true pop variance? How does it compare to the estimated value?
# 9 - so very close.

# 2. Take a sample of size 20 from a uniform distribution in the range [0,20].
#    Calculate mean, sd, variance, median.
#    Calculate empirical quartile values.
samp=runif(20,10,20)
mean(samp);sd(samp);sd(samp);var(samp);median(samp)
quantile(samp)

# 3. Make a vector object consisting of non-NA values in the igf1 column of the juul
#    data frame in the ISwR package.
library(ISwR)
igf1a=na.omit(juul$igf1)
#    Draw the empirical cumulative distribution for the vector.
n=length(igf1a)
plot(sort(igf1a),(1:n)/n,type='s')
abline(h=0.8)
#    What is an approximate ifg1 value for the 80th percentile. Obtain this value from
#    the plot made in (ii). --> 500
#    Obtain the value for (iii) from the ifg1 vector from (i).
quantile(igf1a, probs = 0.8)
sort(igf1a)[round(0.8*n)]

#4. A practice for plotting (1).
#   Plot the point probability histogram of a binomial distribution B(8, 0.75) in blue.
#   Make the range of the x-axis from 0 to 12.
#   Name the x-axis and y-axis labels “X” and “Density”, respectively.
#   Name the title of the plot “Comparison of normal and binomial distributions”.
#   Superimpose the density of a normal distribution N(6, 1.22^2) in red dashed lines with lwd=2.
#   Add an appropriate legend.

x3 = 0:8
plot(x3, dbinom(x3, size=8, prob=0.75), type='h',
     col='blue', xlim=c(0, 12),
     xlab='X', ylab='Density',
     main='Normal and binomial distributions')

# 4(ii)
x2 = seq(0,12,0.02)
points(x2, dnorm(x2,mean=6,sd=1.22),
       type='l', col='red', lwd=2, lty=2)

# 4(iii)
legend('topleft', c('Binomial Distr, n=8, p=0.75','Normal Distr, mu=6, sd=1.22'), 
       col=c('blue','red'), lty=c(1,2), inset=0.02, cex=0.7)

# 5. A practice for plotting (2).
#    Take a sample of size 15 from N(10,2^2), a sample of size 18 from N(11, 4^2),
#    and a sample of size 16 from N(9, 3^2). Make a boxplot in which the distributions
#    of the three samples can be compared. Name the distributions ‘sample 1’, ‘sample 2’,
#    and ‘sample 3’. Hint: one way is to make a list with the samples as components, since
#    the sample sizes are different; another way is to make it in a long-format table and use a formula.

#   Superimpose the actual sample values on the boxplots. Hint: use the stripchart function
#   with the add argument = T, the method argument=‘jitter’, the vertical argument = T, the 
#   pch argument =16, and the col argument =‘red’.

samp1 = rnorm(15, mean=10, sd=2)
samp2 = rnorm(18, mean=11, sd=4)
samp3 = rnorm(16, mean=9, sd=3)

## method (1)
samp.list = list('sample 1'=samp1, 'sample 2'=samp2,
                 'sample 3'=samp3)
boxplot(samp.list)
stripchart(samp.list, add=T, vertical=T,
           method='jitter', pch=16, col='red')

## method (2)
samp.df = data.frame(vals=c(samp1,samp2,samp3),
                     sample=rep(c('sample 1', 'sample 2', 'sample 3'),
                                c(15,18,16)))
boxplot(vals ~ sample, data=samp.df)
stripchart(vals ~ sample, data=samp.df, add=T, vertical=T,
           method='jitter', pch=16, col='red')

# WEEK 2 ####

# You had 50 seeds of blue morning glory and 20 seeds of white morning glory
# mixed well. You kept 10 seeds for yourself and gave away the other seeds.
# Assuming that all 10 seeds grow to make flowers, list up all the possible
# blue:white plant number combinations, and calculate the probability for
# each combination.

seed.count <- paste('Blue',0:10,':White',10:0,sep='')
my.seed <- sample(seed.count,10,replace = T)

choose(70,10)
prob.seed <- choose(50,0:10) * choose(20,10 - 0:10) / choose(70,10)
names(prob.seed) = seed.count

# Computer lab

# 1. There are 40 runners numbered 11 to 50. Randomly select 8 runners.
sample(11:50,8)

# 2. Randomly select two representatives out of a group of people consisting
# of Barbara, Charles, Debbie, Eric, Bob, Nancy, Megan, and Nathan.
sample(c('Barbara','Charles','Debbie','Eric','Bob','Nancy','Megan','Nathan'),2)

# 3. The pond has 300 yellow fish and 100 blue fish. You threw a fishing net
# and caught 50 fish. Assume that both fish kinds are evenly distributed in the pond.

#    i. Simulate your catch in terms of yellow and blue fish.
fish <- c(rep('Y',300),rep('B',100))
caught <- sample(fish,50)

#    ii. Count the number of yellow fish you got. Make R count them.
sum(fish == 'Y')

# 4. Simulate one event consisting of 20 coin tosses. The probabilities for a head and a tail are the same.
sample(c('H','T'),20,replace=T)

# 5. The pond has 300 yellow fish and 100 blue fish. You caught by a fishing pole and released 50 fish one
#    by one and recorded the color of the fish. Assume that all the fish have the same probabilities to be
#    caught, regardless of the color or whether it was caught before.

#    i. Simulate your catch in terms of yellow and blue fish.
catch <- sample(fish,50,replace=T)

# option 2:
caught.fish2 = sample(c('Y','B'), 50, replace=T, prob = c(300, 100))

#    ii. Count the number of yellow fish you got. Make R count them.

sum(catch == 'Y')

# 6. Simulate one event consisting of 20 coin tosses. The probabilities for a head is 0.55.

toss <- sample(c('H','T'),20,replace=T,prob=c(55,45))

# 7. The pond has 300 yellow fish and 100 blue fish. You caught by a fishing pole and 
#    released 50 fish one by one and recorded the color of the fish. Assume that the 
#    probabilities for a yellow fish to be caught is twice as high as that for a blue
#    fish if the same numbers of yellow and blue fishes are in the pond.

#    i. Simulate your catch in terms of yellow and blue fish.

caught.fish3 = sample(c('Y','B'), 50, replace=T, prob = c(3/4 * 2, 1/4 * 2))

#    ii. Count the number of yellow fish you got. Make R count them.

sum(caught.fish3 == 'Y')

# 9.
# Generate the vector x.val that contains numeric values from 0.01 to 0.99 every 0.001.
x.val = seq(0.01,0.99,0.001)

# Use the “length” function to check the number of the elements in x.val.
length(x.val)

# Randomize the order of the elements in the vector x.val and assign the randomized values back to x.val.
x.val = sample(x.val)

# Calculate the vector x1.val, which is the result of applying the function: 
# qnorm(“each element value of x.val”). 
x1.val = qnorm(x.val)

# Visualize the density of x1.val using the histogram.
hist(x1.val,freq=F)

# Make the breaks of the histogram every 0.1. You can use the “breaks” argument to specify the number of bars or the boundary values of each bar.
hist(x1.val,breaks=seq(min(x1.val), max(x1.val)+0.1, 0.1),freq=F)

# Calculate the mean and the standard deviation of x1.val.
mean(x1.val)
sd(x1.val)

# Run: x2.val = seq(-3,3,0.01); lines(x2.val, dnorm(x2.val),col=‘red’) 
# Basically this shows that the histogram is directly comparable with the standard normal distribution’s density. I will explain it more.

x2.val = seq(-3,3,0.01)
lines(x2.val, dnorm(x2.val),col='red')

# 10. Draw a histogram of rbinom(40, 6, 0.4) with the breaks at every integer number.
hist(rbinom(40,6,0.4),breaks=seq(0,6,1))
# OR:
my.binom <- rbinom(40,6,0.4)
hist(my.binom,breaks=seq(0,6,1))

# 11. How many different ways to choose 50 fishes when all of the 400 fishes are numbered?
choose(400,50)

#     i. In the above fishing net case (300 yellow, 100 blue; catching 50 by a single throw of the net),
#         what is the probability of catching 32 yellow fish? Calculate it exactly – not simulating it.
choose(300,32) * choose(100,50-32) / choose(400,50)

#     ii. Calculate the probability of catching 32 yellow fish or fewer.
#         what is the probability of catching 32 yellow fish? Calculate it exactly – not simulating it.
sum(choose(300,32:0) * choose(100,50-32:0) / choose(400,50))

# Part 2, 1: Can you explain what is happening with: sum(head.count==4)?
head.count = c()
sum(head.count==4)
# How many number 4's are there in the head.count vector.

# Part 2, 2: Based on previous studies, this chemical kills 70% of adult mice at the dose
#            of 1.3 ng/g of mouse weight. If you apply the chemical at this dose to 5 adult
#            mice, what are the probabilities you see 0, 1, 2, 3, 4, and 5 mice killed, respectively? 

#            (i) Calculate the values exactly.
choose(5,0:5) * 0.7^(0:5) * 0.3^(5:0)

#            (ii) Calculate the values using 10,000-event simulation for each number of killed mice.
dead.m = c()
for(i in 1:1e4) {
  mice.fate = sample(c('live','dead'),5,replace=T,prob=c(.3,.7))
  dead.m = c(dead.m,sum(mice.fate == 'dead'))
}
table(dead.m)/length(dead.m)

# Part 2, 3: 
lambda = 5
x = 0
if (x==0) {
  div.x = 1
} else {
  div.x = x
}
exp(-lambda) * lambda ^x/prod(1:div.x)

# Part 3, 4:
lambda = 5
x = 0
if (x==0) {
  div.x = 1
} else {
  div.x = x
}
exp(-lambda) * lambda ^x/prod(1:div.x)
dpois(0,5)

head.count = c()
for (event.numb in 1:10000) {
  event = sample(c("H","T"),10,replace=T,prob=c(0.6,0.4)) 
  head.count = c( head.count, sum(event == "H") )}
sum(head.count==4)


# WEEK 1 ####

# Homework:
install.packages("ISwR")
library(ISwR)
summary(tlc)
str(tlc)
?tlc
# female = 1, male = 2

tlc[tlc$sex == 1 & tlc$height > 148.5 & tlc$height < 162.5,"tlc"]

