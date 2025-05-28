fisherexact: Fisher's Exact Test in Go
===========

Fisher's Exact test (FET) is one of the most useful
statistical tests. 

The FET here evaluates a 2x2 contingency table for independence,
returning p-values; typically you'll use the last, two-sided
p-value, by default.

The FET, like the Chi-squared, can be generalized 
to larger tables, but this
repo only does 2x2 tables at the moment.

~~~
// FisherExactTest22 computes Fisher's Exact test
// for independence on a 2x2 contingency table.
//
// n11  n12  | n1_
// n21  n22  | n2_
// ----------+-----
// n_1  n_2  | n
//
// is the layout assumed.
// Three p-values are returned, for each of three
// alternative hypotheses.
func FisherExactTest22(n11, n12, n21, n22 int) (
  less, 
  greater, 
  twoSidedPvalue float64, // you want this, generally.
  )

// ChiSquaredTest22 input is the same.
func ChiSquaredTest22(n11, n12, n21, n22 int, yates bool) (
  pval float64,
  )
~~~

The FET can be used for large and small data. 
For numerical efficiency, the FET is typically 
deployed when small data makes the Chi-squared test's 
asymptotic assumptions unreliable. The nice
thing about the FET is that it works on small data too.
It is appropriate for any size of data. As the 
wikipedia article says,

> [The FET] becomes difficult to calculate with 
> large samples or well-balanced tables, but 
> fortunately these are exactly the conditions 
> where the chi-squared test is appropriate.

https://en.wikipedia.org/wiki/Fisher%27s_exact_test

For comparison, we provide a Chi-squared test implementation
based on gonum/cephes calculations. The small
cephes subpackage (Netlib code by Stephen Mosher) 
required was copied in to avoid depending
on the full gonum library.

https://en.wikipedia.org/wiki/Chi-squared_test

----
author: Jason E. Aten, Ph.D.

License: MIT

FET source material in C++ from https://github.com/samtools/htslib/ (MIT license)

Chi-squared distribution computation copied from gonum's
cephes package. See cephes/ for details/license(3-clause BSD).
https://github.com/gonum/gonum/

https://github.com/gonum/gonum/tree/720fcb9699a9e01862309471af0aac7eb56240bc/mathext/internal/cephes

------
# alternative hypotheses

The R docs for fisher.test explain the hypotheses tested /
corresponding to the returned p-values (less, greater,
two-sided):

> For 2 by 2 tables, the null of conditional independence is
> equivalent to the hypothesis that the odds ratio equals one.
> ‘Exact’ inference can be based on observing that in general, given
> all marginal totals fixed, the first element of the contingency
> table has a non-central hypergeometric distribution with
> non-centrality parameter given by the odds ratio (Fisher, 1935).
> The alternative for a one-sided test is based on the odds ratio,
> so ‘alternative = "greater"’ is a test of the odds ratio being
> bigger than ‘or’ [jea: the odds-ratio, which in this Go package is 
> assumed = 1.0, like the default R value].
>
> Two-sided tests are based on the probabilities of the tables, and
> take as ‘more extreme’ all tables with probabilities less than or
> equal to that of the observed table, the p-value being the sum of
> such probabilities.

The full fisher.test R docs may be helpful here. See 
the literature references at the end for primary sources.

Note that no R code is used in this fisherexact Go package,
and thus the GPL does not apply to its code.

~~~
> ?fisher.test
fisher.test               package:stats                R Documentation

Fisher's Exact Test for Count Data

Description:

     Performs Fisher's exact test for testing the null of independence
     of rows and columns in a contingency table with fixed marginals.

Usage:

     fisher.test(x, y = NULL, workspace = 200000, hybrid = FALSE,
                 hybridPars = c(expect = 5, percent = 80, Emin = 1),
                 control = list(), or = 1, alternative = "two.sided",
                 conf.int = TRUE, conf.level = 0.95,
                 simulate.p.value = FALSE, B = 2000)
     
Arguments:

       x: either a two-dimensional contingency table in matrix form, or
          a factor object.

       y: a factor object; ignored if ‘x’ is a matrix.

workspace: an integer specifying the size of the workspace used in the
          network algorithm.  In units of 4 bytes.  Only used for
          non-simulated p-values larger than 2 by 2 tables.  Since R
          version 3.5.0, this also increases the internal stack size
          which allows larger problems to be solved, however sometimes
          needing hours.  In such cases, ‘simulate.p.values=TRUE’ may
          be more reasonable.

  hybrid: a logical. Only used for larger than 2 by 2 tables, in which
          cases it indicates whether the exact probabilities (default)
          or a hybrid approximation thereof should be computed.

hybridPars: a numeric vector of length 3, by default describing
          “Cochran's conditions” for the validity of the chisquare
          approximation, see ‘Details’.

 control: a list with named components for low level algorithm control.
          At present the only one used is ‘"mult"’, a positive integer
          >= 2 with default 30 used only for larger than 2 by 2 tables.
          This says how many times as much space should be allocated to
          paths as to keys: see file ‘fexact.c’ in the sources of this
          package.

      or: the hypothesized odds ratio.  Only used in the 2 by 2 case.

alternative: indicates the alternative hypothesis and must be one of
          ‘"two.sided"’, ‘"greater"’ or ‘"less"’.  You can specify just
          the initial letter.  Only used in the 2 by 2 case.

conf.int: logical indicating if a confidence interval for the odds
          ratio in a 2 by 2 table should be computed (and returned).

conf.level: confidence level for the returned confidence interval.
          Only used in the 2 by 2 case and if ‘conf.int = TRUE’.

simulate.p.value: a logical indicating whether to compute p-values by
          Monte Carlo simulation, in larger than 2 by 2 tables.

       B: an integer specifying the number of replicates used in the
          Monte Carlo test.

Details:

     If ‘x’ is a matrix, it is taken as a two-dimensional contingency
     table, and hence its entries should be nonnegative integers.
     Otherwise, both ‘x’ and ‘y’ must be vectors or factors of the same
     length.  Incomplete cases are removed, vectors are coerced into
     factor objects, and the contingency table is computed from these.

     For 2 by 2 cases, p-values are obtained directly using the
     (central or non-central) hypergeometric distribution. Otherwise,
     computations are based on a C version of the FORTRAN subroutine
     FEXACT which implements the network developed by Mehta and Patel
     (1983, 1986) and improved by Clarkson, Fan and Joe (1993).  The
     FORTRAN code can be obtained from <https://netlib.org/toms/643>.
     Note this fails (with an error message) when the entries of the
     table are too large.  (It transposes the table if necessary so it
     has no more rows than columns.  One constraint is that the product
     of the row marginals be less than 2^31 - 1.)

     For 2 by 2 tables, the null of conditional independence is
     equivalent to the hypothesis that the odds ratio equals one.
     ‘Exact’ inference can be based on observing that in general, given
     all marginal totals fixed, the first element of the contingency
     table has a non-central hypergeometric distribution with
     non-centrality parameter given by the odds ratio (Fisher, 1935).
     The alternative for a one-sided test is based on the odds ratio,
     so ‘alternative = "greater"’ is a test of the odds ratio being
     bigger than ‘or’.

     Two-sided tests are based on the probabilities of the tables, and
     take as ‘more extreme’ all tables with probabilities less than or
     equal to that of the observed table, the p-value being the sum of
     such probabilities.

     For larger than 2 by 2 tables and ‘hybrid = TRUE’, asymptotic
     chi-squared probabilities are only used if the ‘Cochran
     conditions’ (or modified version thereof) specified by ‘hybridPars
     = c(expect = 5, percent = 80, Emin = 1)’ are satisfied, that is if
     no cell has expected counts less than ‘1’ (‘= Emin’) and more than
     80% (‘= percent’) of the cells have expected counts at least 5 (‘=
     expect’), otherwise the exact calculation is used.  A
     corresponding ‘if()’ decision is made for all sub-tables
     considered.

     Accidentally, R has used ‘180’ instead of ‘80’ as ‘percent’, i.e.,
     ‘hybridPars[2]’ in R versions between 3.0.0 and 3.4.1 (inclusive),
     i.e., the 2nd of the ‘hybridPars’ (all of which used to be
     hard-coded previous to R 3.5.0).  Consequently, in these versions
     of R, ‘hybrid=TRUE’ never made a difference.

     In the r x c case with r > 2 or c > 2, internal tables can get too
     large for the exact test in which case an error is signalled.
     Apart from increasing ‘workspace’ sufficiently, which then may
     lead to very long running times, using ‘simulate.p.value = TRUE’
     may then often be sufficient and hence advisable.

     Simulation is done conditional on the row and column marginals,
     and works only if the marginals are strictly positive.  (A C
     translation of the algorithm of Patefield (1981) is used.)  Note
     that the default number of replicates (‘B = 2000’) implies a
     minimum p-value of about 0.0005 (1/(B+1)).

Value:

     A list with class ‘"htest"’ containing the following components:

 p.value: the p-value of the test.

conf.int: a confidence interval for the odds ratio.  Only present in
          the 2 by 2 case and if argument ‘conf.int = TRUE’.

estimate: an estimate of the odds ratio.  Note that the _conditional_
          Maximum Likelihood Estimate (MLE) rather than the
          unconditional MLE (the sample odds ratio) is used.  Only
          present in the 2 by 2 case.

null.value: the odds ratio under the null, ‘or’.  Only present in the 2
          by 2 case.

alternative: a character string describing the alternative hypothesis.

  method: the character string ‘"Fisher's Exact Test for Count Data"’.

data.name: a character string giving the name(s) of the data.

References:

     Agresti, A. (1990).  _Categorical data analysis_.  New York:
     Wiley.  Pages 59-66.

     Agresti, A. (2002).  _Categorical data analysis_. Second edition.
     New York: Wiley.  Pages 91-101.

     Fisher, R. A. (1935).  The logic of inductive inference.  _Journal
     of the Royal Statistical Society Series A_, *98*, 39-54.
     doi:10.2307/2342435 <https://doi.org/10.2307/2342435>.

     Fisher, R. A. (1962).  Confidence limits for a cross-product
     ratio.  _Australian Journal of Statistics_, *4*, 41.
     doi:10.1111/j.1467-842X.1962.tb00285.x
     <https://doi.org/10.1111/j.1467-842X.1962.tb00285.x>.

     Fisher, R. A. (1970).  _Statistical Methods for Research Workers_.
     Oliver & Boyd.

     Mehta, Cyrus R. and Patel, Nitin R. (1983).  A network algorithm
     for performing Fisher's exact test in r x c contingency tables.
     _Journal of the American Statistical Association_, *78*, 427-434.
     doi:10.1080/01621459.1983.10477989
     <https://doi.org/10.1080/01621459.1983.10477989>.

     Mehta, C. R. and Patel, N. R. (1986).  Algorithm 643: FEXACT, a
     FORTRAN subroutine for Fisher's exact test on unordered r x c
     contingency tables.  _ACM Transactions on Mathematical Software_,
     *12*, 154-161.  doi:10.1145/6497.214326
     <https://doi.org/10.1145/6497.214326>.

     Clarkson, D. B., Fan, Y. and Joe, H. (1993) A Remark on Algorithm
     643: FEXACT: An Algorithm for Performing Fisher's Exact Test in r
     x c Contingency Tables.  _ACM Transactions on Mathematical
     Software_, *19*, 484-488.  doi:10.1145/168173.168412
     <https://doi.org/10.1145/168173.168412>.

     Patefield, W. M. (1981).  Algorithm AS 159: An efficient method of
     generating r x c tables with given row and column totals.
     _Applied Statistics_, *30*, 91-97.  doi:10.2307/2346669
     <https://doi.org/10.2307/2346669>.

See Also:

     ‘chisq.test’

     ‘fisher.exact’ in package ‘exact2x2’ for alternative
     interpretations of two-sided tests and confidence intervals for 2
     by 2 tables.

Examples:

     ## Agresti (1990, p. 61f; 2002, p. 91) Fisher's Tea Drinker
     ## A British woman claimed to be able to distinguish whether milk or
     ##  tea was added to the cup first.  To test, she was given 8 cups of
     ##  tea, in four of which milk was added first.  The null hypothesis
     ##  is that there is no association between the true order of pouring
     ##  and the woman's guess, the alternative that there is a positive
     ##  association (that the odds ratio is greater than 1).
     TeaTasting <-
     matrix(c(3, 1, 1, 3),
            nrow = 2,
            dimnames = list(Guess = c("Milk", "Tea"),
                            Truth = c("Milk", "Tea")))
     fisher.test(TeaTasting, alternative = "greater")
     ## => p = 0.2429, association could not be established
     
     ## Fisher (1962, 1970), Criminal convictions of like-sex twins
     Convictions <- matrix(c(2, 10, 15, 3), nrow = 2,
                           dimnames =
                    list(c("Dizygotic", "Monozygotic"),
                         c("Convicted", "Not convicted")))
     Convictions
     fisher.test(Convictions, alternative = "less")
     fisher.test(Convictions, conf.int = FALSE)
     fisher.test(Convictions, conf.level = 0.95)$conf.int
     fisher.test(Convictions, conf.level = 0.99)$conf.int
     
     ## A r x c table  Agresti (2002, p. 57) Job Satisfaction
     Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
                dimnames = list(income = c("< 15k", "15-25k", "25-40k", "> 40k"),
                          satisfaction = c("VeryD", "LittleD", "ModerateS", "VeryS")))
     fisher.test(Job) # 0.7827
     fisher.test(Job, simulate.p.value = TRUE, B = 1e5) # also close to 0.78
     
     ## 6th example in Mehta & Patel's JASA paper
     MP6 <- rbind(
             c(1,2,2,1,1,0,1),
             c(2,0,0,2,3,0,0),
             c(0,1,1,1,2,7,3),
             c(1,1,2,0,0,0,1),
             c(0,1,1,1,1,0,0))
     fisher.test(MP6)
     # Exactly the same p-value, as Cochran's conditions are never met:
     fisher.test(MP6, hybrid=TRUE)
     
~~~
source: https://github.com/wch/r-source/blob/8329caa0d89a7e036663e1247d4f4ee7a55e756a/src/library/stats/man/fisher.test.Rd#L2

R is licensed under GPL 2 or later.

See https://github.com/wch/r-source/blob/8329caa0d89a7e036663e1247d4f4ee7a55e756a/COPYING

