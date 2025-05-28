fisherexact: Fisher's Exact Test in Go
===========

Fisher's Exact test (FET) is one of the most useful
statistical tests[1]. A fun fact about Sir Ronald A. Fisher
(https://en.wikipedia.org/wiki/Ronald_Fisher) is that
based on citation counts, he is considered to be
the single most influential scientist of all time.
He founded the modern fields of both statistics 
and population genetics.

The FET here evaluates a 2x2 contingency table for independence,
returning p-values; typically you'll use the two-sided
p-value by default. Technically this is a test of
conditional independence, as the marginals of the
table are fixed and conditioned on.

The FET, like the Chi-squared, can be generalized 
to larger tables, but this
repo only does 2x2 tables at the moment.

~~~
// fisherexact.TwoSided22() computes Fisher's Exact test
// for independence in a 2x2 contingency table, and
// returns the p-value for the two-sided null hypothesis
// that the odds-ratio is 1.
//
// n11  n12  | n1_
// n21  n22  | n2_
// ----------+-----
// n_1  n_2  | n
//
// is the layout assumed.
func TwoSided22(n11, n12, n21, n22 int) (twoSidedPvalue float64)

// ChiSquaredTest22 assumes the same layout.
func ChiSquaredTest22(n11, n12, n21, n22 int, yates bool) (
  pval float64,
  )
  
// Also available for testing 1-sided hypotheses.
func FisherExactTest22(n11, n12, n21, n22 int) (
   lessPvalue, 
   greaterPvalue, 
   twoSidedPvalue, 
   probCurrentTable float64,
   )
~~~

The odds-ratio examined here, assuming rows are groups
and columns are outcomes, is defined as follows.

~~~
odds1 = odds of outcome 1 in group 1 = n11/n12

odds2 = odds of outcome 1 in group 2 = n21/n22

OR = odds ratio = odds1/odds2 = (n11 * n22)/(n12 * n21)
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
cephes subpackage (Netlib code by Stephen L. Moshier) 
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

[1] Fisher, R. A. (1935).  The logic of inductive inference.  _Journal
 of the Royal Statistical Society Series A_, *98*, 39-54.
http://csyue.nccu.edu.tw/ch/The%20Logic%20of%20Inductive%20Inference.pdf

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
> assumed = 1.0 under the null-hypothesis, like the default R value].
>
> Two-sided tests are based on the probabilities of the tables, and
> take as ‘more extreme’ all tables with probabilities less than or
> equal to that of the observed table, the p-value being the sum of
> such probabilities.

The full fisher.test R docs may be helpful here; type `?fisher.test`
in R to view them. See the literature references at the end for more.

fisher.test doc source: https://github.com/wch/r-source/blob/8329caa0d89a7e036663e1247d4f4ee7a55e756a/src/library/stats/man/fisher.test.Rd#L2

