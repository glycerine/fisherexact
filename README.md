fisherexact: Fisher's Exact Test in Go
===========

Fisher's Exact test (FET) is one of the most useful
statistical tests. 

The FET here evaluates a 2x2 contingency table for independence,
returning a p-value.

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
func FisherExactTest22(n11, n12, n21, n22 int) (pvalue float64)

// ChiSquaredTest22 input is the same.
func ChiSquaredTest22(n11, n12, n21, n22 int, yates bool) (pval float64)
~~~

The FET can be used for large and small data. 
For numerical efficiency, the FET is typically 
deployed when small data makes the Chi-squared test's 
asymptotic assumptions unreliable, but it
can be used on any size of data or contingency
table. As the wikipedia article says,

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
