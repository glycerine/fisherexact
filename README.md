fisherexact: the Fisher-Exact Test in Go
===========

Fisher's Exact test (FET) is one of the most useful
statistical tests. 

FET evaluates a contingency table for independence.

FET can be used for large and small data. 
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

----
author: Jason E. Aten, Ph.D.

License: MIT
