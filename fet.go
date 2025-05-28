package fisherexact

import (
	"fmt"
	"math"
)

// reference
// https://en.wikipedia.org/wiki/Fisher%27s_exact_test

// based on htslib (MIT licensed) from a decade ago.
// today's htslib can be found https://github.com/samtools/htslib

// Log gamma function
// \log{\Gamma(z)}
// AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
//
// Question: use math.Lgamma / lgamma below instead?
func kf_lgamma(z float64) float64 {

	x := 0.0
	x += 0.1659470187408462e-06 / (z + 7)
	x += 0.9934937113930748e-05 / (z + 6)
	x -= 0.1385710331296526 / (z + 5)
	x += 12.50734324009056 / (z + 4)
	x -= 176.6150291498386 / (z + 3)
	x += 771.3234287757674 / (z + 2)
	x -= 1259.139216722289 / (z + 1)
	x += 676.5203681218835 / z
	x += 0.9999999999995183
	return math.Log(x) - 5.58106146679532777 - z + (z-0.5)*math.Log(z+6.5)
}

// The following computes regularized incomplete gamma functions.
// Formulas are taken from Wiki, with additional input from Numerical
// Recipes in C (for modified Lentz's algorithm) and AS245
// (http://lib.stat.cmu.edu/apstat/245).
//
// A good online calculator is available at:
//
//   http://www.danielsoper.com/statcalc/calc23.aspx
//
// It calculates upper incomplete gamma function, which equals
// kf_gammaq(s,z)*tgamma(s).
///

const KF_GAMMA_EPS = 1e-14
const KF_TINY = 1e-290

// regularized lower incomplete gamma function, by series expansion
func _kf_gammap(s, z float64) float64 {

	sum := 1.0
	x := 1.0
	for k := 1; k < 100; k++ {
		x *= z
		sum += (x / (s + float64(k)))
		if x/sum < KF_GAMMA_EPS {
			break
		}
	}
	return math.Exp(s*math.Log(z) - z - kf_lgamma(s+1) + math.Log(sum))
}

// regularized upper incomplete gamma function, by continued fraction
func _kf_gammaq(s, z float64) float64 {

	var C, D, f float64
	f = 1.0 + z - s
	C = f
	// Modified Lentz's algorithm for computing continued fraction
	// See Numerical Recipes in C, 2nd edition, section 5.2
	for j := 1; j < 100; j++ {
		//a := j * (s - j), b = (j<<1) + 1 + z - s, d;
		a := float64(j) * (s - float64(j))
		b := float64((j<<1)+1) + z - s

		D = b + a*D
		if D < KF_TINY {
			D = KF_TINY
		}
		C = b + a/C
		if C < KF_TINY {
			C = KF_TINY
		}
		D = 1 / D
		d := C * D
		f *= d
		if math.Abs(d-1) < KF_GAMMA_EPS {
			break
		}
	}
	return math.Exp(s*math.Log(z) - z - kf_lgamma(s) - math.Log(f))
}

func kf_gammap(s, z float64) float64 {
	if z <= 1 || z < s {
		return _kf_gammap(s, z)
	}
	return 1 - _kf_gammaq(s, z)
}

func kf_gammaq(s, z float64) float64 {
	if z <= 1 || z < s {
		return 1 - _kf_gammap(s, z)
	}
	return _kf_gammaq(s, z)
}

// Regularized incomplete beta function. The method is taken from
// Numerical Recipe in C, 2nd edition, section 6.4. The following web
// page calculates the incomplete beta function, which equals
// kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
//
//	http://www.danielsoper.com/statcalc/calc36.aspx
func kf_betai_aux(a, b, x float64) float64 {
	if x == 0 {
		return 0
	}
	if x == 1 {
		return 1
	}
	f := 1.0
	C := f
	D := 0.0
	// Modified Lentz's algorithm for computing continued fraction
	for j := 1; j < 200; j++ {
		var aa, d float64
		m := float64(j >> 1)
		if (j & 1) != 0 {
			aa = -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))
		} else {
			aa = m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m))
		}
		D = 1 + aa*D
		if D < KF_TINY {
			D = KF_TINY
		}
		C = 1 + aa/C
		if C < KF_TINY {
			C = KF_TINY
		}
		D = 1 / D
		d = C * D
		f *= d
		if math.Abs(d-1) < KF_GAMMA_EPS {
			break
		}
	}
	return math.Exp(kf_lgamma(a+b)-kf_lgamma(a)-
		kf_lgamma(b)+a*math.Log(x)+b*math.Log(1-x)) / a / f
}

// Regularized incomplete beta function. The method is taken from
// Numerical Recipe in C, 2nd edition, section 6.4. The following web
// page calculates the incomplete beta function, which equals
// kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
//
//	http://www.danielsoper.com/statcalc/calc36.aspx
//
// /
func kf_betai(a, b, x float64) float64 {
	if x < (a+1)/(a+b+2) {
		return kf_betai_aux(a, b, x)
	}
	return 1 - kf_betai_aux(b, a, 1-x)
}

// log\binom{n}{k}
func lbinom(n, k int) float64 {
	if k < 0 || n <= 0 || k > n {
		panic("lbinom error: non-sensical input")
	}
	if k == 0 || n == k {
		return 0
	}
	return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
}

// lgamma in the C standard lib returns "the
// natural logarithm of the absolute value of
// the gamma function of num".
// -- https://en.cppreference.com/w/cpp/numeric/math/lgamma
//
// Special cased here for integer num.
//
// Note: I'm not sure that the Go standard library math.Lgamma
// is the same as math.Log(math.Abs(math.Gamma(num))),
// even assuming infinite precision rather than float64,
// and for numerical precision reasons we don't want to just
// use that expression. I don't _think_ this matters for
// the Fisher Exact test cases needed though since
// contigency table counts are non-negative intergers...
// We panic if it is.
func lgamma(num int) float64 {
	lg, sign := math.Lgamma(float64(num))
	if sign < 0 {
		panic("lgamma error: negative case not handled")
	}
	return lg
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
func hypergeo(n11, n1_, n_1, n int) float64 {

	return math.Exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1))
}

type hgacc_t struct {
	n11 int
	n1_ int
	n_1 int
	n   int
	p   float64
}

// incremental version of hypergenometric distribution
func hypergeo_acc(n11, n1_, n_1, n int, aux *hgacc_t) float64 {

	if n1_ != 0 || n_1 != 0 || n != 0 {
		aux.n11 = n11
		aux.n1_ = n1_
		aux.n_1 = n_1
		aux.n = n
	} else { // then only n11 changed; the rest fixed
		if n11%11 != 0 && n11+aux.n-aux.n1_-aux.n_1 != 0 {
			if n11 == aux.n11+1 { // incremental
				aux.p *= float64(aux.n1_-aux.n11) / float64(n11) *
					float64(aux.n_1-aux.n11) / float64(n11+aux.n-aux.n1_-aux.n_1)

				aux.n11 = n11
				return aux.p
			}
			if n11 == aux.n11-1 { // incremental
				aux.p *= float64(aux.n11) / float64(aux.n1_-n11) *
					float64(aux.n11+aux.n-aux.n1_-aux.n_1) / float64(aux.n_1-n11)
				aux.n11 = n11
				return aux.p
			}
		}
		aux.n11 = n11
	}
	aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n)
	return aux.p
}

// FisherExactTest22 computes Fisher's Exact test
// for independence on a 2x2 contingency table.
//
// n11  n12  | n1_
// n21  n22  | n2_
// ----------+-----
// n_1  n_2  | n
//
// is the layout assumed.
func FisherExactTest22(n11, n12, n21, n22 int) float64 {
	if n11 < 0 {
		panic(fmt.Sprintf("n11 was %v. must be >= 0.", n11))
	}
	if n12 < 0 {
		panic(fmt.Sprintf("n12 was %v. must be >= 0.", n12))
	}
	if n21 < 0 {
		panic(fmt.Sprintf("n21 was %v. must be >= 0.", n21))
	}
	if n22 < 0 {
		panic(fmt.Sprintf("n22 was %v. must be >= 0.", n22))
	}

	var left, right, two float64
	return kt_fisher_exact(n11, n12, n21, n22, &left, &right, &two)
}

func kt_fisher_exact(n11, n12, n21, n22 int, _left, _right, two *float64) float64 {

	var j, max, min int
	var p, q, left, right float64
	var aux hgacc_t
	var n1_, n_1, n int

	n1_ = n11 + n12
	n_1 = n11 + n21
	n = n11 + n12 + n21 + n22 // calculate n1_, n_1 and n
	// max n11, for right tail
	if n_1 < n1_ {
		max = n_1
	} else {
		max = n1_
	}
	min = n1_ + n_1 - n // not sure why n11-n22 is used instead of min(n_1,n1_)
	if min < 0 {        // min n11, for left tail
		min = 0
	}
	*two = 1
	*_left = 1
	*_right = 1
	if min == max {
		return 1 // no need to do test
	}
	q = hypergeo_acc(n11, n1_, n_1, n, &aux) // the probability of the current table

	if q == 0 {
		// https://github.com/samtools/htslib/commit/170047656473a7fadf9f23f3afdfac46b9c52b21
		if n11*(n+2) < (n_1+1)*(n1_+1) {
			// peak to right of n11
			*_left = 0
			*_right = 1
			*two = 0
			return 0
		}
		// peak to left of n11
		*_left = 1
		*_right = 0
		*two = 0
		return 0.0
	}

	// left tail
	p = hypergeo_acc(min, 0, 0, 0, &aux)
	i := min + 1
	for left = 0; p < 0.99999999*q && i <= max; i++ { // loop until underflow
		left += p
		p = hypergeo_acc(i, 0, 0, 0, &aux)
	}
	i--
	if p < 1.00000001*q {
		left += p
	} else {
		i--
	}
	// right tail
	p = hypergeo_acc(max, 0, 0, 0, &aux)
	j = max - 1
	for right = 0; p < 0.99999999*q && j >= 0; j-- { // loop until underflow
		right += p
		p = hypergeo_acc(j, 0, 0, 0, &aux)
	}
	j++
	if p < 1.00000001*q {
		right += p
	} else {
		j++
	}
	// two-tail
	*two = left + right
	if *two > 1 {
		*two = 1
	}
	// adjust left and right
	if abs(i-n11) < abs(j-n11) {
		right = 1 - left + q
	} else {
		left = 1 - right + q
	}
	*_left = left
	*_right = right
	return q
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
