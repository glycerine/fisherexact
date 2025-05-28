package fisherexact

import (
	"fmt"
	"math"
	"testing"
)

func TestFET(t *testing.T) {

	/*
		x := 5.5
		y := 3.0
		fmt.Printf("upper-gamma(%v,%v): %v\n", x, y, kf_gammaq(y, x)*math.Gamma(y))

		a := 2.0
		b := 2.0
		x = 0.5
		fmt.Printf("incomplete-beta(%v,%v,%v): %v\n", a, b, x, kf_betai(a, b, x)/math.Exp(kf_lgamma(a+b)-kf_lgamma(a)-kf_lgamma(b)))
	*/

	nfail := 0
	nfail += test_fisher(2, 1, 0, 31, 1.0, 0.005347593583, 0.005347593583, 0.005347593583)
	nfail += test_fisher(2, 1, 0, 1, 1.0, 0.5, 1.0, 0.5)
	nfail += test_fisher(3, 1, 0, 0, 1.0, 1.0, 1.0, 1.0)
	nfail += test_fisher(3, 15, 37, 45, 0.021479750169, 0.995659202564, 0.033161943699, 0.017138952733)
	nfail += test_fisher(12, 5, 29, 2, 0.044554737835, 0.994525206022, 0.080268552074, 0.039079943857)

	nfail += test_fisher(781, 23171, 4963, 2455001, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(333, 381, 801722, 7664285, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(4155, 4903, 805463, 8507517, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(4455, 4903, 805463, 8507517, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(5455, 4903, 805463, 8507517, 1.0, 0.0, 0.0, 0.0)

	nfail += test_fisher(1, 1, 100000, 1000000, 0.991735477166, 0.173555146661, 0.173555146661, 0.165290623827)
	nfail += test_fisher(1000, 1000, 100000, 1000000, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(1000, 1000, 1000000, 100000, 0.0, 1.0, 0.0, 0.0)

	nfail += test_fisher(49999, 10001, 90001, 49999, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(50000, 10000, 90000, 50000, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(50001, 9999, 89999, 50001, 1.0, 0.0, 0.0, 0.0)
	nfail += test_fisher(10000, 50000, 130000, 10000, 0.0, 1.0, 0.0, 0.0)

	if nfail > 0 {
		t.Fatalf("%v fet test case failures printed above", nfail)
	} else {
		//fmt.Printf("good: all fet tests passed.\n")
	}
}

func differ(obs, expected float64) bool {
	return math.Abs(obs-expected) > 1e-8
}

func fail(testname string, obs, expected float64, n11, n12, n21, n22 int) {
	fmt.Printf("[%d %d | %d %d] %s: %g (expected %g)\n",
		n11, n12, n21, n22, testname, obs, expected)
}

func test_fisher(n11, n12, n21, n22 int, eleft, eright, etwo, eprob float64) (nfail int) {

	var left, right, two float64
	prob := kt_fisher_exact(n11, n12, n21, n22, &left, &right, &two)

	if differ(left, eleft) {
		nfail++
		fail("LEFT", left, eleft, n11, n12, n21, n22)
	}
	if differ(right, eright) {
		nfail++
		fail("RIGHT", right, eright, n11, n12, n21, n22)
	}
	if differ(two, etwo) {
		nfail++
		fail("TWO-TAIL", two, etwo, n11, n12, n21, n22)
	}
	if differ(prob, eprob) {
		nfail++
		fail("RESULT", prob, eprob, n11, n12, n21, n22)
	}
	return
}

func TestChiSquared(t *testing.T) {

	n11 := 10
	n12 := 20
	n21 := 30
	n22 := 40
	yates := true
	pvalChi2 := ChiSquaredTest22(n11, n12, n21, n22, yates)
	//fmt.Printf("yates=%v, pval chi2 = %v\n", yates, pvalChi2)
	if want, got := 0.5040358664525046, pvalChi2; differ(want, got) {
		t.Fatalf("want %v, got %v", want, got)
	}
}
