package fisherexact

import (
	"fmt"
	"math"
	"testing"
)

type AltHypothesis int

const (
	TWO_SIDED AltHypothesis = 0
	LESS      AltHypothesis = 1
	GREATER   AltHypothesis = 2
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

	// compare to R
	// for (i in 1:20) { y=sample(15,4); x=matrix(data=y,byrow=T, nrow=2, ncol=2); t=fisher.test(x); print(paste0("testcase{dat:[]int{",paste0(y,collapse=","),"}, xpval: ",t$p.value, "},"))}
	//  two-side
	cases := []testcase{
		testcase{dat: []int{7, 11, 14, 12}, xpval: 0.372583797871365, alt: TWO_SIDED},
		testcase{dat: []int{15, 3, 4, 7}, xpval: 0.016938783355575, alt: TWO_SIDED},
		testcase{dat: []int{9, 7, 2, 3}, xpval: 0.635117204776647, alt: TWO_SIDED},
		testcase{dat: []int{12, 3, 13, 8}, xpval: 0.295073229514275, alt: TWO_SIDED},
		testcase{dat: []int{9, 11, 14, 13}, xpval: 0.770160207175819, alt: TWO_SIDED},
		testcase{dat: []int{1, 2, 12, 15}, xpval: 1, alt: TWO_SIDED},
		testcase{dat: []int{13, 1, 2, 8}, xpval: 0.000489482250149905, alt: TWO_SIDED},
		testcase{dat: []int{6, 12, 14, 4}, xpval: 0.0175920837235342, alt: TWO_SIDED},
		testcase{dat: []int{6, 5, 10, 7}, xpval: 1, alt: TWO_SIDED},
		testcase{dat: []int{2, 8, 10, 4}, xpval: 0.0360748418360478, alt: TWO_SIDED},
		testcase{dat: []int{7, 6, 14, 1}, xpval: 0.0286231884057971, alt: TWO_SIDED},
		testcase{dat: []int{14, 11, 8, 4}, xpval: 0.723536662357574, alt: TWO_SIDED},
		testcase{dat: []int{9, 8, 5, 11}, xpval: 0.296003611097677, alt: TWO_SIDED},
		testcase{dat: []int{11, 15, 2, 7}, xpval: 0.431106458156121, alt: TWO_SIDED},
		testcase{dat: []int{1, 2, 5, 6}, xpval: 1, alt: TWO_SIDED},
		testcase{dat: []int{8, 1, 5, 14}, xpval: 0.0036231884057971, alt: TWO_SIDED},
		testcase{dat: []int{11, 6, 5, 9}, xpval: 0.155613477924309, alt: TWO_SIDED},
		testcase{dat: []int{15, 13, 8, 3}, xpval: 0.470708489029499, alt: TWO_SIDED},
		testcase{dat: []int{5, 15, 1, 8}, xpval: 0.632815460401667, alt: TWO_SIDED},
		testcase{dat: []int{13, 2, 5, 11}, xpval: 0.00318872821653688, alt: TWO_SIDED},
	}

	// alternative="greater"
	//for (i in 1:20) { y=sample(15,4); x=matrix(data=y,byrow=T, nrow=2, ncol=2); t=fisher.test(x, alternative="greater"); print(paste0("testcase{dat:[]int{",paste0(y,collapse=","),"}, xpval: ",t$p.value, "},"))}
	cases = append(cases, []testcase{
		testcase{dat: []int{5, 4, 3, 9}, xpval: 0.165634674922601, alt: GREATER},
		testcase{dat: []int{6, 10, 4, 2}, xpval: 0.956656346749226, alt: GREATER},
		testcase{dat: []int{5, 15, 6, 10}, xpval: 0.879433938944506, alt: GREATER},
		testcase{dat: []int{14, 5, 1, 11}, xpval: 0.000477180764456481, alt: GREATER},
		testcase{dat: []int{7, 2, 11, 12}, xpval: 0.126822395253986, alt: GREATER},
		testcase{dat: []int{4, 6, 11, 13}, xpval: 0.753365692520309, alt: GREATER},
		testcase{dat: []int{11, 9, 6, 4}, xpval: 0.740596368482425, alt: GREATER},
		testcase{dat: []int{6, 9, 14, 12}, xpval: 0.880919594001881, alt: GREATER},
		testcase{dat: []int{12, 10, 8, 3}, xpval: 0.919030001128468, alt: GREATER},
		testcase{dat: []int{13, 15, 10, 4}, xpval: 0.970413988149141, alt: GREATER},
		testcase{dat: []int{7, 13, 9, 4}, xpval: 0.989407984179953, alt: GREATER},
		testcase{dat: []int{12, 15, 3, 11}, xpval: 0.133232339749066, alt: GREATER},
		testcase{dat: []int{8, 2, 3, 11}, xpval: 0.00693229236774801, alt: GREATER},
		testcase{dat: []int{14, 1, 9, 11}, xpval: 0.00317018909899889, alt: GREATER},
		testcase{dat: []int{12, 3, 6, 10}, xpval: 0.0200246226632015, alt: GREATER},
		testcase{dat: []int{2, 5, 9, 3}, xpval: 0.993728665555291, alt: GREATER},
		testcase{dat: []int{10, 11, 13, 9}, xpval: 0.855338174211227, alt: GREATER},
		testcase{dat: []int{14, 3, 13, 7}, xpval: 0.209242866083801, alt: GREATER},
		testcase{dat: []int{10, 13, 2, 3}, xpval: 0.642735042735043, alt: GREATER},
		testcase{dat: []int{9, 12, 14, 7}, xpval: 0.96919820333098, alt: GREATER},
	}...)

	// alternative="less"
	// > for (i in 1:20) { y=sample(15,4); x=matrix(data=y,byrow=T, nrow=2, ncol=2); t=fisher.test(x, alternative="less"); print(paste0("testcase{dat:[]int{",paste0(y,collapse=","),"}, xpval: ",t$p.value, "},"))}
	cases = append(cases, []testcase{
		testcase{dat: []int{8, 5, 3, 12}, xpval: 0.996300533943555, alt: LESS},
		testcase{dat: []int{8, 6, 12, 4}, xpval: 0.258970514742629, alt: LESS},
		testcase{dat: []int{6, 4, 7, 9}, xpval: 0.886902678691614, alt: LESS},
		testcase{dat: []int{12, 15, 9, 2}, xpval: 0.0384962270390524, alt: LESS},
		testcase{dat: []int{14, 8, 11, 3}, xpval: 0.2853170189099, alt: LESS},
		testcase{dat: []int{2, 15, 14, 12}, xpval: 0.00546833057683645, alt: LESS},
		testcase{dat: []int{14, 1, 6, 2}, xpval: 0.968379446640316, alt: LESS},
		testcase{dat: []int{12, 13, 11, 15}, xpval: 0.754645946510281, alt: LESS},
		testcase{dat: []int{15, 3, 10, 11}, xpval: 0.996768180861618, alt: LESS},
		testcase{dat: []int{14, 4, 10, 9}, xpval: 0.975332093574585, alt: LESS},
		testcase{dat: []int{13, 3, 5, 8}, xpval: 0.997385517767432, alt: LESS},
		testcase{dat: []int{5, 3, 4, 12}, xpval: 0.987261224439849, alt: LESS},
		testcase{dat: []int{1, 14, 15, 3}, xpval: 1.06213292489424e-05, alt: LESS},
		testcase{dat: []int{7, 4, 3, 8}, xpval: 0.98501343857381, alt: LESS},
		testcase{dat: []int{11, 12, 15, 4}, xpval: 0.0390280809581388, alt: LESS},
		testcase{dat: []int{2, 8, 6, 1}, xpval: 0.0133689839572193, alt: LESS},
		testcase{dat: []int{6, 4, 12, 5}, xpval: 0.439345399070799, alt: LESS},
		testcase{dat: []int{15, 3, 5, 2}, xpval: 0.88695652173913, alt: LESS},
		testcase{dat: []int{12, 15, 11, 8}, xpval: 0.274948536192989, alt: LESS},
		testcase{dat: []int{13, 1, 12, 14}, xpval: 0.999807928061543, alt: LESS},
	}...)

	for i, c := range cases {
		// left, right, two float64
		var fisherPval float64
		left, right, two := FisherExactTest22(c.dat[0], c.dat[1], c.dat[2], c.dat[3])
		switch c.alt {
		case LESS:
			fisherPval = left
		case GREATER:
			fisherPval = right
		case TWO_SIDED:
			fisherPval = two
		}
		wantPval := c.xpval
		if differ(wantPval, fisherPval) {
			fmt.Printf("case %v: want %v, got %v\n", i, wantPval, fisherPval)
			//t.Fatalf("case %v: want %v, got %v", i, wantPval, fisherPval)
		}
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
