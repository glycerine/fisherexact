package fisherexact

import (
	"fmt"
	"math"
	"testing"
)

// AltHypothesisFET specifies the null hypothesis
// for the FET. Only used in testing at the moment.
type AltHypothesisFET int

const (
	TWO_SIDED AltHypothesisFET = 0
	LESS      AltHypothesisFET = 1
	GREATER   AltHypothesisFET = 2
)

func (alt AltHypothesisFET) String() string {
	switch alt {
	case TWO_SIDED:
		return "TWO_SIDED"
	case LESS:
		return "LESS"
	case GREATER:
		return "GREATER"
	}
	panic(fmt.Sprintf("unknown AltHypothesisFET: %v", int(alt)))
}

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

	// test cases from
	// https://github.com/samtools/htslib/commit/170047656473a7fadf9f23f3afdfac46b9c52b21

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
	/*
		for (i in 1:20) { y=sample(15,4); x=matrix(data=y,byrow=T, nrow=2, ncol=2); t=fisher.test(x); cat(paste0("testcase{dat:[]int{",paste0(y,collapse=","),"}, xpval: ",t$p.value, ", alt:TWO_SIDED, oddsRatioEstimate:",t$estimate,", nullValue:",t$null.value,"},\n")) }
		for (i in 1:20) { y=sample(15,4); x=matrix(data=y,byrow=T, nrow=2, ncol=2); t=fisher.test(x, alternative="greater"); cat(paste0("testcase{dat:[]int{",paste0(y,collapse=","),"}, xpval: ",t$p.value, ", alt:GREATER, oddsRatioEstimate:",t$estimate,", nullValue:",t$null.value,"},\n")) }
		for (i in 1:20) { y=sample(15,4); x=matrix(data=y,byrow=T, nrow=2, ncol=2); t=fisher.test(x, alternative="less"); cat(paste0("testcase{dat:[]int{",paste0(y,collapse=","),"}, xpval: ",t$p.value, ", alt:LESS, oddsRatioEstimate:",t$estimate,", nullValue:",t$null.value,"},\n")) }
	*/
	cases := []testcase{
		testcase{dat: []int{10, 11, 6, 2}, xpval: 0.237847742795269, alt: TWO_SIDED, oddsRatioEstimate: 0.315539110927868, nullValue: 1},
		testcase{dat: []int{8, 2, 5, 13}, xpval: 0.0163234172387491, alt: TWO_SIDED, oddsRatioEstimate: 9.41594304803736, nullValue: 1},
		testcase{dat: []int{13, 12, 15, 6}, xpval: 0.231825852388689, alt: TWO_SIDED, oddsRatioEstimate: 0.441400756534335, nullValue: 1},
		testcase{dat: []int{7, 4, 13, 3}, xpval: 0.391304347826087, alt: TWO_SIDED, oddsRatioEstimate: 0.41837325843701, nullValue: 1},
		testcase{dat: []int{11, 7, 15, 2}, xpval: 0.121205379714835, alt: TWO_SIDED, oddsRatioEstimate: 0.219084900354817, nullValue: 1},
		testcase{dat: []int{13, 3, 7, 4}, xpval: 0.391304347826087, alt: TWO_SIDED, oddsRatioEstimate: 2.39021012895488, nullValue: 1},
		testcase{dat: []int{10, 2, 9, 5}, xpval: 0.391304347826087, alt: TWO_SIDED, oddsRatioEstimate: 2.67184390683078, nullValue: 1},
		testcase{dat: []int{4, 15, 12, 5}, xpval: 0.00636407114135991, alt: TWO_SIDED, oddsRatioEstimate: 0.119658009316562, nullValue: 1},
		testcase{dat: []int{15, 3, 4, 11}, xpval: 0.00160836004285247, alt: TWO_SIDED, oddsRatioEstimate: 12.3872815442657, nullValue: 1},
		testcase{dat: []int{3, 14, 4, 11}, xpval: 0.678269064392535, alt: TWO_SIDED, oddsRatioEstimate: 0.599213115677309, nullValue: 1},
		testcase{dat: []int{13, 3, 1, 11}, xpval: 0.000341005967604433, alt: TWO_SIDED, oddsRatioEstimate: 38.5534211622107, nullValue: 1},
		testcase{dat: []int{2, 5, 9, 12}, xpval: 0.668338907469342, alt: TWO_SIDED, oddsRatioEstimate: 0.545029610109094, nullValue: 1},
		testcase{dat: []int{10, 9, 3, 12}, xpval: 0.078991149586497, alt: TWO_SIDED, oddsRatioEstimate: 4.24637831104569, nullValue: 1},
		testcase{dat: []int{13, 4, 1, 9}, xpval: 0.00133909653360454, alt: TWO_SIDED, oddsRatioEstimate: 24.8338694584268, nullValue: 1},
		testcase{dat: []int{14, 2, 9, 13}, xpval: 0.00645479609033885, alt: TWO_SIDED, oddsRatioEstimate: 9.47094406080515, nullValue: 1},
		testcase{dat: []int{5, 2, 15, 13}, xpval: 0.672173627262615, alt: TWO_SIDED, oddsRatioEstimate: 2.12118865938244, nullValue: 1},
		testcase{dat: []int{13, 5, 2, 1}, xpval: 1, alt: TWO_SIDED, oddsRatioEstimate: 1.28323666262214, nullValue: 1},
		testcase{dat: []int{9, 1, 13, 3}, xpval: 1, alt: TWO_SIDED, oddsRatioEstimate: 2.02340688404544, nullValue: 1},
		testcase{dat: []int{7, 9, 2, 10}, xpval: 0.223188405797101, alt: TWO_SIDED, oddsRatioEstimate: 3.70666515096608, nullValue: 1},
		testcase{dat: []int{4, 9, 3, 13}, xpval: 0.666500083291688, alt: TWO_SIDED, oddsRatioEstimate: 1.88202746263798, nullValue: 1},
		testcase{dat: []int{3, 12, 13, 15}, xpval: 0.981883029573569, alt: GREATER, oddsRatioEstimate: 0.296690669915481, nullValue: 1},
		testcase{dat: []int{4, 11, 8, 3}, xpval: 0.997287656481357, alt: GREATER, oddsRatioEstimate: 0.149216303259255, nullValue: 1},
		testcase{dat: []int{12, 2, 4, 11}, xpval: 0.00192576570331965, alt: GREATER, oddsRatioEstimate: 14.5399725611462, nullValue: 1},
		testcase{dat: []int{15, 8, 2, 3}, xpval: 0.29010989010989, alt: GREATER, oddsRatioEstimate: 2.70428203532748, nullValue: 1},
		testcase{dat: []int{3, 2, 13, 11}, xpval: 0.603831417624521, alt: GREATER, oddsRatioEstimate: 1.25893261186031, nullValue: 1},
		testcase{dat: []int{6, 12, 3, 2}, xpval: 0.943831911795298, alt: GREATER, oddsRatioEstimate: 0.350852026925064, nullValue: 1},
		testcase{dat: []int{13, 12, 15, 6}, xpval: 0.95129464442598, alt: GREATER, oddsRatioEstimate: 0.441400756534335, nullValue: 1},
		testcase{dat: []int{7, 9, 5, 14}, xpval: 0.234223210975157, alt: GREATER, oddsRatioEstimate: 2.12854504573637, nullValue: 1},
		testcase{dat: []int{5, 12, 1, 13}, xpval: 0.134470399208998, alt: GREATER, oddsRatioEstimate: 5.15574676453379, nullValue: 1},
		testcase{dat: []int{4, 12, 14, 15}, xpval: 0.969437623683502, alt: GREATER, oddsRatioEstimate: 0.365337520344086, nullValue: 1},
		testcase{dat: []int{10, 1, 5, 11}, xpval: 0.00286863792046186, alt: GREATER, oddsRatioEstimate: 19.2437872603241, nullValue: 1},
		testcase{dat: []int{6, 13, 12, 3}, xpval: 0.999419205353335, alt: GREATER, oddsRatioEstimate: 0.124287057961696, nullValue: 1},
		testcase{dat: []int{14, 7, 9, 15}, xpval: 0.0485468092201357, alt: GREATER, oddsRatioEstimate: 3.24016904405364, nullValue: 1},
		testcase{dat: []int{9, 3, 7, 6}, xpval: 0.248135684479741, alt: GREATER, oddsRatioEstimate: 2.47418567997489, nullValue: 1},
		testcase{dat: []int{4, 15, 2, 8}, xpval: 0.669091827712518, alt: GREATER, oddsRatioEstimate: 1.06432110396692, nullValue: 1},
		testcase{dat: []int{4, 2, 12, 5}, xpval: 0.760834893558006, alt: GREATER, oddsRatioEstimate: 0.840097575350736, nullValue: 1},
		testcase{dat: []int{14, 15, 3, 12}, xpval: 0.064839388069586, alt: GREATER, oddsRatioEstimate: 3.62580938021925, nullValue: 1},
		testcase{dat: []int{2, 4, 10, 15}, xpval: 0.773635337595293, alt: GREATER, oddsRatioEstimate: 0.75684388280389, nullValue: 1},
		testcase{dat: []int{10, 11, 12, 9}, xpval: 0.822850937641733, alt: GREATER, oddsRatioEstimate: 0.688117785024529, nullValue: 1},
		testcase{dat: []int{13, 10, 12, 1}, xpval: 0.997749557121081, alt: GREATER, oddsRatioEstimate: 0.114324038136269, nullValue: 1},
		testcase{dat: []int{12, 6, 9, 15}, xpval: 0.98610866489268, alt: LESS, oddsRatioEstimate: 3.23384636554812, nullValue: 1},
		testcase{dat: []int{5, 8, 6, 14}, xpval: 0.811262110880044, alt: LESS, oddsRatioEstimate: 1.4414247884631, nullValue: 1},
		testcase{dat: []int{4, 5, 6, 7}, xpval: 0.639318885448917, alt: LESS, oddsRatioEstimate: 0.936282224569586, nullValue: 1},
		testcase{dat: []int{4, 15, 7, 14}, xpval: 0.30504325665616, alt: LESS, oddsRatioEstimate: 0.541756234023881, nullValue: 1},
		testcase{dat: []int{5, 14, 13, 3}, xpval: 0.00154101834713142, alt: LESS, oddsRatioEstimate: 0.0901376643714314, nullValue: 1},
		testcase{dat: []int{10, 14, 5, 15}, xpval: 0.931750488709851, alt: LESS, oddsRatioEstimate: 2.10595426290551, nullValue: 1},
		testcase{dat: []int{9, 2, 15, 11}, xpval: 0.967259966147619, alt: LESS, oddsRatioEstimate: 3.20224181299334, nullValue: 1},
		testcase{dat: []int{6, 3, 13, 14}, xpval: 0.912579633936697, alt: LESS, oddsRatioEstimate: 2.10877767750628, nullValue: 1},
		testcase{dat: []int{4, 14, 7, 1}, xpval: 0.00327407146629115, alt: LESS, oddsRatioEstimate: 0.0478975493732836, nullValue: 1},
		testcase{dat: []int{7, 5, 6, 10}, xpval: 0.930533494862476, alt: LESS, oddsRatioEstimate: 2.26153897377153, nullValue: 1},
		testcase{dat: []int{6, 2, 10, 8}, xpval: 0.918535469107552, alt: LESS, oddsRatioEstimate: 2.32304100453932, nullValue: 1},
		testcase{dat: []int{8, 5, 9, 4}, xpval: 0.5, alt: LESS, oddsRatioEstimate: 0.720538688839083, nullValue: 1},
		testcase{dat: []int{11, 6, 14, 2}, xpval: 0.131127740137751, alt: LESS, oddsRatioEstimate: 0.272666521958173, nullValue: 1},
		testcase{dat: []int{11, 10, 14, 13}, xpval: 0.627978975549828, alt: LESS, oddsRatioEstimate: 1.0209833267052, nullValue: 1},
		testcase{dat: []int{1, 4, 12, 13}, xpval: 0.260536398467433, alt: LESS, oddsRatioEstimate: 0.281724671612801, nullValue: 1},
		testcase{dat: []int{7, 1, 9, 6}, xpval: 0.973751514335711, alt: LESS, oddsRatioEstimate: 4.39200391118143, nullValue: 1},
		testcase{dat: []int{2, 1, 5, 9}, xpval: 0.948529411764706, alt: LESS, oddsRatioEstimate: 3.3259179969715, nullValue: 1},
		testcase{dat: []int{6, 12, 2, 14}, xpval: 0.969545286008022, alt: LESS, oddsRatioEstimate: 3.37603351887115, nullValue: 1},
		testcase{dat: []int{4, 1, 11, 15}, xpval: 0.982326041280435, alt: LESS, oddsRatioEstimate: 5.1773447503855, nullValue: 1},
		testcase{dat: []int{7, 15, 5, 1}, xpval: 0.0360885491320274, alt: LESS, oddsRatioEstimate: 0.102149685147939, nullValue: 1},
	}

	for i, c := range cases {
		// left, right, two float64
		var fisherPval float64
		left, right, two, probCurrentTable := FisherExactTest22(c.dat[0], c.dat[1], c.dat[2], c.dat[3])
		_ = probCurrentTable
		two2 := TwoSidedTest22(c.dat[0], c.dat[1], c.dat[2], c.dat[3])
		if two2 != two {
			t.Fatalf("i=%v, TwoSided22(%v) and FisherExactTest22.two(%v) disagree", i, two2, two)
		}
		//fmt.Printf("fet: i=%v, input:'%#v', left=%v, right=%v, two=%v, probCurrentTable=%v, oddsRatioEstimate=%v\n", i, c.dat, left, right, two, probCurrentTable, c.oddsRatioEstimate)
		switch c.alt {
		case LESS:
			fisherPval = left
		case GREATER:
			fisherPval = right
		case TWO_SIDED:
			fisherPval = two
		}
		wantPval := c.xpval
		// R computes STATISTIC <- -sum(lfactorial(x)) internally, but not returned.
		// so we cannot compare that. Nor do I see prob current table in its output.
		// wantStat := c.xstat
		if differ(wantPval, fisherPval) {
			fmt.Printf("%v fet case %v: want %v, got %v\n", c.alt, i, wantPval, fisherPval)
			t.Fatalf("%v fet 2-sided pval, case %v: want %v, got %v", c.alt, i, wantPval, fisherPval)
		}
		//if differ(wantProbCurTable, probCurrentTable) {
		//	fmt.Printf("%v case %v: want %v, got %v\n", c.alt, i, wantProbCurTable, probCurrentTable)
		//	t.Fatalf("%v fet test statistic, case %v: want %v, got %v", c.alt, i, wantProbCurTable, probCurrentTable)
		// }
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
