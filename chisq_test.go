package fisherexact

import (
	"testing"
)

func TestChiSquaredTest22(t *testing.T) {

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

	// compare to R
	//> yates=F
	//> yates=T
	//> for (i in 1:10) { y=sample(15,4); x=matrix(data=y,byrow=T, nrow=2, ncol=2); t=chisq.test(x,correct=yates); print(paste0("dat=[]int{",paste0(y,collapse=","),"}; xstat = ",t$statistic,"; xpval = ",t$p.value,"; yates=", yates))}

	cases := []testcase{
		testcase{dat: []int{14, 5, 6, 3}, xstat: 8.15386929811777e-32, xpval: 1, yates: true},
		testcase{dat: []int{2, 5, 14, 8}, xstat: 1.4125639985015, xpval: 0.234631220847672, yates: true},
		testcase{dat: []int{2, 9, 3, 10}, xstat: 1.58875483009546e-31, xpval: 1, yates: true},
		testcase{dat: []int{8, 6, 14, 3}, xstat: 1.30268546812664, xpval: 0.253723274036572, yates: true},
		testcase{dat: []int{15, 8, 13, 4}, xstat: 0.175374497625137, xpval: 0.675378834564654, yates: true},
		testcase{dat: []int{10, 13, 14, 15}, xstat: 0.00417648318697791, xpval: 0.948472008984082, yates: true},
		testcase{dat: []int{12, 13, 2, 10}, xstat: 2.18332556935818, xpval: 0.139512716626223, yates: true},
		testcase{dat: []int{14, 11, 9, 2}, xstat: 1.22984493767102, xpval: 0.267437196119282, yates: true},
		testcase{dat: []int{13, 9, 12, 3}, xstat: 0.953123737373737, xpval: 0.32892544938928, yates: true},
		testcase{dat: []int{11, 6, 4, 7}, xstat: 1.16791443850267, xpval: 0.279830181962734, yates: true},
		testcase{dat: []int{6, 15, 2, 12}, xstat: 0.972222222222222, xpval: 0.324126587757619, yates: false},
		testcase{dat: []int{14, 8, 5, 15}, xstat: 6.3126690243395, xpval: 0.0119878241393436, yates: false},
		testcase{dat: []int{12, 4, 1, 6}, xstat: 7.3043956043956, xpval: 0.00687861300445255, yates: false},
		testcase{dat: []int{3, 7, 11, 2}, xstat: 7.07832722832723, xpval: 0.00780218108505649, yates: false},
		testcase{dat: []int{7, 4, 10, 5}, xstat: 0.0257476728064964, xpval: 0.872518087557924, yates: false},
		testcase{dat: []int{1, 8, 11, 2}, xstat: 11.5891737891738, xpval: 0.000663368781092197, yates: false},
		testcase{dat: []int{13, 8, 7, 5}, xstat: 0.0407967032967033, xpval: 0.839930844519582, yates: false},
		testcase{dat: []int{4, 9, 8, 14}, xstat: 0.11350967872707, xpval: 0.736183257997292, yates: false},
		testcase{dat: []int{13, 8, 9, 6}, xstat: 0.0133580705009276, xpval: 0.907987697529532, yates: false},
		testcase{dat: []int{2, 11, 12, 5}, xstat: 9.01987718164189, xpval: 0.0026705933703964, yates: false},
	}

	for i, c := range cases {
		gotPval := ChiSquaredTest22(c.dat[0], c.dat[1], c.dat[2], c.dat[3], c.yates)
		wantPval := c.xpval
		if differ(wantPval, gotPval) {
			t.Fatalf("case %v: want %v, got %v", i, wantPval, gotPval)
		}
	}
}

type testcase struct {
	dat   []int
	xstat float64
	xpval float64
	yates bool
}
