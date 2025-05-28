package fisherexact

import (
	"fmt"
	"math"

	"github.com/glycerine/fisherexact/cephes"
)

// ChiSquaredTest22 computes a Chi-squared test
// for independene on a 2x2 contingency table.
//
// n11  n12  | n1_
// n21  n22  | n2_
// ----------+-----
// n_1  n_2  | n
//
// is the layout assumed.
func ChiSquaredTest22(n11, n12, n21, n22 int, yates bool) (pval float64) {
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

	margin1 := n11 + n12
	margin2 := n21 + n22
	col1 := n11 + n21
	col2 := n12 + n22
	n := float64(col1 + col2)

	// x__ : expected under independence
	probCol1 := float64(col1) / n
	probCol2 := float64(col2) / n
	probRow1 := float64(margin1) / n
	probRow2 := float64(margin2) / n

	x11 := probRow1 * probCol1 * n
	x12 := probRow1 * probCol2 * n
	x21 := probRow2 * probCol1 * n
	x22 := probRow2 * probCol2 * n

	expected := []float64{x11, x12, x21, x22}
	obs := []float64{float64(n11), float64(n12), float64(n21), float64(n22)}
	//fmt.Printf("obs = '%#v'\n", obs)
	//fmt.Printf("exp = '%#v'\n", expected)
	stat := chiSquaredStatistic(obs, expected, yates)
	pval = 1 - cephes.Igam(0.5, stat/2)
	//fmt.Printf("yates = %v; chi2 stat = %v -> pval %v\n", yates, stat, pval)
	return
}

func chiSquaredStatistic(obs, expected []float64, yates bool) (stat float64) {
	nobs := len(obs)
	nx := len(expected)

	if nobs != nx {
		panic(fmt.Sprintf("chiSquaredStatistic: obs(%v) vs "+
			"expected(%v) length mismatch", nobs, nx))
	}
	for i, a := range obs {
		b := expected[i]
		if a == 0 && b == 0 {
			continue
		}
		d := math.Abs(a - b)
		if yates {
			if d > 0.5 {
				d -= 0.5
			} else {
				d = 0
			}
		}
		stat += d * d / b
	}
	return
}
