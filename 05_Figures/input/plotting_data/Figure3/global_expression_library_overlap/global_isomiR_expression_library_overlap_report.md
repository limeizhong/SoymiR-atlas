# Global isomiR expression-correlation and library-overlap analysis

## Analysis scope

- Mature-arm groups with a detected canonical miRNA: 387
- Mature-arm groups without canonical but with at least two distinct isomiRs: 63
- isomiRs excluded because no canonical and no second distinct isomiR were detected: 81
- Canonical-isomiR pairs excluded because both records had the same Seq-ID: 0
- All within-locus mature-isoform pairs: 877

Pearson correlation was calculated across the 10 tissue-level expression values in S3. Library overlap was calculated as the Jaccard index of detected-library sets from 68282.1588.expression.rawdata.txt.

## Summary

### All within-locus mature-isoform pairs

- Pearson r: median 0.142; IQR -0.121 to 0.471
- Pairs with positive Pearson r: 543/877 (61.9%)
- Pairs with Pearson r >= 0.5: 193/877 (22.0%)
- Library Jaccard: median 0.000; IQR 0.000 to 0.000
- Pairs with no shared detection library: 678/877 (77.3%)

## Overall association

Pearson expression correlation and library Jaccard similarity showed a weak positive association (Spearman rho = 0.197, P = 4.37e-09).
