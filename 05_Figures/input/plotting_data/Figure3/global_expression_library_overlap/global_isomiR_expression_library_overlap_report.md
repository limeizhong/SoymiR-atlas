# Global isomiR expression-correlation and library-overlap analysis

## Analysis scope

- Mature-arm groups with a detected canonical miRNA: 414
- Mature-arm groups without canonical but with at least two distinct isomiRs: 47
- isomiRs excluded because no canonical and no second distinct isomiR were detected: 32
- Canonical-isomiR pairs excluded because both records had the same Seq-ID: 0
- All within-locus mature-isoform pairs: 888

Pearson correlation was calculated across the 10 tissue-level expression values in S3. Library overlap was calculated as the Jaccard index of detected-library sets from 68282.1588.expression.rawdata.txt.

## Summary

### All within-locus mature-isoform pairs

- Pearson r: median 0.141; IQR -0.119 to 0.467
- Pairs with positive Pearson r: 544/888 (61.3%)
- Pairs with Pearson r >= 0.5: 193/888 (21.7%)
- Library Jaccard: median 0.000; IQR 0.000 to 0.000
- Pairs with no shared detection library: 689/888 (77.6%)

## Overall association

Pearson expression correlation and library Jaccard similarity showed a weak positive association (Spearman rho = 0.199, P = 2.21e-09).
