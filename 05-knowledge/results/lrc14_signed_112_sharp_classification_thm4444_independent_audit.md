# Independent audit of THM-4444

**Verdict: PASS.** The sign cones, complete carrier ray, pointwise selector,
continuum integral, finite cutoffs, sharp constants, and \(6/77\) locus all
reproduce without correction.

## Independent path

The referee imports no candidate module. It hash-pins the pre-existing
mixed-parity literal six-sheet engine, independently generates the twelve
signed coefficient placements, and finds exactly three viable sorted cones:
\[
 (2,1,-1),\quad(1,-2,1),\quad(1,2,-1).
\]
For every typed row through height 611 it compares literal sheet geometry
with raw integer carrier enumeration, verifies that the raw set is precisely
the predicted primitive ray, and integrates the normalized profiles using
generic breakpoints rather than the candidate's cone formulas.

## Reproduced census

```text
rows through height 611              28,438
F1/F2/F3                             9,482 / 14,214 / 4,742
strictly above 6/77                  (2,11,20), value 11/140
equal to 6/77                        (1,5,11)
universal continuum integral         9/196
first universal strict-safe height   35
```

The three family leaders and full projection packets agree exactly. The
ordinary run executes 7,484,368 pinned-engine checks and 199,075 referee
checks. The optimized replay has the same transcript.

The proof-line review confirms that norm four alone forces the carrier
covector to zero, the coefficient-two profile is pointwise least in all
three cones, and the strict upper quadrature bound handles the cutoff equality
correctly. This remains a local classification, not chart entry or LRC(14).
