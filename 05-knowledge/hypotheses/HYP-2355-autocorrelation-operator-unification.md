---
id: HYP-2355
name: autocorrelation-operator-unification
status: MIXED — VERIFIED unification (3 problems) + a PROVED refinement of THM-120 (6 homometry classes)
date: 2026-06-07
session: claudebox-2026-06-07-S714
depends_on:
  - THM-120   # |Pf| separates phases at n=6 (prose says "5 triples"; its TABLE lists 6 — see Clarification)
  - THM-174   # det(I+2A)=Pf(S)^2
  - THM-415   # signed-LRC homometry = same Patterson power spectrum |f-hat|^2
  - THM-435   # Pf is the even/odd seam (S713)
  - THM-438   # Paley cluster integrals = Catalan
  - HYP-2350  # unit-distance product anatomy (S712)
provisional_id: true
---

# HYP-2355: The autocorrelation operator MM* unifies the cluster; it refines the Pfaffian classification

## The lens (operator + adjoint + autocorrelation)

Every cluster problem supplies an **operator M** and its **adjoint M\***; the **autocorrelation**
`C = M M*` (positive semidefinite) is the master object, and the answer is an invariant of `C`:
a determinant (Pfaffian), a peak (unit distances), or its flatness (Paley). Verified instances:

| problem | operator M | adjoint M* | autocorrelation MM* | invariant |
|---|---|---|---|---|
| tournament (S713) | skew-adj `S` | `S^T = -S` | `C = -S^2` (PSD) | `|Pf(S)| = sqrt(det C)`; spectrum |
| signed-LRC (THM-415) | `f-hat_eps` (half-system DFT) | conj | `|f-hat|^2` power spectrum | homometry classes |
| unit distance (S711/2) | point masses `sum delta_{p_i}` | reflection | Patterson `A_P` (diff multiset) | peak `A_P(1) = u(n)` |
| Paley (THM-438) | Legendre character `chi` | conj | flat (`|Gauss|^2 = p`) | Catalan moments |
| covering-depth (THM-406) | arc indicators | — | pair overlaps `S_k` | `p_0 = sum(-1)^k S_k` |

Unifying identity: `|Pf(S)| = sqrt(det(MM*))` — the Pfaffian (S713) is the square-root-determinant of
the autocorrelation operator. The "answer" is always a functional of the autocorrelation spectrum.

## RESULT 1 (PROVED, exhaustive): the autocorrelation spectrum refines |Pf| — 6 homometry classes at n=6

`C = -S^2` has eigenvalues `{lambda_k^2}` (each doubled); its spectrum signature is
`(e1,e2,e3) = (sum lambda^2, sum lambda^2 lambda'^2, prod lambda^2)`, with `e1 = n-1` (trace, always 15
at n=6) and `e3 = |Pf|^2` (verified `32768/32768`). Exhaustively over all 32768 tournaments on 6
vertices, there are **exactly 6 distinct spectra**:

```
   (15, 15, 1)  |Pf|=1   3840 tourns      (15, 55, 25) |Pf|=5   4608
   (15, 47, 1)  |Pf|=1  11520 tourns      (15, 63, 49) |Pf|=7   3840
   (15, 39, 9)  |Pf|=3   7680 tourns      (15, 63, 81) |Pf|=9   1280
```

`|Pf| = sqrt(e3)` takes only **5** values `{1,3,5,7,9}` — it MERGES the two `|Pf|=1` classes
`(15,15,1)` and `(15,47,1)`, which the full autocorrelation spectrum (via `e2`) SEPARATES. Likewise the
char-ratio/`e2` value `63` is shared by `|Pf|=7` and `9` (separated by `e3`). So the homometry partition
(equal autocorrelation spectrum) is the join of the `|Pf|` partition and the `e2` partition, and is
strictly finer than `|Pf|` alone. **Clarification of THM-120:** its prose states "exactly 5 skew
eigenvalue triples" but its own table lists 6 rows; this exhaustive recomputation confirms **6** spectra
(the "5" is the count of distinct `|Pf|`). Flagged for a maintainer; not silently editing canon.

Each spectrum class carries a characteristic H (#Hamiltonian-paths) distribution; the densest
autocorrelation class `(15,63,81)` (`|Pf|=9`, S-phase) has the narrowest H-support `{9,15,29}`, while the
flat-ish `(15,47,1)` (`|Pf|=1`) has the widest — `|Pf|` (= det of autocorrelation) anti-tracks H-spread.

## RESULT 2 (VERIFIED): the Patterson autocorrelation FACTORIZES over Minkowski products

For a Minkowski sum `P (+) Q` the difference multiset (Patterson) is the **convolution**
`A_{P(+)Q} = A_P * A_Q` (operator tensor -> autocorrelation convolution). Hence the unit-distance count
is a convolution peak: `u(P(+)Q) = (A_P * A_Q)(radius 1)`. Verified on `K3 [] W7`: direct count `57`,
convolution-formula `|K3| e(W7) + |W7| e(K3) = 3*12 + 7*3 = 57`, Patterson-convolution peak `57`. This is
the autocorrelation-lens form of S712/THM-431, and it re-derives the n=22 obstruction: `22 = 2*11` forces
`A = A_{K2} * A_{(11)}`, whose radius-1 peak is capped at `u(2)*11 + 2*u(11) = 57 < 60`.

## RESULT 3 (VERIFIED): Paley = flat autocorrelation (the perfect operator)

For `p = 3 mod 4` the QR set `Q` is a Paley difference set: its autocorrelation
`c(s) = #{(a,b) in Q^2 : a - b = s}` is **constant** `= (p-3)/4` for all `s != 0` (verified
`p = 7,11,19,23`), i.e. `MM*` is flat off the DC term and `|Gauss sum|^2 = p`.

## The synthesis claim (the hypothesis)

**Flat autocorrelation = quasirandom rigidity = the HARD cases; peaked/convolutional autocorrelation =
extremal structure = the EASY/solved cases.** The LRC transversal core (S642/THM-420), the signed-LRC
homometric collisions (THM-415), and the Paley maximizer are flat-autocorrelation (rigid); the
unit-distance product optima (S712), the AP/lattice LRC tights, and the high-`|Pf|` tournaments are
peaked/convolutional. The recurring multiplicative-vs-additive / quasirandom-core dichotomy is, in one
sentence, **flatness vs peakedness of the autocorrelation operator MM\***.

## Next
- prove RESULT 1's class count for general even n (is the homometry partition = join of `|Pf|` and the
  `e2`/char-ratio partition at all n?);
- the off-diagonal of `-S^2` is the 2-path pair-autocorrelation — relate to H via a spectral formula;
- test the synthesis claim as a *predictive* hardness measure: compute autocorrelation flatness of the
  open LRC frontier configs (n=15,19,21,22) and the unit-distance n=22 candidates.
