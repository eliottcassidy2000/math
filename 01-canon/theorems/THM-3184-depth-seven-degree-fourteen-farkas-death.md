---
id: THM-3184
title: "Depth-seven degree-fourteen even-skeleton Farkas death"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For support (1,3), bank I2, no probability law on the complete 1,820-state
  physical prefix bank through depth seven has nonnegative
  partition-coarsening upset response in every degree through 14.  A readable
  positive rational combination of ten named upsets, using only degrees
  8,10,12,14, is at most -10^-11 on every state.  Thus the depth-seven
  selector cone resurrected through degree 13 in THM-3177 dies at degree 14.
  This is a finite averaged-virtual-prefix obstruction, not a Gaussian-moment
  counterexample or a no-go for deeper selector banks.
source: root/multiscale-newton-flag/product-gamma-width3/2026-08-02
audit: >
  The pure-integer/Fraction companion rebuilds all 1,820 states from two
  hash-pinned transitive sources, independently regenerates partitions and
  bin-packing coarsening, verifies every named upset and row sup-norm, checks
  all state-degree zero-mass identities, clears the readable normalized dual
  to primitive positive integers, and proves the combined coordinate is
  strictly negative on all states.  Normal, optimized, and stored replay agree
  exactly.  An independent immutable referee replayed the immutable candidate,
  independently recomputed all ten minimal antichains, matched both LF hashes,
  and accepted the Farkas logic, strict margin, hard-face anatomy, scope, and
  barcode interpretation.
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3177-depth-seven-degree-thirteen-strict-selector-resurrection
related:
  - THM-3158-depth-five-selector-resurrection-through-degree-twelve
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3169-depth-six-degree-thirteen-nonresurrection-and-mutated-cocycle-certificate
script: 04-computation/gmc_depth_seven_degree_fourteen_farkas_thm3184.py
output: 05-knowledge/results/gmc_depth_seven_degree_fourteen_farkas_thm3184.out
script_sha256: 9d560fc04fd7772a0b92b5973cd85b7ed61a2ec76328badfe5420ae2fcb7d129
output_sha256: d599057d44b12cd222a0d6ecf3755cb5d2ae59d1faa1f7bd6731ad295ab667b3
hash_basis: LF-normalized bytes
---

# THM-3184 -- depth-seven degree-fourteen even-skeleton Farkas death

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3177 proves that the support-`(1,3)`, bank-`I2` cumulative selector cone
is nonempty through degree `13` at physical prefix depth seven.  Its final
section leaves an exact trichotomy at degree `14`: strict continuation,
boundary continuation on named zero facets, or a positive rational Farkas
separator.  The third alternative occurs.  More sharply, only ten even-degree
upset inequalities are needed.

## 1. The cone and its upset coordinates

Retain THM-3177's reduced pole multiset

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1),
```

and let `S_<=7` be its `1,820` multiplicity-valid unordered nonempty
submultisets of sizes at most seven.  For a state `sigma`, degree `N`, and
partition-coarsening upset `U` in `Par(N)`, write

```text
r_(N,U)(sigma)=sum_(mu in U) G_N^sigma(mu).               (1)
```

The response has total mass zero.  Every nonempty proper upset contains
`(N)` and excludes `(1^N)`, and THM-3127 says that a probability law `lambda`
belongs to `C_D^(<=7)` only if

```text
R_(N,U)(lambda)=sum_sigma lambda_sigma r_(N,U)(sigma)>=0  (2)
```

for every such upset and every `5<=N<=D`.

## 2. Ten even-degree facets

For an antichain `A` of partitions of `N`, let `up(A)` be the upset it
generates.  The notation `1^k` means `k` parts equal to one.  The following
table defines ten response rows `r_i`.  The last two columns give a positive
rational multiplier `q_i` and the exact row sup-norm

```text
M_i=max_(sigma in S_<=7) |r_i(sigma)|.                    (3)
```

| `i` | `N` | `|U_i|` | minimal generators of `U_i` | `q_i` | `M_i` |
|---:|---:|---:|:---|:---:|---:|
| 1 | 14 | 130 | `(3,2,1^9)`, `(2,2,2,1^8)` | `429/1354` | 92,609,810,408,824,812,936,364,032 |
| 2 | 8 | 21 | `(2,1^6)` | `101/28748` | 54,841,362,155,328 |
| 3 | 10 | 40 | `(3,1^7)`, `(2,2,1^6)` | `507/24251` | 1,093,237,378,467,812,904 |
| 4 | 12 | 76 | `(2,1^10)` | `203/3189` | 732,937,666,219,205,291,328 |
| 5 | 12 | 74 | `(2,2,1^8)` | `915/24238` | 1,773,893,816,717,680,035,480 |
| 6 | 14 | 128 | `(4,1^10)`, `(3,2,1^9)`, `(2^6,1^2)` | `117/26686` | 11,449,587,449,741,073,916,566,528 |
| 7 | 14 | 121 | `(7,1^7)`, `(3,3,2,1^6)`, `(2^4,1^6)` | `30/29581` | 158,608,091,980,421,055,669,473,280 |
| 8 | 14 | 132 | `(2,2,1^10)` | `913/5536` | 6,483,094,663,875,178,504,128,552 |
| 9 | 14 | 129 | `(5,1^9)`, `(3,3,1^8)`, `(2^3,1^8)` | `5792/19469` | 68,383,863,660,622,477,237,420,032 |
| 10 | 14 | 127 | `(5,1^9)`, `(4,2,1^8)`, `(3,3,1^8)`, `(2^4,1^6)` | `1999/22331` | 136,981,330,444,844,093,510,221,824 |

The independently generated partition order verifies that each displayed set
is a nonempty proper upset and that the generator antichain and cardinality
are exact.

## 3. The pointwise Farkas inequality

Define on the complete state bank

```text
H(sigma)=sum_(i=1)^10 (q_i/M_i) r_i(sigma).               (4)
```

Exact `Fraction` replay gives the uniform pointwise inequality

```text
H(sigma)<=-1/100000000000       for every sigma in S_<=7. (5)
```

The actual largest coordinate is approximately `-3.0334308347e-11`; the
small rational bound in `(5)`, not that decimal, is the promoted margin.

For a fully integral replay, let `L` clear the denominators of all `q_i/M_i`
and divide the resulting ten integers by their gcd.  This gives a primitive
positive integer vector `a=(a_1,...,a_10)`.  The companion prints `a` and
checks independently that

```text
sum_i a_i r_i(sigma)<0             for all 1,820 states. (6)
```

Thus the normalized table is a complete readable dual certificate, while
`(6)` removes any dependence on rational comparison or floating arithmetic.

### Proof of emptiness

Suppose `lambda` were in `C_14^(<=7)`.  Every row in Section 2 is one of its
required upset inequalities, so `(2)` gives `R_i(lambda)>=0`.  Positivity of
`q_i/M_i` would imply

```text
0 <= sum_i (q_i/M_i) R_i(lambda)
   = sum_sigma lambda_sigma H(sigma)
   <= -1/100000000000,                                  (7)
```

because `lambda` is a probability law.  This contradiction proves

```text
C_14^(<=7)=empty.                                        (8)
```

Notice the stronger finite statement: even the relaxation retaining only
the ten displayed inequalities is empty.  Degrees `5,6,7,9,11,13` play no
role in this separator.

## 4. Anatomy of the first hard face

The first degree-fourteen upset has the especially simple complement

```text
Par(14)\U_1={
 (1^14), (2,1^12), (3,1^11), (2,2,1^10), (4,1^10)
}.                                                       (9)
```

Its response row has no zero coordinate: it is positive on exactly `30`
states and negative on the other `1,790`.  The positive counts by prefix
depth `1,...,7` are

```text
(1,2,4,3,3,5,12).                                       (10)
```

Consequently any putative law satisfying just this face must charge its
small positive island.  The other nine rows price that island and produce
the strict global contradiction `(5)`.  This is explanatory anatomy, not a
claim that `U_1` alone is a separator or that the ten-row certificate is the
unique or globally cardinal-minimal one.

## 5. The two-parameter selector barcode

Combining `(8)` with THM-3177 gives the sharp degree transition at depth
seven:

```text
C_13^(<=7) is nonempty,       C_14^(<=7) is empty.        (11)
```

Together with the earlier depth transition,

```text
C_13^(<=6) is empty,          C_13^(<=7) is nonempty,     (12)
```

the feasibility bifiltration has a genuine corner at `(depth,degree)=(7,13)`:
adding one resource layer resurrects a coherent mixed section, while adding
one degree kills every such section.  The degree death is already visible on
the even skeleton `8,10,12,14`.

This makes the holotopy language precise without inventing topology.  Each
nonempty feasible set is convex and hence contractible; the obstruction is
not an internal loop.  It is failure of the nested half-space section to
extend, witnessed by the positive Farkas cocircuit `(4)`.  The preserved
object is the filtered cone together with its changing dual facets.

## 6. Exact verification

Run

```text
python 04-computation/gmc_depth_seven_degree_fourteen_farkas_thm3184.py
python -O 04-computation/gmc_depth_seven_degree_fourteen_farkas_thm3184.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_depth_seven_degree_fourteen_farkas_thm3184.out.
```

The companion uses integer and `Fraction` arithmetic only.  It checks both
transitive LF hashes, reconstructs all `1,820` states and all state-degree
zero-mass identities through degree `14`, independently regenerates the
partition/coarsening order, verifies the ten upsets and all ten sup-norms,
checks `(5)`, clears the normalized dual to a primitive integer vector, and
replays `(6)` coordinatewise.  Floating column generation and support-MILP
compression are discovery history only and are absent from the proof.

## 7. Scope and next boundary

The theorem is an impossibility result for probability averages of the
specific support-`(1,3)`, bank-`I2`, fixed-`Q` virtual prefix currents at
physical depth at most seven.  It does not refute the strict degree-thirteen
law of THM-3177.  It makes no claim about depth eight or any larger bank, no
claim about a response-compatible stopping process or decomposition of the
original product-Gamma response, and no claim of a Gaussian-moment
counterexample, `NC(2)`, `GMC(2)`, or LRC consequence.

The next exact selector question is whether depth eight resurrects degree
fourteen and, if so, whether its new columns cross this even-skeleton dual or
merely mutate it as depth six did at degree thirteen.

QED.
