---
id: THM-3155
title: "Sharp depth-four selector resurrection and degree-twelve death barcode"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  At support (1,3), bank I2, the cumulative selector space through degree 11
  is empty for all physical pole prefixes of depth at most three but becomes
  nonempty at depth four.  An explicit four-state rational law is strictly
  positive on all 23,941 nontrivial coarsening upsets in degrees 5 through 11.
  Adding degree 12 kills the entire 334-state depth-at-most-four bank: five
  upset rows have a primitive positive integer combination strictly negative
  on every state.  Thus depth four and degree 11 form an exact birth/death
  barcode cell.  This is an averaged virtual-prefix theorem, not a sequential
  stopping process or an original-response decomposition.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
audit: >
  Two independent immutable audits reconstructed the depth/horizon
  monotonicity, the 8+33+93+200 state census, the exact four-state law, all
  23,941 strict upset inequalities and their unique minima, and every exact
  max-flow equality through degree 11.  They separately rederived the
  degree-12 singleton/max-flow hostile and the primitive five-positive-row
  Farkas contradiction on all 334 states, including the exact range and
  equality sets.  Both audits confirm the length-singleton compression and
  the fixed-Q averaged-current scope.  Fresh normal and optimized runs
  byte-match the stored 2,128-byte transcript and both declared LF hashes.
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3144-mixed-depth-selector-persistence-death-barcode
  - THM-3149-depth-three-selector-persistence-and-cross-support-wall
related:
  - MISTAKE-354
  - THM-3147-length-singleton-endpoint-jet-facet-observer
  - THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary
script: 04-computation/gmc_depth_four_selector_resurrection_barcode_scout.py
output: 05-knowledge/results/gmc_depth_four_selector_resurrection_barcode_scout.out
script_sha256: ad3d7c14895f10f16622a0918317d20785cba540e48405977019f6706f37d07b
output_sha256: c301fe220893bb5b2becd0c33fd063f54a72168343317c5f28f3603cb3173d42
hash_basis: LF-normalized bytes
---

# THM-3155 -- sharp depth-four selector resurrection and degree-twelve death barcode

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3149 proves that allowing every physical pole prefix through depth three
does not repair the THM-3144 common-section death at degree 11.  Depth four is
different.  Two shallow four-pole states create a strictly positive common
selector through degree 11.  The resurrection is sharp in both directions:
depth three is insufficient, and degree 12 excludes every law on the complete
depth-at-most-four state bank.

## 1. The two-parameter selector spaces

Fix support `(a,b)=(1,3)`, invariant `I=1`, and bank `I2`.  Its reduced pole
multiset is

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                    (1)
```

For `d>=1`, let `S_<=d` be the multiplicity-valid unordered submultisets of
`P` of sizes one through `d`.  For `sigma in S_<=d`, use the fixed-`Q` virtual
prefix current

```text
G_N^sigma(mu)
 =Phi^sigma(h_N)m_mu[Q^sigma]
  -Phi^sigma(m_mu)h_N[Q^sigma]                            (2)
```

from THM-3137, THM-3144, and THM-3149.  If `lambda` is a probability law on
`S_<=d`, write `G_N(lambda)=sum_sigma lambda_sigma G_N^sigma`.

Define the cumulative selector space

```text
C_D^(<=d)={lambda on S_<=d:
           G_N(lambda) is a nonnegative Hasse boundary for 5<=N<=D}. (3)
```

The state counts through depth four are

```text
|S_1|=8, |S_2|=33, |S_3|=93, |S_4|=200,
|S_<=4|=334.                                               (4)
```

THM-3149 gives the first side of the barcode:

```text
C_11^(<=3)=empty.                                          (5)
```

## 2. An explicit strict depth-four resurrection

Consider the law

```text
lambda_(1)       =2291/5000,
lambda_(7)       =  36/5000,
lambda_(1,1,1,1) =1999/5000,
lambda_(1,1,1,2) = 674/5000.                              (6)
```

Its four weights are positive and sum to one.  The last two states are legal:
`P` contains four copies of `1` and three copies of `2`.

The companion enumerates every coarsening upset in degrees 5 through 11 and
evaluates its exact mass under `(6)`.  The only zero rows are the empty and
full upsets, whose masses vanish because every current `(2)` has total mass
zero.  Every nontrivial upset is strictly positive.  The census and exact
minimum are:

| `N` | nontrivial upsets | minimum mass | unique minimizing upset |
|---:|---:|---:|:---|
| 5 | 8 | `11134658556/125` | `{(5)}` |
| 6 | 25 | `3281855316351/125` | `{(6)}` |
| 7 | 45 | `268382847153552/125` | `P_7\{(1^7)}` |
| 8 | 166 | `11849576869236/625` | `P_8\{(1^8)}` |
| 9 | 571 | `23788552125790506/625` | `P_9\{(1^9)}` |
| 10 | 3,586 | `42037071115846938/625` | `P_10\{(1^10)}` |
| 11 | 19,540 | `2271700518167474796/625` | `P_11\{(1^11)}` |

Thus all `23,941` nontrivial upset inequalities are strict.  THM-3127, or
the companion's independent exact max-flow replay in every degree, gives

```text
C_11^(<=4) is nonempty.                                    (7)
```

The law `(6)` was discovered numerically, but neither its statement nor its
proof uses a floating approximation: all displayed checks use exact rational
arithmetic.

## 3. The law itself fails at degree twelve

For orientation, the same law is not a candidate in degree 12.  Its singleton
partition coefficient is

```text
G_12(lambda)(1^12)=872326118251778482488/625>0.             (8)
```

Since the current has total mass zero, the coarsening upset
`P_12\{(1^12)}` has negative mass.  The exact max-flow diagnostic also gives

```text
flow   =2299434408100109034608182728/625,
demand =2299561656443457484826636016/625,
deficit=127248343348450218453288/625
at (9,1,1,1).                                              (9)
```

This only kills `(6)`.  Excluding *every* other depth-four law needs the next
certificate.

## 4. A five-row death certificate for all 334 states

Use the four THM-3149 upsets

```text
U_8       =P_8 \ {(1^8)},
U_10      =P_10 \ {(1^10)},
U_11^(3)  =P_11 \ {(1^11),(2,1^9),(2,2,1^7)},
U_11^(1)  =P_11 \ {(1^11)},                               (10)
```

and add

```text
U_12^(1)=P_12 \ {(1^12)}.                                  (11)
```

Let `R_8,R_10,R_11^(3),R_11^(1),R_12^(1)` be their response
rows on the 334 states.  Take the primitive positive coefficient vector

```text
c_8  =2154825177301047232015910202093994919810226072626362713968741,
c_10 =  23533987357602417932067004154062702538970243225880783736334,
c_3  =      4845881049964897990450105824662373913248601773502462674,
c_1  =    517214042158781933902782482597749952061536320507147656820,
c_12 =      5015113163945262699318503892659649582247080863625246273.
                                                                    (12)
```

The exact combined row

```text
H_4=c_8R_8+c_10R_10+c_3R_11^(3)+c_1R_11^(1)+c_12R_12^(1) (13)
```

is strictly negative on every state in `S_<=4`.  Its range is

```text
-1026703906405262770537884068784375228334079877724481598818095786550552069856
 <=H_4(sigma)<=
-1281646212982665289100583072758234749313467771458224387568535766509728384.
                                                                    (14)
```

The minimum is unique at `(4,5,5,7)`.  The upper equality set is exactly

```text
(1), (3), (1,1,1,1), (1,1,1,2), (2,2,2,3).               (15)
```

Every row in `(13)` is a necessary Hasse upset inequality.  If a probability
law made all degrees 5 through 12 Hasse-positive, averaging `(13)` against it
would be both nonnegative and strictly negative.  Hence

```text
C_12^(<=4)=empty.                                          (16)
```

## 5. The exact birth/death cell

Combining `(5)`, `(7)`, and `(16)` gives

```text
C_11^(<=3)=empty,
C_11^(<=4) is nonempty,
C_12^(<=4)=empty.                                          (17)
```

Therefore depth four is the minimal physical prefix cap that solves the
degree-11 cumulative problem, and degree 11 is the largest cumulative horizon
solvable with prefixes of depth at most four.  By zero-extending `(6)`,
`C_11^(<=d)` remains nonempty for every `d>=4`.  By restriction in degree,
`C_D^(<=4)` remains empty for every `D>=12`.

All five dual rows factor through the THM-3147 length-singleton endpoint jet:
the first four do so by THM-3149, and `U_12^(1)={ell<=11}` is another pure
rank tail.  Thus the exact barcode `(17)` is already visible in that compressed
dual observer, although the positive law `(6)` was checked against the full
partition-lattice cone.

## 6. Exact verification

Run

```text
python 04-computation/gmc_depth_four_selector_resurrection_barcode_scout.py
python -O 04-computation/gmc_depth_four_selector_resurrection_barcode_scout.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_depth_four_selector_resurrection_barcode_scout.out.
```

The companion uses integer and rational arithmetic only.  It reconstructs the
334 physical states, all 23,941 nontrivial upsets through degree 11, every
minimum and max-flow equality, the degree-12 hostile `(8)--(9)`, and all 334
coordinates of `(13)`.  It checks positivity and primitivity of `(12)`, the
exact range `(14)`, and all equality states.  The immutable LF-normalized
hashes are recorded in the frontmatter.

## 7. Scope and next boundary

The theorem concerns probability averages of the derived fixed-`Q` virtual
prefix currents `(2)`.  A previous version said too strongly that an arbitrary
law on these unordered states need not be a sequential stopping distribution.
In the unrestricted finite-state sense every law has such a realization:
sample its terminal prefix, uniformly order the selected labelled poles, and
use the posterior transition kernel.  THM-3163 records the exact construction.
What this theorem does **not** supply is a value-only, prescribed-hazard, or
response-compatible pole process, nor a decomposition of the original
product-Gamma response.  That stronger compatibility is the load-bearing
sidecar; bare terminal-law realizability is automatic.

The result does not determine `C_12^(<=5)`, `C_13^(<=5)`, or any arbitrary
degree/depth limit.  It proves no cross-support portability, arbitrary-radial
NC2, or Gaussian Moment Conjecture.  The next exact question is whether the
two-dimensional barcode continues after adjoining depth-five states; that is
a new theorem, not a consequence of `(17)`.

QED.
