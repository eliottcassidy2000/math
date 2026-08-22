---
id: THM-3158
title: "Sharp depth-five selector resurrection and degree-thirteen death barcode"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  At support (1,3), bank I2, the cumulative selector space through degree 12
  is empty for all physical pole prefixes of depth at most four but becomes
  nonempty at depth five.  An explicit seven-state law of denominator one
  million is strictly positive on all 403,539 nontrivial coarsening upsets in
  degrees 5 through 12.  Adding degree 13 kills the entire 682-state
  depth-at-most-five bank: nine upset rows have a primitive positive integer
  combination strictly negative on every state.  This is an averaged virtual-
  prefix theorem, not a sequential stopping process or original-response
  decomposition.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
audit: >
  An independent immutable audit reconstructed the 682-state census and the
  denominator-one-million law, replayed normal and optimized companions
  byte-for-byte, and checked all 403,539 strict antichain/upset inequalities.
  It separately verified that increasing incomparable generators enumerate
  each upset once, that R5--R8 are lawful upsets with the displayed minimal
  antichains, and that the nine coefficients are primitive and positive.  All
  682 separator coordinates, the exact range and equality states, the N13
  hostile, LF hashes, barcode monotonicity, and fixed-Q averaged-current scope
  passed.
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3147-length-singleton-endpoint-jet-facet-observer
  - THM-3155-sharp-depth-four-selector-resurrection-through-degree-eleven
related:
  - THM-3144-mixed-depth-selector-persistence-death-barcode
  - THM-3149-depth-three-selector-persistence-and-cross-support-wall
script: 04-computation/gmc_depth_five_selector_resurrection_barcode_scout.py
output: 05-knowledge/results/gmc_depth_five_selector_resurrection_barcode_scout.out
script_sha256: 9aa84269691b7e84d193914180681bcc2d6551d478242cd06e51cb97f7005d2c
output_sha256: 0b1955f6d34ac88da2be3fccaab38e6da89d06a03577005fc68e92c19ee2fe1e
hash_basis: LF-normalized bytes
---

# THM-3158 -- sharp depth-five selector resurrection and degree-thirteen death barcode

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3155 proves that depth four resurrects the support-`(1,3)`, bank-`I2`
selector through degree 11, but that degree 12 kills every law on the complete
depth-at-most-four state bank.  Depth five repeats the phenomenon one step
later.  A sparse seven-state law survives every degree through 12, whereas an
exact nine-facet certificate kills every depth-five law at degree 13.

## 1. The cumulative selector spaces

Fix invariant `I=1`, support `(a,b)=(1,3)`, and bank `I2`.  The reduced pole
multiset is

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                    (1)
```

For `d>=1`, let `S_<=d` be the multiplicity-valid unordered nonempty
submultisets of `P` of sizes at most `d`.  For `sigma in S_<=d`, use the
fixed-`Q` virtual prefix current

```text
G_N^sigma(mu)
 =Phi^sigma(h_N)m_mu[Q^sigma]
  -Phi^sigma(m_mu)h_N[Q^sigma].                           (2)
```

For a probability law `lambda` on `S_<=d`, put
`G_N(lambda)=sum_sigma lambda_sigma G_N^sigma` and define

```text
C_D^(<=d)={lambda: G_N(lambda) is a nonnegative Hasse boundary
                   for every 5<=N<=D}.                   (3)
```

The depth layers have sizes

```text
(|S_1|,|S_2|,|S_3|,|S_4|,|S_5|)=(8,33,93,200,348),
|S_<=5|=682.                                              (4)
```

The inherited first side of the next barcode cell is

```text
C_12^(<=4)=empty                                          (5)
```

by THM-3155.

## 2. A seven-state strict resurrection through degree twelve

Let `lambda_sigma=n_sigma/10^6`, with the only nonzero numerators

| state `sigma` | `n_sigma` |
|:---|---:|
| `(1)` | 97,856 |
| `(1,1,1)` | 56,951 |
| `(1,1,1,1)` | 140,643 |
| `(1,1,2,2,2)` | 194,498 |
| `(1,2,2,3,3)` | 7,398 |
| `(2,2,2,3,3)` | 118,572 |
| `(5,5,6,7,8)` | 384,082 |

The numerators are positive, primitive, and sum to `10^6`; every displayed
state is a legal submultiset of the pole multiset in `(1)`.

The companion streams every antichain of every partition-coarsening poset in
degrees 5 through 12.  By the antichain/upset bijection, this is an exhaustive
check of every Hasse-dual inequality without materializing the degree-12 bank.
Only the empty and full upsets have zero mass.  Every nontrivial upset has
strictly positive mass.  The exact census and minima are:

| `N` | all upsets | nontrivial | minimum numerator over `10^6` | unique minimizer |
|---:|---:|---:|---:|:---|
| 5 | 10 | 8 | 53,951,073,740,640 | `{(5)}` |
| 6 | 27 | 25 | 14,713,182,072,397,560 | `{(6)}` |
| 7 | 47 | 45 | 447,905,708,691,282,144 | `P_7\{(1^7)}` |
| 8 | 168 | 166 | 175,399,501,380,224,640 | `P_8\{(1^8)}` |
| 9 | 573 | 571 | 2,519,721,209,269,962,000 | `P_9\{(1^9)}` |
| 10 | 3,588 | 3,586 | 4,663,094,943,523,464,240 | `P_10\{(1^10)}` |
| 11 | 19,542 | 19,540 | 303,959,323,177,696,749,696 | `P_11\{(1^11)}` |
| 12 | 379,600 | 379,598 | 29,594,579,030,356,898,846,112 | `P_12\{(1^12)}` |

Thus all `403,539` nontrivial inequalities are strict.  The companion also
replays exact max-flow equality in every degree.  THM-3127 gives

```text
C_12^(<=5) is nonempty.                                   (6)
```

The numerical search that found the law is not evidence for `(6)`; the proof
uses only the displayed rational law and exhaustive exact arithmetic.

## 3. The displayed law fails at degree thirteen

For orientation, this same law is not a degree-13 solution.  Its singleton
coefficient is positive:

```text
G_13(lambda)(1^13)=299630262562066957334637/62500.         (7)
```

Equivalently the singleton-complement upset is negative.  The exact transport
diagnostic gives

```text
flow   =2890357753645205037202525717623/25000,
demand =3613065094746902954107599900417/31250,
deficit=471610761586630417771013553/125000
at (10,1,1,1).                                            (8)
```

This failure alone says nothing about other laws on `S_<=5`.  Their universal
exclusion is the content of the next certificate.

## 4. Nine exact upset facets kill all 682 states

Write `up(A)` for the upset generated by a set `A` of partitions.  Use the
nine necessary Hasse inequalities

```text
R0=P_8 \{(1^8)},
R1=P_10\{(1^10)},
R2=P_11\{(1^11)},
R3=P_12\{(1^12)},
R4=P_13\{(1^13)},

R5=up{(2,2,1^7)} in P_11,

R6=up{(4,2,1^7),(3,3,1^7),(3,2,2,1^6),(2,2,2,2,1^5)},
R7=up{(4,2,1^7),(3,3,3,1^4),(3,2,2,2,1^4),
      (2,2,2,2,2,1^3)},
R8=up{(4,2,1^7),(3,3,1^7),(3,2,2,1^6),
      (2,2,2,2,2,1^3)}                                   (9)
```

where `R6,R7,R8` lie in `P_13`.  Regard each `Ri` as its 682-coordinate
response row.  In the order `(R0,...,R8)`, take the primitive positive integer
vector

```text
79966346203432495238210892467836051822404864748631537035246352913301789115999332193913950905464025191074830103726371948258801337873186194518,
4385616886032821959435837177631104375265963397917293088514990249982327757237457317872391016149793462027134777492436937280742150390259271553,
323165012303885417971182559526221438425147324616172737787579961666908713240511160571617207173787027117769558974661863206105407364575577287,
7744440757659416200591308784108711415468181697333643558792830072051161414809397905536549772356510263760667669854566940871988150246747864,
64198555842605394277362557914437771970209622650006876621643073243158046817590826310245420481922857881786832819412787764831035747707567,
1238919946181093630207021295990664443370069758642989075525265110060651909301010930905141007031425295949856900214021045838758589520118799,
57173834483208673247111573734889531022707614332675275902916730146155429982624862452119240877347978257760901011795049436871545379049,
3220669545142804044646293595442679616378197077725758366899484010296430542402560592047500094580557804747642540554956255672701117951,
43322631946675983735399838865716397620757342948083795752004906378365718688259231686082777317064579418881565313974450622142635728100.
                                                                    (10)
```

The resulting row `H_5=sum_i c_i R_i` is strictly negative on every one of
the 682 states.  Its exact range is

```text
-131646038360152726651793725918294945734506800899964245935473472676004061252003844475399955296629305866663476692762840778970142279919601748623581179536276608
 <=H_5(sigma)<=
-7905909118176927518212656514112818602784521715342679413460937962315990351772853589181871947048820698679626976418226870569912439470329653112885762708992.
                                                                    (11)
```

The minimum is unique at `(8)`.  The upper equality set is exactly

```text
(1), (1,1,2), (1,1,1,1), (1,1,2,2,2),
(1,2,2,3,3), (2,2,2,3,3),
(4,5,5,6,7), (4,5,6,7,8), (5,5,6,7,8).                  (12)
```

If a probability law satisfied all degree-5-through-13 Hasse inequalities,
its average against `(10)` would be both nonnegative and strictly negative.
Therefore

```text
C_13^(<=5)=empty.                                         (13)
```

The separation oracle discovered `(9)`, but the proof of `(13)` is the direct
exact 682-coordinate verification of `(10)--(12)`.

## 5. The next exact barcode cell

Combining `(5)`, `(6)`, and `(13)` gives

```text
C_12^(<=4)=empty,
C_12^(<=5) is nonempty,
C_13^(<=5)=empty.                                         (14)
```

Hence depth five is the minimal physical prefix cap that solves the cumulative
degree-12 problem, and degree 12 is the largest cumulative horizon solvable at
depth at most five.  Zero extension makes `C_12^(<=d)` nonempty for every
`d>=5`; degree restriction makes `C_D^(<=5)` empty for every `D>=13`.

The new dual carrier is still tiny compared with the 682-state primal bank,
but unlike the earlier barcode cells it includes three genuinely noncentral
degree-13 upsets.  The wall is therefore no longer visible from singleton
complements alone.  This is the first exact signal that the compressed dual
observer must grow as the staircase advances.

## 6. Exact verification

Run

```text
python 04-computation/gmc_depth_five_selector_resurrection_barcode_scout.py
python -O 04-computation/gmc_depth_five_selector_resurrection_barcode_scout.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_depth_five_selector_resurrection_barcode_scout.out.
```

The companion uses integer and `Fraction` arithmetic only.  It reconstructs
all 682 states, streams all 403,539 nontrivial upsets through degree 12,
checks every strict minimum and max-flow equality, verifies the degree-13
hostile `(7)--(8)`, and evaluates the nine-row separator on every state.  It
also checks legality, positivity, normalization, primitivity, the exact range,
and every equality state.  LF-normalized hashes are in the frontmatter.

## 7. Scope and next boundary

The theorem concerns probability averages of the derived fixed-`Q` virtual
prefix currents `(2)`.  Multiplicity-valid unordered states do not by
themselves define a sequential pole-removal process, and no such stopping
realization is proved here.  Nor is there yet a decomposition of the original
product-Gamma response through these selector currents.

The theorem does not determine `C_13^(<=6)`, the eventual full-bank behavior,
arbitrary-radial NC2, or the Gaussian Moment Conjecture.  The immediate finite
question is whether depth-six states cross the nine-facet wall; the conceptual
question is whether the successive sparse dual carriers arise from a single
pole-subtraction recurrence rather than unrelated LP certificates.

QED.
