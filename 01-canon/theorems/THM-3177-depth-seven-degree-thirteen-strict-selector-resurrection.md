---
id: THM-3177
title: "Depth-seven degree-thirteen strict selector resurrection"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For support (1,3), bank I2, an explicit primitive thirteen-state law of
  denominator 100,000,000 is strictly positive on every nontrivial
  partition-coarsening upset in every degree 5 through 13.  Hence the
  cumulative degree-13 selector cone, empty through physical prefix depth
  six, resurrects at depth seven.  The exact verifier rebuilds all 1,820
  states from two hash-pinned transitive sources and uses an independent
  partition generator, bin-packing coarsening test, and integer Dinic cut.
  This is an averaged virtual-prefix theorem, not a response-compatible
  stopping process, original-response decomposition, or GMC extension.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
audit: >
  Fresh normal and optimized immutable replays match the stored 18-line
  transcript byte-for-byte after LF normalization.  A separately written
  exact auditor used the pinned response vectors, an independent partition
  generator and bin-packing coarsening test, and NetworkX integer min-cut; it
  reproduced all nine strict minima, cut sizes, and complete generator
  antichains.  Static hostile audit also accepted the nontrivial-upset cut,
  infinity capacity, posterior reach/lumping, transitive hashes, scope, and
  time-only positive-hazard no-go.  No theorem defect was found.
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3158-depth-five-selector-resurrection-through-degree-twelve
  - THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary
  - THM-3169-depth-six-degree-thirteen-nonresurrection-and-mutated-cocycle-certificate
related:
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
script: 04-computation/gmc_depth_seven_selector_resurrection_thm3177.py
output: 05-knowledge/results/gmc_depth_seven_selector_resurrection_thm3177.out
pole_prefix_dependency_sha256: 151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7
signed_bank_dependency_sha256: 15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f
script_sha256: b10529ca9c8c1d3a248bd7c8cd89ef0a8a67d46584d0f3d70bb2fa9781da137f
output_sha256: c893600be19715abc3ade31b7471bb0986c0417d81378b7f85893418c18a9134
hash_basis: LF-normalized bytes
---

# THM-3177 -- depth-seven degree-thirteen strict selector resurrection

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3158 proves that the support-`(1,3)`, bank-`I2` selector resurrects
through degree 12 at physical prefix depth five.  THM-3169 then shows that
degree 13 remains impossible after adjoining every depth-six state, even
though 24 of those new states cross the old depth-five separating wall.
The next layer changes the primal status: depth seven admits a law which is
strict on every nontrivial Hasse-dual upset through degree 13.

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

The physical depth layers through seven have sizes

```text
(8,33,93,200,348,507,631),              |S_<=7|=1820.     (4)
```

THM-3169 supplies the inherited obstruction

```text
C_13^(<=6)=empty.                                         (5)
```

## 2. A primitive thirteen-state law

Let `lambda_sigma=n_sigma/100000000`, with the only nonzero numerators

| state `sigma` | `n_sigma` |
|:---|---:|
| `(1)` | 329,341 |
| `(2)` | 506,132 |
| `(1,1,2)` | 353,515 |
| `(1,1,1,1)` | 589,431 |
| `(5,5,6,7,8)` | 13,469,020 |
| `(1,1,1,1,2,2)` | 448,483 |
| `(1,1,1,1,2,2,2)` | 510,392 |
| `(1,1,2,2,2,3,3)` | 446,099 |
| `(1,2,2,2,3,3,4)` | 313,313 |
| `(1,4,5,5,6,7,8)` | 12,884,975 |
| `(2,3,3,5,6,7,8)` | 17,513,294 |
| `(3,3,4,4,6,7,8)` | 2,642,438 |
| `(4,4,5,5,6,7,8)` | 49,993,567 |

Every displayed state is a legal physical submultiset of `(1)`.  The
numerators are positive, primitive, and sum to `100000000`.

## 3. Exact minimization over every upset

For each `5<=N<=13`, the companion first verifies zero total mass for every
one of the `1,820` state responses.  Empty and full upsets therefore have
mass zero.  Every nonempty proper upset in the partition-coarsening lattice
contains the top `(N)` and excludes the bottom `(1^N)`.

The verifier independently generates the partition lattice, tests refinement
by exact bin-packing, and agrees with the hash-pinned coefficient source on
all `361` partitions and all `22,809` ordered coarsening pairs.  It then
minimizes the raw numerator weight over all nontrivial upsets by an
integer-capacity Dinic cut.  Closure edges have capacity one larger than the
total finite capacity, while source/top and bottom/sink edges force the two
anchors.  Thus a positive minimum proves every nontrivial upset inequality,
without enumerating antichains or trusting a floating oracle.

The exact minima are as follows.  The displayed mass is the raw numerator;
the normalized Hasse mass is that integer divided by `100000000`.  The last
column lists the complete minimal-generator antichain of the minimizing
upset returned by the deterministic exact cut.

| `N` | raw minimum | upset size | minimal generators |
|---:|---:|---:|:---|
| 5 | 109,628,640 | 1 | `(5)` |
| 6 | 72,639,032,313,522,840 | 1 | `(6)` |
| 7 | 2,990,897,801,796,238,560 | 14 | `(2,1^5)` |
| 8 | 318,185,095,107,504 | 21 | `(2,1^6)` |
| 9 | 9,346,991,684,459,633,328 | 29 | `(2,1^7)` |
| 10 | 3,288,235,113,012,380,640 | 41 | `(2,1^8)` |
| 11 | 2,177,378,045,882,508,776,256 | 55 | `(2,1^9)` |
| 12 | 20,818,394,379,371,373,887,136 | 75 | `(3,1^9)`, `(2,2,1^8)` |
| 13 | 1,303,441,821,531,365,341,151,616 | 96 | `(3,1^10)`, `(2,2,2,2,2,1,1,1)` |

All nine minima are strictly positive.  THM-3127's finite
Strassen--Hasse duality therefore gives

```text
lambda in C_13^(<=7),             so C_13^(<=7) is nonempty. (6)
```

The floating column-generation search which discovered the thirteen states
is not evidence for `(6)`.  The proof is the displayed rational law together
with the exact all-upset minima above.

## 4. The sharpened death--persistence--resurrection barcode

Combining `(5)--(6)` with THM-3158 gives

```text
C_12^(<=5) is nonempty,
C_13^(<=5) is empty,
C_13^(<=6) is empty,
C_13^(<=7) is nonempty.                                  (7)
```

Thus seven is the minimal physical prefix depth which solves the cumulative
degree-13 selector problem.  The transition is not monotone wall crossing:
depth six crosses the depth-five wall but remains separated by a mutated
eleven-facet cocycle, while depth seven crosses the full degree-13 cone.
Zero extension makes `C_13^(<=d)` nonempty for every `d>=7`, and degree
restriction gives nonemptiness for every horizon `D<=13` at depth seven.

This is a concrete holotopy lesson.  Resource layers can attach columns on
both sides of one projected wall without changing feasibility, and a later
layer can resurrect feasibility only after the dual observer has mutated.
The invariant object is the filtered cone and its changing upset witnesses,
not one persistent scalar separator.

## 5. Abstract sequential realization is automatic, not sufficient

THM-3163 applies to every law on a finite prefix bank, hence to the law in
Section 2.  Its exchangeable labelled lift has, in the displayed state order,

```text
(4,3,18,1,1,3,1,6,8,8,6,1,1)                           (8)
```

terminal lifts: `61` labelled terminal atoms in total.  The exact posterior
chain reaches `1,620` labelled states with `6,692` positive directed edges.
Strong lumping by equal pole values gives `289` multiplicity states and `908`
positive edges.  The companion checks every stochastic identity and
propagates the terminal law exactly.

This is only THM-3163's **abstract** sequential realization.  Its transitions
use the full future terminal law.  They are not value-only, prescribed-hazard,
selector-current-compatible, or compatible with a positive decomposition of
the original product-Gamma response.  The Markov chain therefore supplies no
missing response sidecar and no GMC consequence.

There is also a sharp elementary obstruction to the simplest prescribed-
hazard model.  Suppose stopping depends only on time and every available pole
value has a positive value/time removal hazard.  The positive singleton stops
at `(1)` and `(2)` force the depth-one stopping probability to be positive.
Values `3,4,5,6,7,8` occur in positive terminal states, so their positive
initial hazards would then force positive singleton terminal masses at each of
those values.  The displayed law assigns all six masses zero, a contradiction.
This excludes only the stated time-only-stop/positive-hazard model; it does not
exclude THM-3163's state-dependent posterior construction.

## 6. Exact verification

Run

```text
python 04-computation/gmc_depth_seven_selector_resurrection_thm3177.py
python -O 04-computation/gmc_depth_seven_selector_resurrection_thm3177.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_depth_seven_selector_resurrection_thm3177.out.
```

The companion uses integer and `Fraction` arithmetic only.  It enforces and
prints both transitive LF hashes, reconstructs all `1,820` states, checks all
`16,380` state-degree balances, independently audits the partition order,
performs every exact minimum cut, and verifies the posterior-chain census.
The two dependency hashes and the companion/output hashes are recorded in the
frontmatter.

## 7. Scope and next boundary

The theorem concerns probability averages of derived fixed-`Q` virtual
prefix currents.  It proves no response-compatible pole-local transition,
no decomposition of the original response, no arbitrary-radial NC2 theorem,
no extension of GMC, and no LRC consequence.

It also makes **no degree-14 claim**.  A discovery-only floating solve comes
arbitrarily close to a degree-14 boundary, but its rationalization violates a
named degree-14 upset.  That is unresolved boundary evidence, not an
infeasibility certificate.  The next exact question is whether
`C_14^(<=7)` contains a strict rational law, only a boundary law on named zero
facets, or admits a positive rational Farkas separator.

QED.
