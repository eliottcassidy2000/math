---
id: THM-3163
title: "Universal finite-prefix Markov realization and physical sidecar boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT / INDEPENDENT AUDIT REQUIRED.
  Every probability law on the subsets of a finite labelled pole multiset is
  the terminal law of an explicit state-dependent Markov chain that removes
  one remaining label at a time.  A law on unlabeled submultisets has an
  exchangeable labelled lift, so the chain descends to multiplicity states.
  In particular THM-3158's seven-state law has an exact 29-terminal labelled
  realization.  Bare sequential realizability is therefore automatic and is
  not the missing GMC sidecar; value-only, response-compatible, or original-
  current-compatible transport remains unproved.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
depends_on: []
related:
  - MISTAKE-354
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3155-sharp-depth-four-selector-resurrection-through-degree-eleven
  - THM-3158-depth-five-selector-resurrection-through-degree-twelve
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
script: 04-computation/gmc_finite_prefix_markov_realization_thm3163.py
output: 05-knowledge/results/gmc_finite_prefix_markov_realization_thm3163.out
script_sha256: 35ee2d40987768c23f755fc53201cb13e794ad77be52e055be565e0b6e3832db
output_sha256: f58aa203072330c18ccec25be7448fa15ad4f8011d3fcee8da3e1445b66fc4cf
hash_basis: LF-normalized bytes
---

# THM-3163 -- universal finite-prefix Markov realization and physical sidecar boundary

**PROVED CANDIDATE + VERIFIED-EXACT / INDEPENDENT AUDIT REQUIRED.**

The product-Gamma selector staircase produces probability laws on finite
banks of legal pole-prefix states.  It is tempting to call the absence of a
sequential stopping process the remaining difficulty.  In the unrestricted
finite-state sense that difficulty does not exist: every terminal law has an
explicit one-letter-at-a-time Markov realization.  What remains substantive
is compatibility with the algebra being transported.

## 1. The posterior prefix kernel

Let `E` be a finite labelled set and let `pi` be any probability law on
`2^E`, possibly with mass at the empty set.  Think of `S subset E` as the set
of labels already removed.  For `s=|S|`, define the reach weight

```text
r(S)=sum_(T superset S) pi(T)/binom(|T|,s).                 (1)
```

Whenever `r(S)>0`, put

```text
q(S)=pi(S)/r(S),                                           (2)

p(S,e)=1/r(S) sum_(T superset S union {e})
        pi(T)/[binom(|T|,s)(|T|-s)]       (e notin S).      (3)
```

At state `S`, stop with probability `q(S)` and otherwise move to
`S union {e}` with probability `p(S,e)`.  States with `r(S)=0` are
unreachable and may be completed arbitrarily.

## 2. Exact stochasticity and terminal law

All terms in `(1)--(3)` are nonnegative.  For each strict superset `T` of
`S`, its contribution to the sum over all possible next labels is

```text
sum_(e in T\S) pi(T)/[binom(|T|,s)(|T|-s)]
 =pi(T)/binom(|T|,s).                                      (4)
```

Adding the stop contribution `pi(S)` proves

```text
q(S)+sum_(e notin S)p(S,e)=1.                              (5)
```

The probability of reaching `S` is exactly `r(S)`.  This is clear
probabilistically: first sample `T` from `pi`, choose a uniformly random order
of `T`, and expose its initial segments.  Conditional on `T`, the set of the
first `s` labels is uniform among the `s`-subsets of `T`, giving `(1)`.

There is also a direct Markov induction.  If `|S|=s>0`, then the total arrival
mass from its immediate predecessors is

```text
sum_(e in S) r(S\{e}) p(S\{e},e)
 =sum_(T superset S) pi(T)
   sum_(e in S) 1/[binom(|T|,s-1)(|T|-s+1)]
 =sum_(T superset S) pi(T)/binom(|T|,s)
 =r(S),                                                     (6)
```

because

```text
binom(t,s-1)(t-s+1)=s binom(t,s).                           (7)
```

Therefore the probability of stopping at `S` is

```text
r(S)q(S)=pi(S).                                            (8)
```

This proves that every finite terminal law is realized by a time-homogeneous
state-dependent Markov chain on the prefix graph.

## 3. Passage to a physical pole multiset

Let a physical pole multiset contain `m_v` labelled copies of each value `v`.
For an unlabeled submultiset `sigma`, with multiplicities `k_v`, the number of
labelled lifts is

```text
L(sigma)=product_v binom(m_v,k_v).                          (9)
```

Lift an unlabeled law `lambda` by assigning mass

```text
pi(T)=lambda(sigma)/L(sigma)                                (10)
```

to every labelled lift `T` of `sigma`.  This law is invariant under the
product of the permutation groups of equal-valued labels.  Formulas
`(1)--(3)` are equivariant under the same group, so the labelled Markov chain
descends to a Markov chain on multiplicity-valid unordered prefix states.
Each step removes one remaining physical pole copy.

Thus **every** law in every finite selector simplex `Delta(S_<=d)` has a
formal sequential stopping realization.  No convex or Hasse-positivity
hypothesis is needed.

## 4. Exact realization of the THM-3158 law

For the support-`(1,3)`, bank-`I2` pole multiset

```text
(1,1,1,1,2,2,2,3,3,4,4,5,5,6,7,8),                       (11)
```

THM-3158's denominator-one-million law is supported on

```text
(1), (1,1,1), (1,1,1,1), (1,1,2,2,2),
(1,2,2,3,3), (2,2,2,3,3), (5,5,6,7,8).                   (12)
```

The seven states have respectively

```text
4,4,1,6,12,1,1                                             (13)
```

labelled lifts, for 29 labelled terminal states in total.  Applying
`(1)--(3)` reaches 232 labelled prefix states and has 728 positive directed
transitions.  The exact companion checks stochastic normalization at every
reachable state, propagates the chain independently, and recovers the 29
labelled terminal masses and all seven unlabeled masses exactly.

Consequently the THM-3158 law does have a sequential pole-removal realization
in this broad sense.  This adds no Hasse positivity beyond THM-3158 and no
response identity beyond the prescribed terminal distribution.

## 5. Why this does not solve the selector problem

The construction uses the entire terminal law in every posterior transition.
The probability of removing a label can depend on the complete current state
and on all future terminal supersets.  It need not be determined only by the
pole value and depth, need not arise from the scalar row polynomial, and need
not intertwine the original response with the fixed-`Q` virtual-prefix
currents.

The realization is also nonunique.  For the point mass at a two-label set,
either deterministic ordering gives the same terminal law and a different
transition kernel.  Terminal data therefore do not supply the missing
transport canonically.

This separates four notions that must not be conflated:

```text
finite terminal-law realizability                       automatic;
value-only or prescribed-hazard realization             additional;
selector-current-compatible depth transport             additional;
positive decomposition of the original response         additional. (14)
```

THM-3160 identifies the algebraic reason for the third gap: same-degree
selector currents lose cross-degree Pluecker minors.  The posterior chain
above carries no such minors and cannot repair that projection loss.

## 6. Exact verification

Run

```text
python 04-computation/gmc_finite_prefix_markov_realization_thm3163.py
python -O 04-computation/gmc_finite_prefix_markov_realization_thm3163.py
```

and compare with

```text
05-knowledge/results/gmc_finite_prefix_markov_realization_thm3163.out.
```

The companion uses `Fraction` arithmetic.  It tests `(1)--(8)` on a dense
four-label probability law, then constructs and independently propagates the
29-terminal labelled lift of the exact THM-3158 law.  The general proof is the
finite identity `(4)--(8)`; the companion supplies exact positive and
repeated-value controls.

## 7. Scope

This theorem proves only abstract finite Markov realization.  It does not
make a selector law independent-deletion, value-only, depth-hazard, or
response-compatible.  It gives no original-response decomposition,
arbitrary-radial NC2 theorem, or Gaussian Moment Conjecture consequence.

The corrected research question is not whether the sparse selector laws can
be stopping distributions.  They all can.  It is whether one can choose a
lawful transition kernel that simultaneously carries the original response,
the required cross-degree endpoint data, and Hasse positivity.

QED (candidate pending independent audit).
