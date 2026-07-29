---
id: THM-2895
title: "Singleton-complement parity descent and four-root j6 closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.  The
  singleton-complement flag lemma descends cover size by two, with grammars
  5->3->1 and 6->4->2 and a sharp common-denominator wall at seven slots.
  Its exact marked-suffix application closes four seven-body/six-slot roots
  through 25 branches, 784 H4-pair residuals, and five recursive H2 edges.
source: root-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2893-complement-cap-finite-core-flag-lemma
related:
  - THM-2894-unmarked-residual-semilattice-order-and-group-clutch-no-go
  - THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas
verification:
  - 04-computation/lrc14_j6_h4_pair_residual_exact_kernel_codex_20260729.py
  - 04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895.py
  - 05-knowledge/results/lrc14_j6_suffix_parity_flag_closure_thm2895.out
  - 04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.py
  - 05-knowledge/results/lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.out
---

# THM-2895 -- singleton-complement parity descent and four-root j6 closure

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.**

## 1. Abstract theorem: the parity descent

Let `C` be a nonempty finite union of rational intervals, with

```text
h=|C|,                 c(w)=|C intersect D_w|,
D_w={t in R/Z : ||wt||<=1/14},
```

and let `V` be a cofinite allowed set of positive integer labels.  Fix
`3<=p<=6`, and suppose `B_1` is a global singleton cap

```text
c(w)<=B_1                    for every w in V.             (1)
```

If

```text
B_1 < (8-p)h/7,                                           (2)
```

put

```text
theta=h-B_1,
H_p={w in V : c(w)>=theta/(p-1)}.                         (3)
```

Then `H_p` is finite, every hypothetical `p`-cover contains at least two
vertices of `H_p`, and:

- if `p>=4`, it is enough to prove that the literal residual behind every
  unordered pair in `H_p` has no `(p-2)`-cover;
- if `p=3`, it is enough to check only those pairs in `H_p` whose union
  coverage is at least `theta`, and prove that no remaining singleton covers
  their literal residual.

Thus repeated fresh applications have the two parity grammars

```text
5 -> 3 -> 1,                  6 -> 4 -> 2.                (4)
```

Each arrow requires a newly reconstructed nonempty residual, allowed suffix,
singleton cap, and finite core.  Formula `(4)` is a recursion grammar, not a
claim that one cap propagates automatically.

### Proof

Apply THM-2893 with

```text
(k,s,ell)=(p,p-1,2).
```

The complement cap is exactly `(1)`.  Its finite-core inequality is

```text
B_1 < (7-s)h/7 = (8-p)h/7,
```

which is `(2)`.  THM-2893 therefore gives

```text
|K intersect H_p| >= k-s+1=2
```

for every hypothetical `p`-cover `K`.  When `p>=4`, `ell=2<s=p-1`, so the
heavy-flag condition on a selected core pair is vacuous.  The literal
residual has `p-2` labels left.  When `p=3`, `ell=s=2`, so the selected pair
itself must be heavy.  This proves the two clauses and `(4)`.

## 2. Why the ladder stops sharply at seven slots

The strict boundary in `(2)` is structural.  Choose a common denominator
`D` for all endpoints of `C`.  If `w` is a multiple of `D`, each component
of `C` is a union of complete periods of the `w`-comb.  Exactly one seventh
of every period belongs to `D_w`, up to endpoints, and hence

```text
c(w)=h/7.                                                (5)
```

A finite prefix exclusion leaves infinitely many allowed multiples of `D`.
At `p=7`, condition `(2)` would demand `B_1<h/7`, contradicting `(5)`.
At equality `B_1=h/7`, the high threshold in `(3)` is also `h/7`, so all
those multiples lie in the high core: it is genuinely infinite.  This is
the exact singleton-complement wall between the six-slot rung and the
seven-slot sector, not an artifact of a chosen scan horizon.

The two labels extracted at every live step explain the parity in `(4)`.
They do not define a tournament.  The pair relation is symmetric, and the
theorem-bearing state is the marked suffix together with the literal
residual, as required by THM-2894.

## 3. Four-root conclusion

Let `E` be any of

```text
(2,8,9,10,11,13,14),
(1,3,9,10,11,12,14),
(2,5,9,11,12,13,14),
(2,3,4,5,6,7,8).                                      (6)
```

For every six-element set `Q subset {15,16,...}`, the exact conclusion is

```text
|G_(E union Q)|>0.                                      (7)
```

These are four of the `3432` seven-body roots.  The first three are the
rank-one hostile cases from the earlier complement-cap battery; the fourth
is a low-stratum control.

## 4. Ordered suffix reduction

The globally sealed root gates have least sizes

```text
(19,20,21,13),
```

giving `73` unique-earliest-apex branches.  THM-2893's marked suffix
refinement and a global residual top-five scalar test close `48` branches.
Exactly `25` marked branches remain.

On every remaining branch, the exact singleton maximum satisfies

```text
B_1<3h/7.
```

The `p=5` instance of Section 1 therefore applies.  Its `H_4` cores have
sizes between `3` and `16`; their complete unordered-pair universe has

```text
784
```

literal residuals.

For a residual `R`, let `q_1>=q_2>=q_3` be its globally sealed singleton
coverages and let `B_2` be its exact globally sealed pair-union cap.  The
two sufficient three-cover obstructions

```text
q_1+q_2+q_3 < |R|,             B_2+q_1 < |R|             (8)
```

are incomparable on this universe.  The first closes `771/784`, the second
`773/784`, and their adaptive union closes `779/784`.  The computation pays
`7551` finite pair-union evaluations.

The rank stratification is informative:

```text
rank one:      pairs 361; top3 348; B2+q1 350; union 356; failures 5;
later ranks:   pairs 423; top3 423; B2+q1 423; union 423; failures 0.  (9)
```

Thus each cheap certificate separately closes every non-rank-one row.  The
only recursive work is concentrated in the two hostile rank-one carriers.

## 5. Recursive closure of the five hard residuals

The five survivors in `(8)` occur only on the rank-one branches:

```text
E=(2,8,9,10,11,13,14), a=19:
  H_4-pairs (37,108), (37,125), (108,125);

E=(2,5,9,11,12,13,14), a=16:
  H_4-pairs (23,25), (25,34).                           (10)
```

For each residual in `(10)`, a fresh cap satisfies `B_1<5|R|/7`.  The
`p=3` instance of Section 1 gives high cores and heavy edges

```text
H_2=(17,23,46), edge (17,23),        on the first three;
H_2=(20,37),    edge (20,37),        on (23,25);
H_2=(20,23),    edge (20,23),        on (25,34).         (11)
```

There is exactly one heavy edge in each core.  Behind those edges, the
longest residual-component lengths and geometric singleton horizons are

```text
ell=121/34776,  13/3500,  125/28728,  33/7252,  19/4508;
W  =41,         38,       32,         31,       33.      (12)
```

THM-2893's longest-component seal excludes every `w>W`.  Exact
noncontainment checks for all allowed `15<=w<=W` pay

```text
23+20+15+13+15=86
```

tests and find zero singleton covers.  Literal and direct reconstruction
agree on all five terminal residuals.  Hence all five rows in `(10)`, all
`784` pair obligations, all `25` nonscalar suffix branches, and finally all
four roots in `(6)` close.

## 6. Verification ledger

The locked verifier, exact kernel, and transcript are

```text
04-computation/lrc14_j6_h4_pair_residual_exact_kernel_codex_20260729.py
SHA-256 b82f318bf89ffd3ab4c918c87736461d068e03f25941aa25a0961d0f74b4d70a

04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895.py
SHA-256 970d77503f8d56d737e223dabb3c3562d7b19cd018ca75398e3deb054715e5f6

05-knowledge/results/lrc14_j6_suffix_parity_flag_closure_thm2895.out
SHA-256 c11260f6544a319e1cc1862c9221b188a4314860422470e465b82e7ce492b1b4
```

The primary verifier pins its inherited exact kernels, the ranked-suffix
source and transcript, every aggregate above, the five rows in
`(10)`--`(12)`, and a canonical ledger over all `25` branches and `784`
residual caps.  Ordinary and optimized replays are byte-identical.  The
canonical ledger digest is

```text
c541e41f688dd9873e57b5316a8c9b28c496e9bd808a65c62f7d20b7a3b87d4e.
```

A second implementation imports no repository LRC module.  It reconstructs
the interval unions, periodic antiderivative, all `73` suffix branches, all
`784` pair residuals, and all five recursive rows:

```text
04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.py
SHA-256 9023a4042dc8def3f8781668e721549972fb1458d07f2fab89bf7ac3a6f745cc

05-knowledge/results/lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.out
SHA-256 e559bef5f6a3c2934eaac463985e641705b92eef4ba22c0f0b6e9a5e66c0a63d
```

Its ordinary and optimized replays are also byte-identical, and its
independent ledger digest is

```text
eae3f0a54ea86dfb61a63e7340f5199689bbf8714f63da28c9276f289611bac4.
```

## 7. Exact scope

The abstract parity descent is general under `(1)`--`(3)`.  The finite-exact
application proves only the four roots in `(6)`.  It does not prove the full
seven-body/six-slot rung and does not prove unrestricted LRC(14).  THM-2896
supplies the all-root first-apex gates but no apex-carrier closure.
