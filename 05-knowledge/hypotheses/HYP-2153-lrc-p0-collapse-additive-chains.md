---
id: HYP-2153
status: S599 covering-depth master object plus exact S602 bounded audit; classification open
source: opus-2026-06-03-S599; codex-2026-06-03-S602
related:
  - HYP-1802
  - HYP-1810
  - HYP-2084
  - HYP-2135
  - HYP-2138
  - HYP-2146
  - HYP-2141
  - HYP-2114
  - HYP-2151
  - HYP-2152
  - THM-358
  - THM-401
---

# HYP-2153: the `p0` collapse family is an additive-chain shell family, not just AP

## Covering-Depth Definition

Opus-S599 formalizes the user's `p_0` notation as the zero-depth cell of the
LRC circular-arc covering problem.  For a speed set `V` with `m=|V|` and gap
`delta`, define

```text
depth_delta(t) = #{v in V : ||v t|| < delta},
p_k(delta)    = meas{t : depth_delta(t)=k}.
```

Then the lonely-measure is `p_0(delta)`.  At the reduced LRC threshold
`delta=1/(m+1)`, the worry-set is the critical stratum where `p_0=0` but
boundary lonely points still exist.

S602 uses the equivalent endpoint-critical predicate with `n=m+1`: for a
primitive row `V` with `|V|=n-1` and threshold `1/n`, call `V` `p0`-collapsed
when:

```text
1. the safe set has no open interval;
2. boundary lonely witnesses exist;
3. the lcm of boundary-witness denominators is n;
4. the boundary witnesses are exactly {a/n : gcd(a,n)=1}.
```

This is stronger than "there is a witness at `1/n`."  It says the whole
visible boundary quotient collapses to the unit skeleton of `Z/nZ`, while the
full endpoint-protection quotient may still live at `n*lcm(V)`.

## Claim

The `p0`-collapsed family is larger than the arithmetic progression family,
and the visible non-AP branch is additive:

```text
(1,3,4,7)      with 4=1+3 and 7=3+4,
(1,3,4,5,9)    with 4=1+3, 5=1+4, and 9=4+5.
```

The working classification problem is:

```text
p0-collapse = unit-boundary skeleton
              + two-seed additive-chain generation
              + THM-401 shell constraints.
```

Prime `C=2n-1` should leave only AP plus the small `{2}` transversal flips
visible in the bounded floor branch.  Composite `C` can add non-transversal
sporadics through nonunit shell ramification, as in the known `n=8` rows.

The precise converse remains open.  S599 verifies AP and the named chains, plus
a non-chain control `(1,2,4,7)` with `p_0=6/35>0`.  S602 finds that every
collapsed row in its targeted boxes is a two-seed addition chain, but also
records many two-seed chains that are not collapsed.  The missing ingredient is
the shell/endpoint filter.

## S599 Covering-Depth Evidence

`04-computation/lrc_covering_depth_distribution_s599.py` makes `p_0` literal:
it computes the full depth distribution `p_k`.  It verifies:

```text
E[depth] = sum_k k p_k = 2 m delta,
```

so the first moment is configuration-independent.  Collapse is therefore a
higher-moment/additive-overlap phenomenon, not a first-moment measure bound.

At the critical threshold, S599 records `p_0=0` for AP rows and for the named
chains `(1,3,4,7)` and `(1,3,4,5,9)`, while the non-chain control
`(1,2,4,7)` has `p_0=6/35>0`.

S599 also approaches the collapse from below and fits

```text
p_0(delta_c - epsilon) ~ epsilon^alpha.
```

For AP rows and `(1,3,4,7)`, the fitted exponent is `alpha=1`: a singleton-wall
or one-Helly-stage collapse, tying this subproblem to HYP-2146/HYP-2151/HYP-2152.

## S602 Evidence

`04-computation/lrc_p0_collapse_additive_chains_s602.py` reruns exact S356
interval arithmetic in targeted primitive boxes:

```text
k=3, max_speed=12
k=4, max_speed=16
k=5, max_speed=13
k=6, max_speed=11
k=7, max_speed=14
```

It finds `9` `p0`-collapsed primitive rows:

```text
n=4: (1,2,3)
n=5: (1,2,3,4), (1,3,4,7)
n=6: (1,2,3,4,5), (1,3,4,5,9)
n=7: (1,2,3,4,5,6)
n=8: (1,2,3,4,5,6,7),
     (1,2,3,4,5,7,12),
     (1,4,5,6,7,11,13)
```

Every one is a two-seed addition chain.  AP accounts for `5/9`; the two named
sporadics are the unique primitive canonical `C`-transversal floor rows with
flip-set `{2}` for `n=5,6`; the two non-AP `n=8` rows are not perfect
transversals and have only nonunit missing shells in `C=15`.

For the user's named rows, S602 records:

```text
(1,3,4,7):
  first witness = 1/5, witnesses = 4, witness_Q = 5,
  C = 9, perfect_transversal = True, flipset = (2,),
  recipe = 1 seed; 3 seed; 4=1+3; 7=3+4.

(1,3,4,5,9):
  first witness = 1/6, witnesses = 2, witness_Q = 6,
  C = 11, perfect_transversal = True, flipset = (2,),
  recipe = 1 seed; 3 seed; 4=1+3; 5=1+4; 9=4+5.
```

## Tournament-Analysis Reading

S602 uses proof lenses as tournament vertices:

```text
AP_only,
unit_boundary_skeleton,
two_seed_addition_chain,
top_has_sum_pair,
top_previous_two_sum,
perfect_C_transversal,
flipset_{2},
nonunit_C_hole,
C_prime.
```

The observable is `(collapse_hits, -false_positives, lens_name)` over the
scanned universe.  The tournament is transitive:

```text
score_hist = {0:1, 1:1, ..., 8:1}
directed_3_cycles = 0
singleton SCCs = 9
Hamiltonian path =
unit_boundary_skeleton >
two_seed_addition_chain >
top_has_sum_pair >
perfect_C_transversal >
AP_only >
top_previous_two_sum >
C_prime >
flipset_{2} >
nonunit_C_hole.
```

The important assumption challenge is that runners or raw arcs are not the
right vertices here.  The quotient preserves the collapse predicate and the
additive/shell labels, but it destroys endpoint phase order inside `Q(V)`.
That is acceptable for this subproblem because the task is classification of
the collapse family, not proof of endpoint protection itself.

## Proof Program

1. Prove the equivalence between S599's covering-depth condition `p_0=0` at
   `delta=1/n` and S602's endpoint-critical/unit-boundary predicate, including
   the boundary-witness caveat.
2. Use HYP-1810 to pin the unit-boundary skeleton for collapsed rows.
3. Use THM-401/HYP-2135: a missed unit shell gives a `2/(2n-1)` witness, so a
   collapsed row must be a perfect `C`-transversal or live in nonunit shells.
4. Classify perfect-transversal floor rows by flip-set.  The bounded data says
   `{2}` is the only non-AP primitive flip that stays at the floor for `n=5,6`;
   all other one-flips jump above the second-gap edge in S553.
5. Classify composite-`C` nonunit branches by additive-chain generation and
   gcd-shell ramification.  This is the branch that matters for `n=14`, where
   `C=27=3^3`.

## See

`04-computation/lrc_p0_collapse_additive_chains_s602.py`,
`05-knowledge/results/lrc_p0_collapse_additive_chains_s602.out`,
`07-reflections/lrc-p0-collapse-additive-chains-s602.md`,
`04-computation/lrc_covering_depth_distribution_s599.py`,
`05-knowledge/results/lrc_covering_depth_distribution_s599.out`,
`07-reflections/lrc-covering-depth-distribution-the-master-object-s599.md`,
HYP-1802, HYP-1810, HYP-2084, HYP-2135, HYP-2138, HYP-2114,
HYP-2146, HYP-2151, HYP-2152, THM-401.
