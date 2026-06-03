---
id: HYP-2153
status: NEW SUBPROBLEM; exact S602 bounded audit, classification open
source: codex-2026-06-03-S602
related:
  - HYP-1802
  - HYP-1810
  - HYP-2084
  - HYP-2135
  - HYP-2138
  - HYP-2141
  - HYP-2114
  - THM-358
  - THM-401
---

# HYP-2153: the `p0` collapse family is an additive-chain shell family, not just AP

## Operational Definition

The user's `p_0=0` notation is not yet a stable symbol in the LRC files.  Until
that notation is formalized, S602 uses the exact endpoint-collapse predicate
from S357/S359.

For a primitive speed row `V` with `|V|=n-1` and threshold `1/n`, call `V`
`p0`-collapsed when:

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

The `p0`-collapsed family is larger than the arithmetic progression family.
The first non-AP branch is additive:

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

1. Formalize the `p_0=0` notation and prove it equals the endpoint-collapse
   predicate above, or record the exact difference.
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
HYP-1802, HYP-1810, HYP-2084, HYP-2135, HYP-2138, HYP-2114, THM-401.
