---
id: HYP-2135
status: OPEN support-calculus formalization plus exact S591 named-row audit
source: codex-2026-06-03-S591
related:
  - THM-401
  - HYP-2138
  - HYP-2141
  - HYP-2137
  - HYP-2136
  - HYP-2132
  - HYP-2134
  - HYP-2133
  - HYP-2122
  - HYP-2118
  - HYP-2088
  - HYP-2084
  - HYP-2083
  - HYP-2059
  - THM-400
  - THM-397
---

# HYP-2135: after THM-401, LRC should be checked by a labelled sumset-support calculus

## Claim

THM-401 pins the natural odd modulus

```text
C = 2n - 1
```

but the remaining proof object is not the raw support of `V+V`.  It is a
labelled support calculus with at least five layers:

```text
S0 speed-shell support       shells hit by v_i mod C up to +/-.
S1 pair-sum shell support    shells hit by v_i+v_j mod C up to +/-.
S2 denominator support       actual pinch denominators D=v_i+v_j.
S3 shield support            which appearing D are killed by some D|v_k.
S4 visibility/lift support   unit-visible missing shells vs nonunit or lifted fibers.
```

The slogan is:

```text
speed support says which C-shells are occupied;
pair support says which pinches exist;
shield support says which pinches are already blocked;
unit support says which holes invert to immediate 2/C clocks;
lift support says when the witness has left the C-ledger and needs endpoint/DUN data.
```

This separates three things that were easy to conflate after HYP-2132:

```text
perfect antipodal transversal,
large or complete pair-sum support,
proof-relevant lonely witness.
```

They are not equivalent.

## Formal Definitions

Let `V={v_1,...,v_{n-1}}`, `C=2n-1`, and let

```text
sh_C(r) = 0                    if r = 0 mod C,
        = min(r mod C, -r mod C) otherwise.
```

For shells `A={1,...,n-1}` define:

```text
speed_C(V)(a) = #{i : sh_C(v_i)=a}
pair_C(V)(a)  = #{i<j : sh_C(v_i+v_j)=a}
Den(V)        = {v_i+v_j : i<j}
Low_n(V)      = Den(V) cap {1,...,n}
Shield(D,V)   = {v in V : D|v}
U_C           = {a in A : gcd(a,C)=1}.
```

Then:

```text
missing speed shells = {a : speed_C(V)(a)=0}
doubled speed shells = {a : speed_C(V)(a)>1}
unit holes           = missing cap U_C
nonunit holes        = missing - U_C
unshielded low D     = {D in Low_n(V) : Shield(D,V)=empty}
```

The extra lift coordinate records exact witness denominators from the pinch
lemma (HYP-2059).  If the exact maximizer has denominator outside the small
`n`/`C` clocks, the support ledger has projected away information and must be
lifted.

## Elementary Lemmas

**Lemma 1 (AP support identity).** For `AP_n={1,...,n-1}`:

```text
speed_C(AP_n)(a)=1 for every shell a,
pair_C(AP_n) misses exactly shell 1,
Low_n(AP_n)={3,4,...,n},
the only unshielded low denominator is D=n.
```

Proof sketch.  The first line is immediate because every lower shell is
occupied once.  Pair sums of distinct AP speeds range from `3` to `2n-3`; up to
antipodes mod `2n-1`, these hit shells `2,...,n-1` and cannot hit shell `1`.
For `D=3,...,n-1`, the speed `D` itself shields `D`; `D=n` appears as a pair
sum but speed `n` is absent.

S591 verifies this identity for `n=4..20` with no failures.

**Lemma 2 (unit-hole inverse witness).** Suppose a speed row has no zero residue
mod `C` and misses a unit shell `a`.  Let `k=a^{-1} mod C`.  Then at

```text
t = k/C
```

no speed is at distance `1/C`, so the minimum distance is at least `2/C`.

Proof sketch.  A speed would land at residue `+/-1` after multiplication by
`k` exactly when its original residue lay in shell `a`.  That shell is missing.
The zero-residue caveat is necessary because zero maps to zero.

**Lemma 3 (raw pair support is only a shadow).** Complete or large pair-shell
support does not decide the LRC floor.  S591 gives:

```text
open_gap_n7_S573: speed shells perfect, pair shells complete, M=5/33.
Vstar_n14:        pair shells complete, nonunit speed hole, M=1/14.
unit_shift_AP:    pair shells complete, unit speed hole, but D=n is shielded and M=1/8.
```

So `V+V` support must be labelled by speed shells, denominator shields, and
visibility/lift data before it becomes proof-relevant.

## S591 Evidence

`04-computation/lrc_sumset_support_calculus_s591.py` audits named rows from
S571-S573 and n=14.

Key outputs:

```text
n=11 C=21: unit shells 6, nonunit shells 4.
n=12 C=23: all 11 shells are unit-visible.
n=13 C=25: unit shells 10, nonunit shells 2.
n=14 C=27: unit shells 9, nonunit shells 4.
```

This rotates the finite-checking frontier: total `n=12` is clean for the
`C=2n-1` support ledger because `C=23` is prime, while `n=11,13,14` have
nonunit residual strata.

Named-row lessons:

```text
AP_n14:
  perfect speed-shell transversal; pair-shell support misses 1;
  D=n=14 is the only unshielded low denominator; M=1/14.

Vstar_n14:
  misses only nonunit shell 12 and doubles nonunit shell 3;
  pair-shell support is complete; D=14 remains unshielded; M=1/14.

doubled_apex_edge_n14:
  misses unit shell 13 and doubles shell 1;
  inverse clock t=25/27 gives M=2/27 exactly.

unit_shift_AP_n14:
  misses unit shell 1, but speed 14 shields D=14;
  inverse shell gives only the 2/27 edge lower bound, while exact M=1/8.

open_gap_n7_S573:
  perfect speed-shell transversal and complete pair-shell support;
  exact witness denominator is 33, so the C-ledger alone is a projection.

nonunit_hole_n8 rows:
  miss nonunit shell 6 and double nonunit shell 3;
  exact witness denominator is 23, again a lifted pair denominator.
```

The Tournament Analysis in S591 uses proof lenses as vertices, not runners or
arcs.  The proof ranking is transitive:

```text
low_denominator_shields
> unit_visible_holes
> witness_denominator
> speed_shell_transversal
> nonunit_holes
> gcd_strata
> pair_sum_shell_holes
> raw_sumset_size
```

It has score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`, zero directed
3-cycles, singleton SCCs, one Hamiltonian path, and 18 edge flips versus a
cost-only ranking.  The flips are useful: cheap raw support is not the right
proof lens.

## New Falsifiable Sub-Hypotheses

**H11 Support-lift descent.** Every strict sub-edge row with no unit holes is
either a floor row controlled by an unshielded `D=n` pincer, or it has a lifted
exact witness denominator that can be discharged by HYP-2088's `D/U/N` ledger
plus endpoint-owner labels.

**H12 Nonunit ramification.** Nonunit speed-shell holes are ramified fibers of
the `Z/C` support sheaf.  They are invisible to unit inverse clocks, but should
become visible after one of:

```text
gcd-stratum descent,
CRT lift,
endpoint-owner peeling,
pair-denominator shield accounting.
```

For `n=14`, the only nonunit shell strata are gcd `3` shells `{3,6,12}` and
gcd `9` shell `{9}`.

**H13 Unshielded-n floor gate.** In the additive/fold branch, the AP-like floor
rows are characterized not by maximal pair support but by the survival of an
unshielded `D=n` denominator after all low `D<n` denominators are shielded.
Shielding `D=n` pushes the row into the loose branch unless another endpoint
or Phi obstruction replaces the delta clock.

**H14 Support sheaf gluing.** The right n=14 certificate object is a sheaf over
antipodal speed shells:

```text
base:    shells in Z/27 up to +/-;
stalk:   pair denominators, shields, endpoint owners, D/U/N obligations;
gluing:  unit action, reflection, CRT, and lifted pair denominators;
defect:  nonunit ramification or monodromy in a forgotten endpoint fiber.
```

A proof quotient is allowed only if it preserves whether a fiber contains a
unit-visible hole, an unshielded low denominator, or a lifted exact witness.

**H15 Additive-chain collapse labels.** S602 adds a new label to the support
calculus: whether the row is generated by a two-seed addition chain, and
whether the top speed is a sum of lower chain terms.  In the targeted
`p0`-collapse boxes, all collapsed rows are two-seed addition chains.  The
prime/transversal sporadics `(1,3,4,7)` and `(1,3,4,5,9)` are the `{2}` flip
rows with top equal to the previous two chain terms; the composite `n=8`
non-transversal rows have only nonunit missing shells and still have top-sum
pairs.  Thus raw `V+V` support should be refined by additive-generation labels
before it is used as a floor-row classifier.  See HYP-2153 and
`lrc_p0_collapse_additive_chains_s602.py`.

## Relation To THM-401 And HYP-2132

THM-401 proves the modulus identity.  HYP-2135 proposes the next formal layer:
the residual finite arithmetic problem should be expressed as a labelled
support calculus at `Z/(2n-1)`, not as raw `S+S` size or unlabelled pair-shell
coverage.

This also clarifies addition versus multiplication:

```text
addition:
  creates pair denominators and summand shells;

multiplication by units:
  inverts unit shell holes into explicit 2/C clocks;

nonunit multiplication:
  creates ramified fibers that require lift/descent labels;

divisibility:
  shields pair denominators and turns visible folds into blocked pinches.
```

The proof route for `n=11,12,13` should be cheaper than the paper's exact
prime-fiber search precisely when this calculus closes before the fallback
finite sieve.  The `n=12, C=23` row is the clearest test: with no nonunit shell
strata, any missing speed shell is unit-visible.

## Relation To HYP-2141

HYP-2141 supplies the tournament-side restriction: LRC-accessible comparisons
are round, so the regular tight beat is the interval circulant rather than a
Paley or QR beat.  In this support-sheaf language, the base shell order is not
arbitrary support data; it is the round additive order.  Unit action remains a
witness symmetry, while the beat itself is the interval structure that supplies
the AP row and the `D=n` pincer.

## Assumption Challenge

Do not make tournament vertices be runners by default.  In this hypothesis,
the vertices are support lenses, shell fibers, denominators, missing-shell
obligations, and proof labels.

The quotient preserves:

```text
unit-visible missing shells,
nonunit holes,
unshielded low pair denominators,
lifted exact witness denominators.
```

It destroys exact time order, unmarked tournament isomorphism class, and raw
phase geometry.  If that destruction mixes floor/open-gap/loose rows, the
quotient must be lifted with endpoint-owner, pincer, or D/U/N labels.

## See

`04-computation/lrc_sumset_support_calculus_s591.py`,
`05-knowledge/results/lrc_sumset_support_calculus_s591.out`,
`04-computation/lrc_p0_collapse_additive_chains_s602.py`,
`05-knowledge/results/lrc_p0_collapse_additive_chains_s602.out`,
`07-reflections/lrc-sumset-support-calculus-s591.md`,
THM-401, HYP-2153, HYP-2138, HYP-2141, HYP-2137, HYP-2136, HYP-2132, HYP-2134, HYP-2133, HYP-2122, HYP-2118,
HYP-2088, HYP-2084, HYP-2083, HYP-2059, THM-400, THM-397.
