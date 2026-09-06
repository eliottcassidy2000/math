# Independent audit of the joint-position grid supplier

**Status: INDEPENDENT ANALYTIC/SOURCE AUDIT PASS + EXACT REPLAY.** This
audit accepts the two theorems in
[the primary note](long_frontier2_sep06_lrc.md), with their stated finite
scale universes and actual-entry hypotheses. It does not prove universal
LRC(14). The auditor read the complete proof and standalone source, then
replayed normal and optimized execution against the frozen output.

## 1. All translates, open endpoints, and common-factor transport

The literal arc intersection construction is complete in units
`1/(14pq)`: it intersects ordered danger arcs, advances both cursors at
simultaneous ends, and rejoins the artificial split at zero. The independent
clipped-length multiset is consistent with every retained component.
For the two declared pairs and `620<=t<=790`, each component has length
strictly less than `1/t`; consequently its projected phase arc contributes
either zero or one grid point. Disjoint original components cannot count
the same grid point.

The sweep correctly treats danger sets as open. At each phase endpoint it
removes all ending arcs, evaluates the endpoint while starting arcs are
absent, and then adds those starts to obtain the next open-cell count.
Wrapped arcs give the initial count just before zero, and cyclic closure
is checked. Endpoint counts are no larger than their right-cell counts;
the resulting minimum covers every real phase. The separate literal path
checks every endpoint and every open-cell midpoint for both pairs at
`t=621,633,645,687,761`. It does not replace the complete sweep over the
84 admissible scales.

For the physical triple `h*(101*113,101*128,113*128)`, the first and third
members reduce to `(101,128)` after division by `113h`; the first and
second reduce to `(113,128)` after division by `101h`. Those factors are
units modulo `t` under the stated hypotheses. Multiplication therefore
permutes any translated grid into another translated grid, so each
all-translates minimum survives the reduction. A swapped pair reference
in the draft prose was repaired before this acceptance; the source and
finite certificate were already correctly typed.

Subtracting one pair intersection from the sum of seven danger counts is
a valid pointwise union bound. Taking the larger of two available pair
credits is valid; summing them is not used. The exact inequality
`max(kappa_1,kappa_2)>=7ceil(t/7)-t+5` gives at least five weak-safe
points on every translate. The primary-pair failures and secondary repairs
are genuine literal half-phase controls, with the declared common factors
odd at those five actual scales.

## 2. Complete actual entries and strict safety

The full graph reconstruction, primitive/disjoint shape checks, and strict
physical sum bound cover all 81 actual scales. The internal atlas trees
are connected with coefficient heights at most 356. Cross pairs are
coprime and have reduced sums greater than 356, so there are exactly the
stated components of sizes six and seven.

The mixed-support exclusion covers every nonzero distinguished coefficient,
not only its minimal multiple. For `tv,tw,u`, that coefficient is divisible
by `t gcd(v,w)>Q`. For `u,w,tv`, division by `gcd(u,w)` forces an integer
multiple of `tv`, whose magnitude is at least `t min(V)>394Q`, whereas
the two bounded internal terms have magnitude at most `394Q`. These also
exclude mixed support-two degenerations. Internal weighted tree relations
span dimension eleven, giving both containments in `W_(Q,3)=V_dec`.
The native signed-box controls independently verify the two orientations
on all 231 mixed supports at five named scales.

Every physical subset of size seven through twelve has gcd one: a mixed
subset contains a coprime cross pair, and the only possible pure subset
is all seven primitive `U` labels. Its complementary word is therefore
all ones. The source uses the pinned complete profile bank and checks
all 4095 such subsets at each of the five full controls, rather than
substituting only numerical ceilings.

At the half-grid, the six odd `V` labels have clearance `1/2` after the
scale factor is restored. Among `U`, only one label is odd. An even `U`
label has odd reduced phase denominator because `t` is odd, so it cannot
have clearance exactly `1/14`. At any weak-safe point, at most one label
can lie on the boundary. A sufficiently small motion into its good side
preserves every other strict inequality. This establishes strict safety;
the literal first-safe-point controls independently attain strict clearance
in all 81 rows.

## 3. Consumer boundary, replay, and scope

All 21 primitive `U` pairs have larger coefficient at least 113. Since
`t<=790<7*113`, every separate open-component credit is zero. This
verifies the stated failure of bounds which sum only those credits;
it does not assert failure of all phase methods. The listed forward,
endpoint, dual, whole-arc and nested-origin inequalities have the proper
six-plus-seven normalization, and their explicit source checks match the
limited sufficient-consumer claims.

The generic theorem ranges over exactly 84 grid sizes, with arbitrary
remaining tail heights subject to its hypotheses. The actual theorem
ranges over exactly 81 boxed rows at the displayed fixed shape heights.
Neither is an unbounded scale theorem, an unsafe family, or an inference
from sampled phases.

Replay from the repository root:

```bash
python3 -B 04-computation/long_frontier2_sep06_lrc.py > /tmp/long_frontier2_lrc_audit_normal.out
python3 -B -O 04-computation/long_frontier2_sep06_lrc.py > /tmp/long_frontier2_lrc_audit_optimized.out
cmp /tmp/long_frontier2_lrc_audit_normal.out /tmp/long_frontier2_lrc_audit_optimized.out
cmp /tmp/long_frontier2_lrc_audit_normal.out 05-knowledge/results/long_frontier2_sep06_lrc.out
```

All three outputs are byte-identical and pass **49,880 always-active exact
gates**. The source imports no research producer. Its inherited profile
bank is hash checked, with read-only `git show HEAD:path` fallback for
sparse omission. No Git mutation is needed by the verifier.

```
source SHA256
83ddf56931ae05e8d467ee130f11f862279433df3a649960ae1ae12019f8b614
frozen and both replay outputs SHA256
90382eb5d6aad90b3f0988e9055764ca1238a7cb62107e79a0fe7b841d379dc5
semantic SHA256
54f257677f163e86e1f59dc5cfdaac6d65a4d00756474fb88b0481ceddcf76ed
inherited profile-bank SHA256
935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f
```

No mathematical correction remains. The primary source and frozen output
were not edited by this audit.
