---
id: THM-4507
title: "Finite polynomial-valuation corrections cannot make logarithmic Collatz height globally nonincreasing"
status: >
  PROVED elementary + INDEPENDENTLY AUDITED; FINITE-EXACT controls.
  For any finite prime set S, finitely many nonzero P_i in Z[X], c>0,
  arbitrary real function f of their valuation vector, and bounded g,
  c log n + f((v_p(P_i(n)))_(p,i)) + g(n) cannot be nonincreasing at
  every sufficiently large shortcut Collatz step, or every sufficiently
  large accelerated odd step. Positive integer polynomial roots may be
  excluded by a finite core. The same obstruction holds for the maximal-
  rise reset map C(n)=U^(v_2(n+1))(n), since each shadow period is one C step.
  The proof constructs arbitrarily long high
  positive segments with identical endpoint valuation vectors and
  unbounded endpoint height ratio. It also excludes any correction bounded
  on each fixed valuation fibre, including arbitrary interactions with
  finite colours. Positive witnesses change with requested block length;
  no divergent positive orbit or Collatz proof is claimed.
depends_on: []
related:
  - 01-canon/theorems/THM-4483-forced-charges-two-place-ranks.md
  - 05-knowledge/results/collatz_blueprint_20260921_energy.md
note: 05-knowledge/results/crossroads_crossing_20260926_resource.md
script: 04-computation/experiments/crossroads_crossing_20260926_resource.py
output: 05-knowledge/results/crossroads_crossing_20260926_resource.out
audit_note: 05-knowledge/results/crossroads_crossing_20260926_phase_audit.md
audit_script: 04-computation/experiments/crossroads_crossing_20260926_phase_audit.py
audit_output: 05-knowledge/results/crossroads_crossing_20260926_phase_audit.out
manifest: 05-knowledge/results/crossroads_crossing_20260926_manifest.json
audit: >
  Root derived the proof; the colour lane independently checked rational
  word legality, denominator units, avoidance of all polynomial roots,
  precision at every selected prime, complete positive-core avoidance,
  and the odd-map extension. Its independent code imports no source under
  audit and checks sixteen actual blocks, nine polynomial forms and four
  prime sets. The crossing lane separately checked a 273-step witness
  with prime set {2,3,11}. Root replayed all scripts normally and under
  Python -O; outputs agree. The fibre-bounded corollary uses the same
  fixed-vector endpoint argument, with a fibre-dependent finite bound.
source: crossroads-crossing-20260926; reserved and checkpointed before promotion
---

# THM-4507 — finite polynomial-valuation potential obstruction

**PROVED + INDEPENDENTLY AUDITED.** Full proof and quantifiers are in the
[resource note](../../05-knowledge/results/crossroads_crossing_20260926_resource.md).

Let T be shortcut Collatz. Fix a finite set of primes S and nonzero
integer polynomials P_i. No potential

    V(n)=c log n+f((v_p(P_i(n)))_(p,i))+g(n), c>0, g bounded,

with arbitrary real f, decreases weakly at every sufficiently large T
step. The same holds for the odd accelerated map and the maximal-rise
reset map C(n)=U^(v_2(n+1))(n). More generally the
correction may be any function bounded on each fixed valuation fibre.

## Mechanism

Choose B divisible by ord_p(2), and L divisible by ord_p(2),ord_p(3), for
all p>=5 in S. For sufficiently large a=1 mod L, the parity word
w=1^a0^B has expanding multiplier M=3^a/2^(a+B)>1 and negative rational
fixed point

    r=-(3^a-2^a)/(3^a-2^(a+B)).

Its denominator is a unit at every selected prime. The distinct r tend
to -1, so choose one outside the finite set of polynomial roots.
Simultaneous CRT approximation produces positive integers following w^k,
with enough initial 2-adic precision to pay for all k(a+B) divisions,
and enough odd-prime precision to fix every polynomial valuation.
Their endpoints n,y have identical feature vectors and

    y=r+M^k(n-r),       y/n>M^k.

Every intermediate state can be kept above any prescribed finite core.
The correction is bounded on this one fixed fibre; c log(y/n) is
unbounded. Summed nonincrease inequalities give a contradiction.

## Boundary and inherited scope

THM-4483 already uses expanding rational cycles to obstruct linear
2-adic counter banks, including a class of infinite banks. This theorem's
additional scope is arbitrary nonlinear dependence on finitely many
polynomial valuations at any fixed finite set of primes.

It excludes neither selected-return potentials nor unbounded information
inside each fixed feature fibre. The positive source changes with k.
Confusing these arbitrarily long finite shadows with one infinite
positive shadow would reverse the result's quantifiers.
