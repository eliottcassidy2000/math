---
id: THM-3267
title: "Norm-phase factorization ladder and projective-determinant blindness"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED in the chosen
  F_169 Singer model. The THM-3255 C_12 norm phase is exactly recoverable
  from a full nonzero increment q, has exactly its C_2 parity quotient after
  projectivizing q, and has no nontrivial quotient through the determinant
  coordinate v=det(q,R). The pair ([q],v) still retains exactly parity.
  Affine translations destroy even that parity, while the THM-3246 negative
  seam is a one-per-phase, twelve-distinct-direction transversal in every
  primitive Singer gauge. The pinned nearest physical contracts do not
  supply the owner-to-endpoint intertwiner needed to realize these abstract
  coordinates on one literal ancestry sheet. No physical current, row
  exclusion, safety certificate, or LRC(14) decrement is supplied.
source: root/creative-synthesis-cont/2026-08-03
audit: >
  The exact companion pins eighteen inherited artifacts, replays THM-3246,
  reconstructs F_169 and every primitive Singer gauge, and exhausts all
  projective, determinant, mixed, affine-translation and seam fibres. Normal
  and optimized executions byte-match the frozen transcript; the source has
  no assertion node or floating literal. An independent implementation used
  integer pair arithmetic and geometric scalar orbits rather than the
  candidate's exponent-defined fibres; it recovered the field order, norm
  phase, all fibre populations, both translation hostiles, the complete
  translation-difference census, and the seam matching. A separate contract
  audit confirmed the stated boundary in each pinned nearest source. It also
  repaired the transcript label distinguishing 168 translations from 28,056
  nonzero-to-nonzero transitions.
depends_on:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word
  - THM-3247-heisenberg-central-fourier-decomposition-and-canonical-current-cyclicity
  - THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-phase-marker-boundary
related:
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-3252-singer-compactified-owner-hodge-word-universal-charged-cyclicity
  - THM-3253-positive-owner-mass-newton-cyclicity-and-maximal-common-heisenberg-module
script: 04-computation/lrc_norm_phase_physical_factorization_no_go_scout_20260803.py
output: 05-knowledge/results/lrc_norm_phase_physical_factorization_no_go_scout_20260803.out
script_sha256: a16658e5b2ad41d9b32ebd8b7b5bd7ffedc26a913b7a0bff75d9fad2f2a617e3
output_sha256: 392d903d8451cb9df895050bab54ead3f515580950b6a9dd243061a3409a857c
hash_basis: LF-normalized bytes
---

# THM-3267 -- full phase, projective parity, determinant blindness

**PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED in one chosen
Singer model.**

THM-3255 finds a missing eleven-dimensional multiplicative sector indexed by
the twelve values of the field norm. This theorem determines exactly which
abstract endpoint coordinates retain that phase, and where the information
is destroyed. The result is a factorization boundary, not the missing
physical LRC bridge.

## 1. Chosen norm phase

Work in

~~~text
F_169=F_13[u]/(u^2-2),       alpha=1+2u.
~~~

The element alpha has order 168 and

~~~text
Norm(alpha)=6 in F_13^*,     ord(6)=12.                (1)
~~~

For q=alpha^j, define

~~~text
phi(q)=log_6 Norm(q)=j mod 12 in C_12.                 (2)
~~~

Every norm fibre has fourteen points. Thus a full nonzero q, together with
the chosen field/norm gauge, recovers the complete C_12 phase exactly.

## 2. Projectivization retains exactly parity

A projective direction [q] is an F_13^*-orbit of size twelve. In the
canonical exponent gauge its exponents are

~~~text
r+14k,      0<=k<12.                                  (3)
~~~

Their phases are r+2k modulo 12. Every direction therefore contains exactly
the six phases of one parity, each twice. Seven of the fourteen directions
carry the even phases and seven carry the odd phases.

Consequently [q] determines phi(q) modulo 2. It determines no finer quotient
of C_12: within one direction all six values of the fixed parity occur. The
sharp surviving quotient is C_2, as also predicted by gcd(12,14)=2.

The companion verifies this statement geometrically, using scalar orbits,
over all 8,064 primitive Singer gauges and all

~~~text
8064*14=112896                                           (4)
~~~

projective fibres. Each fibre has the same six-phases-twice census.

## 3. Determinant blindness

Retain the THM-3247 endpoint-plane coordinate

~~~text
v=det(q,R),       q!=0,       R in F_13^2.             (5)
~~~

For fixed q and v there are thirteen origins R. Hence each of the thirteen
v-fibres across all q has

~~~text
168*13=2184 states, with every phase appearing 182 times. (6)
~~~

Thus v alone carries no nontrivial quotient of the norm phase.

Nor does v repair the information lost by projectivization. There are 182
fibres of ([q],v); every one contains 156 states and the six phases of the
parity selected by [q], each exactly 26 times. Therefore

~~~text
full q + chosen norm  -> exact C_12 phase,
[q]                   -> exact C_2 parity,
v                      -> trivial phase quotient,
([q],v)                -> exact C_2 parity.            (7)
~~~

## 4. Affine translations destroy parity

Although phi is a function of the full affine point q, it is not invariant
under the standard H_13 point action. Translation by (1,0) gives

~~~text
alpha^0=(1,0) -> (2,0)=alpha^70,    phase 0 -> 10,
alpha^1=(1,2) -> (2,2)=alpha^40,    phase 1 -> 4.      (8)
~~~

The second transition changes parity. Among the 167 nonzero-to-nonzero
transitions for this translation, 13 preserve phase and each nonzero phase
difference occurs 14 times.

Across all 168 nonzero translations there are 28,056 nonzero-to-nonzero
transitions. Phase difference zero occurs 2,184 times, while each of the
eleven nonzero differences occurs 2,352 times. Thus full affine covariance
leaves no nonconstant invariant quotient of phi. Any physical use of the
phase must carry it as a transforming sidecar or supply additional structure.

## 5. The negative seam is a transversal, not a marker

THM-3246's negative seam has exponent set

~~~text
N={0,1,2,3,4,5,162,163,164,165,166,167}.              (9)
~~~

It contains exactly one exponent in every C_12 phase and meets twelve
distinct projective directions in every one of the 8,064 primitive gauges.
In the canonical gauge, its phase-to-direction matching is

~~~text
0->0, 1->1, 2->2, 3->3, 4->4, 5->5,
6->8, 7->9, 8->10, 9->11, 10->12, 11->13,             (10)
~~~

with directions 6 and 7 missing. This is a complete phase transversal, not
a selected phase and therefore not the marker missing in THM-3255.

## 6. Physical type boundary

The exact factorization ladder above lives in the chosen abstract Singer
model. The closest pinned physical constructions do not yet realize it:

- THM-2791 retains a same-ancestry partial rail-sheet germ but stops before
  endpoint allocation;
- THM-3234 records that Singer compactification forgets q, interval geometry,
  phase endpoints and physical ancestry;
- THM-3247 supplies the abstract current J_q(R)=P_(R+q)Q_R but leaves the
  physical descent/intertwiner open; and
- THM-3252 and THM-3253 use external or abstract plane relocations rather than
  a physical owner-to-endpoint map.

This is a statement about those pinned nearest contracts, not a universal
nonexistence theorem. The first missing physical record is an allocated

~~~text
(q,R) on the same literal ancestry sheet,                     (11)
~~~

together with a compatible F_169 Singer structure or equivalent absolute
norm gauge. The cheapest decisive next test is to construct that record on
the two cylinders of the THM-2791 partial germ and check the owner/current
square while retaining the literal contributor labels.

No such square is constructed here. Nothing in this theorem is a physical
owner current, cellwise safety statement, row exclusion, lonely-runner
witness, or decrement for LRC(14).

## 7. Exact reproduction

Run

~~~text
python3 04-computation/lrc_norm_phase_physical_factorization_no_go_scout_20260803.py
python3 -O 04-computation/lrc_norm_phase_physical_factorization_no_go_scout_20260803.py
~~~

and compare LF-normalized bytes with the declared transcript. The companion
uses exact finite-field arithmetic, exhausts every stated fibre and gauge,
replays the inherited seam companion, and pins every selected connection
contract.

QED.
