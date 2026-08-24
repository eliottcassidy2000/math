---
id: THM-3903
title: "Strict Arithmetic-Kakeya cell optima: two-by-two is 7/4 and saturated empty-wildcard two-by-three is 11/6"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the public strict
  equal-suffix paid-row forcing semantics, the exact optimum on dimensions
  [2,2] is 7/4.  The witness has m=4, r=3, n=4, t=0 and firing trace 1+3.
  A symbolic 60-family audit excludes its only smaller cyclic possibility.
  In the saturated [2,3] ladder with empty wildcard set and all seven edge
  slots active, four seeds are necessary and sufficient, so the exact optimum
  in that restricted universe is 11/6.  This is a local forcing-certificate
  theorem, not an arbitrary-[2,3] theorem, an Arithmetic-Kakeya bound, a
  private-verifier audit, or a 1.675 result.
source: root / strict k=2 cyclic-cell and saturated-ladder exact searches, 2026-08-23
audit: >
  INDEPENDENT POSITIVE-CONTROL AND TRACE-COVERAGE AUDIT PASS.  The primary
  path replays the seven strict
  paid rows in the maintained forcing engine, checks the explicit first
  firing combination, both rank drops, every physical-row deletion, and the
  loose-rule hostile.  A separate SymPy implementation works in raw (a,b)
  coordinates and solves all 10 seed-position multisets times 6 prospective
  first firing pairs.  Its generic-rank census contains no required (6,4)
  case.  For saturated [2,3], two independent positive controls agree; exact
  rank-family and Groebner audits eliminate every two-round and multiround
  three-seed trace, while a separate combinatorial enumerator confirms that
  the three canonical trace families are exhaustive.  The rational Groebner
  elimination is one exact CAS path, not independently reimplemented.  Normal,
  optimized, and frozen streams byte-match; finite-field exhaustions are
  hostile evidence only.
depends_on:
  - THM-2850-paid-distinguished-coloop-round-budget-and-slope-grammar-rank-defect
  - THM-3875-arithmetic-kakeya-k1-forest-one-round-score-floor
script: 04-computation/ak_strict_c4_score_seven_four_thm3903.py
output: 05-knowledge/results/ak_strict_c4_score_seven_four_thm3903.out
script_sha256: 71d0ef85a0fb023c8a01524efbf8bb7c93b120178105b6b31634ffa379502ba0
output_sha256: e6fb8600aad22fed5a84809e9135b526807cf5e3c44deb670b820af75e5432cc
semantic_sha256: 1359b2f2924587642590d7fc2756546c0912cc6a4c62ed323e8aa63d40421544
independent_audit_script: 04-computation/ak_strict_c4_score_seven_four_independent_audit_thm3903.py
independent_audit_output: 05-knowledge/results/ak_strict_c4_score_seven_four_independent_audit_thm3903.out
independent_audit_script_sha256: 13a8d295499d22c9a998ed13dcd98e0a23fb8fef87b0ad6bd57ac88d2c8e5997
independent_audit_output_sha256: 93b2804560df95c10690a76ee245db6148bbad4c21d2c1c4292ae55aed0e9760
independent_audit_semantic_sha256: d1e21eca692789bcbe00f3a232ef0e1bb7cee890f24123a9145f6fc9e4248888
saturated_2x3_positive_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_positive_control.py
saturated_2x3_positive_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_positive_control.out
saturated_2x3_positive_script_sha256: 9fe8b9f4b408e747d22ae4f669b4b95a56502853672764d9871c9be7981e9e17
saturated_2x3_positive_output_sha256: ed382cd9d43502d76ad1698d62c01bb667ae93f06deee04cb5c53c131be38a44
saturated_2x3_positive_semantic_sha256: 4a285d760f3510cf222a1022d3bdcdfd109cda528b28408a1328171aa3f723d9
saturated_2x3_independent_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_positive_control_independent.py
saturated_2x3_independent_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_positive_control_independent.out
saturated_2x3_independent_script_sha256: 943847cb82a486071a6971707dcd9e93c02cc50e15ac95b1d50e0ba394893008
saturated_2x3_independent_output_sha256: 8b9a0c95af6be3d62d8baf290ea9e6b56b1ccc25a78ae191dd7f0c48656faf37
saturated_2x3_independent_semantic_sha256: fd12df1760197bb5a63c124e3bd3351bf4cc600cd214aa38b26be6c91fdbdd2a
saturated_2x3_helper_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_all_cost10_kernel_search.py
saturated_2x3_helper_script_sha256: 2534eba7916af3261a56c26a22196273b4ca35da936f4c11e03616188344be94
saturated_2x3_f4_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f4_symbolic_audit.py
saturated_2x3_f4_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_full_f4_symbolic_audit.out
saturated_2x3_f4_script_sha256: 10bbc29434f79d5a4781911b8198f0fccde0f842c56beb4356d28a5e62b8d1a8
saturated_2x3_f4_output_sha256: dfa0af90f707aed3f76d61077ef7087d7076f9a996981d51ed17b4e72d7bbb77
saturated_2x3_zero_cover_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f3_zero_cover_audit.py
saturated_2x3_zero_cover_script_sha256: e293c3e8ab39fb02c9e0dd723ad59790de25df0633dc3d7e4b86c29c80d1f22b
saturated_2x3_f3_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f3_exact_audit.py
saturated_2x3_f3_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_full_f3_exact_audit.out
saturated_2x3_f3_script_sha256: ebc272ce627a755b966cdee473b626ee71f95cfc28eea72f0fdc57483310c12a
saturated_2x3_f3_output_sha256: 4747d5a4afcadb7a5bc9c1f66913ad95c42fae0453412d051c17604259f3c2d0
saturated_2x3_f2_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f2_exact_audit.py
saturated_2x3_f2_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_full_f2_exact_audit.out
saturated_2x3_f2_script_sha256: 878f77240a512e7f0aa9947595655a67f261ea3eeca051944713ade066b99732
saturated_2x3_f2_output_sha256: dcb128dc9501865f9edd757a163bcba18848bb719fe71d00507687a9caa8a27c
saturated_2x3_rank9_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_rank9_f3_audit.py
saturated_2x3_rank9_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_full_rank9_f3_audit.out
saturated_2x3_rank9_script_sha256: 4c81650c92bbdcab8515ff06c9b9078157f4a2e9c0dc7b856d15e29df43080ed
saturated_2x3_rank9_output_sha256: 6880cd09710128176e0cb02ae3eb15129e86a66e2e0b9fc0768f96346b4b172f
saturated_2x3_multiround_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_multiround_rank_variety.py
saturated_2x3_multiround_script_sha256: dd56fec396f77629e7d991c18752bf45c8354de190c7612a78b7497508fe42a3
saturated_2x3_rung3_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_multiround_rung3.out
saturated_2x3_rung3_output_sha256: 33b7052b35c045d9e170345ae172981bc8d7d27669f451c7fff9ac653d545dfa
saturated_2x3_corner3_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_multiround_corner3.out
saturated_2x3_corner3_output_sha256: c79ff38fd3a4c9ef86a080dd746b61eabe7415fc683fc47ec42a67e183f7e089
saturated_2x3_corner4_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_multiround_corner4.out
saturated_2x3_corner4_output_sha256: 96882d07db9990a9a0593336e60a23a35617dc6d0c3441ef32e73e22748801de
saturated_2x3_trace_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_trace_coverage_independent.py
saturated_2x3_trace_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_trace_coverage_independent.out
saturated_2x3_trace_script_sha256: c68c2d54fe55851f3226102bafd7b2d8c6df071fd24e143cc66c18967d88ed3b
saturated_2x3_trace_output_sha256: f338d2a26ddc6ef03ec104a1ebc88430459bec4e54978866955ec074c838296d
saturated_2x3_trace_semantic_sha256: 0657c9b3951350c3d8bae697d08c9b170a68be1be1b29ba62194c106d9732adc
saturated_2x3_finite_field_script: 04-computation/ak_strict_saturated_2x3_thm3903/ladder_modular_exhaust.cpp
saturated_2x3_finite_field_output: 05-knowledge/results/ak_strict_saturated_2x3_thm3903/ladder_modular_gauge_exhaust.out
saturated_2x3_finite_field_script_sha256: 84f3093e3574b2ebd7051aa659aab2d767b541223a25e7c6ab8d30f95e90ea43
saturated_2x3_finite_field_output_sha256: 4eb532f60896edaa8dc826b27836416709b8bbc0c2d448304182e08ccc5f9eae
hash_basis: raw LF bytes
---

# THM-3903 -- strict `[2,2]` and saturated empty-wildcard `[2,3]`

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  This theorem concerns
only the repository's public strict/equal-suffix paid-row forcing semantics.
It neither imports a certificate-to-`AK(s)` theorem nor audits a private
benchmark verifier.  In particular, neither exact cell result below proves an
Arithmetic--Kakeya bound.

## 1. Exact witness

Name the vertices of the strict `[2,2]` cell cyclically

```text
A=(1,1),       B=(1,2),       C=(2,2),       D=(2,1).
```

Use the following raw labels.  Their normalized slopes are
`rho=(a-b)/(a+b)`:

```text
f_1(1)   = (5,-3),       rho=4,       creates A-D and B-C;
f_2(1,1) = (0,2),        rho=-1,      creates A-B;
f_2(2,1) = (4,-2),       rho=3,       creates D-C;

seed A   = (-2,4),       rho=-3;
seed B   = (2,0),        rho=1;
seed C   = (47,-27),     rho=37/10.
```

Every label has `a+b!=0`.  The first coordinate entry pays for its two strict
matching rows, so

```text
m=4,       r=3,       n=4,       t=0,
score=(m+r)/(n-t)=7/4.                                  (1)
```

Order the seven physical rows as

```text
[A-D, B-C, A-B, D-C, seed A, seed B, seed C].            (2)
```

The integer coefficient vector

```text
(-6,-14,21,-6,-15,35,-2)                                (3)
```

is supported only at `D` and equals `6(1,-1)` there in raw coordinates,
or `12(0,1)` after `(a,b)->(s,d)=(a+b,a-b)`.  Thus `D` fires.  After `D`
becomes a wildcard, the restriction of the remaining rows to `A,B,C` has
rank six, so all three vertices fire simultaneously.  The strict trace is

```text
(D), (A,B,C),                rank profile 7 -> 6 -> 0.    (4)
```

The three cycle constraints in normalized coordinates have rank three and
kernel generated by `(-1/10,1/10,1,0)`, so `D` is the unique first-round
firing.  Changing only the final seed slope from `37/10` to `47/10` destroys
that first cofactor and stalls, showing that the designed cancellation is a
real codimension-one condition.

Deleting any one of the seven physical rows stalls before the first firing.
Deleting the shared `f_1(1)` slot removes both horizontal rows, so every legal
paid-slot deletion also fails.  The displayed strict witness is therefore
deletion-minimal for its labels.

## 2. Lower bound on the complete strict cell

Consider any successful strict `[2,2]` cell with nonvertical labels.  If an
initial wildcard is present, or if one of the four cycle edges is absent,
the live simple graph after wildcard projection is a forest.  THM-3875 then
forces score at least `2`.  Every candidate below two therefore has no
initial wildcards and contains all four cell edges, hence has `m=4`.

After any nonempty first firing, the remaining live graph is a forest.
THM-3875's one-round barrier at that state permits at most one further
nonempty round.  THM-2850 then gives `g/u>=3/2`, where `g=m+r` and `u=n-t`.
Thus fewer than two seeds are impossible.

It remains to exclude two seeds.  Such a cell would have

```text
g=6,       u=4,       score=3/2.                          (5)
```

Equality in THM-2850 forces two firing rounds of sizes `(2,2)`, initial
generator rank six, and terminal restriction rank four on the two surviving
vertices.

Normalize every nonzero row to `(1,rho)`.  Let `alpha` be the shared
horizontal slope, `beta,gamma` the vertical slopes, and `s_0,s_1` the two
seed slopes.  For each seed-position multiset and prospective two-vertex
first firing set `F`, the cell has two cycle constraints `W` on the four
distinguished coordinates.  Firing exactly `F` requires

```text
rowspan(W)=span{e_v : v in F};                            (6)
```

in particular both rows of `W` vanish away from `F`.  These are linear
equations in `(alpha,beta,gamma,s_0,s_1)`.

The independent audit solves all

```text
10 unordered seed-position multisets x 6 firing pairs = 60 families.
```

On each resulting linear family it computes the generic pair

```text
(initial generator rank, terminal restriction rank).
```

The exact census is

```text
(4,2):  4 families
(4,3):  8 families
(5,2): 30 families
(5,3):  4 families
(6,2): 14 families.                                      (7)
```

The required pair `(6,4)` never occurs.  Specialization can only lower a
generic rank, so no exceptional rational slope choice escapes the census.
Therefore two seeds are impossible, three are necessary, and Section 1
attains the exact optimum `7/4`.

## 3. Loose-rule hostile and scope

Under the literal suffix-unconstrained rule-(1), the same paid first-axis
entry creates a complete-bipartite transporter bank rather than the two
strict matching rows.  The paid data then acquire an eighth uncharged basis
row, changing the trace to the single round `(A,B,C,D)`.  No part of the
proof above uses that row; it is retained solely as a hostile control for the
known unsound loose semantics.

The `[2,2]` result improves the prior strict local record `13/7` and reaches
the Katz--Tao `7/4` rung inside this one forcing cell.  It does **not** reach
`1.675`, establish the strong tournament inequality, validate the private
verifier, or prove the certificate-to-Arithmetic--Kakeya implication.

## 4. Saturated empty-wildcard `[2,3]` has exact optimum `11/6`

This section moves one rectangle beyond the cyclic cell, but in a deliberately
fixed universe.  Name the ladder vertices

```text
A=(1,1), B=(1,2), C=(1,3), D=(2,1), E=(2,2), F=(2,3).
```

The first-axis slot supplies the three tied rungs `AD,BE,CF`; all four
second-axis slots `AB,BC,DE,EF` are active.  Thus `m=7`.  The initial wildcard
set is empty.  Seeds may repeat at a vertex, and every label is rational and
nonvertical.  Within exactly this universe, four seeds are necessary and
sufficient.

### 4.1 Exact four-seed witness

Use the normalized slopes `rho=(a-b)/(a+b)` and raw representatives

```text
AD,BE,CF: rho=4,     raw (5,-3)
AB:       rho=-1,    raw (0,2)
BC,EF:    rho=-8,    raw (-7,9)
DE:       rho=3,     raw (4,-2)

seed A:   rho=-3,    raw (-2,4)
seed B:   rho=1,     raw (2,0)
seed E:   rho=37/10, raw (47,-27)
seed C:   rho=-7,    raw (-6,8).
```

Then

```text
m=7, r=4, n=6, t=0, score=(m+r)/(n-t)=11/6.
```

In row order

```text
AD, BE, CF, AB, BC, DE, EF, seed A, seed B, seed E, seed C,
```

the integral combination

```text
(-6,-14,0,21,0,-6,0,-15,35,-2,0)
```

equals `6(1,-1)` at `D` and zero elsewhere.  Hence `D` fires.  The remaining
ten columns on `A,B,C,E,F` have rank ten, so those five vertices fire together.
The raw rank profile is `11 -> 10 -> 0`, with trace

```text
(D), (A,B,C,E,F).
```

The maintained-engine replay and a separate raw-coordinate SymPy replay agree
on the combination, ranks and trace.  Each legal paid-slot or seed deletion
stalls before all six vertices fire.  The same data under the unsound loose
rule create thirteen rows and fire all six vertices at once; this is a hostile
control, not part of the strict proof.

### 4.2 Exact exclusion of every two-round three-seed trace

Suppose instead that `r=3`.  Then `g=m+r=10`, `u=6`, and the candidate score is
`5/3`.  There are

```text
C(6+3-1,3)=56
```

unordered seed-position multisets, including repetitions.  The three rung
rows share one slope `alpha`; subtracting `alpha` from every slope is an
invertible shear fixing the distinguished direction, so the audit gauges
`alpha=0`.

Let `A` be the normalized `10 x 12` generator and let `W` be its four-row
cycle-slope matrix.  Write `lambda=10-rank(A)`.  THM-2850's two-round balance
window leaves exactly

```text
lambda=0: first-wave sizes 2,3,4;
lambda=1: first-wave size 3;
lambda>=2: impossible.
```

Every surviving branch is eliminated over `Q`:

```text
rank 10, waves 4+2:  840 linear families; required ranks (10,4) absent;
rank 10, waves 3+3: 2576 branches = 2296 generic failures + 280 unit ideals;
rank 10, waves 2+4:  840 cases    =  342 generic failures + 498 unit ideals;
rank  9, waves 3+3: 1120 cases, 1078 families; required ranks (9,6) absent.
```

For `4+2` and rank-nine `3+3`, exact linear-family generic ranks suffice:
specialization cannot raise the terminal rank.  For full-rank `3+3`, the two
live square-cycle rows force a finite coordinate-zero cover; a separate exact
cover audit proves completeness.  Each cover is then saturated by real Gram
sums certifying the required initial and terminal ranks.  For `2+4`, all
rank-at-most-two minors are imposed and every admissible pair of nonzero
initial and terminal minors is saturated.  Over `Q`, a Gram sum of squares is
nonzero exactly when one of its minors is nonzero.  Thus a unit Groebner basis
really excludes the corresponding rational rank stratum rather than merely a
generic point.

### 4.3 Multiround closure

After a first firing, every vertex deletion except one outer corner or one
whole outer rung leaves a forest.  THM-3875 then permits only one further
nonempty round, already covered above.  Deleting an outer rung leaves a `C4`
and permits at most three total rounds.  Deleting one outer corner leaves a
`C4` with a leaf; it permits three rounds, or four only when the other endpoint
of that outer rung fires alone next.  Reflection reduces the first wave to
`{A,D}` or `{A}`.

For each live set `L`, any successful continuation must have every component
pinned; otherwise THM-2850's unpinned-component kernel prevents that component
from beginning to fire.  On a pinned state form

```text
W_L = C_L diag(rho) B_L,
```

where `C_L` is an exact cycle basis of the restricted incidence matrix.  The
headroom is `H_L=rank(W_L)`.  If `F` is the next firing wave, its distinguished
columns being coloops forces

```text
rank(W_L restricted to L\F)=H_L-|F|.
```

The audit enumerates every necessary nonincreasing headroom profile, imposes
all upper-rank minors, and uses separate Rabinowitsch saturations of Gram sums
for every exact lower-rank condition.  Across all 56 seed multisets the exact
census is

```text
outer rung,   3 rounds: 3360 = 742 generic failures + 2618 unit ideals;
outer corner, 3 rounds: 3920 = 782 generic failures + 3138 unit ideals;
outer corner, 4 rounds: 4480 = 982 generic failures + 3498 unit ideals;
potential rational survivors: 0.
```

A structurally separate combinatorial enumerator finds 148 forest-compatible
ordered three-round traces, 56 four-round traces, and none with five or six
rounds; after symmetry and headroom expansion it independently reproduces the
three Groebner universes `3360,3920,4480`.  This audits trace coverage, not the
Groebner elimination itself.

One round is excluded by THM-2850.  If a successful ladder had fewer than
three seeds, adding legal seed rows would preserve every existing singleton
row-space witness and produce a successful three-seed ladder, now excluded.
Therefore `r<=3` is impossible, while Section 4.1 attains `r=4` and `11/6`.

Complete common-slope-gauged searches over `F_2,F_3,F_5` also find zero
three-seed successes, including all 4,375,000 assignments over `F_5`.  This is
independent hostile evidence only; it is not the rational elimination.

The theorem does **not** cover inactive edge slots, a nonempty initial
wildcard set, arbitrary `[2,3]` instances, larger rectangles, the private
verifier, or any certificate-to-Arithmetic--Kakeya import.  Those frontiers
remain **OPEN**.  The restricted score `11/6` is larger, not smaller, than the
strict `[2,2]` record `7/4`.

## 5. Exact replay

Run from the repository root:

```bash
python3 -B 04-computation/ak_strict_c4_score_seven_four_thm3903.py
python3 -B -O 04-computation/ak_strict_c4_score_seven_four_thm3903.py
python3 -B 04-computation/ak_strict_c4_score_seven_four_independent_audit_thm3903.py
python3 -B -O 04-computation/ak_strict_c4_score_seven_four_independent_audit_thm3903.py
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_positive_control.py
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_positive_control.py
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_positive_control_independent.py
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_positive_control_independent.py
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f4_symbolic_audit.py
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f4_symbolic_audit.py
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f3_exact_audit.py
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f3_exact_audit.py
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f2_exact_audit.py
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_f2_exact_audit.py
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_rank9_f3_audit.py
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_full_rank9_f3_audit.py
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_multiround_rank_variety.py rung3
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_multiround_rank_variety.py rung3
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_multiround_rank_variety.py corner3
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_multiround_rank_variety.py corner3
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_multiround_rank_variety.py corner4
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_multiround_rank_variety.py corner4
python3 -B 04-computation/ak_strict_saturated_2x3_thm3903/ladder_trace_coverage_independent.py
python3 -B -O 04-computation/ak_strict_saturated_2x3_thm3903/ladder_trace_coverage_independent.py
g++ -O2 -std=c++17 04-computation/ak_strict_saturated_2x3_thm3903/ladder_modular_exhaust.cpp -o .scratch/ladder_modular_exhaust_thm3903
./.scratch/ladder_modular_exhaust_thm3903 2 gauge
./.scratch/ladder_modular_exhaust_thm3903 3 gauge
./.scratch/ladder_modular_exhaust_thm3903 5 gauge
```

The first two companions verify the `[2,2]` witness and its independent
60-family no-go.  The saturated-ladder companions verify the `[2,3]` witness,
two-round strata, multiround varieties and independent trace coverage.  Every
listed normal and optimized stream byte-matches its raw-LF frozen output; the
three long multiround commands accept `--start/--stop` for reproducible slices.
The three finite-field invocations concatenate, in the displayed order, to the
pinned hostile transcript; they are not dependencies of the rational proof.
**QED.**
