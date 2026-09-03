---
id: THM-4382
title: "LRC14 signed one-four-one comb exact measure and sharp maximum"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED; LRC(14) OPEN. Every primitive distinct positive odd three-unit
  triple satisfying a signed coefficient pattern (1,4,1) has scale-three
  failure-comb measure at most 12/301, uniquely attained by {1,11,43}, with
  coefficient-four carrier 11. No universal nonresonant, seam-entry, ledger,
  or LRC(14) consequence is asserted.
source: root + lrc_nonresonant_next + clean-room referee / next-sharp session, 2026-09-03
depends_on:
  - THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound
related:
  - THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow
  - THM-4370-lrc14-septimal-valuation-pigeonhole-seven-forced-lower-tails
  - THM-4372-lrc14-septimal-depth-transport-and-off-valuation-four-sevenths-rebate
primary_script: 04-computation/lrc14_signed_one_four_one_comb_exact_measure_thm4382.py
primary_output: 05-knowledge/results/lrc14_signed_one_four_one_comb_exact_measure_thm4382.out
primary_script_sha256: f4411755a8abd6ed7f025f3e2f3c1d26dad3a7c572daf0cfc41c9ebafcaf87ff
primary_output_sha256: 53b86fbda179b34313e96697975161696cbd6ebea3aafe6115f6a78624f07126
independent_referee_script: 04-computation/lrc14_signed_one_four_one_comb_exact_measure_independent_referee_thm4382.py
independent_referee_output: 05-knowledge/results/lrc14_signed_one_four_one_comb_exact_measure_independent_referee_thm4382.out
independent_referee_script_sha256: 7915d60b4a244278ac7813b9c37b00fe8e0413905013eebd1c405d78c186c915
independent_referee_output_sha256: 54dcf24e58816a7633fc073348df8c02667f8e2257f29abdf990f34462e45072
hash_basis: raw LF bytes
audit: >
  PASS. The 116,699-check clean-room verifier independently reconstructs the
  relation classification, labeled-owner defect elimination, exact component
  and quadrature formulas, all 31 small-product cases, unique equality, and
  disjointness from the signed-(1,2,1) sector. Its full-x-circle route agrees
  with the structurally different primary verifier; normal, optimized, and
  hash-seeded replays match their frozen streams.
---

# THM-4382 -- Signed-(1,4,1) scale-three comb exact measure and sharp maximum

**PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THE SHARP SECTOR MAXIMUM IS `12/301` AT `{1,11,43}`. THE
UNIVERSAL NONRESONANT SECTOR, SEAM ENTRY, AND LRC(14) REMAIN OPEN.**

For a positive three-unit `w` and `j in Z/3Z`, retain THM-4373's
physical danger sheet

```text
D_(w,j)={x in R/Z: ||w(x+j/3)||<1/14},
F_(a,b,c)=disjoint union over pi in S_3 of
          D_(a,pi(0)) intersect D_(b,pi(1)) intersect D_(c,pi(2)),
```

where the disjoint-union wording ignores only measure-zero endpoints.

The exact statement is:

> For every primitive triple of distinct positive odd integers, each prime to
> three, which after a permutation and sign choices satisfies a signed
> coefficient pattern `(1,4,1)`, the scale-three failure comb has Haar measure
> at most `12/301`.  Equality occurs only for the unordered triple
> `{1,11,43}`, with coefficient-four speed `11`.

This is a proper subfamily of the signed-`(1,2,1)`-nonresonant sector.  It is
not the universal nonresonant bound, a seam-entry theorem, or LRC(14).

## 1. Inheritance and live concepts

- Closest proved mechanism: THM-4373,
  `01-canon/theorems/THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound.md`,
  supplies the scale-three owner quotient, determinant components, and
  period-three quadrature.
- Canonical hostile/near miss: THM-4373 records `(1,11,43)` only as the finite
  nonresonant maximum through height `199`; that finite census is explicitly
  not a universal theorem.
- Least-used relevant sidecar: THM-4348,
  `01-canon/theorems/THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow.md`,
  warns that a short relation or wall shadow does not provide physical entry,
  global owner transport, or an LRC closure.
- Concept board: the integer defect `delta`; the determinant `K`; the
  period-three error `E`; and the separation from the `(1,2,1)` sector.

No matching correction to THM-4373/4382, the failure-comb statement, or the
constant `12/301` appears in `01-canon/MISTAKES.md` as of this audit.

## 2. Exact relation classification

Designate the coefficient-four speed by `b` and sort the coefficient-one
speeds as `p<q`.  Positivity leaves exactly

```text
4b=q+epsilon p,        epsilon=+1 (mean) or -1 (difference).       (1)
```

All odd three-units square to one modulo `24`.  In the mean branch, oddness of
`b` is `p+q=4 mod 8`, while `3` not dividing `b` forces `q=p mod 3`.
Multiplication by `p` and CRT therefore give `pq=19 mod 24`, equivalently

```text
q=19p mod 24,          b=(p+q)/4.                              (2a)
```

In the difference branch, the corresponding conditions are `q-p=4 mod 8`
and `q=-p mod 3`, hence `pq=5 mod 24`, equivalently

```text
q=5p mod 24,           b=(q-p)/4.                              (2b)
```

The converses are immediate.  The two congruences cannot both hold for a
three-unit `p`.  Also

```text
gcd(p,b,q)=gcd(p,q),                                      (3)
```

because the common endpoint gcd is odd and divides `4b`.  Thus primitivity is
exactly `gcd(p,q)=1`.  Distinctness only removes the difference presentation
`(p,q)=(1,5)` from the primitive parametrization: generally `b=p` there is
equivalent to `q=5p`; `b=q` is impossible for `p<q`.  A second coefficient-four
carrier is also impossible under the distinct three-unit hypotheses.  Indeed,
if `p` were a second carrier, `4p=q+tau b` for `tau=+-1`; combining with `(1)`
gives

```text
(16-tau epsilon)p=(4+tau)q,
```

whose four cases force a repeated speed or one of `p,b,q` to be divisible by
three.  The carrier `q` is excluded because `q` is strictly the largest speed.

Consequently `(2a)--(3)`, `(2b)--(3)`, and distinctness give an exact,
nonduplicated classification of the primitive signed-`(1,4,1)` triples.

## 3. The owner congruence kills the integer defect

As in THM-4373, pass from the physical phase `x` to `y=3x mod 1` and put

```text
r=3/14.
```

For an eligible speed `w`, let `n_w` be the unique nearest integer and
`e_w=wy-n_w`, so `|e_w|<r`.  From `(1)` define

```text
delta=n_q+epsilon n_p-4n_b,
4e_b=e_q+epsilon e_p+delta.                              (4)
```

If all three speeds are eligible, strict inequalities give

```text
|delta|<6r=9/7,
```

so `delta` is one of `-1,0,1`.  Define the exact owner colour

```text
o_w=-w^(-1)n_w mod 3.                                   (5)
```

If the three colours are distinct, relative to `o_b` the other two colours
differ by opposite nonzero residues.  Reducing `(4)` modulo three and using
`b=q+epsilon p mod 3` gives

```text
delta = q(o_b-o_q)+epsilon p(o_b-o_p) mod 3.
```

But `(2)` gives `q=epsilon p mod 3`, so the right side vanishes.  The analytic
three-value range therefore forces

```text
delta=0.                                                (6)
```

Now `n_q+epsilon n_p` is divisible by four.  Since `(1)` also gives
`q=-epsilon p mod 4`, the endpoint determinant

```text
K=q n_p-p n_q                                           (7)
```

satisfies

```text
K=-p(n_q+epsilon n_p)=0 mod 4,
K=4k.                                                   (8)
```

This implication is exact in the other direction.  Endpoint eligibility and
`4|K` are equivalent to `4|(n_q+epsilon n_p)`.  Set that quotient equal to
`n_b`; then

```text
e_b=(e_q+epsilon e_p)/4,
```

so `|e_b|<r/2<r`, making `n_b` the middle nearest integer.  Furthermore

```text
o_q-o_p=K/(pq) mod 3.                                   (9)
```

Under `(6)`, the speed and nearest-integer relations also give
`o_p+o_b+o_q=0 mod 3`.  Thus the endpoint owners differ if and only if all
three owners are distinct, and `(8)--(9)` yield the exact owner gate

```text
K=4k,                  3 does not divide k.             (10)
```

This is the load-bearing new step.  Merely bounding `delta` does not kill it;
the owner permutation does.

## 4. Exact component and quadrature formulas

Primitivity gives `gcd(p,q)=1`.  For each integer `k`, all solutions of

```text
q n_p-p n_q=4k                                         (11)
```

form one orbit under `(n_p,n_q)->(n_p+p,n_q+q)`, exactly translation of `y`
by one.  Hence each `k` produces at most one circle component.  The endpoint
intervals have radii `r/p,r/q`, while their centre separation is
`4|k|/(pq)`.  Therefore

```text
lambda_k=max(0,min(2r/q, r/p+r/q-4|k|/(pq))).           (12)
```

Put

```text
A=3(q-p)/56,                 B=3(q+p)/56,
f(t)=4/(pq) ((B-t)_+-(A-t)_+).                          (13)
```

Then `(12)` is exactly `lambda_k=f(|k|)`.  The class `k=0` fails the owner
condition, and positive/negative determinants have equal lengths.  Thus the
exact measure is

```text
mu(F_(p,b,q))=2 sum_(k>=1, 3 does not divide k) f(k).    (14)
```

For `t>=0`, retain THM-4373's elementary period-three error

```text
E(t)=sum_(k>=1,3 does not divide k)(t-k)_+ - t^2/3.
```

Writing `rho=t mod 3`, `0<=rho<3`, direct summation gives

```text
E(t)= -rho^2/3                    (0<=rho<=1),
      rho-1-rho^2/3               (1<=rho<=2),
      2rho-3-rho^2/3              (2<=rho<3).           (15)
```

It is period three and `-1/3<=E<=0`.  Since
`B^2-A^2=9pq/784`, substituting `(15)` into `(14)` gives

```text
mu(F_(p,b,q))
 =3/98 + 8/(pq)(E(B)-E(A))
 <=3/98 + 8/(3pq).                                     (16)
```

The integral check is `integral_0^infinity f(t)dt=9/392`.

## 5. Large products and every small case

At `pq=289`, the ceiling in `(16)` is

```text
3385/84966 < 12/301
```

with gap `101/3653538`; the ceiling decreases with `pq`.  Thus only `pq<289`
requires enumeration.  Here `p<17`.  Applying `(2)--(3)` and distinctness
leaves exactly the following 31 presentations; every measure is from `(14)`
and independently from a definition-level rational circle-cell union.

| `p` | `b` | `q` | branch | `mu(F)` |
|---:|---:|---:|:---|---:|
| 1 | 5 | 19 | mean | `4/133` |
| 1 | 7 | 29 | difference | `6/203` |
| 1 | 11 | 43 | mean | `12/301` |
| 1 | 13 | 53 | difference | `12/371` |
| 1 | 17 | 67 | mean | `12/469` |
| 1 | 19 | 77 | difference | `18/539` |
| 1 | 23 | 91 | mean | `18/637` |
| 1 | 25 | 101 | difference | `24/707` |
| 1 | 29 | 115 | mean | `24/805` |
| 1 | 31 | 125 | difference | `24/875` |
| 1 | 35 | 139 | mean | `30/973` |
| 1 | 37 | 149 | difference | `32/1043` |
| 1 | 41 | 163 | mean | `36/1141` |
| 1 | 43 | 173 | difference | `36/1211` |
| 1 | 47 | 187 | mean | `40/1309` |
| 1 | 49 | 197 | difference | `6/197` |
| 1 | 53 | 211 | mean | `48/1477` |
| 1 | 55 | 221 | difference | `48/1547` |
| 1 | 59 | 235 | mean | `48/1645` |
| 1 | 61 | 245 | difference | `54/1715` |
| 1 | 65 | 259 | mean | `54/1813` |
| 1 | 67 | 269 | difference | `60/1883` |
| 1 | 71 | 283 | mean | `60/1981` |
| 5 | 7 | 23 | mean | `4/115` |
| 5 | 13 | 47 | mean | `12/329` |
| 5 | 11 | 49 | difference | `12/343` |
| 7 | 1 | 11 | difference | `0` |
| 7 | 5 | 13 | mean | `4/637` |
| 7 | 11 | 37 | mean | `62/1813` |
| 11 | 7 | 17 | mean | `4/187` |
| 13 | 1 | 17 | difference | `2/91` |

The unique maximum is `12/301` at `(p,b,q)=(1,11,43)`.  Together with the
large-product ceiling, this proves the infinite sharp statement.

## 6. Disjointness from signed `(1,2,1)`

Assume `(1)` and suppose a signed `(1,2,1)` relation also exists.  Since
`b<q` and `p<q`, the coefficient-two carrier cannot be `q`.

If it is `b`, write `2b=q+tau p`, `tau=+-1`.  Comparing twice this equation
with `(1)` gives

```text
q=(epsilon-2tau)p.
```

The four cases make `q` negative, equal to `p`, or divisible by three.

If the coefficient-two carrier is `p`, write `2p=q+tau b`.  Substitution of
`4b=q+epsilon p` gives

```text
(8-tau epsilon)p=(4+tau)q.
```

For `tau=+1`, the ratios are `q/p=7/5` or `9/5`, forcing respectively
`3|b` or `3|q`.  For `tau=-1`, the ratios are `3` or `7/3`, forcing
`3|q` or `3|p`.  All contradict the three-unit hypothesis.  Hence the signed
`(1,4,1)` and signed `(1,2,1)` sectors are disjoint.  This disjointness does
not make the former equal to the full nonresonant sector.

## 7. Exact LRC scope

The theorem controls only the scale-three failure comb for three tails in the
classified relation sector.  A valid conditional body corollary is:

> If a finite positive core `C` has nonempty safe set
> `G_C={y: ||cy||>=1/14 for all c in C}` and
> `mu(G_C)>=12/301`, then `3C union {p,b,q}` is `1/14`-lonely for every
> classified signed-`(1,4,1)` triple.

Indeed, the body-safe preimage under `x->3x` is nonempty compact with the same
measure, while the strict failure comb is a proper open subset of the circle.
Strict measure separation forbids containment; at equality, a nonempty
compact subset cannot equal a proper open subset of the connected circle, and
the nonempty open difference has positive measure.

What does **not** follow:

- no lower bound `mu(G_C)>=12/301` is supplied for every relevant core;
- no arbitrary nonresonant triple is classified;
- no hypothetical LRC(14) counterexample is forced to contain this relation;
- no THM-4348 wall-shadow, renewal, or entry obligation is discharged;
- no universal nonresonant bound, seam-entry theorem, ledger decrement, or
  LRC(14) proof follows.

Common dilation by an odd three-unit preserves Haar measure and permutes sheet
labels, so a nonprimitive extension is available, with equality exactly at
dilates of `{1,11,43}`.  The audited theorem itself needs only the stated
primitive scope.

## 8. Independent verification and frozen result

The independent referee uses four mutually checking exact routes:

1. the determinant series `(14)`;
2. the period-three quadrature `(16)`;
3. owner-labelled exact `y`-circle cells, including the full `delta/K/k` iff;
4. the original `x`-circle definition, unioned over all six sheet assignments.

It exhausts all `47,499` primitive distinct odd three-unit triples through
speed `199`, finding `515` signed-`(1,4,1)` triples, `1,023`
signed-`(1,2,1)` triples, and zero overlap.  It also checks `9,900` endpoint
congruence iff statements, `13,381` owner cells, all 31 small-product cases,
and eight larger controls.  Total explicit checks: `116,699`.

The primary verifier separately checks 782 formula/direct/permutation rows,
50 complete wall-cell controls, 7,512 literal owner/determinant cells, the
31-case table, common-dilation and same-endpoint hostiles, and all 45,961
remaining primitive triples through height 199 after the two closed relation
sectors. That last residual census is finite evidence only.

Reproduction from the repository root:

```powershell
python -B 04-computation/lrc14_signed_one_four_one_comb_exact_measure_thm4382.py
python -O -B 04-computation/lrc14_signed_one_four_one_comb_exact_measure_thm4382.py
python -B 04-computation/lrc14_signed_one_four_one_comb_exact_measure_independent_referee_thm4382.py
python -O -B 04-computation/lrc14_signed_one_four_one_comb_exact_measure_independent_referee_thm4382.py
```

Both modes match their committed output streams exactly. The raw-LF SHA-256
values are recorded in the front matter. The referee's proof, script, and
output were frozen before the primary artifact was opened; post-freeze
comparison found full agreement and no scope discrepancy.

**QED.**
