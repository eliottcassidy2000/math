---
id: THM-4383
title: "LRC14 signed two-five-one comb exact measure and sharp maximum"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED; LRC(14) OPEN. Every primitive distinct positive odd three-unit
  triple satisfying a signed coefficient pattern (2,5,1) has scale-three
  failure-comb measure at most 12/371, uniquely attained by {1,11,53}, with
  coefficient-five carrier 11. The gcd-five endpoint seam requires a retained
  middle-integrality sidecar. No universal nonresonant, seam-entry, ledger,
  all-tail, or LRC(14) consequence is asserted.
source: root + lrc251_probe + clean-room referee / next-sharp session, 2026-09-03
depends_on:
  - THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound
related:
  - THM-4382-lrc14-signed-one-four-one-comb-exact-measure-and-sharp-maximum
  - THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow
  - THM-4370-lrc14-septimal-valuation-pigeonhole-seven-forced-lower-tails
  - THM-4372-lrc14-septimal-depth-transport-and-off-valuation-four-sevenths-rebate
primary_script: 04-computation/lrc14_signed_two_five_one_comb_exact_measure_thm4383.py
primary_output: 05-knowledge/results/lrc14_signed_two_five_one_comb_exact_measure_thm4383.out
primary_script_sha256: 4edcec46a605098184c9582bf7789cb97d6a5880e2f3b6a4cf6c952d30ee6d9d
primary_output_sha256: 7b12211f2b00a8f5223c8773a67e8c89d0d81078461f403b784186c2ab5dfd0e
independent_referee_script: 04-computation/lrc14_signed_two_five_one_comb_exact_measure_independent_referee_thm4383.py
independent_referee_output: 05-knowledge/results/lrc14_signed_two_five_one_comb_exact_measure_independent_referee_thm4383.out
independent_referee_script_sha256: 9278e8ebf724ba5b4bb9acde8c9eb112c9a3e9430913a36ca8e01fd7310ab60b
independent_referee_output_sha256: 51392fcb89fa5f3e41c0762614f352dc8c2b7c9dba575bb7536b6dd1e207f3a2
hash_basis: raw LF bytes
audit: >
  PASS. The 971,507-check clean-room verifier independently reconstructs the
  sign and gcd classification, labeled-owner defect elimination, cyclic
  quotient with its gcd-five middle-integrality sidecar, exact component and
  quadrature formulas, all 75 small-product presentations, unique equality,
  and scale-independent coefficient-sector disjointness. It checks 97,809
  exact full-x wall cells. Its explicit checks remain live under optimized
  Python; normal, optimized, and hash-seeded replays match the frozen stream.
  The primary's optimized replay is determinism evidence only because that
  producer uses Python assert.
---

# THM-4383 -- Signed-(2,5,1) scale-three comb exact measure and sharp maximum

**PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THE SHARP SECTOR MAXIMUM IS `12/371` AT `{1,11,53}`. THE
UNIVERSAL NONRESONANT SECTOR, SEAM ENTRY, ALL-TAIL TRANSFER, AND LRC(14)
REMAIN OPEN.**

## Outcome

The class suggested by the height-199 survivor `(1,11,53)` does **not** remain
rank two.  The three labelled physical owners force the only possible middle
nearest-integer defect to be zero.  After this forcing, the integer lift
lattice modulo physical translation is infinite cyclic, with one determinant
coordinate.  The exact Haar formula is a period-three quadrature, and its
sharp primitive maximum is

```text
mu(F_(p,b,q)) <= 12/371,
```

with equality only at the coefficient-labelled presentation

```text
(p,b,q)=(1,11,53),       5b=2p+q.
```

Equivalently, after forgetting the coefficient labels, equality occurs only
for the speed set `{1,11,53}`.  Without primitive normalization, the equality
sets are its common odd three-unit dilates.

More precisely:

> **Theorem.**  Let `p,b,q` be distinct positive odd
> integers prime to three, with `gcd(p,b,q)=1`.  Suppose that for some
> `sigma,tau in {+1,-1}`,
>
> ```text
> 5b=2 sigma p+tau q > 0.                              (1)
> ```
>
> For the physical scale-three failure comb `F_(p,b,q)` defined in THM-4373,
>
> ```text
> mu(F_(p,b,q)) <= 12/371.                             (2)
> ```
>
> Equality holds exactly for `(p,b,q)=(1,11,53)` with
> `(sigma,tau)=(+1,+1)`.

Every signed relation `+/-2p +/-5b +/-q=0` has the form `(1)` after an
overall sign change.  A finite coefficient-pair audit also proves that a
distinct odd three-unit triple has at most one signed-(2,5,1) presentation,
and that this sector is disjoint from the signed-(1,2,1) and signed-(1,4,1)
sectors of THM-4373 and THM-4382.

## Inheritance pass and concept board

Closest proved mechanism: THM-4373,
`01-canon/theorems/THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound.md`,
where labelled owners turn the signed-(1,2,1) middle relation into one
determinant and period-three quadrature.

Canonical hostile: `(1,11,53)`, the inherited exact height-199 maximum after
removing the signed-(1,2,1) and signed-(1,4,1) sectors, with measure `12/371`.

Corrected near miss: relation plus triple eligibility does not by itself make
the middle integer defect zero.  At `y=1/53`, all three speeds `(1,11,53)` are
eligible, but `(n_1,n_11,n_53)=(0,0,1)`, the defect is `-1`, and the owners are
`(0,0,1)`.  It is the **distinct labelled-owner** condition that kills the
defect.

Least-used relevant sidecar: the integral middle defect before owner
scalarization.  The endpoint determinant is introduced only after that defect
has been forced to zero.

The live concept board was:

| object | preserved predicate / operation | lost coordinate / hostile |
|---|---|---|
| physical lifts | three distinct labelled owners under `y=3x` | a common sheet shift is harmless, individual owner labels are not |
| signed `(2,5,1)` relation | exact middle defect `delta` | interval geometry alone permits `delta=+/-1` |
| defect-zero lift lattice | quotient by `y -> y+1` | endpoint gcd may be five, so quotienting endpoints alone is unsafe |
| determinant classes | residue `k mod 3` and exact overlap length | `k=0 mod 3` has owner collision |
| coefficient patterns | exact signed vectors | scalar Haar measure forgets which speed owns coefficient 2, 5, or 1 |

Meta-patterns used: preserve the selected side rather than only the walls;
attack a proposed bound before extending it; use redundant paths as detectors;
and test structured adversaries.  The source is the physical comb, the target
is a cyclic determinant lattice, the preserved predicate is literal owner
distinctness, and the discarded coordinate is the common physical sheet.

## 1. Physical owners and the exact middle defect

Use THM-4373's definitions

```text
D_(w,j)={x in R/Z: ||w(x+j/3)||<1/14},
F_(a,b,c)=disjoint union over pi in S_3 of
          D_(a,pi(0)) intersect D_(b,pi(1)) intersect D_(c,pi(2)).
```

Put

```text
y=3x mod 1,                  r=3/14.
```

For an eligible speed `w`, let `n_w` be the unique nearest integer to `wy`,
and write

```text
e_w=wy-n_w,                 |e_w|<r,
o_w=-w^(-1)n_w mod 3.                                  (3)
```

The `o_w` are the actual lift owners up to the same common sheet shift for all
three speeds.  Thus physical failure is exactly eligibility of all three
speeds plus

```text
{o_p,o_b,o_q}=Z/3Z.                                    (4)
```

For `(1)`, define

```text
delta=5n_b-2 sigma n_p-tau n_q.                        (5)
```

The speed relation gives

```text
5e_b=2 sigma e_p+tau e_q-delta,
|delta|<8r=12/7.
```

Therefore

```text
delta in {-1,0,1}.                                     (6)
```

Let `p=rho mod 3`, where `rho in {+1,-1}`.  Reducing `(1)` modulo three and
using that all speeds are units forces

```text
b=-sigma rho mod 3,          q=-sigma tau rho mod 3.   (7)
```

The sum of the three owners is

```text
S=rho(-n_p+sigma n_b+sigma tau n_q) mod 3,
delta=-sigma rho S mod 3.                              (8)
```

Three distinct residues modulo three have sum zero.  Hence `(4)` implies
`delta=0 mod 3`; the strict range `(6)` upgrades this to the integer identity

```text
delta=0.                                               (9)
```

This is the decisive owner-defect forcing.  Conversely, if the two endpoint
errors satisfy `|e_p|,|e_q|<r` and `(9)` holds, then

```text
e_b=(2 sigma e_p+tau e_q)/5,
|e_b|<3r/5<r,                                          (10)
```

so the middle tail is automatically eligible.

## 2. The owner-defect lattice is cyclic

On the defect-zero lattice define

```text
k=b n_p-p n_b.                                         (11)
```

Because `q=tau(5b-2 sigma p)`, one has

```text
gcd(p,b)=gcd(p,b,q)=1.                                 (12)
```

After using `(9)` to recover

```text
n_q=tau(5n_b-2 sigma n_p),
```

the integer lattice is just `Z^2` with coordinates `(n_p,n_b)`.  The map
`(n_p,n_b) -> k` is surjective by `(12)`, and its kernel is exactly
`Z(p,b)`.  Consequently

```text
{delta=0 integer lifts}/Z(p,b,q)  is isomorphic to Z,  (13)
```

with no torsion.  Thus the class collapses to one determinant even in the
exceptional endpoint stratum `gcd(p,q)=5`.

Indeed, the endpoint determinant is

```text
q n_p-p n_q=5 tau k.                                   (14)
```

From `(7)`, the difference `o_b-o_p` is nonzero exactly when `k` is nonzero
modulo three.  Under `(9)`, the owner sum is already zero, so two unequal
owners force the third residue to be the missing one.  Hence the complete
physical criterion is

```text
delta=0,       |e_p|,|e_q|<r,       3 does not divide k. (15)
```

The endpoint gcd has only the values one and five: every common divisor of
`p,q` divides `5b`, while primitivity excludes a common divisor with `b`.
The smallest primitive `gcd(p,q)=5` example is `(p,b,q)=(5,7,25)`.  It shows
why an endpoint-only Bezout argument would be false.  Formula `(13)`, which
uses the primitive pair `(p,b)`, handles both gcd strata uniformly.  Literal
physical-wall checks at `(5,7,25)`, `(5,13,55)`, and `(55,23,5)` agree with
the formula below.

At the sharp tuple the lattice is especially explicit:

```text
(n_1,n_11,n_53)
  =t(1,11,53)+k(0,-1,-5),              t,k in Z,       (16)
```

and failure retains exactly the classes `3 does not divide k`.

## 3. Exact Haar formula

For fixed `k`, the endpoint interval centres are separated by

```text
5|k|/(pq).                                             (17)
```

Their radii are `r/p` and `r/q`.  Define

```text
A=3|q-p|/70,                 B=3(p+q)/70,
f(t)=5/(pq) [(B-t)_+-(A-t)_+].                         (18)
```

Then the exact intersection length of the determinant-`k` component is
`f(|k|)`.  For each `k`, `(13)` gives one orbit under translation by one in
`y`.  Different `k` cannot overlap because the nearest-integer vector is
unique.  Positive and negative determinant classes have equal length, while
`k=0 mod 3` fails the owner condition.  Therefore

```text
mu(F_(p,b,q))
 =2 sum_(k>=1, 3 does not divide k) f(k).              (19)
```

Let the period-three quadrature error be

```text
E(t)=sum_(k>=1, 3 does not divide k)(t-k)_+ - t^2/3.   (20)
```

Exact summation in the two retained residue classes gives `E(t+3)=E(t)` and

```text
-1/3 <= E(t) <= 0.                                    (21)
```

Substituting `(18)` into `(19)` yields the exact closed formula

```text
mu(F_(p,b,q))
 =6/245 + 10/(pq) [E(B)-E(A)].                        (22)
```

This formula depends only on the coefficient-two and coefficient-one
endpoints `p,q`; all three sign branches have the same formula.  It is not an
unlabelled quotient statement: the proof used the physical middle speed and
all three owners to obtain `(9)` and `(15)` first.

For `(p,b,q)=(1,11,53)`,

```text
A=78/35,       B=81/35.
```

Only `k=+/-1,+/-2` contribute.  Each component is the entire narrower
endpoint interval and has length `3/371`, so

```text
mu(F_(1,11,53))=4*(3/371)=12/371.                     (23)
```

## 4. Sharp global bound inside the signed-(2,5,1) sector

From `(21)--(22)`,

```text
mu(F_(p,b,q)) <= 6/245+10/(3pq).                      (24)
```

For `pq>=425`, this is strictly below `12/371`; at the boundary the gap is

```text
12/371-(6/245+10/(3*425))=8/662235.                   (25)
```

It remains to enumerate `pq<425`.  The verifier exhausts every ordered odd
three-unit endpoint pair in that range, forms every positive integral

```text
b=(2p+q)/5       or       b=|2p-q|/5,
```

and applies distinctness and primitive filters.  Exactly `75` presentations
remain.  Exact evaluation of `(19)` and `(22)` gives the unique maximum
`(1,11,53)` with `12/371`.  The next two are the endpoint-gcd-five controls

```text
(55,23,5), (5,13,55),       each with 12/385.
```

Thus `(2)` is sharp and its equality case is unique in primitive coefficient
labels.

## 5. Duplicate and earlier-sector collision audit

A signed relation is a coefficient vector obtained by permuting `(5,2,1)`
and choosing signs, modulo overall sign.  There are exactly `24` such vectors.
If a positive speed triple admitted two presentations, its primitive ray
would be the cross product of two distinct coefficient vectors.  The exact
audit finds `45` positive primitive rays.  Every one has a repeated entry, an
even entry, or an entry divisible by three.  Hence none is a distinct odd
three-unit triple.

The identical cross-product audit between signed `(5,2,1)` vectors and the
signed `(2,1,1)` and `(4,1,1)` vectors finds respectively `39` and `51`
positive primitive rays.  Again none is a distinct odd three-unit ray.  This
proves:

1. the signed-(2,5,1) coefficient presentation is unique in the stated
   universe;
2. the class is disjoint from the signed-(1,2,1) sector of THM-4373; and
3. it is disjoint from the signed-(1,4,1) sector of THM-4382.

This is an exhaustive finite coefficient-pattern proof, not a bounded speed
census.

## 6. Hostiles, failure boundary, and scope

Two minimal structural hostiles identify the necessary sidecars:

1. **Owner forcing is necessary.**  At the first positive `q`-tooth centre
   `y=1/53` for `(1,11,53)`, all three tails are eligible, but
   `(n_1,n_11,n_53)=(0,0,1)`, `delta=-1`, and owners `(0,0,1)` collide.
   Relation plus eligibility alone does not force the defect to vanish.
2. **Endpoint coprimality is false.**  `(5,7,25)` is the lexicographically
   first primitive positive signed-(2,5,1) three-unit presentation with
   `gcd(p,q)=5`.  The correct primitive coordinates are `(p,b)`, not `(p,q)`.

There is no counterexample to the theorem in its stated universe.
Its sharp equality case is a positive control, and the two hostiles above are
counterexamples only to tempting proof shortcuts.

The result closes one infinite short-relation sector of the scale-three
triple-comb problem.  It does not classify arbitrary nonresonant triples,
where the nearest-integer lattice modulo `Z(a,b,c)` is generally rank two.  It
also scalarizes away the body-component address, prefix, current, first-exit,
and renewal data.  Consequently it supplies no scale-three all-tail transfer
and does not prove LRC(14).

The next genuinely open comb lane is the first remaining short coefficient
pattern after removing signed `(1,2,1)`, `(1,4,1)`, and `(2,5,1)`; it should
be selected by a fresh exact census rather than inferred from scalar endpoint
data.

## Reproduction

From the repository root:

```powershell
python -B 04-computation/lrc14_signed_two_five_one_comb_exact_measure_thm4383.py
python -B -O 04-computation/lrc14_signed_two_five_one_comb_exact_measure_thm4383.py
python -B 04-computation/lrc14_signed_two_five_one_comb_exact_measure_independent_referee_thm4383.py
python -B -O 04-computation/lrc14_signed_two_five_one_comb_exact_measure_independent_referee_thm4383.py
```

The primary uses exact rational arithmetic and performs:

- a literal physical-wall reconstruction of `(1,11,53)` across `385` cells,
  checking physical sheet labels against nearest-integer owners;
- the determinant sum and independent quadrature formula;
- all `75` small-product presentations needed after the analytic cutoff;
- `21` direct physical-wall controls through speed `45`, covering the sum,
  `2p-q`, and `q-2p` sign branches;
- three explicit endpoint-gcd-five physical controls;
- the two hostile examples above; and
- the complete signed-coefficient collision audits (`45`, `39`, and `51`
  positive rays).

The independent referee performs 971,507 explicit runtime checks, including
97,809 full-`x` wall cells and all 75 small-product presentations. It retains
the gcd-five hostile

```text
(p,b,q,x)=(5,7,25,13/1050),
K=q n_p-p n_q=-5,          2n_p+n_q=1,
```

which disproves an endpoint-only determinant gate while leaving the theorem's
coprime-`(p,b)` quotient intact. Normal and optimized referee modes both keep
all checks live and reproduce the frozen stream. The primary uses bare Python
assertions, so its optimized equality is only a determinism check.

Raw-LF artifact hashes are

```text
primary script:   4edcec46a605098184c9582bf7789cb97d6a5880e2f3b6a4cf6c952d30ee6d9d
primary output:   7b12211f2b00a8f5223c8773a67e8c89d0d81078461f403b784186c2ab5dfd0e
referee script:   9278e8ebf724ba5b4bb9acde8c9eb112c9a3e9430913a36ca8e01fd7310ab60b
referee output:   51392fcb89fa5f3e41c0762614f352dc8c2b7c9dba575bb7536b6dd1e207f3a2
```

This proves only the signed `(2,5,1)` scale-three three-tail comb theorem.
The universal nonresonant triple, seam-entry, body-safe-set, ledger, all-tail,
and LRC(14) problems remain **OPEN**.

QED.
