---
id: THM-3743
title: "Lonely-runner polyhedron Khinchin-flatness relation reduction"
status: >
  CITED + PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED;
  NUMERICAL CAP SUPERSEDED BY THM-4009. This theorem's l1<=356 implication
  remains valid, but THM-4009 now gives ||a||_2<14 and ||a||_1<=50.
  A hypothetical primitive LRC(14) counterexample has a nonzero integer speed
  relation a with ||a||_1<=356.  An l1-minimal such relation is a Graver
  element.  Its support-two branch lies in an exact 19,314-ratio atlas; its
  support-at-least-three branch is a genuinely multiway bounded partition
  identity.  THM-2052's rank-twelve case is already terminal; in its
  rank-eleven case the new row either raises the rank and gives an explicit
  135-digit speed cap, or is already short inside the unresolved space.
  The reduction is necessary, never sufficient, and does not prove LRC(14).
source: root + overnight-jc-lrc-khinchin + lrc14-cover-defect-bridge / 2026-08-23
audit: >
  PASS.  The primary companion checks quotient-lattice duality and projected
  widths on 15 primitive controls and 9,228 projected vertices, computes the
  exact flatness cap, Graver split, pair atlas, ambient coefficient universe,
  rank-eleven Hadamard terminal, and AP hostile in 18,076 gates.  A separately
  implemented referee rechecks every integer constant, all 19,314 reduced
  pair ratios, the labelled count, the triple hostile, the AP boundary time,
  and the 135-digit terminal.  Normal and optimized replays agree and
  line-normalize exactly to the frozen outputs.  Literature use is limited to
  the cited
  lonely-runner zonotope equivalence and the explicit general flatness bound.
depends_on:
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
related:
  - THM-4009-euclidean-covering-transference-short-relation-compression
  - THM-778-centered-christoffel-endpoint-skew-product
  - THM-2051-fejer-bv-small-relation-alternative-for-lrc14
  - THM-2169-bounded-relation-on-every-lrc-deletion
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
script: 04-computation/lrc14_khinchin_flatness_relation_audit_20260823.py
output: 05-knowledge/results/lrc14_khinchin_flatness_relation_audit_20260823.out
script_sha256: 16358a45d3ee6fc6c4ad6a6fb5780e1cf37a245e17ea78a051794abea1a68397
output_sha256: 7e657f35b943704a14b3557120482f8637fcac0b4e6d192db337bde9c33ba8bb
independent_script: 04-computation/lrc14_khinchin_flatness_relation_independent_audit_thm3743.py
independent_output: 05-knowledge/results/lrc14_khinchin_flatness_relation_independent_audit_thm3743.out
independent_script_sha256: aa766e031f2b608fdd91605372dfed99ab35da83e313b88e6156736e451463e8
independent_output_sha256: 3f249fa31348220da03c59865696da27f293c09e67dd08ed4a53f81af939f51b
hash_basis: raw working-tree bytes
---

# THM-3743 -- flatness yields one short exact speed relation

**CITED + PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  The external inputs are the lonely-runner projected
zonotope equivalence of
[Beck--Hosten--Schymura](https://arxiv.org/abs/1606.01783) and the explicit
general bound of Barvinok, as recorded by
[Averkov--Hofscheier--Nill](https://arxiv.org/abs/1911.03511).  Everything
after those two inputs is proved below.  No optimality or literature-priority
claim is made.

> **Numerical-cap update (2026-08-24).** The flatness argument below remains
> correct, but its `l1<=356` output is no longer the sharp proved reduction in
> this repository. [THM-4009](THM-4009-euclidean-covering-transference-short-relation-compression.md)
> uses the same projected zonotope's Euclidean radius-`3/7` inball and
> Banaszczyk covering transference to prove `||a||_2<14`, hence
> `sum a_i^2<=195`, `||a||_1<=50`, and `|a_i|<=13`.

Let `n=(n_1,...,n_13)` be a primitive vector of distinct positive integer
speeds.  Put

```text
V=n^perp,
pi:R^13 -> V                    orthogonal projection,
Lambda=pi(Z^13),
Z(n)=pi([1/14,13/14]^13).                              (1)
```

If `n` were a counterexample to LRC(14), then `Z(n)` would be disjoint from
`Lambda`.  Consequently there is a nonzero `a in Z^13` such that

```text
a dot n=0,                         ||a||_1<=356.         (2)
```

An `l1`-minimal nonzero relation may be chosen conformally indecomposable,
hence is a Graver element of the one-row speed matrix.

## 1. Exact quotient dual and width

The projected lattice has the native dual

```text
Lambda^*=Z^13 intersect n^perp.                         (3)
```

Indeed, for `a in V` one has `<a,pi(z)>=<a,z>`.  Integral pairing with
every `z in Z^13` is therefore equivalent, by the standard basis, to
`a in Z^13`; membership in `V` is exactly `a dot n=0`.

Projection also preserves the pairing with the cube.  Independent endpoint
choice in every coordinate gives, for `a in Lambda^*`,

```text
width_a Z(n)
 =(13/14-1/14) sum_i |a_i|
 =(6/7)||a||_1.                                        (4)
```

This is the exact algebraic core: the abstract flat covector is an ordinary
integer speed relation, and its width is proportional to its `l1` norm.

## 2. Cited flatness bound and the integer cap

The Beck--Hosten--Schymura equivalence makes a hypothetical counterexample's
closed zonotope lattice-point-free, which is stronger than the interior
condition required by flatness.  Khinchin's Flatness Theorem therefore
supplies a nonzero direction in (3) with

```text
(6/7)||a||_1 <= Flt(12).                               (5)
```

The explicit dimension-`d` estimate of Barvinok, as recorded by
Averkov--Hofscheier--Nill, is

```text
Flt(d) <= sqrt((d+1)(2d+1)/6) d^(3/2).                 (6)
```

At `d=12`,

```text
Flt(12)^2 <= (13*25/6)12^3=93600,
Flt(12) <=60 sqrt(26),
||a||_1 <=70 sqrt(26)<357.                             (7)
```

Because `||a||_1` is an integer, (2) follows.

The zonotope's translation and its closed-boundary convention cause no loss:
flatness width is translation invariant, and complete disjointness from the
lattice is more restrictive than absence of interior lattice points.

## 3. The minimum direction is a Graver relation

Choose `a` of minimum `l1` norm among the nonzero vectors in (3).  If

```text
a=b+c
```

with nonzero kernel vectors `b,c` conformal to `a` and coordinatewise
dominated by it, then

```text
||a||_1=||b||_1+||c||_1,
```

contradicting minimality of `a` because each summand is a shorter relation.
Thus `a` is conformally indecomposable.

There are two semantically different branches.

### Support two

For speeds `n_i,n_j` and `g=gcd(n_i,n_j)`, primitiveness forces, up to sign,

```text
(a_i,a_j)=(n_j/g,-n_i/g),
(n_i+n_j)/g=||a||_1<=356.                              (8)
```

There are exactly `19,314` unordered coprime ratios `p<q` with
`p+q<=356`, or `1,506,492` choices after selecting the labelled coordinate
pair.  The `19,314` count is a ratio count, not a placement count.  Each
ratio can be routed through THM-778's exact centered-Christoffel word only
after retaining its gcd scale and the owner/root/phase sidecars.

### Support at least three

This is a primitive signed partition identity of total coefficient mass at
most `356`.  It cannot be replaced by pair ratios without losing its
multiway cancellation.  The exact hostile

```text
3-2(4)+5=0
```

has relation `(1,-2,1)` of norm four on speeds `(3,4,5)`, while every pair
relation is longer.  This branch belongs to relation-code/Fourier carriers,
not a manufactured tournament or a scalar continued fraction.

## 4. Exact join with THM-2052

Let `W` be THM-2052's rational span of support-at-most-three relations of
height `Q=91^6`.  This section inherits THM-2052's conditional dependence on
the settled twelve-speed LRC.  If `dim W=12`, that theorem is already in its
finite terminal.  Suppose `dim W=11` and choose `a` from (2).  Then

```text
a notin W  => dim span(W,a)=12,
a in W     => W already contains a nonzero vector of l1 norm <=356. (9)
```

In the first branch choose eleven independent rows of `W`.  Each has
Euclidean norm at most `sqrt(3)Q`, whereas `||a||_2<=356`.  Cofactors and
Hadamard give

```text
max_i n_i <= floor(356*3^(11/2)*Q^11)                  (10)
```

and the exact terminal

```text
296721347184071259951513500572385227832063299530012809732642802018052529730741866571256016315918699080945038248181476620474893833141605.
```

It has `135` digits.  The second branch is the genuine remaining task:
classify the rank-eleven star spaces that already contain a short dense
Graver vector.

There is a sharp support qualification.  If `a` has support two or three,
then its coefficient height is at most `356<91^6`, so it is itself one of
THM-2052's generators and automatically lies in `W`.  Such a row can refine
the rank-eleven code but cannot be the outside-`W` rank-twelve increment.
Only support at least four can enter either branch of (9), and span membership
must then be tested rather than inferred from support.

The square-triangular Pell/conductor compiler supplies an exact hostile to a
selector shortcut.  Its support-two image meets the cap in only
`(1,5),(5,29),(29,169)`, but all three rows are inside `W`; a safe thirteen-
speed Pell prefix contains them and has loneliness maximum `99/338>1/14`.
Thus Pell alignment is an address refinement, not a badness predicate.

The ambient candidate universe is finite but enormous:

```text
sum_(s=1)^13 2^s C(13,s) C(356,s)
 =1978967793896659449022201064.                        (11)
```

This is an atlas-pruning theorem, not an instruction to enumerate (11)
blindly.

## 5. Hostile, loss ledger, and stopping boundary

The arithmetic progression `(1,2,...,13)` has the short relation
`(1,-2,1,0,...)` and width `24/7`, yet `t=1/14` is a valid loneliness
boundary time.  A short relation is therefore necessary for a hypothetical
counterexample, never sufficient.

The typed connection is

```text
source:      lattice-point-free projected lonely-runner zonotope
target:      integer speed-relation lattice
map:         minimum lattice-width covector
preserved:   counterexample => one bounded exact resonance
destroyed:   support, sign partition, endpoint owner, root/word, arrival
sidecars:    THM-778 for support two; THM-2052/2169 for higher support
cheapest next test:
             short-vector intersection in each rank-eleven star space.      (12)
```

Flatness may be iterated on lattice-free affine slices, but after the first
pull the simple formula (4) is gone: later covectors live modulo previous
directions and the sliced zonotope no longer has independent cube endpoints.
Twelve bounded independent relations cannot be inferred without a
quotient-gauge sidecar controlling lifted ambient coefficients.

## 6. Reproduction

```bash
python3 -B 04-computation/lrc14_khinchin_flatness_relation_audit_20260823.py
python3 -B -O 04-computation/lrc14_khinchin_flatness_relation_audit_20260823.py
python3 -B 04-computation/lrc14_khinchin_flatness_relation_independent_audit_thm3743.py
python3 -B -O 04-computation/lrc14_khinchin_flatness_relation_independent_audit_thm3743.py
```

The computations audit exact constants and finite universes; Sections 1--5
carry the implications.  The theorem does not recover owners, semantic
arrival, simultaneous loneliness, or a contradiction.  LRC(14) remains
open.  **QED.**
