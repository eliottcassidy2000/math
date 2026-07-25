---
id: THM-2309
title: "Owner-aligned pivot packets and the visible-height separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every chosen
  blocker owner and every omitted unit have an exact rank-six relation packet
  with that owner plus the other five units as a thirteen-unit pivot. Two
  exact target grafts make all nine columns bright modulo thirteen; the row
  heights are at most the scalar l1 height. The full mod-thirteen relation
  code then splits as this packet plus the two-dimensional target-blocker
  plane. In an arbitrary rank-six packet, one-blocker plus five-unit pivots
  are classified exactly by the cocircuit outside the unit hyperplane, and a
  bright representable packet can miss the selected owner completely. If at
  least two scalar coordinates are seven-units, a CRT lift gives an
  owner-aligned packet bright at both primes and hence a large all-91-unit
  exact-address bank, but at row-dependent Bezout height. No theorem
  identifies the algebraic target-plane support with the analytic terminal
  word or forces these packets into the degree-526 visible span, so no scalar
  row is excluded and LRC(14) remains open.
source: codex-2026-07-25-owner-aligned-pivot
depends_on:
  - THM-2301-essential-affine-arrangement-and-visible-rank-six-address-bank
related:
  - THM-2287-anchored-scalar-rank-six-plucker-flag-and-finite-label-atlas
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2307-dual-rank-six-reconstruction-spectrum-and-selector-no-go
  - THM-2305-canonical-blocker-word-handoff-hypergraph
script: 04-computation/lrc14_owner_aligned_pivot_thm2309.py
output: 05-knowledge/results/lrc14_owner_aligned_pivot_thm2309.out
script_sha256: a17a015f9187b59a36cef541100ae29a83395690294148ddee6ec5a45e0ea889
output_sha256: f487850adfe0e0c5cc1f71efabf85f5e8c150bc42453b08f6d87afc41a322b4b
hash_basis: working-tree bytes (LF)
---

# THM-2309 -- owner-aligned pivots exist, but not at visible height

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2301 proves a coefficient-independent rank-six packet, but its pivot
atlas does not prescribe which blocker occurs in a unit minor. THM-2305
canonically selects a source owner `j` and retains both remaining blockers
as its possible terminal word. The apparently missing matroid incidence is

```text
pivot = selected source j + five guard/unit columns,
complement = both target blockers + one omitted unit.             (1)
```

There are only six such pivots for a fixed source. This theorem separates
two questions which the bounded atlas had merged:

```text
Does the full exact relation lattice contain (1)?          yes, explicitly;
does the degree-526 visible packet have to contain (1)?    no.             (2)
```

The positive construction can also be made bright at both seven and
thirteen, subject to the sharp scalar septimal-support condition, but its
height is row-dependent. It therefore supplies genuine Kakeya/Fourier
needles without landing them on the prescribed handoff current.

## 1. Scalar notation

Let `E=U disjoint_union B`, where

```text
U={u_0,u_1,...,u_5},             |U|=6,
B={j,a,b},                       |B|=3.              (3)
```

The six labels in `U` are the guard and five unit speeds. The three labels
in `B` are the blocker speeds; `j` is the selected source owner and `a,b`
are its two possible targets. Let

```text
w=(w_e)_(e in E) in Z^E,
Lambda(w)={x in Z^E:x.w=0},                          (4)
```

where every `w_e` is nonzero and

```text
13 does not divide w_u,             u in U,
13 divides w_t,                     t in B.          (5)
```

Put

```text
S=sum_(e in E)|w_e|.                                 (6)
```

These are exactly the mod-thirteen scalar types on every remaining
first-depth-one row. No cover hypothesis is needed for the algebra below.

## 2. An explicit exact owner pivot, with both targets bright

Fix any omitted unit `u_0`. Put

```text
P={j} union (U minus {u_0}).                         (7)
```

For every `k in P`, define the exact pair relation

```text
r_k=w_k e_(u_0)-w_(u_0)e_k in Lambda(w).             (8)
```

Choose two distinct unit labels

```text
k_a,k_b in U minus {u_0}
```

and form the two target relations

```text
t_a=w_a e_(u_0)-w_(u_0)e_a,
t_b=w_b e_(u_0)-w_(u_0)e_b.                         (9)
```

Replace `r_(k_a)` by `r_(k_a)+t_a` and `r_(k_b)` by
`r_(k_b)+t_b`. Call the resulting six-row integer matrix `R_(j,u_0)`.
Every row is still in `Lambda(w)`.

On the six columns `P`, target grafting changes nothing. The matrix is
diagonal with diagonal entry `-w_(u_0)`, so

```text
det R_(j,u_0)[P]=(-w_(u_0))^6,
13 does not divide det R_(j,u_0)[P].                 (10)
```

Thus (7) is the desired owner-aligned pivot and its complement is exactly

```text
E minus P={a,b,u_0}.                                 (11)
```

All nine columns of `R_(j,u_0)` are nonzero modulo thirteen:

- every pivot column has its nonzero diagonal entry;
- columns `a` and `b` have entry `-w_(u_0)` in their graft rows;
- column `u_0` has entry `w_k+w_a` or `w_k+w_b` in the
  two graft rows and entry `w_k` in the other unit rows, all congruent to
  the corresponding nonzero `w_k` modulo thirteen.

The source row can have zero in column `u_0` modulo thirteen, but the
column already has five nonzero unit-row entries. Hence the packet is
thirteen-bright at every label, including both terminal targets.

Every row has infinity norm at most `S`. Indeed an ungrafted row has only
the two coefficients `w_k,-w_(u_0)`, while a grafted row has coefficients

```text
w_k+w_a, -w_(u_0), -w_(u_0)
```

or the analogous triple with `b`; the triangle inequality and (6) apply.
Consequently

```text
sum_(rows r of R_(j,u_0)) ||r||_infinity <=6S.       (12)
```

This works for every source `j` and every one of the six choices of
omitted unit `u_0`.

## 3. The owner-aligned thirteen-adic address bank

The six rows have rank six modulo thirteen by (10), and no column is dark.
Apply THM-2301's sharp essential-arrangement bound. At depth `n>=1`, the
packet contains at least

```text
9*12^5*13^(6n-6)                                    (13)
```

distinct unanchored all-thirteen-unit relation addresses. At any prescribed
coordinate, including the source `j` or either target, at least

```text
9*12^4*13^(5n-5)                                    (14)
```

addresses can be normalized to coefficient `1 modulo 13^n`.

Take centered representatives of the six row-combination coefficients.
By (12), every address in (13)--(14) has an exact integer representative
of height at most

```text
3S(13^n-1).                                         (15)
```

Unlike THM-2301's height-`1453(13^n-1)/2` bank, (15) depends on the scalar
row. In exchange, the selected source and both target labels occur in the
prescribed pivot/complement configuration.

## 3a. The exact two-target quotient

The owner pivot has a sharper dual meaning. Reduce modulo thirteen and put

```text
K=(w modulo 13)^perp subset F_13^9,
L=rowspace(R_(j,u_0) modulo 13).                    (15a)
```

Then `dim K=8` and `dim L=6`. Because `L[P]=F_13^P`, every `x in K` has a
unique `ell in L` agreeing with it on the six pivot coordinates. The
residual `x-ell` is supported on

```text
C={a,b,u_0}.
```

On `C`, the scalar word has only the nonzero coordinate `w_(u_0)`.
Orthogonality forces the residual's `u_0` coordinate to vanish. Therefore

```text
K=L direct_sum span(e_a,e_b),
K/L canonically identified with F_13^{ {a,b} }.     (15b)
```

Every mod-thirteen relation class has a unique owner-gauge representative
on the two target blockers. The `169` quotient classes split as

```text
zero:                         1,
pure a-axis:                 12,
pure b-axis:                 12,
mixed a,b support:          144.                    (15c)
```

Projectively, the `14` nonzero directions are two pure axes and twelve
mixed directions. This is the same **type decomposition** as THM-2305's
terminal words `{a}`, `{b}`, and `{a,b}`. It is not yet an identification
of the objects: THM-2305 stratifies points of a current-service set, while
(15b) stratifies residue classes of exact relation addresses. The exact
missing landing statement is

```text
terminal word sigma
  -> a nonzero contributing relation address x
     with support(pi_(a,b)(x))=sigma,                (15d)
```

where `pi_(a,b)` is projection along (15b). A fork belongs naturally to
the mixed locus, not to a cosmetically oriented edge.

## 4. What an arbitrary rank-six packet can force

The failure at bounded visible height is matroidal and exact. Let `F` be
any field and let `R` be a rank-six matrix over `F` with six unit columns

```text
V_1,...,V_6
```

and three blocker columns. Suppose

```text
alpha_1 V_1+...+alpha_6 V_6=0,
alpha_i!=0 for every i.                             (16)
```

This is the reduction of the scalar relation modulo thirteen. Put

```text
H=span(V_1,...,V_6).                                 (17)
```

Then `rank(H)<=5`.

If `rank(H)<=4`, no set consisting of one blocker and five units can be a
basis. If `rank(H)=5`, (16) spans the one-dimensional relation space among
the units. Since every `alpha_i` is nonzero, every five of the six unit
columns are independent. Therefore, for a blocker `t`,

```text
{V_t} union {five unit columns} is a basis
iff V_t notin H.                                    (18)
```

The answer is all-or-none: `t` occurs in all six such pivots or in none.
Moreover,

```text
C^*=B minus H                                       (19)
```

is a nonempty cocircuit of the represented matroid, and the exact number
of one-blocker plus five-unit bases is

```text
6|C^*| in {6,12,18}.                                (20)
```

With the convention

```text
Delta_(t,r)
 =det(V_t,V_1,...,omit V_r,...,V_6),
```

there is, for every `t in C^*`, a nonzero scalar `gamma_t` such that

```text
Delta_(t,r)=(-1)^r gamma_t alpha_r.                  (21)
```

For `t notin C^*`, all six minors vanish. Thus the Pluecker data reconstruct
the cocircuit perfectly, but it does not select the analytically chosen
owner. The intrinsic object is the unary incidence `t in C^*`, not a
pairwise orientation. Turning the tied blockers in `H` into a tournament
would add information not present in the packet.

There is also an exact link to THM-2307. When `rank(H)=5`, the annihilator
of `H` in the row-coefficient space is one-dimensional. Its corresponding
row in `R` is supported exactly on `C^*`. Hence

```text
e_t in rowspace(R) iff C^*={t}.                     (21a)
```

Thus THM-2307's literal-dark branch detects a singleton cocircuit, not the
selected source. Its no-literal branch still permits `C^*={a,b}`, giving
all twelve target pivots but none of the six source pivots. The quotient
blocker gains in `F^6/H` form a projective three-label vector with zero
entries off `C^*`; this unary support and its gain ratios are legitimate,
but they do not define a tournament orientation.

## 5. A sharp bright hostile packet

The selected-owner failure survives rank six, complete column brightness,
the scalar annihilator, and an explicit all-unit address. Work over
`F_13` with basis `e_1,...,e_6`. Set

```text
V_(u_i)=e_i,                         1<=i<=5,
V_(u_6)=-(e_1+...+e_5),

V_j=e_1+2e_2,
V_a=e_6,
V_b=e_3+2e_4.                                      (22)
```

Then

```text
sum_(i=1)^6 V_(u_i)=0,                              (23)
```

so the scalar annihilator has blocker weights zero and all six unit weights
one. Every column is nonzero and the whole packet has rank six. The
coefficient vector

```text
lambda=(1,1,1,1,1,1)
```

evaluates on `(j,a,b,u_1,...,u_6)` as

```text
(3,1,3,1,1,1,1,1,8),                               (24)
```

an all-unit word.

Nevertheless, the unit hyperplane is

```text
H=span(e_1,...,e_5).
```

Both `V_j` and `V_b` lie in `H`; only `V_a` lies outside it. Hence the
selected source `j` occurs in none of the six desired pivots, while `a`
occurs in all six.

This is not merely an abstract scalar-type fantasy. For example, take the
positive primitive integer word

```text
(w_j,w_a,w_b;w_(u_1),...,w_(u_6))
 =(13,13^3,2*13^5;1,14,27,40,53,66).                (25)
```

It has live strict valuation profile `(1,3,5)` and reduces to the scalar
word used in (23). Since `w_(u_1)=1`, every row of (22) can be lifted to an
exact integer relation by choosing integer representatives away from
`u_1` and then solving the `u_1` coefficient exactly. Thus the hostile
matroid occurs inside the exact relation lattice of a correctly typed
scalar word.

Equation (25) is not asserted to be a scalar cover, and its lifted rows are
not asserted to lie in the degree-526 survivor span. The witness proves the
logical no-go:

```text
rank six + no dark columns + all-unit Kakeya addresses
does not force the selected-owner pivot.                         (26)
```

Here `C^*={a}`, so `e_a` is the literal-dark row predicted by (21a). The
packet sees a pure **target** cocircuit while missing the source owner.

## 6. A simultaneous septimal/thirteen-adic packet

There is a stronger full-lattice construction when the scalar word has at
least two septimally visible coordinates. Assume now that `w` is primitive
and put

```text
supp_7(w)={e in E:7 does not divide w_e}.             (27)
```

Suppose

```text
|supp_7(w)|>=2.                                      (28)
```

Choose the omitted unit `u_0` so that the complement

```text
C={a,b,u_0}
```

meets `supp_7(w)`. This is always possible: if a target is in the support,
any omitted unit works; otherwise (28) forces some unit into the support,
and that unit can be omitted.

We use the following elementary field-section lemma.

### Field-section lemma

Let `q` be a prime power, let `v in F_q^9` have support of size at least two,
and let `P` be a six-set whose three-element complement `C` meets that
support. There is a `6 by 9` matrix `R_q` such that

```text
(R_q)[P]=I_6,
R_q v=0,
every column of R_q is nonzero.                     (29)
```

To prove this, write the three complement columns as `X_c`. Equation
`R_qv=0` is

```text
sum_(c in C) v_c X_c=-v_P.                          (30)
```

If exactly one `v_c` is nonzero, the corresponding column is forced to
`-v_P/v_c`, which is nonzero because the total support has size at least
two. Give the zero-weight complement labels arbitrary nonzero columns. If
at least two `v_c` are nonzero, choose all but the last corresponding
columns nonzero while avoiding the single choice that would make the last
one zero, and then solve for the last column. This proves (29).

Apply the lemma over `F_7` to `w modulo 7`. Apply it over `F_13` to
`w modulo 13`: its support is exactly `U`, and `u_0 in C`. In both cases
use the same owner pivot `P`. Combine the two matrices entrywise by the
Chinese remainder theorem to obtain a matrix `Y` modulo `91`, centered in
`[-45,45]`.

Choose a Bezout vector

```text
z in Z^9,                 z.w=1,
B(w)=||z||_infinity.                                  (31)
```

For each centered row `y` of `Y`, the integer `y.w` is divisible by `91`.
Set

```text
r=y-(y.w)z.                                           (32)
```

Then

```text
r in Lambda(w),
r= y modulo 91,
||r||_infinity<=45(1+S B(w)).                         (33)
```

The resulting exact six-row packet has the owner pivot `P` invertible
modulo both seven and thirteen, and every column is nonzero modulo both
primes.

Applying the sharp arrangement theorem separately at seven and thirteen
and then using CRT gives, at every depth `n>=1`, at least

```text
(3*6^5)(9*12^5)91^(6n-6)
 =52242776064*91^(6n-6)                              (34)
```

distinct exact relation addresses whose nine coordinates are all units
modulo `91`. At any prescribed coordinate, at least

```text
(3*6^4)(9*12^4)91^(5n-5)
 =725594112*91^(5n-5)                               (35)
```

can be normalized to `1 modulo 91^n`. Centered row-combination
coefficients and (33) give the exact height bound

```text
135(1+S B(w))(91^n-1).                               (36)
```

Every coordinate of an address in (34)--(35) is a seven-unit. Therefore
all nine interval Fourier factors of lengths `5/7` or `6/7` are nonzero.
These are genuine full-lattice Fourier needles with the source owner in
the prescribed pivot and both targets in its complement.

Condition (28) is sharp for an all-seven-unit bank. If
`supp_7(w)={e_0}`, every exact relation `m.w=0` satisfies

```text
m_(e_0)=0 modulo 7.                                  (37)
```

No all-seven-unit exact relation exists at any height. A Fourier term can
still survive by having `m_(e_0)=0` literally; a nonzero multiple of seven
at that coordinate kills the interval factor. Thus (37) is a congruence
obstruction in the full lattice, not the literal support deletion which
THM-2301 obtains only inside its bounded visible carrier.

## 7. Composition with the canonical handoff, and the exact loss

The connection ledger is:

```text
source:
  the exact relation lattice Lambda(w), its reductions modulo seven and
  thirteen, THM-2301's affine-arrangement count, and THM-2305's selected
  source owner j with target word among {a},{b},{a,b};

map:
  choose the owner-complement pivot, build a star section, graft both
  target columns, and optionally combine independent seven- and
  thirteen-adic sections by CRT before exact Bezout lifting;

preserved:
  all nine scalar labels, the selected source, both target blockers in
  the pivot complement, the exact target-plane quotient (15b), exact
  integer relation equations, all-depth injectivity, and in (34)--(35)
  nonzero interval Fourier support;

destroyed:
  the coefficient-independent height 526, membership in W_vis, the actual
  Jackson survivor that forced rank six, the source root character,
  endpoint component, signed amplitude, base phase, and terminal-word
  mass;

sharp hostile controls:
  (22)--(26) for selected-owner incidence in a bounded bright packet,
  (37) for a fully active septimal bank, and THM-2299/2301's phase
  cancellation even after an anchored rank-six minor;

needed sidecar:
  prove that THM-2305's bounded signed handoff coefficient lands on one
  owner-aligned exact address in the visible carrier with terminal word
  equal to its target-plane support as in (15d), or retain a
  terminal-component phase current while replacing the visible packet by
  the row-dependent one.                                            (38)
```

Thus the composition resolves the matroid existence question: a suitable
owner pivot is always present in the full exact lattice, and under (28) it
supports many exact nonzero Fourier needles. The remaining obstruction is
short **actual-Fourier visibility plus phase landing**, not basis
existence. No scalar profile is excluded and LRC(14) remains open.

## 8. Exact companion

The companion verifies the direct star and target graft for every omitted
unit on the typed row (25), the hostile packet's ranks, scalar relation,
all-unit word, and all eighteen owner/unit minors, the field-section
construction over every nontrivial septimal support pattern, the CRT
exact-lift identities, and all address constants and height ledgers. Run

```bash
python3 04-computation/lrc14_owner_aligned_pivot_thm2309.py
python3 -O 04-computation/lrc14_owner_aligned_pivot_thm2309.py
```

Both executions must match the stored output byte-for-byte after LF
normalization. QED.
