---
id: THM-3341
title: "U-spine square-hypotenuse transplant, Markov adjacency, and norm-17 triangular-plane torsors"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the U-spine
  hypotenuse C_t=2t^2+2t+1 and plane scalar Q_t=2C_t+1, the complete positive
  intersection table is: C_t is square on one negative-Pell orbit and never
  triangular; Q_t is never square and is triangular on exactly two norm-17
  Pell orbits; the C- and Q-value sets are disjoint.  Consecutive C-square
  roots are precisely the fixed-two Markov coordinates compiled by each
  square-triangular row.  Gaussian squaring transplants the positive Berggren
  middle ray into the sparse C-square depths of the U-ray, with state-dependent
  drift rather than a branch homomorphism.  Equal-hypotenuse Boolean fibre
  ranks are unbounded even inside this square selector.  The first nontrivial
  square-selector fibre is C_696=985^2, with two explicit ancestry words.
  Triangular Q_t has tournament-size residue 3 mod 16 but does not select
  Boolean grade or an orientation.  No LRC, Jacobian, or tournament existence
  claim follows.
source: codex-2026-08-12-u-spine-square-triangular-intersections
audit: >
  independent read-only Pell/Markov/branch audit (complete intersection
  table, negative-Pell and norm-17 descents, Gaussian-square transplant,
  unbounded Boolean-fibre proof, first collision, and scope boundaries:
  ACCEPT); fresh hostile audit found and repaired the Markov orientation/
  boundary and shared-field wording (MISTAKE-369), then accepted every
  infinite classification, ancestry word, fibre claim, replay, and hash;
  dependency-free normal/-O companion replay: ACCEPT
depends_on:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
related:
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
script: 04-computation/u_spine_square_triangular_torsors_thm3341.py
output: 05-knowledge/results/u_spine_square_triangular_torsors_thm3341.out
script_sha256: 4d67cc02b8929b257c7aaca4bec2edefa4d23533191be92042af963fb1d3ffd7
output_sha256: 9739622dad5e57108d966d0aa41cb1949974bfe78e4169bba0f4beba562fcbc2
hash_basis: working-tree bytes (LF)
---

# THM-3341 -- square and triangular intersections on the U-spine

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is a repository synthesis and proof interface; no literature-priority
claim is made.  It continues THM-3334's parabolic branch and uses THM-3335's
proved square-triangular/Markov compiler.  THM-3339's golden ancestry ray is a
different discriminant-`5` carrier and is used only as a scope comparison.

## 1. The two scalar sequences and the complete intersection table

Extend the physical branch `t>=1` to its algebraic boundary `t=0` and put

```text
C_t=t^2+(t+1)^2=2t^2+2t+1,
Q_t=2C_t+1=(2t+1)^2+2.                                  (1)
```

THM-3334 identifies

```text
S_t=(2t+1,2t(t+1),C_t)                                  (2)
```

as the consecutive-parameter primitive Pythagorean U-spine.  Its Berggren
depth is `t-1` for `t>=1`.  The complete positive intersection table is:

| condition, `t>=1` | complete answer |
|---|---|
| `C_t` square | infinitely many: `t=3,20,119,696,4059,...` |
| `C_t` triangular | none |
| `C_t` square-triangular | none |
| `Q_t` square | none |
| `Q_t` triangular | infinitely many: `t=6,23,221,798,7524,...` |
| `Q_t` square-triangular | none |
| `C_t=Q_s` for nonnegative `s,t` | none |

At the degenerate boundary, `C_0=1` is square-triangular and
`Q_0=3=T_2`.  These are not positive-parent descendant planes.

The five negative rows have elementary proofs.  Write `X=2t+1`.  If
`C_t=T_n`, then with `Y=2n+1`,

```text
Y^2-4X^2=5,
(Y-2X)(Y+2X)=5.                                         (3)
```

Positive factor pairs force `(Y,X)=(3,1)`, hence `t=0`.  Also `X` is odd,
so

```text
Q_t=X^2+2=3 mod 8,                                      (4)
```

which is never a square.  Finally

```text
C_t=1 mod 4,                 Q_s=3 mod 4,                (5)
```

so the two value sets are disjoint.  The square-triangular exclusions follow
from the corresponding square or triangular exclusions.

## 2. Every square `C_t` is one negative-Pell row

The square equation is

```text
C_t=m^2
iff X^2-2m^2=-1,                 X=2t+1.                 (6)
```

All positive solutions of (6) are

```text
X_k+m_k sqrt(2)=(1+sqrt(2))^(2k+1),             k>=0.    (7)
```

For completeness, a direct descent proves exhaustion.  When `m>1`, the
smallest possibility is `m>=5`, and multiplication by `3-2sqrt(2)` gives

```text
(X,m) -> (3X-4m,3m-2X).                                 (8)
```

Both coordinates are positive, the second is strictly smaller than `m`, and
(6) is preserved.  Repetition reaches `(1,1)`; reversing (8) generates (7)
uniquely.

Let the Pell numbers be

```text
pell_0=0,       pell_1=1,       pell_(r+2)=2pell_(r+1)+pell_r. (9)
```

Expanding (7) yields the exact selector

```text
m_k=pell_(2k+1),
t_k=(pell_(2k)+pell_(2k+1)-1)/2,                        (10)
```

with recurrences

```text
m_(k+2)=6m_(k+1)-m_k,
t_(k+2)=6t_(k+1)-t_k+2,                                 (11)
```

and seeds `(t_0,m_0)=(0,1)`, `(t_1,m_1)=(3,5)`.  Thus

```text
t_k: 0,3,20,119,696,4059,...,
m_k: 1,5,29,169,985,5741,... .                           (12)
```

The square roots are exactly the odd-index Pell numbers.  This negative-Pell
selector is not THM-3335's square **even-leg** selector: the latter selects
`t=1,8,49,288,...` where `2t(t+1)` is square.  Section 3 gives their exact
adjacency relation.

## 3. Square-triangular rows connect consecutive square hypotenuses

Take one positive THM-3335 row

```text
R^2=T_N=N(N+1)/2,                     A=2N+1.            (13)
```

Define

```text
M_-=A-2R,                 M_+=A+2R,
t_-=2R-N-1,               t_+=2R+N.                     (14)
```

Direct substitution gives

```text
C_(t_-)=M_-^2,             C_(t_+)=M_+^2,               (15)
M_-^2+M_+^2+4=6M_-M_+,                                   (16)
M_-M_+-(2R)^2=1.                                        (17)
```

Thus `(2,M_-,M_+)` is a Markov triple and the doubled square-triangular root
is the Pell--Cassini carry between its two nontrivial coordinates.  Conversely,
modulo swapping the two non-two coordinates, orient a positive fixed-two
Markov solution uniquely as

```text
1<=M_-<=M_+.
```

The equality row `(M_-,M_+)=(1,1)` is the algebraic boundary `R=N=0`.
THM-3335 proves that the strict rows `M_-<M_+` correspond exactly and uniquely
to positive square-triangular rows.  Their inverse is

```text
R=(M_+-M_-)/4,              N=(M_-+M_+-2)/4.             (18)
```

Equations (10) and (14) identify the two square hits as consecutive rows:

```text
(t_-,M_-)=(t_(k-1),m_(k-1)),
(t_+,M_+)=(t_k,m_k).                                     (19)
```

The first rungs are:

| `R^2=T_N` | `(M_-,M_+)` | `(t_-,t_+)` |
|---|---|---|
| `1=T_1` | `(1,5)` | `(0,3)` |
| `36=T_8` | `(5,29)` | `(3,20)` |
| `1225=T_49` | `(29,169)` | `(20,119)` |
| `41616=T_288` | `(169,985)` | `(119,696)` |

This makes the previously visible identity

```text
29*169-70^2=1                                             (20)
```

the `N=49,R=35` instance of (17).  The cannonball coincidence is unique
inside this alignment.  If `N=2h+1` and

```text
1^2+...+h^2=(2R)^2,
```

then using (13) reduces the equation to `h/6=4`, hence

```text
(h,2R,N)=(24,70,49).                                     (21)
```

This uniqueness is scoped to the displayed square-triangular alignment; it
is not an independent classification of every square-pyramidal equation.

## 4. Gaussian squaring is a state-dependent branch transplant

Let

```text
A_+=[1 2 2;2 1 2;2 2 3]                                 (22)
```

be the positive/middle Berggren child and let `U` be THM-3334's parabolic
child.  Define Gaussian squaring on an ordered primitive triple by

```text
Gamma(a,b,c)=(|a^2-b^2|,2ab,c^2).                        (23)
```

The `A_+` ray is exactly the consecutive-leg sequence

```text
(3,4,5), (21,20,29), (119,120,169),
(697,696,985), ...                                      (24)
```

with the odd leg kept first.  Its hypotenuses are `m_1,m_2,...` from (12),
and the unordered legs are `(t_k,t_k+1)`.  Therefore (23) gives, for every
`j>=0`,

```text
Gamma(A_+^j(3,4,5))
  =S_(t_(j+1))
  =U^(t_(j+1)-1)(3,4,5).                                (25)
```

For example,

```text
(3,4,5)       -> (7,24,25)=U^2(3,4,5),
(21,20,29)    -> (41,840,841)=U^19(3,4,5).               (26)
```

The drift is not constant.  From (14) and (19),

```text
t_k-t_(k-1)=2N_k+1,                                     (27)
```

so successive `A_+` edges map to `17,99,577,...` U-edges after the first
image.  Equation (25) is a state-dependent branch transplant, not a
conjugacy, constant-length substitution, or semigroup homomorphism.

## 5. Boolean fibres remain unbounded inside the square selector

The square-hypotenuse subfamily already contains unbounded Gaussian
factor-choice ranks.  Pell numbers obey the strong-divisibility law

```text
gcd(pell_a,pell_b)=pell_gcd(a,b),
pell_a divides pell_b whenever a divides b.             (28)
```

This follows from the Pell addition identities and the ordinary Euclidean
algorithm on indices; the companion checks all pairs below `50` independently.

Choose distinct odd primes `ell_1,...,ell_r`, set `L=product ell_i`, and put

```text
m=pell_L.                                                 (29)
```

Because `L` is odd, (10) places `m^2` on the `C_t` square selector.  The
pairwise-coprime integers `pell_(ell_i)>1` all divide `m`, so

```text
omega(m)>=r.                                              (30)
```

Every prime divisor of `m` is `1 mod 4`: from
`m^2=t^2+(t+1)^2` and `gcd(t,t+1)=1`, reduction modulo an odd divisor of `m`
makes `-1` a square.  THM-3334's proved fixed-hypotenuse torsor therefore has

```text
|X_(m^2)|=2^(omega(m)-1),                                 (31)
```

which is unbounded as `r` grows.  The selector does not tame ancestry
multiplicity.

The first nontrivial square-selector fibre is

```text
t=696,             C_t=985^2=970225,             985=5*197. (32)
```

Its two primitive parents are

```text
(1393,970224,970225),
(522767,817344,970225).                                  (33)
```

Their root-to-child ancestry words are respectively

```text
U^695,                    UUUUUDADUDDU.                   (34)
```

The scalar square value in (32) loses a depth-`695` parabolic ancestry and a
depth-`12` wild ancestry at once.

## 6. Triangular `Q_t` is exactly two norm-17 orbits

Write `Q_t=T_nu`, `X=2t+1`, and `Y=2nu+1`.  Then

```text
Y^2-8X^2=17.                                             (35)
```

All positive solutions of (35) descend under multiplication by
`3-sqrt(8)`:

```text
(Y,X)->(3Y-8X,3X-Y).                                     (36)
```

For `X>=5`, both coordinates are positive and the second decreases.  The
small cases leave exactly the reduced seeds `(5,1)` and `(7,2)`.  Under the
forward unit `3+sqrt(8)`, parity of `X` alternates.  Keeping the required odd
`X` therefore leaves exactly two step-two orbits:

```text
Y+X sqrt(8)=(5+sqrt(8))(17+6sqrt(8))^j,                  (37)
```

or

```text
Y+X sqrt(8)=(37+13sqrt(8))(17+6sqrt(8))^j,               (38)
```

for `j>=0`.  In `(t,nu)` coordinates, both obey

```text
(t,nu)->(17t+6nu+11,48t+17nu+32),                        (39)
```

with seeds `(0,2)` and `(6,18)`.  Thus every nonnegative hit is generated,
and the first positive hits are

```text
(t,nu)=(6,18),(23,66),(221,626),(798,2258),(7524,21282),... . (40)
```

Equation (39) preserves

```text
nu=2 mod 16.                                             (41)
```

If `Q_t=T_nu` is read only as a tournament edge count, its tournament size is
therefore `nu+1=3 mod 16`.  This is a cardinality label, not an orientation,
owner, or skew-EW construction.

Nor does triangularity select Boolean grade.  The first three positive rows
give

```text
t=6:      C_t=85=5*17,             |X_C|=2,
t=23:     C_t=1105=5*13*17,        |X_C|=4,
t=221:    C_t=98125=5^4*157,       |X_C|=2.              (42)
```

The first size-two and size-four fibres happen to be triangular-Q rows, but
the grade immediately falls again.

## 7. Three quadratic carriers and the stopping boundary

The nearby discriminants have different jobs:

| carrier | equation | arithmetic role |
|---|---|---|
| full U-spine | `C_t=2t^2+2t+1`, polynomial discriminant `-4` | Gaussian split primes and CRT Boolean fibres |
| square `C_t` | `X^2-2m^2=-1` | negative-Pell sparse depth selector |
| triangular `Q_t` | `Y^2-8X^2=17` | two norm-17 depth selectors |
| THM-3335 clock | `x^2-8R^2=1` | square even leg / fixed-two Markov compiler |
| THM-3339 golden ray | `|n^2-mn-m^2|=1` | three exact ancestry rays |

The visible intersections are typed maps, not equality of quadratic fields or
actions.  The three Pell lanes in the middle of the table do share one field
and fundamental positive unit:

```text
Q(sqrt(8))=Q(sqrt(2)),
3+sqrt(8)=3+2sqrt(2),
17+6sqrt(8)=(3+2sqrt(2))^2.                              (43)
```

The negative-Pell square selector and THM-3335 clock advance under the same
unit on different norm/parity cosets; the odd-`X` norm-17 branches advance
under its square.  Their typed states and selected predicates differ even
though their quadratic field does not.  In particular:

- THM-3335 joins adjacent negative-Pell roots by (14)--(19);
- Gaussian squaring gives the variable-length branch map (25);
- the `-4` in THM-1310's cubic discriminant remains only a square factor and
  does not identify these Pell or Gaussian fibres with Keller monodromy;
- triangular counts in (40) and the four/six cardinalities of Boolean fibres
  manufacture no tournament orientation;
- no fixed-cusp LRC owner, phase, endpoint word, or global first-exit state is
  supplied.

No LRC row is removed, no Jacobian or Dixmier conjecture is settled, and no
skew-EW or tournament existence theorem follows.

## 8. Exact reproduction

Run

```bash
python 04-computation/u_spine_square_triangular_torsors_thm3341.py
python -O 04-computation/u_spine_square_triangular_torsors_thm3341.py
```

The dependency-free companion checks the negative-Pell formulas and
recurrences; every `t<=10^6` against independent square/triangular tests;
fifteen square-triangular/Markov adjacency rows; the cannonball alignment;
nine exact Gaussian branch images and U-depths; Pell strong divisibility;
the first square-selector Boolean collision and both ancestry words; both
norm-17 branches; and the first Boolean-grade controls.  Normal and optimized
runs byte-match the stored transcript and end in `ALL CHECKS PASSED`.

QED.
