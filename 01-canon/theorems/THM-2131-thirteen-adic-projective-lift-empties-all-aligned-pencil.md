---
id: THM-2131
title: "Thirteen-adic projective lifting empties the all-aligned guard pencil"
status: >
  PROVED. In the rank-two LRC torus model, suppose the guard is 13 times a
  character, all eight terminals are nonzero modulo 13, and their reductions
  lie in one projective line. Four simultaneous singleton phases force the
  divided guard into the exact integer span of every four terminals. If the
  terminals share a projective line modulo 13^k, an explicitly periodic
  169-point next-digit transversal restricts their bands to eight double
  strips over F_13^2 while fixing the guard. THM-2124 gives an (8) or (7,1)
  lift. The split lift would make seven terminal bands cover the guard-safe
  torus, contrary to THM-2097 in rational rank two and THM-2128 in rational
  rank one. Alignment therefore propagates through every power of thirteen.
  Infinite determinant divisibility and the exact span relation collapse the
  guard and all terminals to rational rank one, a contradiction. This empties
  THM-2124's all-aligned guard-blocker/no-terminal-blocker branch; it does not
  close the complementary prime-invoice branch or prove LRC(14).
source: codex-2026-07-22-LRC-thirteen-adic-projective-lift
depends_on:
  - THM-2097
  - THM-2124
  - THM-2128
related:
  - THM-2123
  - THM-2125
  - THM-2130
---

# THM-2131 -- thirteen-adic projective lift

Let `Gamma` be a rank-two character lattice and put

```text
K=Hom(Gamma,R/Z),                    epsilon=1/14.     (1)
```

Take nonzero characters

```text
g=13u,                    c_1,...,c_8 in Gamma.        (2)
```

Assume:

1. `g,c_1,...,c_8` span `Gamma tensor Q`;
2. every `c_i mod 13` is nonzero, and all eight reductions lie in one
   projective line in `Gamma/13Gamma`;
3. the inherited LRC specialization is present: some primitive integral
   cocharacter `d` makes `g.d` positive odd and the `c_i.d` positive and
   pairwise distinct;
4. outside a null set `E`, the terminal danger bands cover the guard-safe
   region:

```text
C={X in K:||g.X||>1/7}
 subset union_(i=1)^8 {X:||c_i.X||<epsilon}.          (3)
```

We prove that these assumptions are inconsistent. Together with THM-2128,
this empties both direction partitions in THM-2124's guard-`13`-blocker branch
when no terminal is itself a `13`-blocker.

## 1. Four singleton phases capture the divided guard integrally

Fix four terminal labels `I subset {1,...,8}`. We first prove the exact phase
containment

```text
{Y:||c_i.Y||<epsilon for every i in I}
 subset {Y:||u.Y||<=1/7}.                             (4)
```

Suppose instead that the strict reverse inequalities hold on a nonempty open
set `O`. Multiplication by thirteen is a finite covering of `K`. The image
`[13](E)` is null: on finitely many local covering charts, `[13]` is a smooth
local diffeomorphism and hence takes null sets to null sets. Choose

```text
Y in O minus [13](E).                                 (5)
```

Every one of the 169 roots `X` of `13X=Y` is outside `E`, and

```text
g.X=13u.X=u.Y.                                       (6)
```

Thus the guard is safe on the entire affine root plane

```text
{X:13X=Y}=X_0+K[13] isomorphic to F_13^2.             (7)
```

Every terminal reduction has the same nonzero projective normal on this
plane. The exact thirteen-grid count used in THM-2116 and THM-2124 says that
terminal `i` selects one affine line when

```text
||c_i.Y||<=epsilon,                                   (8)
```

and otherwise selects at most two parallel affine lines. The four labels in
`I` therefore select four lines in total, while the other four select at most
eight. Twelve parallel lines cannot cover the thirteen-line affine pencil.
This contradicts (3) at one of the roots in (7), and proves (4).

Put

```text
K_I=intersection_(i in I) ker(c_i).                   (9)
```

Equation (4) confines the compact subgroup `u(K_I)` of the circle to the
closed radius-`1/7` arc. That subgroup must be trivial. Indeed, the full
circle contains `1/2`, while a nontrivial finite circle subgroup of order `m`
contains an element of norm

```text
floor(m/2)/m >= 1/3.                                  (10)
```

Hence `u` kills `K_I`. The exact annihilator identity, equivalently Smith
normal form for the character map

```text
(c_i)_(i in I):K -> (R/Z)^I,                          (11)
```

is

```text
K_I^perp=image(Z^I -> Gamma)=sum_(i in I) Z c_i.      (12)
```

It is the integer span, not merely its saturation. Consequently

```text
u in sum_(i in I) Z c_i                  for every |I|=4. (13)
```

For the induction below, fix once and for all the instance

```text
u in Zc_1+Zc_2+Zc_3+Zc_4.                            (14)
```

## 2. Seven terminal bands are globally irredundant

We will repeatedly use the following consequence of the inherited LRC
specialization.

> **Seven-band lemma.** No seven of the terminal bands can cover `C` almost
> everywhere.

Fix a seven-label set `J`. If

```text
rank_Q{g,c_j:j in J}=2,                               (15)
```

then `g` and those seven terminals satisfy THM-2097 with the same primitive
specialization `d`. That theorem gives a nonempty open set in `C` outside all
seven danger bands, contradicting an almost-everywhere seven-band cover.

If the rank in (15) is one, let `alpha` be the primitive generator of the
common rational character line. Since `g=13u` and the saturated lattice line
is `Z alpha`, there are integers `k,r_j` with

```text
u=k alpha,                    c_j=r_j alpha.           (16)
```

Every `r_j` is a unit modulo thirteen because `c_j mod 13` is nonzero. The
primitive character `alpha:K->T` is surjective and pushes Haar measure to Haar
measure. Since every predicate depends only on `alpha.X`, an almost-everywhere
seven-band cover would therefore descend to

```text
{t:||13k t||>1/7}
 subset union_(j in J) {t:||r_j t||<1/14}.            (17)
```

This is excluded by THM-2128's one-dimensional seven-comb theorem. The two
rank cases prove the seven-band lemma.

## 3. The next-digit transversal

Let

```text
N=13^k,                         k>=1,                 (18)
```

and suppose the eight terminal reductions lie in one primitive projective
line modulo `N`; here this means a primitive cyclic direct summand of
`Gamma/NGamma`. Because every `c_i` remains nonzero modulo thirteen, each is
a unit multiple of a primitive generator of this summand. Lifting that
generator and making an integral basis change gives a basis `a,b` of `Gamma`
and integers

```text
c_i=R_i a+N h_i b,             13 not|R_i.            (19)
```

Equation (14) puts `u mod N` in the same cyclic module. Absorbing its first
coordinate into an integer `A` gives

```text
u=Aa+N h_0b.                                           (20)
```

For integer representatives `p,q in {0,...,12}`, define a point `z_pq in K`
by

```text
a.z_pq=p/13,                  b.z_pq=q/(13N).          (21)
```

The 169 points in (21) are distinct. They are an indexed next-digit
**transversal**, not a subgroup: wrapping `q` by thirteen changes `b.z` by
`1/N`. This loss is harmless because every relevant character evaluation is
periodic in the indices. Equations (19)--(21) give

```text
c_i.(X+z_pq)=c_i.X+(R_i p+h_i q)/13       in R/Z,
g.(X+z_pq)=g.X+A p+h_0 q=g.X              in R/Z.      (22)
```

Replacing `p` by `p+13` changes the first terminal expression by the integer
`R_i`; replacing `q` by `q+13` changes it by the integer `h_i`. The guard
changes integrally as well. Thus (22) defines honest functions on the index
plane `F_13^2`, even though (21) is not a homomorphism from that plane to `K`.

Clean the exceptional set by all 169 translations:

```text
E_N^#=union_(p,q in {0,...,12})(E-z_pq).              (23)
```

This finite union is null. Since the nonzero character `g` is surjective,
`C` has Haar measure `5/7`, so `C minus E_N^#` is nonempty and has positive
measure.

For `X in C minus E_N^#`, every point `X+z_pq` remains guard-safe by (22),
lies outside `E`, and is therefore covered by a terminal band. On the index
plane the mask of terminal `i` is

```text
A_i(X)={(p,q) in F_13^2:
 ||c_i.X+(R_i p+h_i q)/13||<epsilon}.                 (24)
```

The linear form

```text
ell_i(p,q)=R_i p+h_i q                               (25)
```

is nonzero because `R_i` is a unit. Its values are a translated
thirteen-grid on the circle. An interval of length `1/7<2/13` contains at
most two grid values, so (24) is a union of at most two parallel affine lines
with projective normal

```text
[R_i:h_i] in P^1(F_13).                               (26)
```

The eight masks cover the indexed plane. THM-2124 applies exactly and says
that the multiplicity partition of the normals (26) is

```text
(8) or (7,1).                                         (27)
```

The affine offsets `c_i.X` in (24) have not been discarded; THM-2124 is
uniform in those offsets.

## 4. The split lift would make one terminal redundant

Suppose the second pattern in (27) holds, and call its seven aligned labels
`J`. Their masks are unions of affine lines in one common pencil. If their
union missed one line of that pencil, the exceptional transverse double
strip would meet the missed thirteen-point line in at most two points. It
could not finish the cover. Hence the seven aligned masks already cover all
of `F_13^2`.

This holds for every `X in C minus E_N^#`. Taking the index `(p,q)=(0,0)` in
(24) shows

```text
C subset union_(j in J) {X:||c_j.X||<epsilon}
                                             almost everywhere. (28)
```

Equation (28) contradicts the seven-band lemma. Therefore the split pattern
is impossible, and all eight next-digit normals in (26) are projectively
equal.

Since every `R_i` is a unit modulo thirteen, there is a `delta in F_13` such
that

```text
h_i=delta R_i mod 13                  for all i.       (29)
```

Choose an integral lift `d_0` of `delta` and put

```text
a'=a+N d_0 b.                                          (30)
```

The pair `(a',b)` is again an integral basis (the change from `(a,b)` is
unipotent), and (19), (29) give

```text
c_i-R_i a'=N(h_i-d_0R_i)b in 13N Gamma.               (31)
```

Thus all eight terminals lie in the one primitive projective line generated
by `a'` modulo `13N`. In determinant form, up to the basis unit
`det(a,b)=+/-1`, the same lift is

```text
det(c_i,c_j)=+/-N(R_i h_j-R_j h_i)=0 mod 13N.         (32)
```

We have proved the induction step

```text
common projective line mod 13^k
  implies common projective line mod 13^(k+1).        (33)
```

## 5. Infinite lifting forces rational rank one

Hypothesis 2 supplies the base `k=1` of (33). Induction gives, for every
pair `i,j` and every `k>=1`,

```text
13^k divides det(c_i,c_j).                            (34)
```

An integer divisible by every power of thirteen is zero. Hence

```text
det(c_i,c_j)=0                         for all i,j,    (35)
```

so the eight nonzero terminals lie on one rational character line. The exact
integer relation (14) puts `u`, and therefore `g=13u`, on that same line.
This contradicts the rational rank-two assumption in item 1 and proves the
theorem. QED.

## 6. Scope, assumption challenge, and Tournament Analysis

The proof closes only THM-2124's `(8)` direction pattern with a guard
`13`-blocker and eight terminal nonblockers. THM-2128 already closes `(7,1)`.
This theorem does not by itself exclude simultaneous guard and terminal
blockers, the complementary nonblocker-guard prime invoice, ranks nine
through twelve, the remaining finite rank-seven boxes, or LRC(14).

The challenged assumption was that a next-digit Kakeya needle had to be a
finite subgroup. The points (21) are instead a section of one base-13 digit
box. Passing to their indices preserves the guard and every terminal danger
predicate exactly by the periodic laws (22), while it forgets the geometric
wrap by `1/N`. The integer coefficients `N h_i` are precisely the sidecar
that makes that forgotten wrap invisible to all relevant characters. This is
why the affine offset drift survives the projective lift.

Natural candidate vertices were runners, terminal characters, projective
directions at each `13`-adic depth, digit-grid points, affine strips, and proof
obligations. The intrinsic pairwise observable on terminals is

```text
v_13(det(c_i,c_j)),                                   (36)
```

a symmetric ultrametric similarity, not an orientation. Artificially
orienting a pair by which lift coefficient is larger produces ties and loses
the common cyclic module needed in (20). Score histograms, directed cycles,
SCCs, edge flips, and Hamiltonian paths therefore do not preserve the cover
predicate. The faithful carrier is the `13`-adic dendrogram of projective
directions, together with the affine masks (24), the exact span sidecar (14),
and the null-set cleaning (23). The operative self-similarity is the lift

```text
projective pencil mod N -> double-strip pencil on one digit
                         -> projective pencil mod 13N, (37)
```

not a tournament on the eight runner labels.
