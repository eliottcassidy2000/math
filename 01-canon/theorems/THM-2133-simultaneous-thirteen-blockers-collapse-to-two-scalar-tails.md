---
id: THM-2133
title: "Simultaneous thirteen-blockers collapse to two scalar flood tails"
status: >
  PROVED. In the rank-two LRC guard/eight-terminal torus cover, suppose the
  guard and at least one terminal are divisible by thirteen. After common-
  content descent, one through six terminal blockers are impossible. The
  unique-blocker branch first forces the seven residual terminals into a
  minimal mod-13 pencil; a compact-subgroup cross lemma and a digit-section
  lift then propagate the guard, blocker, and pencil through every 13-adic
  level, contradicting rational rank two. With seven terminal blockers, the
  divided guard and blockers must instead lie on one rational line and the
  cover is exactly a scalar seven-comb containment. Its normalized valuation
  profile has guard a 13-unit and either six unit teeth plus one deeper tooth,
  or five unit teeth plus two deeper teeth. Necessary fourwise and pairwise gcd
  invoices are forced in those two tails. This leaves two explicit one-dimensional
  flood tails; it does not prove LRC(14).
source: codex-2026-07-22-LRC-simultaneous-thirteen-blockers
depends_on:
  - THM-2080
  - THM-2097
  - THM-2124
  - THM-2128
related:
  - THM-2122
  - THM-2125
  - THM-2130
  - THM-2131
external: Settled Lonely Runner Conjecture for at most eight integer speeds (nine total runners).
---

# THM-2133 -- simultaneous thirteen-blockers

Let `Gamma` be a rank-two character lattice and put

```text
K=Hom(Gamma,R/Z),                    epsilon=1/14.     (1)
```

Let a nonzero guard `g` and eight nonzero terminals `c_1,...,c_8` span
`Gamma tensor Q`. Assume the inherited LRC specialization: some primitive
integral cocharacter `d` makes

```text
g.d positive and odd,
c_i.d positive and pairwise distinct.                 (2)
```

Suppose, outside a null set `E`, that

```text
{X:||g.X||>1/7}
 subset union_(i=1)^8 {X:||c_i.X||<epsilon}.          (3)
```

Assume the guard is a thirteen-blocker and write

```text
g=13u,
B={i:c_i in 13Gamma},                  b=|B|,
c_i=13v_i                              for i in B.     (4)
```

If all nine characters in (3) are divisible by thirteen, substitution
`Y=13X` gives the same cover for all nine divided characters. Perform the
maximal common-content descent first. If the divided guard becomes nonzero
modulo thirteen, the instance exits this theorem for the THM-2123/2125
nonblocker-guard lane. We treat the complementary normalized branch in which
the guard remains a blocker and not all terminals are blockers. Thus

```text
1<=b<=7.                                               (5)
```

The omitted value `b=8` is precisely another common-content descent step, not
an additional normalized branch.

We prove:

1. `1<=b<=6` is impossible;
2. if `b=7`, the divided characters `u,v_i` lie on one rational character
   line, and (3) is equivalent to a scalar seven-comb cover;
3. after common scalar descent, that cover has guard coefficient prime to
   thirteen and exactly five or six terminal coefficients prime to thirteen;
4. the two surviving valuation profiles obey necessary gcd invoices recorded in
   Section 7.

Together with THM-2131's exclusion of `b=0`, this classifies the normalized
guard-`13`-blocker branch down to two scalar tails.

## 1. Root capacity turns multiple blockers into a smaller cover

Assume `b>=2`. If a point `Y` satisfies

```text
||u.Y||>1/7,
||v_i.Y||>epsilon                         for i in B, (6)
```

then all blockers are safe on every root `X` of `13X=Y`, while the guard is
safe there because

```text
g.X=u.Y,                    c_i.X=v_i.Y.              (7)
```

Choose `Y` off `[13](E)`; this is possible in any nonempty open set of phases
satisfying (6), because a finite smooth covering maps null sets to null sets.
The 169 roots form an affine plane over `F_13`. Every nonblocker terminal cuts
that plane in at most two affine lines, hence in at most 26 points. Since

```text
26(8-b)<=156<169,                                    (8)
```

the nonblockers cannot cover the root plane. Therefore (3) forces the exact
phase containment

```text
{Y:||u.Y||>1/7}
 subset union_(i in B){Y:||v_i.Y||<epsilon}           (9)
```

up to a null set. A strict counterexample would have an open neighborhood,
so the distinction between exact and almost-everywhere containment is
harmless below.

Put

```text
H=u.d=g.d/13,
q_i=v_i.d=c_i.d/13.                                  (10)
```

The integer `H` is positive odd, and the `q_i` are positive and pairwise
distinct. If `2<=b<=6`, THM-2080 supplies a positive-measure set of `t` with

```text
||Ht||>1/7,                    ||q_i t||>epsilon.     (11)
```

The geodesic `Y=td` contradicts (9). Hence

```text
2<=b<=6 is impossible.                               (12)
```

## 2. Seven blockers are exactly a scalar seven-comb

Let `b=7`, and let `c_0` be the sole nonblocker. If

```text
rank_Q{u,v_i:i in B}=2,                               (13)
```

THM-2097 applies with guard `u`, the seven terminals `v_i`, and the
specialization (10). It supplies a nonempty open strict escape from (9), a
contradiction. Thus all eight divided characters lie on one saturated
rational line. For a primitive generator `alpha`, write

```text
u=H alpha,                    v_i=r_i alpha.          (14)
```

After reversing `alpha` if needed, the coefficients are positive and
distinct on the terminal side. The ambient rank-two assumption makes `c_0`
rationally transverse to `alpha`. Haar pushforward under the primitive
character `alpha:K->T` turns (9) into

```text
C_H subset union_(i=1)^7 D_(r_i)        almost everywhere, (15)
C_H={t:||Ht||>1/7},
D_r={t:||rt||<1/14}.
```

This scalar containment is also sufficient for the original cover: if
`Y=13X` and the original guard is safe, (15) makes some `v_i.Y=c_i.X`
dangerous. Hence the transverse nonblocker `c_0` is inert, and (3) is
equivalent to (15).

It remains to eliminate the unique-blocker branch and then sharpen (15).

## 3. A compact-subgroup cross lemma

We use the following exact consequence of two small phases.

> **Cross lemma.** Let `a,b,p,q` be characters of a compact torus and suppose
> ```text
> {Y:||a.Y||<epsilon, ||b.Y||<epsilon}
>  subset {Y:||p.Y||<=1/7} union {Y:||q.Y||<=epsilon}. (16)
> ```
> Then
> ```text
> p in Za+Zb                    or q in Za+Zb.         (17)
> ```

Put `L=ker(a) intersect ker(b)`. If neither `p` nor `q` kills `L`, push Haar
measure on `L` through their nontrivial image subgroups in the circle. In any
nontrivial compact circle subgroup, a closed radius-`1/7` arc and a closed
radius-`1/14` arc each have relative measure at most `1/2`. Equality can occur
only for the order-two image, where the preimage of the arc is an index-two
kernel. Thus (16) could cover `L` only if both images had order two and two
proper index-two subgroups covered `L`. Equal kernels cover only half of `L`;
distinct kernels cover three quarters. This is impossible.

Therefore one of `p,q` kills `L`. Exact annihilator duality gives

```text
L^perp=Za+Zb,                                         (18)
```

which proves (17). This is an integer span, not its saturation.

There is a coefficient refinement when `a,b` are rationally independent. If
only `p` kills `L`, write `p=ma+nb`. If `|m|+|n|>=3`, prescribe small phases
with the signs of `m,n` so that `||p.Y||>1/7`; on their common fiber the
nontrivial image of `q` has a point outside its radius-`1/14` arc. This
contradicts (16). Hence

```text
p=ma+nb,                       |m|+|n|<=2.            (19)
```

If only `q` kills `L`, the same argument uses the smaller threshold for `q`
and a guard-safe point in the nontrivial `p`-image. It gives

```text
q=+/-a                         or q=+/-b.             (20)
```

Every finite coset used here has a point of norm at least `1/4`, so the strict
thresholds `1/7` and `1/14` have margin.

## 4. The unique blocker forces a common residual pencil

Assume now

```text
g=13u,                  c_*=13v,
c_1,...,c_7 notin 13Gamma.                            (21)
```

Choose a clean phase `Y` on which `u` is guard-safe and `v` is terminal-safe.
Such phases have positive measure: the two excluded bands have total measure
at most `3/7`. Remove the null set `[13](E)` before choosing `Y`, so all 169
roots are clean. On those roots the guard and blocker are safe, so the seven
residual terminals must cover the affine root plane. They are seven double
strips. THM-2124's sharp minimal-layer theorem says their projective
normals are all equal:

```text
c_1,...,c_7 mod 13 lie in one primitive line L_1.    (22)
```

If two residual labels `i!=j` have singleton phases at `Y`, their two masks
contain 13 points each and the other five contain at most 26 each. Their total
capacity is

```text
2*13+5*26=156<169.                                   (23)
```

Consequently every pair obeys the exact containment

```text
{Y:||c_i.Y||<epsilon, ||c_j.Y||<epsilon}
 subset {Y:||u.Y||<=1/7} union {Y:||v.Y||<=epsilon}. (24)
```

The cross lemma yields the permanent exact-span alternative

```text
u in Zc_i+Zc_j                    or
v in Zc_i+Zc_j                    for every i!=j.     (25)
```

## 5. One-sided projective alignment is impossible

Suppose, at some modulus `N=13^k`, all residuals lie in one primitive cyclic
summand `L_N` of `Gamma/NGamma`. Reducing (25) modulo `N` shows that at least
one of `u,v` belongs to `L_N`.

Assume first that `u mod N` belongs to `L_N` but `v mod N` does not. Then the
second alternative in (25) is impossible for every pair. Equations (19) and
(25) imply that at most two residual characters lie off the rational `u`-
line: after fixing one off-line label, every other off-line label is, up to
sign, `c_1+u` or `c_1-u`, and both possibilities cannot coexist because their
pair has no integral coefficient-`l1`-norm-two expression for `u`. Thus at
least five residuals lie in `Q u`.

Write `u=H alpha` on its primitive line and the `r>=5` proportional residuals
as `q_i alpha`. Their coefficients `q_i` are units modulo thirteen, hence
units modulo `N`; any one of them identifies `L_N` with the reduction of
`Z alpha`. Therefore `v mod N` outside `L_N` implies `v notin Q u`, so the
blocker is genuinely transverse to the `alpha`-fiber. Also (2) makes the
scalar guard coefficient `13H` odd. If
`r<=6`, THM-2080 gives a scalar phase safe for these residuals and safe for
the guard `13H`; if `r=7`, THM-2128's seven-comb theorem gives the same strict
escape. On the fiber of `alpha` over this phase, the blocker and the at most
`7-r` remaining residuals are transverse. Their at most three danger bands
have total fiber measure at most `3/7`, so a point avoids all of them. This
contradicts (3).

Assume instead that `v mod N` belongs to `L_N` but `u mod N` does not. Now the
first alternative in (25) is impossible. Equation (20) shows that, after one
off-line residual is fixed, all other six equal `+/-v`; hence at least six
residuals lie on `Qv`. A proportional residual has unit coefficient modulo
`N`, so it identifies `L_N` with the primitive rational `v`-line. Hence
`u mod N` outside `L_N` implies `u notin Qv`; the guard is genuinely
transverse. Settled LRC for at most eight speeds supplies a scalar
phase safe for the blocker `13v` and all six or seven proportional residuals.
On its fiber, guard danger has measure `2/7` and the at most one remaining
terminal danger has measure `1/7`. A point is guard-safe and avoids every
terminal, again contradicting (3).

Thus every residual projective line that contains one of `u,v` contains both:

```text
u mod N in L_N                 and v mod N in L_N.    (26)
```

## 6. A digit section lifts the sevenfold pencil forever

The base line `L_1` is supplied by (22). Equation (25) puts at least one of
`u,v` in it, and Section 5 puts both in it. Inductively suppose (26) holds for
`N=13^k`. Choose a lattice basis `a,b` and write

```text
c_i=R_i a+N h_i b,                  13 not|R_i,
u=Aa+N h_0b,
v=Ba+N h_*b.                                         (27)
```

For `p,q in {0,...,12}`, define the digit section

```text
a.z_pq=p/13,                       b.z_pq=q/(13N).    (28)
```

It is not a subgroup, but every relevant evaluation is periodic in the
indices:

```text
c_i.(X+z_pq)=c_i.X+(R_i p+h_i q)/13,
13u.(X+z_pq)=13u.X,
13v.(X+z_pq)=13v.X.                                  (29)
```

Clean `E` by the 169 translations in (28), and choose `X` where the guard is
safe and the blocker is safe. Their joint safe set has measure at least
`4/7`. The seven residual masks cover the index plane `F_13^2`. By the
minimal seven-strip theorem in THM-2124, all next-digit normals

```text
[R_i:h_i] in P^1(F_13)                               (30)
```

are equal. As in THM-2131, a unipotent basis change lifts the common residual
line from `N` to `13N`.

Reducing (25) modulo `13N` puts at least one of `u,v` in the lifted line;
Section 5 puts both there. Induction therefore proves common alignment through
every power of thirteen. Every determinant among

```text
u,v,c_1,...,c_7                                      (31)
```

is divisible by all `13^k`, hence is zero. All nine characters in (21) have
rational rank one, contradicting the ambient rank-two assumption. Thus

```text
b=1 is impossible.                                   (32)
```

## 7. The two scalar tails and their gcd invoices

Return to the scalar containment (15) from the `b=7` branch. Divide a common
power of thirteen from `H,r_1,...,r_7`. Then:

* all seven `r_i` cannot be units, by THM-2128;
* `13` cannot divide `H`: suppose instead that exactly `k` of the `r_i` are
  units, where `1<=k<=6` after common descent. If a phase in `C_(H/13)` is
  safe for every divided nonunit tooth, then all thirteen of its roots are
  guard-safe. The `k` unit masks cover at most `2k` roots, so they cannot cover
  all thirteen. Hence every such phase must be dangerous for one of the
  `7-k<=6` divided nonunit teeth. This forces a smaller comb containment with
  at most six teeth, contrary to THM-2080;
* at least five `r_i` are units. If at most four were units, settled LRC would
  make every divided nonunit tooth safe; on thirteen roots the unit guard
  forbids at most four points and the unit teeth at most eight, leaving an
  escape.

Thus exactly one of the following holds:

```text
(I)  13 not|H; six unit coefficients q_1,...,q_6;
     one coefficient 13v;

(II) 13 not|H; five unit coefficients q_1,...,q_5;
     two coefficients 13v_1,13v_2.                   (33)
```

In case (I), every four-subset `I subset {1,...,6}` obeys

```text
gcd(q_i:i in I) divides v.                            (34)
```

Indeed, otherwise the common kernel in the scalar circle of those four unit
characters contains a phase on which `v` is strictly safe. Perturb within the
open singleton/deep-safe neighborhood to avoid the inherited null exceptional
set, every interval boundary, and all thirteen translated null sets. On its
thirteen roots, the guard forbids at most four points, the four chosen unit
teeth are singletons, and the other two unit teeth forbid at most four points.
Total capacity is twelve, a contradiction.

In case (II), every pair `i!=j` obeys the two-colour alternative

```text
gcd(q_i,q_j) divides v_1
 or gcd(q_i,q_j) divides v_2.                         (35)
```

If neither divisibility held, the cross lemma on the common scalar kernel of
`q_i,q_j` would give a phase where both `v_1,v_2` are safe. The guard uses at
most four root points, the two chosen teeth are singletons, and the other
three unit teeth use at most six points, again totaling twelve. A small
perturbation inside this open singleton/deep-safe neighborhood avoids the
inherited null exceptional set, all interval boundaries, and all thirteen
translated null sets while preserving those strict capacities.

Equations (34)--(35) are necessary arithmetic obligations extracted here: a
fourwise gcd ledger in case (I), and a two-coloured pair-gcd graph in case
(II). They are not sufficient for the scalar containment. For example,

```text
H=7,       (r_i)=(1,2,3,4,5,6,13)                   (36)
```

pays every fourwise invoice, but `t=3/8` is uncovered. Likewise

```text
H=7,       (r_i)=(1,2,3,4,5,13,26)                  (37)
```

pays every pair-colour invoice, but `t=1/6` is uncovered.

## 8. Scope, flood tails, and Tournament Analysis

Together, THM-2131 and this theorem leave no unconstrained normalized guard-
blocker pencil. The sole simultaneous-blocker survivors are the scalar profiles (33),
which still require an actual impossibility proof. In particular, the theorem
does not close the non-guard-blocker fivefold pencil, finite rank-seven boxes,
ranks nine through twelve, or LRC(14).

The challenged assumption was that deeper terminal content only reduced root
capacity. Seven deep blockers instead create a scalar flood tail: the
transverse eighth terminal becomes irrelevant, while the remaining unit teeth
must pay a higher-order gcd ledger. The faithful carrier in case (II) is a
two-coloured graph on the five unit labels, with colour `s` recording whether
`gcd(q_i,q_j)|v_s`. This is not a tournament: an edge may have both colours,
and orienting it discards the divisibility predicate. In case (I), even the
graph loses information because the obligations live on four-subsets. The
correct objects are a labelled pair cover and a four-uniform divisibility
hypergraph, together with the scalar comb phases and strict root-capacity
sidecar.
