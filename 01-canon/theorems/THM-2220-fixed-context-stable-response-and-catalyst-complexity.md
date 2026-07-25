---
id: THM-2220
title: "Fixed-context stable response and catalyst complexity"
status: >
  PROVED (abstract commutative nonexpansive metric-monoid theorem) + CITED
  APPLICATION (Gordian distance and Brittenham--Hermiller). Every fixed
  context response rho_(nx)(y) has the same asymptotic slope ell_hash(x),
  with finite-scale error at most 2ell(y). Catalytic capacity along powers
  is sublinear, and an eta-optimal catalyst saving kappa(nx)-eta has length
  at least half that saving. Absolute-cost response ideals obey a lax
  monoidal law; on any fixed finite prime-knot alphabet they have finite
  antichain bases. Exact diagonal convolution fails already for T(2,7) and
  its mirror. The full min-plus continuation kernel is the required exact
  sidecar. No positive knot catalyst is produced.
source: klein-2026-07-24-fixed-context-stable-response
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
---

# THM-2220 -- fixed-context stable response and catalyst complexity

## 1. Setup

Let `(M,+,0)` be a commutative monoid with a metric `d` satisfying joint
nonexpansivity

```text
d(a+c,b+d) <= d(a,b)+d(c,d).                         (1)
```

Put

```text
ell(x)=d(x,0),
rho_x(y)=d(x+y,y),
ell_cat(x)=inf_(y in M) rho_x(y),
ell_hash(x)=lim_(n->infinity) ell(nx)/n.             (2)
```

The last limit exists by subadditivity.  A context `y` therefore defines a
subadditive response pseudolength

```text
ell^[y](x)=rho_x(y).                                 (3)
```

## 2. Fixed-context stable-response theorem

> **Theorem.** For every `x,x',y in M` and every positive integer `n`:
>
> ```text
> rho_(x+x')(y) <= rho_x(y)+rho_(x')(y),             (4)
>
> ell_cat(nx) <= rho_(nx)(y) <= ell(nx),             (5)
>
> 0 <= ell(nx)-rho_(nx)(y) <= 2ell(y).               (6)
> ```
>
> Consequently the sequence `n |-> rho_(nx)(y)` is subadditive and
>
> ```text
> lim_(n->infinity) rho_(nx)(y)/n
>   =inf_(n>=1) rho_(nx)(y)/n
>   =ell_hash(x),                                    (7)
> ```
>
> independently of the fixed context `y`.

### Proof

Insert the intermediate point `x'+y`.  The first leg is obtained by adding
the same summand `x'` to the two endpoints defining `rho_x(y)`:

```text
rho_(x+x')(y)
 =d(x+x'+y,y)
 <=d(x+x'+y,x'+y)+d(x'+y,y)
 <=rho_x(y)+rho_(x')(y).                             (8)
```

This proves (4), and taking `x=mx`, `x'=nx` proves subadditivity in the copy
number.

The lower inequality in (5) is the definition of `ell_cat`.  Adding `y` to
both endpoints defining `ell(nx)` gives the upper inequality:

```text
rho_(nx)(y)=d(nx+y,y)<=d(nx,0)=ell(nx).              (9)
```

For the other side, use the three-leg path

```text
nx -> nx+y -> y -> 0.
```

The two exterior legs have cost at most `ell(y)`, so

```text
ell(nx)
 <=d(nx,nx+y)+d(nx+y,y)+d(y,0)
 <=2ell(y)+rho_(nx)(y).                              (10)
```

Equations (9)--(10) prove (6).  Divide (6) by `n` and use
`ell(nx)/n -> ell_hash(x)`.  This gives the limit in (7); Fekete's lemma
applied to (4) gives the infimum formula.  QED.

Thus all fixed-context lengths `ell^[y]` can differ at finite scale, but
their homogeneous shadow is the same.  A fixed catalyst changes the
`n`-copy cost by at most an additive endpoint charge `2ell(y)`, never by a
new linear slope.

## 3. Catalytic capacity is sublinear along powers

Define

```text
kappa(x)=ell(x)-ell_cat(x)>=0.                       (11)
```

Then

```text
kappa(nx)/n -> 0.                                    (12)
```

Here is a short proof which does not assume (12) from the
localization--homogenization commutation theorem.  Fix positive `k` and a
context `c`.  Iterate (4) and then apply the endpoint estimate (10) to get

```text
ell(nkx) <= n rho_(kx)(c)+2ell(c).                   (13)
```

After division by `n`, use the subsequence homogeneity

```text
lim_(n->infinity) ell(nkx)/n=k ell_hash(x)
```

and pass to the limit:

```text
k ell_hash(x)<=rho_(kx)(c).                          (14)
```

Infimize over `c`:

```text
k ell_hash(x)<=ell_cat(kx)<=ell(kx).                 (15)
```

Divide by `k` and let `k` tend to infinity.  The two exterior terms converge
to `ell_hash(x)`, hence so does the middle term.  Subtracting
`ell_cat(kx)/k` from `ell(kx)/k` proves (12).

There is also an exact complexity invoice.  If a context `y_n` saves

```text
s_n=ell(nx)-rho_(nx)(y_n),                           (16)
```

then (6) gives

```text
0<=s_n<=2ell(y_n).                                   (17)
```

For an exact optimal context, when one exists,

```text
ell(y_n)>=kappa(nx)/2.                               (18)
```

For an `eta_n`-optimal context
`rho_(nx)(y_n)<=ell_cat(nx)+eta_n`, the repaired form is

```text
ell(y_n)>=max(0,kappa(nx)-eta_n)/2.                  (19)
```

Thus a growing catalytic saving requires a growing context, even though
(12) says that the saving is always sublinear in `n`.

## 4. Cost ideals form a lax monoidal invariant

For a cost threshold `m`, define the absolute-cost context ideal

```text
A_m(x)={y in M:rho_x(y)<=m}.                         (20)
```

It is an additive upper ideal: if `y in A_m(x)`, then simultaneous
translation gives `y+z in A_m(x)` for every `z`.

For subsets of `M`, use Minkowski addition.  The response inequality with
two independently chosen contexts is

```text
rho_(x+x')(y+z)<=rho_x(y)+rho_(x')(z).               (21)
```

Indeed, pass through `x'+y+z`; the first leg is a translate of the pair
defining `rho_x(y)`, and the second is at most `rho_(x')(z)` after adding
`y`.  Therefore

```text
A_m(x)+A_n(x') subset A_(m+n)(x+x').                (22)
```

Equivalently, the epigraph

```text
E(x)={(y,m):m>=rho_x(y)} subset M x R_(>=0)          (23)
```

is an upper ideal and

```text
E(x)+E(x') subset E(x+x').                          (24)
```

This is the corrected monoid-valued statement: the response epigraph is
**lax monoidal**, not exactly monoidal.

If `d` is integer-valued, let

```text
I_s(x)={y:rho_x(y)<=ell(x)-s},
sigma(x,x')=ell(x)+ell(x')-ell(x+x').                (25)
```

If `s+t>=sigma(x,x')`, (22) gives the defect-tax law

```text
I_s(x)+I_t(x')
 subset I_max(0,s+t-sigma(x,x'))(x+x').             (26)
```

If `s+t<sigma(x,x')`, the same displayed formula follows from the
universal bound `rho_(x+x')(w)<=ell(x+x')`, which gives membership in
`I_0`.  The same split proof works with the Minkowski sum on the left
replaced by `I_s(x) intersection I_t(x')`, using (4) in the first case.
The `sigma` tax is forced by converting absolute-cost bounds into saving
relative to `ell(x+x')`; no sharpness of this tax is claimed here.

On a fixed finite prime alphabet `P_1,...,P_r` for knots, every ideal in
(20) or (25), after intersection with that alphabet submonoid, has a finite
minimal antichain in `N^r` by Dickson's lemma.
Equation (22) says that pairwise sums of source antichain generators lie in
the target ideal; equivalently each such sum dominates some target-minimal
generator.  This is a finite antichain calculus, but only in the inclusion
direction.

## 5. Gordian specialization

For knots under connected sum, (7) becomes the context-independent formula

```text
u_hash(K)
 =lim_(n->infinity) d_G(nK#J,J)/n
 =inf_(n>=1) d_G(nK#J,J)/n                          (27)
```

for **every fixed knot `J`**.  Moreover

```text
u_cat(nK)<=d_G(nK#J,J)<=u(nK),                      (28)

0<=u(nK)-d_G(nK#J,J)<=2u(J),                        (29)

[u(nK)-u_cat(nK)]/n -> 0.                           (30)
```

Integer-valuedness makes optimal contexts exist, and (18) becomes

```text
u(J_n)>= [u(nK)-u_cat(nK)]/2                         (31)
```

for every optimal catalyst `J_n`.

If all `J_n` use a fixed nonempty prime alphabet and
`J_n=#_i a_(n,i)P_i`, put `U=max_i u(P_i)`.  Since

```text
u(J_n)<=sum_i a_(n,i)u(P_i)<=U ||a_n||_1,           (32)
```

with

```text
s_n=u(nK)-d_G(nK#J_n,J_n),
```

an extra saving `s_n` forces the exponent-height escape

```text
||a_n||_1>=s_n/(2U).                                 (33)
```

So an unbounded saving cannot remain inside any bounded prime-exponent box.

## 6. Brittenham--Hermiller hostile audit: why equality is false

Let

```text
K=T(2,7),              Kbar=mirror(K),              (34)
```

so `u(K)=u(Kbar)=3`, while Brittenham--Hermiller give

```text
u(K#Kbar)<=5.                                        (35)
```

Signed half-signature is additive, and its absolute difference is
Gordian-1-Lipschitz.  For `nK#J` versus `J`, and separately for the mirror
with the opposite sign, it gives the lower bound `3n`; unknotting the `n`
displayed summands gives the upper bound `3n`.  Hence for every context
`J` and every `n`,

```text
d_G(nK#J,J)=3n,
d_G(nKbar#J,J)=3n.                                  (36)
```

Thus the stable-response theorem is exact on each summand: the context
cannot contract either calibrated side.

Nevertheless, at the identity context,

```text
rho_(K#Kbar)(U)=u(K#Kbar)<=5
  <6=rho_K(U)+rho_Kbar(U).                           (37)
```

This is precisely the pure geodesic-bypass mechanism of THM-2176.  It shows
that (4), (22), and (24) cannot be promoted to equality.

More sharply, define the min-plus convolution of diagonal responses by

```text
(rho_x square rho_z)(w)
 =inf_(a+b=w) [rho_x(a)+rho_z(b)].                   (38)
```

Equation (21) always gives

```text
rho_(x+z)(w)<=(rho_x square rho_z)(w).               (39)
```

The knot monoid is conical by Schubert prime decomposition, as imported in
THM-2191, so the only decomposition of `U` is `U#U`.
Consequently

```text
(rho_K square rho_Kbar)(U)=3+3=6,
rho_(K#Kbar)(U)<=5.                                 (40)
```

The response diagonal therefore has no exact connected-sum convolution
law.  The strongest survivor is the lax epigraph law (24); exact
composition requires restoring the off-diagonal intermediate and using the
full min-plus kernel

```text
P_K(A,B)=d_G(K#A,B),
P_K tensor P_L=P_(K#L).                              (41)
```

This audit also prevents a false conclusion from (7): all fixed-context
responses having the same homogeneous shadow does **not** make
`u_hash` additive.  Brittenham--Hermiller give

```text
u_hash(K)=u_hash(Kbar)=3,
u_hash(K#Kbar)<=5.                                  (42)
```

The last inequality follows by repeating the at-most-five crossing-change
certificate independently on each of `n` copies and dividing by `n`.

## 7. Equality and failure boundaries

1. The `2ell(y)` term is an endpoint invoice, not a claim of optimality.
   Improving the constant for knots would require a knot-specific
   cancellation or concordance argument.
2. Equation (30) is only `o(n)`.  It does not imply that
   `u(nK)-u_cat(nK)` is bounded.
3. A common stable slope does not exclude finite catalysis.  Conversely, a
   homogenization gap does not produce a catalyst; THM-2191's discrete
   metric example remains the hostile converse.
4. The finite-antichain statement is per fixed prime alphabet and is not an
   effective global bound over all prime knots.
5. The theorem neither computes `u(K#Kbar)` nor produces a positive knot
   catalyst.  It identifies the exact scale law any future catalyst must
   obey.

QED.
