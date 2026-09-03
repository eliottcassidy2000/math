# Historical THM-4401 draft -- Casimir leaves of the rank-two Poisson counterexample

**Status: SUPERSEDED RESEARCH DRAFT.** This outline was independently audited
and promoted to canonical THM-4401. Use that theorem, rather than this
historical draft, as the current truth source. The exact primary companion is
`04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py`.

Primary input: the fully read local source
`/tmp/arxiv2608.23777.6Wrvip/explicit_rank_two_poisson_counterexample_updated.tex`,
SHA-256
`bcadd413c2086dbcb2f8a2339f6b19ac7e82926016d68681dd06f86c914ed2e2`.
The derivation below starts again from its displayed core polynomials, rather
than importing a computer transcript from the paper.

Closest inherited mechanisms are
`THM-1300-jacobian-counterexample-dixmier-A3-explicit.md` (the underlying
three-dimensional Keller core),
`THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension.md`
(the earlier exact symplectic-suspension mechanism),
`THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md`
(the planar centralizer firewall), and
`THM-3438-weighted-lift-keller-degree-spectrum.md` (degree-three and `S_3`
monodromy under stabilization).  The least-used decisive sidecars here are
`THM-3554-punctured-kummer-collision-surface-normal-form.md` and
`THM-3555-catalan-thickening-universal-cubic-root-cover.md`, together with
`THM-3561-rational-keller-danielewski-polynomial-completion.md` for the
complementary Lagrangian slice. THM-4397 was subsequently proved and is used
by the promoted theorem only to identify the two symplectic gauges.

The field `k` has characteristic zero.  Write the preprint's independent
source coordinates as `(x,y,beta,D)` and put

```text
w = 1+xy,
alpha = 2-3xy-x^2 beta,
R = 2x-3x^2y-x^3 beta = x alpha,
S = y+3xw^2 beta+3xy^2(4+3xy),
T = -(w^3 beta+y^2w(4+3xy))/2.
```

The preprint proves that these are polynomial source coordinates, that
`R,T,S` are independent of `D`, and that

```text
{D,R}=1,       {S,T}=1,
{R,S}={R,T}={D,S}={D,T}=0.
```

Thus `R` is a Casimir of the three-variable Jacobian-Poisson core and
`{R,-}=-partial_D` on the full coordinate ring.

## Candidate theorem

For `rho in k`, Hamiltonian reduction of the four-dimensional map at `R=rho`
is the surface map

```text
f_rho : X_rho=Spec A_rho -> A^2_(S,T),
A_rho = k[x,y,beta]/(R-rho).
```

Equivalently, one may fix a global translation cross-section `D=d`; the
two-dimensional slice `R=rho,D=d` induces the same map `(S,T)`.

### 1. Nonzero leaves

If `rho!=0`, then `x` is a unit because `x alpha=rho`.  Define new leaf
coordinates

```text
t=x^(-1),                    u=t(1+xy)=t+y.                 (1)
```

Here `u` is a new marked-root coordinate, not the preprint's temporary
notation `u=xy`.  The inverse formulas are

```text
x=t^(-1),       y=u-t,       beta=5t^2-3tu-rho t^3.        (2)
```

Direct substitution gives mutually inverse ring maps

```text
A_rho ~= k[t,t^(-1),u].                                  (3)
```

On this Laurent plane the reduced outputs are

```text
S=2t+4u-3rho u^2,
T=(rho u^3-u^2-ut)/2,                                    (4)
det partial(S,T)/partial(t,u)=-t.                         (5)
```

Since `t` is a unit, `f_rho` is everywhere etale.  Its source is not an
affine plane: `t` is a nonconstant unit.

### 2. Marked-root incidence and the exact deleted divisor

Let

```text
P_rho(U)=rho U^3-2U^2+S U+4T.                            (6)
```

Equations `(4)` imply

```text
P_rho(u)=0,       P_rho'(u)=2t,                          (7)
t=(S-4u+3rho u^2)/2.                                    (8)
```

Hence the same formulas on the affine completion `A^2_(t,u)` present it as

```text
k[t,u] ~= k[S,T,U]/(P_rho(U)).                           (9)
```

For `rho!=0`, `(9)` is finite flat of degree three, with basis
`1,U,U^2`.  Its Jacobian is `-t`, so the completion ramifies precisely on
the boundary divisor `t=0`; the original Casimir leaf is exactly the
simple-marked-root locus `P_rho'(u)!=0`.

For an explicit comparison with the universal marked-root cover of
`THM-3555-catalan-thickening-universal-cubic-root-cover.md`, put

```text
v = u-2/(3rho),
p = S/rho-4/(3rho^2) = 2t/rho-3v^2,
q = 4T/rho+2S/(3rho^2)-16/(27rho^3).                    (10)
```

The source change `(t,u)<->(v,p)` is triangular polynomial, the target
change `(S,T)<->(p,q)` is affine, and `(6)` becomes

```text
v^3+pv+q=0,             q=-v^3-pv.                      (11)
```

Moreover

```text
p+3v^2=2t/rho.                                          (12)
```

Thus every nonzero Casimir leaf is exactly the THM-3555 universal cubic
root cover after deleting its ramification parabola.  Over an algebraically
closed field, its set-theoretic image is the target plane minus the one
triple-root/cusp point

```text
(S,T)=(4/(3rho),-2/(27rho^2)).                           (13)
```

Off the cubic discriminant it has three simple points per fiber; over the
smooth discriminant it retains the one simple marked root; at `(13)` every
root is multiple and the punctured source has no point.

The exact discriminant is

```text
Disc_U(P_rho)
 =4(S^2-rho S^3+32T-36rho ST-108rho^2T^2),             (14)
```

and its pullback is

```text
4t^2(9rho^2u^2-8rho t-12rho u+4).                       (15)
```

The second factor in `(15)` is where the two *unmarked* roots collide.  It
must not be mistaken for ramification at the marked source point: that
ramification is exactly `t=0` by `(5)` and `(7)`.

### 3. The zero leaf and the THM-3554 boundary

At `rho=0`, the factorization `R=x alpha` is a comaximal factorization:

```text
1=alpha/2+x(3y+x beta)/2.                               (16)
```

Therefore Chinese remainders give the scheme-level disjoint union

```text
A_0 ~= k[y,beta] x k[t,t^(-1),u],
X_0 ~= A^2 disjoint-union (G_m x A^1).                  (17)
```

On the plane component `x=0`,

```text
(y,beta) |-> (S,T)=(y,-beta/2-2y^2)                    (18)
```

is a polynomial automorphism, with inverse
`y=S, beta=-2T-4S^2`.

On the curved component, `(4)` specializes to

```text
S=2t+4u,       T=-(u^2+ut)/2,
S^2+32T=4t^2.                                           (19)
```

The affine completion is now finite flat of degree two because

```text
P_0(U)=-2U^2+S U+4T.                                    (20)
```

In the normalized outputs of the preprint, `S=F_2` and `T=-F_1/2`, so

```text
S^2+32T=F_2^2-16F_1.                                   (21)
```

Equations `(19)--(21)` are exactly the punctured Kummer normal form of
`THM-3554-punctured-kummer-collision-surface-normal-form.md`: after using
`(t,S)` as Laurent source coordinates, the map is `(t,S)->(S,4t^2)`.
The known triple fiber at `(R,T,S)=(0,1/8,0)` splits as one point from the
automorphic plane component and the two Kummer points `t=+/-1` from the
Laurent component.

This is why the zero reduction keeps the `1+2` collision without producing
a planar counterexample: the whole source is disconnected, its only
`A^2` component is invertible, and the noninjective component has the
nonconstant unit `t`.

## A complementary connected two-dimensional slice

The same cubic identity gives a second useful near-descent.  Fix
`(T,S)=(tau,sigma)` with `tau!=0`.  On the core curve, `w=1+xy` is a unit
because the displayed formula for `T` writes `-2tau` as `w` times a
polynomial.  Put

```text
h=x/w,                  d=1-sigma h-6tau h^2.           (22)
```

There are exact identities

```text
y=S+6Th,                wd=1,
R=2h-Sh^2-4Th^3.                                          (23)
```

Conversely, setting `y=sigma+6tau h`, one inverse parameterization is

```text
x=h/d,
beta=-2tau d^3-4y^2d^2-3hy^3d.                         (24)
```

It follows that

```text
{T=tau,S=sigma} ~= Spec k[h,d^(-1)] x A^1_D,           (25)
(h,D) |-> (R,D)=(2h-sigma h^2-4tau h^3,D),
Jac_(h,D)(R,D)=2d.                                      (26)
```

Thus `(25)` is a connected punctured-plane symplectic slice and `(26)` is
etale.  At the distinguished values `(tau,sigma)=(1/8,0)`,

```text
d=1-3h^2/4,
R=2h-h^3/2=-h(h-2)(h+2)/2.                              (27)
```

The three collision points are the three retained roots `h=0,+/-2`; the
deleted points are `h=+/-2/sqrt(3)`.  Completing `(25)` to `A^2_(h,D)`
restores the Jacobian factor `2d`, which vanishes on precisely the two
deleted vertical lines.  This is a connected punctured planar triple-cover
counterpart to the disconnected `1+2` zero-Casimir slice, but it is still
not a map `A^2->A^2` with constant Jacobian.

## Exact descent obstructions and quantifier firewall

There are three distinct questions, and the result answers only the first
two.

1. **Literal Poisson descent preserving all four outputs is impossible.**
   The required four-by-four generator bracket matrix has determinant one
   and rank four.  For four functions on a planar symplectic algebra, their
   bracket matrix has the form `J Omega_2 J^T` and rank at most two over the
   fraction field.  Equivalently, if planar polynomials `r,d,s,t` obeyed
   `{d,r}=1`, `{r,s}={r,t}=0`, and `{s,t}=1`, then
   `THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md`
   would put `s,t` in `k[r]`, contradicting `{s,t}=1`.

2. **The legitimate symplectic reduction exists, but it is punctured.**
   For `rho!=0` its source has a nonconstant unit; at `rho=0` it is a
   disjoint plane-plus-Laurent union.  Consequently `R-rho` is not a
   polynomial coordinate hypersurface, so
   `THM-3546-invariant-graph-keller-descent-criterion.md` cannot turn this
   natural invariant slice into `A^2`.  Filling the missing divisor makes
   the Jacobian `-t`, not a constant.  This is the first exact obstruction
   after performing the correct reduction.

3. **The collision-carrying Lagrangian base change is also punctured.**
   Merely pulling back the two-form to a Lagrangian surface gives `0=0`, but
   the natural slice is stronger than an arbitrary restriction: the source
   surface `{D=0,S=0}` is the scheme-theoretic inverse image of the target
   Lagrangian plane with coordinates `(R,T)`.  Since the ambient map is
   etale, this base change is etale and contains the three-point collision.
   Under `(R,T,S)=(F_3,-F_1/2,F_2)`, it is exactly the middle-coordinate leaf
   analyzed in
   `THM-3561-rational-keller-danielewski-polynomial-completion.md`.  That
   theorem proves its source is `A^2` minus a nonempty curve, not `A^2`; its
   rational planar Keller pair has the triple collision, while polynomial
   filling changes the target to a non-plane Danielewski surface.  Thus this
   Lagrangian route reaches a genuine punctured planar near-counterexample,
   not `JC(2)`.

4. **Arbitrary destabilization remains open.**  The preprint itself gives
   the exact ordinary-map factorization
   `Phi_pt=sigma o (G x id) o Psi`; the four-dimensional map is a stable
   copy of the three-dimensional Keller core after polynomial source and
   target automorphisms.  If `G` were polynomially left-right equivalent to
   `F x id_(A^1)` for a planar map `F` (equivalently, if `Phi_pt` were
   equivalent to `F x id_(A^2)`), then `F` would be a noninvertible planar
   Keller map of generic degree three and `S_3` monodromy.  Those invariants
   are compatible with, not contradictory to, a hypothetical planar
   counterexample.  No cancellation theorem or such stable equivalence is
   proved here.

One further firewall is useful.  The original-variable coordinate
`R=x(2-3xq)` has no planar polynomial Jacobian mate by
`THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate.md`.
The two natural symplectic slices analyzed above and the canonical
Lagrangian base change are all punctured models.  This does not exclude an
accidental Keller restriction on some other, noncanonical surface; without
the inverse-image/base-change condition, Lagrangianity alone supplies no
planar Jacobian determinant.

## Strongest transferable lemma

The rank-two Poisson counterexample does not cancel the universal cubic
ramification inside an affine plane.  It **turns the derivative into a
unit by deleting its divisor**:

```text
marked-root derivative = 2t,
source leaf = {t!=0},
affine completion Jacobian = -t.                        (28)
```

This gives a precise construction target for planar JC: export the marked
cubic's ramification divisor to nonproper escape at infinity while keeping
the source ring a polynomial plane (hence with no nonconstant units).  The
current four-dimensional symplectic suspension achieves the first half only
because the reduced source is Laurent.  This is the exact connection to the
THM-3554/THM-3555 punctured-Kummer/universal-cubic program, not a proof or
counterexample for `JC(2)`.

## Replay

```bash
python3 -B 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py
python3 -B -O 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py
diff -u \
  <(python3 -B 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py) \
  <(python3 -B -O 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py)
```

The normal and optimized outputs are byte-identical. Final primary and
independent hashes are pinned in THM-4401.
