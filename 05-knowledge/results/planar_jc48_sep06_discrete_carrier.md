# Polynomial carrier descent is rigid, even in one direction

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
September 6, 2026. This classifies polynomial source automorphisms that
preserve the full specified carrier. It does not classify fixed individual
source polynomials, arbitrary polynomial maps, or Keller pairs. JC(2) is OPEN.

## 1. Inheritance and the new operation

Work over an algebraically closed characteristic-zero field `k`. Let

```text
R=k[x,t],       p=t(1+x^2t),       y=xtp,       u=x^2t,
A=k[p,y] subset R,                Delta=p^3-y^2.
```

The closest proved mechanism is the exact carrier theorem in
[the Hamiltonian note, §§3–4](planar_jc48_sep06_hamiltonian.md), inheriting
[THM-4308, source-normal bracket and Hasse truncation](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md).
For `S in A`, that note proves

```text
{-u/2+H,S} in A for every H in A
iff S in k+p^2 Delta A,
```

and excludes local nilpotence for every nonconstant fixed-source carrier.
The present operation is an arbitrary **polynomial discrete automorphism**
of the source, with no assumed flow, Hamiltonian generator, or isotopy.
It extends the obstruction beyond that excluded additive-group class.

The live board is: one-way carrier descent; the actual exceptional fibre;
rational canonical maps; completed formal flows; and a single fixed source
versus the entire carrier. The inherited hostile is a finite Hamiltonian
repair whose first all-order source variation has a cusp pole. The corrected
near miss is to equate a formal automorphism with a polynomial one. The
underused sidecar is the **whole fibre** of the carrier map, including its
two nonisomorphic components, rather than just its rational inverse.

The incoming [continuing7 synthesis](continuing7_20260906_synthesis.md)
provides no source-carrier theorem to import. Its singular-fibre closure
retains whole fibres when a rational chart denominator vanishes. That is a
procedural connection here: inspecting the exceptional fibre strengthens a
two-way normalizer argument to one-way descent. No mathematical map between
its Laurent response and this source ring is claimed. Targeted tracked
canon/results searches found no exact prior carrier-normalizer statement;
there is no literature priority claim.

## 2. Complete one-way polynomial descent theorem

**Theorem.** For a polynomial automorphism `Phi` of `A2_(x,t)`, the following
are equivalent:

1. `Phi^*(A) subset A`.
2. `Phi^*(A)=A`.
3. For a unique `lambda in k*`,
   `Phi^*(x)=lambda^-1 x`, `Phi^*(t)=lambda^2 t`.

In this case

```text
Phi^*(p)=lambda^2 p,    Phi^*(y)=lambda^3 y,
Phi^*(u)=u,            J_(x,t)(Phi)=lambda.              (1)
```

Consequently the only Jacobian-one polynomial source automorphism with
even one-way carrier descent is the identity. The result requires a
polynomial **automorphism**, not merely a polynomial endomorphism having
constant Jacobian.

**Proof.** Consider the actual polynomial carrier map

```text
pi:A2_(x,t) -> A2_(p,y),       pi=(p,y).
```

Its rational inverse and image are

```text
t=Delta/p^2,       x=yp/Delta,
image(pi) = {(0,0)} union {p Delta !=0}.                 (2)
```

Indeed `Delta=t p^2`; if `p!=0`, then `t!=0` and `Delta!=0`, and the
displayed inverse reconstructs the unique source point. If `p=0`, then
`y=0`. Its only positive-dimensional fibre is therefore

```text
E=pi^-1(0,0)=L disjoint_union M,
L={t=0} isomorphic to A1,
M={1+x^2t=0} isomorphic to G_m.                          (3)
```

Every other target fibre is empty or a singleton. Notice in particular
the missing nonzero cusp and line points; `pi` is not surjective.

Assume one-way descent. Since `p,y` are algebraically independent by (2),
there is a polynomial map `f:A2_(p,y)->A2_(p,y)` with
`pi Phi=f pi`. It need not initially be assumed invertible. The closed
curve `Phi(L)` lies in the fibre over `f(0)`, so (2)–(3) force `f(0)=0`.
Thus `Phi(E) subset E`. Each image component is a closed irreducible curve
contained in the union of two irreducible curves, so is an entire component:
a proper closed subset of an irreducible curve has dimension zero. The two
images are distinct because `Phi` is an automorphism. They cannot be
swapped: `A1` and `G_m` are not isomorphic, as their rings have respectively
only constant units and nonconstant units. Hence `Phi(L)=L`, `Phi(M)=M`.

The principal prime ideals of these two divisors are preserved. For
constants `c,d in k*`, therefore,

```text
Phi^*(t)=ct,       Phi^*(1+x^2t)=d(1+x^2t).
```

Restricting the second equality to `t=0` gives `d=1`. Subtract one and
cancel `t` to obtain `c Phi^*(x)^2=x^2`. Unique factorization gives
`Phi^*(x)=b x`, `c=b^-2`, with `b in k*`. Set `lambda=b^-1`. This is (1).
The converse is immediate and gives equality of carriers automatically.

The first proof attempt used Poisson divisors and assumed two-way carrier
equality. Root supplied the exceptional-fibre argument above; it removes
that assumption and requires no Poisson calculation. It has been checked
independently against the empty/singleton fibres and the no-swap step.

## 3. Exact source-form and completion consequences

Let `G_H=-u/2+H(p,y)`. A polynomial source automorphism satisfies

```text
Phi^*(G_H) in -u/2+A for every H in A                  (4)
```

if and only if it is one of the scalings (1). In fact the three tests
`H=0,p,y` already suffice: subtraction implies `Phi^*(p),Phi^*(y) in A`,
so §2 applies. Conversely (1) fixes `u` and gives

```text
Phi^*(G_H)=-u/2+H(lambda^2 p,lambda^3 y).
```

Thus a nonzero symplectic source correction cannot be completed to a
polynomial source automorphism satisfying (4), regardless of its generator
or discrete construction. This is an obstruction to the **whole-carrier
completion strategy**. A map paying one particular source polynomial may
fail (4), and is not excluded. The theorem also says nothing about whether
an arbitrary Keller endomorphism is invertible; invertibility of `Phi` is
an explicit hypothesis.

For comparison, the rational Poisson bracket on the carrier field is

```text
{F,G}=-(Delta/p) J_(p,y)(F,G).                           (5)
```

If one assumes two-way carrier equality, its polynomial automorphism has
constant carrier Jacobian `j`. Covariance under a source automorphism of
constant Jacobian `kappa` gives
`j Delta Phi^*(p)=kappa p Phi^*(Delta)`.
The distinct irreducible pole `p` and zero `Delta` must be preserved up to
units. Thus `Phi^*(p)=alpha p`; an automorphism over this coordinate has
`Phi^*(y)=beta y+h(p)`. Preserving `Delta` forces `h=0` and
`beta^2=alpha^3`. Formula (5) then gives
`kappa=alpha^2/beta`, recovering (1). This is a consistency route for
the two-way subcase, not the proof of one-way descent.

## 4. Positive formal and rational controls, with their exact losses

The inherited universal-carrier theorem already permits nontrivial
symplectic automorphisms of the completed carrier. The following explicit
control records the boundary of the new discrete obstruction. For

```text
S=p^2 Delta=t^5(1+x^2t)^4,             L={-,S},
```

literal differentiation gives

```text
L(p)=2py Delta,             L(y)=Delta(5p^3-2y^2),
L(Delta)=-4y Delta^2,       L(u)=10p^3y,
L(x)=t^4(1+x^2t)^3(5+9x^2t),
L(t)=-8xt^6(1+x^2t)^3.                                 (6)
```

Thus `L(A) subset A` and `L(-u/2) in A`. In `R[[epsilon]]`,
`exp(epsilon L)` is well defined coefficient by coefficient, is inverted
by `exp(-epsilon L)`, preserves the bracket, and takes every `G_H` into
`-u/2+A[[epsilon]]`. This proves a genuine all-order formal operation.
It is nontrivial by (6). There is no conflict with §2: formal units such
as `1+epsilon F` need not be scalar units, and an infinite time series is
not a polynomial source automorphism at a nonzero scalar time. No such
specialization is presumed to exist. The inherited factorial-closure
theorem already proves that this `L` is not locally nilpotent.

There is also a different, explicitly rational canonical operation. For
a scalar or formal parameter `w`, set

```text
Psi_w(x)=x/(1-wx),           Psi_w(t)=t(1-wx)^2.          (7)
```

It has Jacobian one, inverse `Psi_-w`, composition law
`Psi_v Psi_w=Psi_(v+w)`, and fixes `u`. It therefore fixes `G_0=-u/2`
exactly. This is a rational one-parameter action generated by `u`, which
is outside `A`; it is not a polynomial additive-group action on the
whole source. The control `H=0` is not the genuine fixed low source jet
of THM-4308.

For `w!=0`, (7) fails polynomial source regularity and carrier descent:

```text
Psi_w(p)=p(Delta-wyp)^2/Delta^2,
Psi_w(y)=y(Delta-wyp)^3/Delta^3.                         (8)
```

The numerators have nonzero residues modulo `Delta`, respectively
`w^2 p^6` and `-w^3 p^9`, so these have genuine cusp poles of orders two
and three. Fixing one source `G_0` does not imply preserving all `G_H`.
This is the first failed implication if one tries to use (7) as a repair
of the normalized whole carrier. Its surviving use is on the rational
source chart, or for the explicitly fixed control polynomial only.

A second hostile is the polynomial symplectic translation `(x,t)->(x+1,t)`.
It does not descend through `pi`: on the collapsed component
`M`, parametrized by `x=z,t=-z^-2`, its transformed `p` is
`(2z+1)/z^4`, which is nonconstant. Preserving the polynomial source and
the bracket is insufficient without the carrier condition.

## 5. Contract, stopping point, and reproduction

Source: polynomial automorphisms of the original `(x,t)` plane. Target:
endomorphisms of the specified polynomial carrier `A`. Map: polynomial
descent through the actual carrier map `pi`. Preserved data: the whole
collapsed fibre, its component isomorphism types, and, when requested,
the source Jacobian. Lost data upon passage to the fraction field: missing
points, pole regularity and the distinction between `A1` and `G_m` as
entire components. The required sidecar is (3), not just the birational
inverse (2). Root's strengthening turns this sidecar into an exact
classification rather than a sufficient obstruction.

The strongest positive survivor is the formal operation (6). The
rational operation (7) is an exact fixed-control alternative with explicitly
lost carrier regularity. A productive next question is to characterize
polynomial source automorphisms paying one **genuine fixed** `G_H` while
not preserving `A`, or to change the carrier itself. Retrying nonconstant
LNDs inside `k+Delta A`, or adding more finite rows to an all-carrier
polynomial symplectic automorphism, cannot evade the proved obstruction.

The [standalone exact source](../../04-computation/planar_jc48_sep06_discrete_carrier.py)
and [frozen output](planar_jc48_sep06_discrete_carrier.out) use no inherited
producer. Their finite universe comprises the symbolic carrier and both
inverse directions; an indeterminate nonzero scaling; the two same-fibre
translation controls; two independent rational action parameters and the
exact cusp residues; and the universal formal carrier's ordinary versus
carrier-ring iterates through order three. Divergence and the second
symplectic exponential coefficient are checked separately. This is a
fixed exact identity/control set, not a sampled automorphism census or a
finite-order substitute for the proof in §2.

```sh
python3 -B 04-computation/planar_jc48_sep06_discrete_carrier.py
python3 -B -O 04-computation/planar_jc48_sep06_discrete_carrier.py
```

Both modes pass **67 always-active gates** with byte-identical output.

```text
source SHA-256:
f748d1b7f1ce09b2a1ab7f72bf9fb397b2bb4bc3ae1a24951041ad1b90ca8dcd
output SHA-256:
528bb8a3a7ca013d17e80b32cea8cfd0b7c31a4ebc3ddd692a151f881808884e
semantic SHA-256:
250f154a1a40184013f186e918a3f2026fdf4f4a35b1f0a7e9b0c5b8ff88c7ec
```

The [independent audit](planar_jc48_sep06_discrete_carrier_audit.md) passes
the one-way descent proof, actual empty/singleton fibres, no-swap step,
source-form equivalence, independent Poisson route, formal completion
scope, rational pole orders and translation hostile. Full source review
and normal/optimized/frozen 67-gate replay also pass. Audit SHA-256:
`4ca743ca25723c93aad4e809865ea7c241e595c231508549247b6a574d536ce6`.
All three owned files are frozen. No old frozen artifact is changed,
and no new theorem ID is reserved.
