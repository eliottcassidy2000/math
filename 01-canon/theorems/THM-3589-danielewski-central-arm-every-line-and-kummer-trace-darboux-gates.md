---
id: THM-3589
title: "Danielewski central-arm every-line and Kummer-trace Darboux gates"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On every nontrivial critical-value compiler target, a hypothetical
  polynomial Darboux pair restricts to immersed noninjective polynomial
  curves on both central arms.  Each coordinate has degree at least two on
  each arm, so its weight support reaches from at most -2n to at least 2.
  The transverse central-arm jet forces a third weight, 1 or -n, in each
  coordinate; hence every two-weight coordinate is impossible, and an
  exact three-by-three pair must use opposite middle weights.  In particular
  no affine-linear target coordinate can have a polynomial Darboux mate.
  On every smooth nonconstant squarefree surface c^n e=Sigma(b), no
  nonconstant element of C[b,ce] has a polynomial Darboux mate.  The
  nonconstant-arm branch has no rational mate; the constant-arm branch is
  sharply only polynomial, as f=b has an explicit rational mate.  No
  Darboux pair and no counterexample to JC(2) is constructed.
source: root / planar-Jacobian counterexample hostile session, 2026-08-21
audit: >
  PASS.  An independent hostile audit rederived the every-line degree
  invoice, homogeneous arm survival, transverse central jet, complete
  C[b,ce] split, Kummer-residue obstruction, and sharp rational and
  constant-Sigma hostiles.  Ordinary and optimized runs are byte-identical
  to the stored 711-gate transcript; the AST has no assertion gates.
depends_on:
  - THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier
related:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-3586-nodal-cylinder-cap38-width-period-and-second-conductor-keller-gates
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
script: 04-computation/jc2_danielewski_central_arm_kummer_trace_thm3589.py
output: 05-knowledge/results/jc2_danielewski_central_arm_kummer_trace_thm3589.out
script_sha256: 3ec1c2c10381edabbe2867541509fe010b65de728bd5383e9502ff69b7b4b38b
output_sha256: 3672e567af811d53122c7cb557511e378ba98419f3f9fb5b8a50943ba892357a
hash_basis: raw LF bytes
---

# THM-3589 -- Danielewski central-arm every-line and Kummer-trace Darboux gates

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
gives necessary conditions for a polynomial Darboux
pair on the smooth targets of THM-3581.  It neither supplies such a pair nor
settles `JC(2)`.

Work over `C`.  For `n>=2` and a squarefree nonconstant polynomial
`Sigma in C[b]`, put

```text
Y=Spec A,                 A=C[b,c,e]/(c^n e-Sigma(b)),       (1)
```

with the symplectic Poisson bracket

```text
{b,c}=c^n,              {c,e}=-Sigma'(b),
{b,e}=-n c^(n-1)e.                                      (2)
```

The two results below use different hypotheses.  The every-line theorem
uses the nontrivial compiler map of THM-3581.  The stable-subalgebra theorem
uses only `(1)--(2)`.

## 1. The two central arms of a critical-value compiler

Use the notation of THM-3581.  Thus `S` is nonconstant and squarefree,
`S(0)!=0`, all its critical values are nonzero, and

```text
P+n tP'=S^(n-1),       B=tP^n,       E(0)!=0,
W=E(B)/S^n,            t=x^nq,                              (3)

Phi(x,q)=(b,c,e)=(B(t),xP(t)S(t),qW(t)).                    (4)
```

The central fibre `b=0` in `Y` is the union of the two affine lines

```text
L_c={b=e=0},                     L_e={b=c=0}.              (5)
```

Since `P(0)=S(0)^(n-1)` and `W(0)=E(0)/S(0)^n`, the two source coordinate
lines map isomorphically onto them:

```text
Phi(x,0)=(0,S(0)^n x,0),
Phi(0,q)=(0,0,E(0)q/S(0)^n).                              (6)
```

The nonconstant `S` makes the generic degree of `Phi` exceed one; more
concretely, every nonzero point of `L_e` has the full collision fibre from
THM-3581.  Thus `Phi` is etale and noninjective.

## 2. Every-line theorem and nonlinear boundary debt

Suppose

```text
F,G in A,                         {F,G}=1.                (7)
```

Then

```text
H=(F,G) o Phi:A2_(x,q)->A2                              (8)
```

is a planar Keller map and is noninjective because `Phi` is.  Gwozdziewicz's
cited theorem says that a complex planar Keller map which is injective on
one affine source line is an automorphism.  Therefore `H` is noninjective on
every affine source line, in particular on `q=0` and `x=0`.  Equations `(6)`
give

```text
(F,G)|L_c and (F,G)|L_e are noninjective polynomial maps A1->A2.  (9)
```

They are also immersed.  Indeed `(7)` makes `(F,G):Y->A2` etale, so its
differential cannot kill the nonzero tangent of either smooth line.

We use the elementary sharp lemma:

> If `(f,g):A1_C->A2_C` is immersed and noninjective, then
> `deg f>=2`, `deg g>=2`, and `max(deg f,deg g)>=3`.

If, say, `f` has degree one, it alone is injective.  If it is constant,
immersion makes `g'` nowhere zero; over `C` this forces `g` to be affine
linear and again injective.  The roles may be exchanged.  Consequently
both degrees are at least two.  If both were quadratic, equality at two
distinct parameters would force the same midpoint for both quadratics, where
both derivatives vanish.  That contradicts immersion.  Consequently

```text
deg_c F(0,c,0), deg_c G(0,c,0) >=2,
deg_e F(0,0,e), deg_e G(0,0,e) >=2,
max of each displayed pair >=3.                        (10)
```

The familiar nodal parameterization

```text
s |-> (s^2-1,s(s^2-1))                                (11)
```

is immersed, noninjective, and has degrees `(2,3)`, so the degree-two floor
is sharp as a statement about polynomial curves.

### 2.1 Weight-support consequence

Give `A` the grading

```text
wt(b)=0,                       wt(c)=1,       wt(e)=-n. (12)
```

A homogeneous piece of nonnegative weight `r` is `c^r f(b)`.  A piece of
weight `-s<0` has the unique regular form

```text
c^(nm-s)e^m h(b),             m=ceil(s/n),
f(b)=Sigma(b)^m h(b).                                  (13)
```

On `L_c`, only nonnegative pieces survive.  On `L_e`, a negative piece
survives only when `s=nm`, in which case it restricts to `e^m h(0)`.
Equation `(10)` therefore forces, in **each** of `F` and `G`,

```text
a surviving weight r>=2,
a surviving weight -nm with m>=2, hence -nm<=-2n.      (14)
```

In particular each support has span at least `2n+2`.  This is a necessary
location gate, not a count-only Darboux nonentry theorem.  Across the pair,
one output has a surviving positive weight at least three, and one (not
necessarily the same output) has a surviving weight `-nm` with `m>=3`.

### 2.2 The transverse central-arm jet

The two arm invoices meet at

```text
p_0=(b,c,e)=(0,0,0).                                    (15)
```

Here the compiler target has `Sigma(b)=bE(b)` with `E(0)!=0`, so
`Sigma'(0)!=0`.  Linearizing `c^n e=Sigma(b)` at `p_0` gives `db=0` on
the tangent plane; thus `(c,e)` are local tangent coordinates.  Evaluating
`{F,G}=1` at `p_0` gives

```text
-Sigma'(0)(F_c G_e-F_e G_c)(p_0)=1.                     (16)
```

Consequently each row of the `(c,e)`-Jacobian is nonzero and the two rows
are independent.  The high positive piece forced by `(14)` has weight at
least two and contributes no linear `c`-jet.  The low negative piece has
weight `-nm`, `m>=2`, and contributes no linear `e`-jet.  A nonzero linear
`c`-jet can come only from a surviving weight-one piece; a nonzero linear
`e`-jet can come only from a surviving weight `-n` piece.  Therefore, in
**each** of `F,G`, the high and low pieces from `(14)` must be joined by at
least one further piece of weight `1` or `-n`:

```text
each output has at least three nonzero weights.          (17)
```

This kills every `2 x m` and `m x 2` compiler Darboux cell, regardless of
degree.  If both outputs have exactly three weights, each has exactly one
middle jet piece, and `(16)` forces their middle weights to be opposite:
one is `1` and the other is `-n`.  This is a local second-conductor gate,
not a proof that the remaining `3 x 3` or wider cells are nonempty.

Every affine-linear target polynomial

```text
alpha b+beta c+gamma e+delta                            (18)
```

has degree at most one on both central arms.  Hence no nonconstant
affine-linear target polynomial can be either coordinate of a pair `(7)`.
This conclusion holds for every nontrivial compiler carrier, not only the
degree-thirteen row of THM-3581.

## 3. The stable subalgebra C[b,ce]

Return to arbitrary squarefree nonconstant `Sigma` in `(1)`, and put

```text
z=ce.                                                    (19)
```

The exact bracket is

```text
{z,b}=(n-1)Sigma(b),                  c^(n-1)=Sigma/z.   (20)
```

We prove:

```text
For every nonconstant f in C[b,z], there is no P in A with {P,f}=1. (21)
```

The proof has two branches, and their different mate scopes are
load-bearing.

### 3.1 A nonconstant arm restriction: no rational mate

Suppose some root `beta` of `Sigma` has

```text
f(beta,z) nonconstant in z.                             (22)
```

On the generic curve `f(b,z)=w`, take `b` as a separating parameter.  At a
generic point over `b=beta`, one has `z!=0` and `f_z!=0`.  If a rational
`P in Frac(A)` obeyed `{P,f}=1`, then along this curve

```text
dP=-db/[(n-1)Sigma(b)f_z(b,z)].                         (23)
```

The right side has a nonzero simple residue at that point.  Passing from
the `(b,z)` curve to the full surface fibre adjoins the finite separable
Kummer sidecar in `(20)`.  Pullback multiplies the nonzero residue by the
ramification index; equivalently, field trace would make a nonzero multiple
of `(23)` exact downstairs.  Both are impossible because an exact rational
differential has zero residues.  Thus `(22)` excludes even a rational mate.

### 3.2 Constant arm restrictions: no regular mate

Suppose instead that `f(beta,z)` is constant in `z` for every root `beta`
of `Sigma`.  Put `h(b)=f(b,0)`.  Squarefreeness gives the exact decomposition

```text
f(b,z)=h(b)+Sigma(b)g(b,z),             g in C[b,z].    (24)
```

For a regular `P in A`, the derivation `{P,b}` vanishes on every affine arm

```text
D_beta={b=beta,c=0}.                                    (25)
```

This is immediate from `(2)`: the images of the generators under `{- ,b}`
are divisible by `c^(n-1)`.  On `D_beta`, both `{P,b}` and `Sigma` vanish,
so `(24)` gives

```text
{P,f}|D_beta
 =(h'+Sigma'g){P,b}+Sigma{P,g}=0.                       (26)
```

It cannot equal one.  This proves `(21)`.

### 3.3 Sharp rational hostile

The word `regular` in the second branch cannot be removed.  The element
`f=b` lies in `C[b,z]`, is constant on every arm, and has the rational mate

```text
P=c^(1-n)/(n-1),                       {P,b}=1.          (27)
```

Thus the Kummer sidecar can repair the constant-arm branch only by acquiring
a pole on the arms.  It cannot repair the residue debt in `(22)`.

## 4. Preservation and loss ledger

The two gates have complementary sources:

```text
every-line gate:
  source       A2 --Phi--> Y --(F,G)--> A2
  preserved    etaleness and the inherited collision
  observed     two honest affine source lines / two central target arms
  lost         the other inverse sheets and their boundary approach
  sidecar      the immersed noninjective conductor parameterization

stable-subalgebra gate:
  source       A with the Kummer coordinate c
  target       C[b,z], z=ce
  preserved    the Hamiltonian time differential on f=w
  lost         the (n-1)-st-root choice c^(n-1)=Sigma/z
  sidecar      c, whose trace preserves residues but whose poles give (27)
```

The first says any genuine compression must already encode a nodal-type
self-intersection on both simplest boundary lines.  The second says it must
also leave the natural two-variable stable subalgebra `C[b,ce]`.  Neither
condition is sufficient for a Darboux pair.

## 5. Hypothesis and scope boundaries

1. `S` nonconstant is essential in Sections 1--2.  At the degree-one plane
   boundary the compiler is injective and the every-line contradiction
   disappears.
2. The every-line theorem is over `C` and is CITED.  It is not inferred from
   the finite collision alone.
3. Section 3 needs a smooth **nonconstant** squarefree target.  Multiple
   roots produce singular arms where the symplectic and residue statements
   change type.  If `Sigma` is constant, there are no arms and the hostile
   `(27)` is regular, so the polynomial conclusion is false.
4. Section 3 excludes **polynomial** mates for all `f in C[b,ce]`; only its
   nonconstant-arm branch excludes rational mates.  Equation `(27)` is the
   mandatory hostile to any stronger wording.
5. A polynomial outside `C[b,ce]`, with nonlinear restrictions on both
   central arms and mixed Kummer dependence, remains open.

## 6. Exact verification contract

The companion checks, with explicit gates and no Python `assert`:

- the brackets `(2)`, `(20)`, and the rational hostile `(27)` for `2<=n<=8`;
- the two central-line formulas `(6)` on cyclic compiler rows
  `2<=n<=5`, `1<=r<=3`;
- exact arm decompositions `(24)` and the vanishing `(26)` for a polynomial
  basis on several squarefree targets;
- nonzero residue controls for both linear and nonlinear arm restrictions;
- the sharp immersed nodal curve `(11)`, the weight-survival invoice
  `(13)--(14)`, and the transverse jet/count gates `(15)--(17)`.

The use of Gwozdziewicz's theorem, the generic residue/trace argument, and
the universal arm decomposition are proof-driven.  The finite rows are
hostile controls, not extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_danielewski_central_arm_kummer_trace_thm3589.py
python3 -O 04-computation/jc2_danielewski_central_arm_kummer_trace_thm3589.py
```
