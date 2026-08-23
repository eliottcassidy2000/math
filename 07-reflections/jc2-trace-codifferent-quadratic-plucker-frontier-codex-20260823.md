# Trace-codifferent quadratic Pluecker frontier

**Status: FINITE-EXACT bounded obstruction and positive exterior signal;
the full quadratic rank-two intersection is OPEN.**  This reflection extends
the cheap next test in proved THM-3784 using the polynomial target observables
of proved THM-3779.  It constructs neither a polynomial Keller pair nor a
planar Jacobian counterexample.

## Universe and rank convention

For

```text
t=x^2(1+xy),
R_r=x t^r,
U=t/x+t^(m+1),
P=1/x+((m+1)/m)t^m,
V=UP,
E=P(V-1),
```

the exact probe uses

```text
W1(m)=span_Q{R_0,...,R_m,U,V,E},                    m=2,3,
W2=span_Q of generator-degree 1 and 2 monomials in W1(2).
```

`W2` is the span of the resulting source polynomials, not a free word
space.  For a source space `W`, write

```text
beta: Lambda^2 W -> Q[x,y],             F wedge G |-> J(F,G).
```

A constant in `Im(beta)` is only a linear sum of brackets.  One Darboux
pair requires a decomposable preimage, equivalently an alternating matrix
of rank two.  This distinction is load-bearing below.

## The complete linear census

The exterior ranks are

```text
             dim   wedges   rank(nonconstant)   rank(+constant)
m=2           6      15             15                 15
m=3           7      21             20                 20
```

Thus no two arbitrary linear combinations in either `W1(2)` or `W1(3)`
have a nonzero constant bracket.  The full `m=3` kernel is one-dimensional:

```text
J(R1,R2)-(1/3)J(R0,R3)=0,
```

because `J(R_a,R_b)=(b-a)x^4t^(a+b-1)`.  Its constant is zero.

The nearest target-field pair remains

```text
J(U,E)=2V-1.
```

It is axis-local only: on the three components `(x=0,A=0,B=0)` of `U=0`,
the bracket values are `(1,-1,-1)`.  For `m=2`, both functions are in the
maximal target observable, their field traces are `(3U,3E)`, and their axis
valuations are `(1,0)`.

## A genuine quadratic component equalizer

At `m=2`, `xP=1+(3/2)R2` and `E=P(V-1)` give

```text
1=V-(3/2)R2-R0 E+(3/2)R2 V.                         (1)
```

This is the cheapest positive signal from mixing codifferent rungs with
THM-3779 observables.  It repairs the component spectrum exactly:

```text
term             x=0        A=0        B=0
V                 1          0          0
R2                0          0         -1
R0 E              0         -1          1/2
R2 V              0          0          0
right side         1          1          1
```

It also passes the trace sidecar:

```text
Tr(V)=3V,  Tr(R2)=-2,  Tr(R0 E)=0,  Tr(R2 V)=-2V,
```

so the right side has trace `3=Tr(1)`.  The four axis valuations are
`0,5,1,5`.  Component spectrum, trace, and axis valuation detect different
parts of the cancellation.

The 27 raw quadratic words have source rank 24, and (1) shows that this
space already contains the constants.  Therefore the correct Pluecker
ambient space is

```text
W2/Q,                         dim=23.                 (2)
```

Exact reduction through the cubic cover gives

```text
W2 intersect k(U,P)
 =span{1,U,V,E,U^2,UV,UE,VE,E^2},             dim=9. (3)
```

Thus (1) is a mixed cancellation back into the already-proved maximal
target observable; it does not enlarge THM-3779.

## The constant exterior tensor is not a pair

First, the mixed-leg image has

```text
rank_nonconstant(W1 tensor (W2/Q))
 =rank_with_constant=56.                             (4)
```

Hence no constant-bracket pair in `W2` can have either leg in `W1`; both
would need genuine quadratic content modulo `W1`.

On the full exterior space,

```text
dim Lambda^2(W2/Q)=253,
coefficient monomials=621 through total degree 91,
rank(nonconstant)=111,
rank(+constant)=112,
dim beta^(-1)(1)=141.                                (5)
```

A sparse exact point of this affine fibre is

```text
 J(R0,E)
+(9/2) J(R0,R1 U)
+18    J(R0,R2 E)
-30    J(R1,R1 E)
+(45/8)J(R1,U^2)
+(9/2) J(R0 R1,U V)
=1.                                                   (6)
```

Its alternating coefficient matrix on `W2/Q` has rank **6**, and its
four-Pfaffian on `(R0,R1,E,R1 E)` is `30`.  It is therefore a sum of at
least three decomposable brackets, not a Darboux pair.  Formula (6) is the
hostile against the invalid implication

```text
constant in Im(beta)  ==>  decomposable constant preimage.
```

## Exact Pluecker slices and the open stopping point

Every constant-bracket two-plane has an invertible linear jet at the
origin.  Modulo constants it has a unique normalized ordered basis

```text
F=x+h,                       G=E+k,
```

with `h,k` in the 21-dimensional zero-linear-jet kernel of `W2/Q`.  Thus
the full rank-two problem is a 42-variable bilinear chart, or equivalently
the intersection of (5) with the rank-two cone cut out by
`binom(23,4)=8855` four-Pfaffians.

The exact primary probe exhausts all 441 one-correction slices

```text
F=x+a h_i,                    G=E+b h_j.
```

It linearizes with `c=ab`, solves the coefficient equations over `Q`, and
then imposes the Segre equation `c=ab`.  No slice survives.  Adjoining the
missing source coordinate `y` is the positive control:

```text
J(x,y)=1,                       skew rank=2.
```

An independent subaudit reproduced (1)--(6) and supplied these additional
stopping diagnostics:

1. **FINITE-EXACT:** the fixed-`x` derivative map `G |-> G_y` has rank 21,
   while adjoining the constant target `1` raises it to 22.  Hence no
   `G in W2` satisfies `J(x,G)=1`.
2. **FINITE-MODULAR diagnostic, not a theorem:** over `F_1009`, the span of
   `c0 wedge K` and `K wedge K` has rank `8178` in the 8855-dimensional
   fourth exterior power, and `c0 wedge c0` lies in that span.  Therefore
   the lowest single quadratic/Pfaffian separator fails at this prime.
3. **FINITE-MODULAR diagnostic, not a theorem:** a degree-three XL system
   over `F_101` produced no unit certificate (`9090` generated rows, final
   pivot rank `4956`).

The modular failures do not exhibit a rank-two point and do not prove that
one exists.  They show why another exterior-rank statistic is unlikely to
settle the problem.  The full 42-variable rank-two intersection remains
**OPEN**.  The next honest attack is saturated elimination in the normalized
chart, stratified by highest degree, component restrictions, and axis
valuation.

## Trace-Poisson boundary

For a finite separable Poisson extension `L/K` and `b in K`, the derivation
`delta_b={-,b}` is the unique extension of its restriction to `K`.  Trace
commutes with that derivation:

```text
Tr({a,b})={Tr(a),b}.                                  (7)
```

If `Tr(a)=0` and `{a,b}=c in K`, then

```text
0=Tr(c)=[L:K]c,
```

so `c=0` in characteristic zero.  Every trace-zero rung `R_r`, `r<m`, and
every trace-zero linear combination therefore has no nonzero base-field
constant mate.

This does **not** force both legs outside `K` for an arbitrary mixed element
of nonzero trace.  In that case (7) only gives the necessary condition
`{Tr(a),b}=[L:K]c`, which need not contradict anything.  The broader
quantifier must remain open.

## Exact scope and reproduction

```text
no pair in W1(m), m=2,3                         FINITE-EXACT
no m=2 W2 pair with either leg in W1            FINITE-EXACT
no normalized one-correction pair in m=2 W2     FINITE-EXACT
constant exterior tensor (6), skew rank six     FINITE-EXACT
full m=2 W2 rank-two intersection                OPEN
arbitrary Darboux pair / JC(2)                   OPEN
```

Run from the repository root:

```powershell
$env:PYTHONHASHSEED='0'
python 04-computation/jc2_trace_codifferent_quadratic_plucker_probe_20260823.py
python -O 04-computation/jc2_trace_codifferent_quadratic_plucker_probe_20260823.py
```

Python `3.13.14`, SymPy `1.14.0`; normal and optimized outputs are identical
and their LF-normalized transcript byte-matches the frozen result.

```text
primary script
  1ad8ee283130b608253ad5fcd85f0f50d1b14b23feb2ca6c90fe269ab9b26722
primary frozen output
  9fefce0cecbf6e37a18b2c41f9beff09522007059273e2b7defccc5fbbf985c8
primary active gates
  CHECKS=23

independent subaudit script
  ebb9a18337fd1d2bf9a44d27527ca86842994c0f3acf76dfdf1ee47a6807a22c
independent subaudit output
  b4d280181a8b3d083d6c3928fd425cffa54233263fce765f5b79b25d4775bfb7
independent subaudit report
  c872cf2a3b122de63f8473a054245e54db69b5490bb559f1b22ef511fb091a69
```
