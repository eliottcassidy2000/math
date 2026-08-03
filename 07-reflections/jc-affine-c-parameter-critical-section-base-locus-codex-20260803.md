# Affine-c critical-section base locus and square discriminant

**Status:** FINITE-EXACT + INDEPENDENTLY AUDITED PARTIAL SCOUT on the one-parameter slices
`C=c+x`, `E'=1` in the two THM-3212 accessory fields.  The exact
`gcd(a,H)` locus and its first normal structure are determined.  Only a
bounded rational census is made for the remaining repeated-`H` locus.  This
note is not canon and supplies no inverse cover or conclusion about `JC(2)`.

## Inheritance pass

The closest proved mechanism is
[THM-3289](../01-canon/theorems/THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go.md):
after the finite clutch passes, the coupled affine lane has a degree-53
saturated resultant `H`.  The closest positive control is the previous
[linear-subresultant scout](jc-affine-c0-e0-linear-subresultant-critical-section-codex-20260803.md),
where `c=1` makes the penultimate row `a(x)y+b(x)` a global graph because
`gcd(a,H)=1`.

The hostile near miss was the suggestion that a homogeneous pair `[-b:a]`
might automatically repair a zero of `a`.  It does not: the exceptional
locus below has `a=b=0`, so this is a base locus, not a second affine chart.
The least-used relevant sidecar is the first-normal-map discipline of
[THM-2985](../01-canon/theorems/THM-2985-multiparameter-normal-map-and-arc-factor-separation.md):
the lawful next object at a transverse base point is its normal direction or
blow-up, not division by a coefficient seen on one chart.

The live concept board was:

| Object | Preserved datum | Lost datum / next test |
|---|---|---|
| `a(x)y+b(c,x)` | primitive linear subresultant ideal | fails to choose a root when both coefficients vanish |
| `K_i[x]/(a)` | exact norm and parameter multiplication | forgets the other 17 residual roots |
| `Res_x(ST,H)` | owner-boundary incidence | mixes finite-gate and live `S` contact until factored |
| `Res_c(H,H_x)` | all repeated-`H` incidences over `x` | its degree-119 residual still needs projection to parameter space |

## The quotient-algebra square

Let `K_i` be either cubic accessory field and specialize

```text
C=c+x,             d=C'=1,             k=E'=1.
```

After removing the same degree-44 boundary as in THM-3289, exact symbolic
arithmetic gives

```text
H(c,x)=H_0(x)+c H_1(x)+c^2 H_2(x),
deg_x(H_0,H_1,H_2)=(53,52,36),

S_1=a(x)y+b_0(x)+c b_1(x),
deg a=36,          deg b_0=44,          deg b_1=36.       (1)
```

The polynomial `a` is squarefree and irreducible over each `K_i`.  In the
field

```text
A_i=K_i[x]/(a),
```

both `b_1` and `H_2` are units.  Put

```text
c_*=-b_0/b_1 in A_i.                                     (2)
```

Its reduced polynomial representative has degree 24.  The load-bearing
identities are

```text
b = b_1(c-c_*)                         mod a,
H = H_2(c-c_*)^2                       mod a.             (3)
```

Thus the apparently large parameter resultant is a quotient-algebra norm,
not a new Sylvester computation.  The characteristic polynomial of
multiplication by `c_*` is a monic polynomial `D_i(c)` satisfying

```text
deg D_i=36,       D_i squarefree and irreducible over K_i,

Res_x(a,b)=unit * D_i,
Res_x(a,H)=unit * D_i^2.                                 (4)
```

Only zero loci and multiplicities in `(4)` are invariant; the displayed
monic `D_i` is the chosen normalization.  Raw PRS rows are never identified
across normalizations, in accordance with MISTAKE-360.

Consequently, over an algebraic closure,

```text
gcd_x(a,H(c,x)) != 1       iff       D_i(c)=0.            (5)
```

Because `D_i` is irreducible of degree 36 even over `K_i`, it has no rational
parameter root.  This upgrades the isolated `c=1` graph to every rational
`c` in this slice as far as the `a`-unit gate is concerned.

## What happens on `D_i`

At a `D_i` point, equations `(2)--(3)` give

```text
a=b=0.
```

Hence `[-b:a]=[0:0]`: a homogeneous coefficient pair does **not** provide a
two-chart selector.  The preceding degree-two subresultant row survives; all
three of its coefficients are nonzero in `A_i`.  Therefore the two localized
cubics have gcd of degree exactly two in `y`, rather than a hidden cubic gcd.

The base locus is nevertheless clean.  The Jacobian of `(a,b)` with respect
to `(x,c)` is triangular at a base point, with determinant

```text
a'(x)b_1(x).
```

Both factors are units modulo `a`.  Thus the coefficient-map base locus is
transverse.  Blowing up the ideal `(a,b)` resolves the rational projective
map `[-b:a]` and records a first normal direction, exactly the kind of sidecar
suggested by THM-2985.  It does not by itself select either of the two common
`y` roots, construct a second target coordinate, or produce an inverse.

The same locus is genuinely discriminantal for `H`.  Holding `c` fixed and
reducing at `(a,c-c_*)` gives

```text
H_x=0,                 H_xx is a unit.                   (6)
```

So the distinguished residual root has exact `x`-multiplicity two.  Since
`H` and `H_x` are quadratic in `c`, eliminate `c` by writing

```text
P=H_2 H_0'-H_2'H_0,
Q=H_2 H_1'-H_2'H_1,
R=H_1 H_0'-H_1'H_0,
E_x=P^2-QR=Res_c(H,H_x).                                 (7)
```

In both fields,

```text
deg_x(P,Q,R,E_x)=(88,87,104,191),
E_x=a^2 E_res,             deg E_res=119.                (8)
```

The `a^2` factor is exact.  Factoring and projecting the degree-119 residual
back to `c` was deliberately stopped: an exploratory exact field gcd/factor
step on the degree-191 eliminant exceeded 90 seconds before producing a row,
and a direct, lower-degree multivariate parameter resultant had already shown
severe coefficient swell.  No claim is made about the full residual
parameter discriminant.

## Boundary and finite-clutch walls

Let `g=ST=rad(V)`.  The parameter resultants are much smaller:

```text
B_i(c)=Res_x(g,H),             deg B_i=6,
F_i(c)=Res_x(g,1-A'(c+x)),     deg F_i=5.                (9)
```

Both are squarefree.  Their factor-degree profiles over `K_i` are

```text
B_i: (1,1,1,1,2),       F_i: (1,1,1,2),
B_i = unit * F_i * L_i,
gcd(D_i,B_iF_i)=1,                                        (10)
```

where the linear factor `L_i` is exactly the live `q_3` wall at `S`.  Thus
the degree-36 subresultant degeneracy is disjoint from both the finite clutch
and all owner-boundary contacts.

Write `C_S=c+x_S` at the simple `S` root and
`V=v_1t+v_2t^2+...`.  The three named local values are

```text
C_S^(finite)=1/2,
C_S^(live)=-1/2-2v_2/(3v_1^2),
C_S^(exceptional)=-v_2/(3v_1^2).                         (11)
```

The apparent `q_4` exceptional-slope wall is their exact midpoint:

```text
C_S^(exceptional)
 = (C_S^(finite)+C_S^(live))/2.                          (12)
```

It is neither a root of `B_i` nor of `F_i`, and `D_i` avoids all three named
values.  This explains geometrically why the exceptional `q_4` factor is not
a third boundary component: it lies between the two actual `S` contacts.
The symbol `c` in THM-3289's local formulas means `C_S`, not the global
intercept used in this scout.

## Rational hostile census

The exact universe was

```text
c=p/q in Q,     gcd(p,q)=1,     max(|p|,q)<=12,
```

which has 183 parameters and includes `c=1` and `c=1/2`.  For each parameter,
squarefreeness of the degree-53 `H(c,x)` was certified after reduction through
a simple residue root of the accessory cubic:

```text
4111: p=1013, alpha -> 29,
3211: p=1013, alpha -> 121.                              (13)
```

The degree remains 53 and `gcd(H,H')=1` in each residue field.  A nonzero
modular discriminant proves the characteristic-zero discriminant is nonzero;
this is an exact certificate, not probabilistic evidence.  Two further good
residue embeddings per field are frozen as hostile fallbacks, although the
first embedding certifies all 183 cases.

There are no repeated-`H` parameters in this rational universe.  The only
wall hits are

```text
c=-2, -3/2,
```

and in both accessory fields they lie simultaneously on `B_i` and `F_i`.
They are finite-clutch failures, hence already supply the explicit owner
critical point of THM-3289; they are not counterexamples to the live graph
lane.  No rational value hits `D_i`.

## Connection contract and remaining frontier

The exact map is

```text
(H(c,x), a(x)y+b(c,x))
   -> multiplication by c_* on K_i[x]/(a)
   -> characteristic polynomial D_i(c).                 (14)
```

It preserves the `a`-degeneracy locus, its norm multiplicity two, and the
degree-two common-`y` fibre.  It destroys the other residual roots and the
choice between the two `y` roots on `D_i`.  The first restoration sidecar is
the transverse blow-up of `(a,b)`; the cheapest next test is to pull the
surviving quadratic subresultant to that exceptional divisor and ask whether
its two roots split canonically or have nontrivial monodromy.

What remains open is substantial:

- factor and project the degree-119 residual repeated-`H` eliminant;
- audit whether the blow-up quadratic splits over the exceptional divisor;
- deform `d`, `k`, `V`, or `A` rather than fixing `d=k=1`;
- construct a second-coordinate cofactor or inverse cover.

Nothing here proves an inverse, a Keller mate, `JC(2)`, or `DC(2)`.

## Reproduction

Run

```text
python3 04-computation/jc_affine_c_parameter_section_degeneracy_partial_scout_20260803.py
python3 -O 04-computation/jc_affine_c_parameter_section_degeneracy_partial_scout_20260803.py
```

and compare LF-normalized output with
`05-knowledge/results/jc_affine_c_parameter_section_degeneracy_partial_scout_20260803.out`.

The frozen SHA-256 digests are

```text
script  6448f0a8d8238358adca613610cede3fca72dac210ba0487d126b6d466849697
output  28e0395a4e90be88ebd413f27ebdaf82ee4ba51233d72fca42a423448bd6faa5
```

The independent companion
`04-computation/jc_affine_c_parameter_section_degeneracy_independent_audit_20260803.py`
imports or executes no primary scout.  It separately reconstructs the generic
subresultants, both quotient algebras, the two characteristic polynomials,
norm identities, transverse base ideal, quadratic fibre, exact double root,
and named wall factorization.  Normal and optimized runs byte-match
`05-knowledge/results/jc_affine_c_parameter_section_degeneracy_independent_audit_20260803.out`.
Its frozen digests are

```text
independent script  4afadde9302723c5af1a9a209525733c453f928332f21d82e1125f3ddb5662f8
independent output  d72ed2fb78bd3e237ea5893a75573d54f2c58fa9080799a8633e73547aaba81d
```
