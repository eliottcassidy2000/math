---
id: THM-3595
title: "Danielewski separated transverse time-form rational nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On every smooth c^N e=Sigma(b) with N>=1 and squarefree deg(Sigma)>=2,
  a coordinate f(c)+g(e) with both f and g nonconstant has no rational
  constant-bracket mate.  On the quadratic exponent-two target, the same
  is true for every affine-linear lambda b+mu c+nu e with nu nonzero.
  The forced generic-fibre time form is nonzero holomorphic except in the
  explicit logarithmic genus-zero boundary cases, where nonzero residues
  still forbid exactness.  No Darboux
  pair and no JC(2) consequence is constructed.
source: root / planar-Jacobian construction-hostile session, 2026-08-21
audit: >
  Two independent hostile audits rederived the arbitrary-polynomial
  separated fibre, geometric connectedness, finite cancellation, the full
  infinity/Riemann--Hurwitz invoice, the unique logarithmic boundary, the
  A13 first-kind row, every affine-linear quadratic row with nu nonzero,
  and all sharp hostiles.  Normal and optimized companions are
  byte-identical to the stored 114679-gate output.
depends_on: []
related:
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-3552-two-face-cyclic-fiber-holomorphic-exactness-obstruction
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
  - THM-3591-danielewski-arm-blind-separated-darboux-nonentry
script: 04-computation/jc2_danielewski_separated_timeform_rational_nonentry_thm3595.py
output: 05-knowledge/results/jc2_danielewski_separated_timeform_rational_nonentry_thm3595.out
script_sha256: 6486357dca3b51c3787dcb1f5443f2f47b18dcf69cdf202dc63cf0cf79a5a200
output_sha256: 6ea913519271bf96b2bfd3c5e8fd0897b789487f544a36aa14dffeda4617786e
hash_basis: raw LF bytes
---

# THM-3595 -- Danielewski separated transverse time-form rational nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is a
full-function-field hostile to planar Jacobian
counterexample construction through a smooth Danielewski target.  It proves
neither the nonexistence of arbitrary Darboux pairs on that target nor any
case of `JC(2)`.

Work over `C`.  Let `N>=1`, let `Sigma in C[b]` be squarefree of degree
`D>=2`, and put

```text
Y_(N,Sigma): c^N e=Sigma(b),                            (1)

{b,c}=c^N,       {c,e}=-Sigma'(b),
{b,e}=-N c^(N-1)e.                                      (2)
```

## 1. Separated transverse coordinates have no rational mate

Let `f,g in C[T]` be nonconstant, with

```text
A=deg f>=1,                    B=deg g>=1,
Q=f(c)+g(e).                                             (3)
```

Then there are no `P in C(Y_(N,Sigma))` and `kappa in C*` satisfying

```text
{P,Q}=kappa.                                             (4)
```

This is stronger than polynomial nonentry: denominators do not help.

### 1.1 The forced time form

Let `w` be transcendental.  On the generic plane curve

```text
C_w: f(c)+g(e)=w                                        (5)
```

the standard residue form is

```text
omega=-dc/g'(e)=de/f'(c).                               (6)
```

For generic `w`, `(5)` is smooth.  Pull it back along

```text
Sigma(b)=y,                     y=c^N e.                 (7)
```

to obtain the generic `Q`-fibre in `Y`.  The Hamiltonian vector field
`D_Q={-,Q}` satisfies

```text
D_Q(c)=-Sigma'(b)g'(e),          D_Q(e)=Sigma'(b)f'(c).
```

Hence

```text
eta=omega/Sigma'(b)                                      (8)
```

satisfies `eta(D_Q)=1`.  If `(4)` held, then on the generic fibre

```text
dP=kappa eta.                                            (9)
```

The rest of the proof is the exact divisor invoice for `eta`.

There is no hidden component choice.  The polynomial `f(c)+g(e)` is
noncomposite: a hypothetical
outer derivative would divide both nonzero polynomials `f'(c)` and `g'(e)`,
whose gcd in `C[c,e]` is one.  The characteristic-zero
Bertini--Krull/Stein criterion therefore makes `(5)` geometrically integral.

The degree-`D` lift in `(7)` is geometrically connected as well.  The finite
branch values of `Sigma` are fixed and nonzero.  By `(10)`, the `y`-map of
`(5)` is unramified over each of them for transcendental `w`.  Any branch
of the `y`-map over zero is harmless, because squarefreeness makes the
`Sigma`-cover unramified over zero.  Thus the two Galois
closures have no common finite branch value.  A nontrivial characteristic-zero
cover of `P1` cannot branch only over infinity, so the covers are linearly
disjoint.  The pullback consequently has degree exactly `D` and is
geometrically integral.

### 1.2 Finite places are regular

At an ordinary affine point, `(6)` is a nowhere-vanishing local generator:
when `g'(e)!=0`, `c` is a local parameter and the first expression is a
unit; when `g'(e)=0`, smoothness of `(5)` gives `f'(c)!=0`, and the second
expression is a unit.

Now let `alpha` be a critical point of `Sigma`, put

```text
r=ord_alpha Sigma'>=1,             v=Sigma(alpha).
```

Squarefreeness implies `v!=0`.  On the curve `c^N e=v`, the restriction of
`Q` is the rational function

```text
h_v(c)=f(c)+g(v/c^N).                                   (10)
```

Its critical values are finite unless it is constant; in the constant case
the generic fibre misses this level altogether.  Thus for generic `w`, every
intersection above `y=v` is transverse.  With local parameters

```text
y-v=t,                       t=u^(r+1),
Sigma'(b)=unit*u^r,
```

the pullback of a regular differential `unit*dt` has order `r`, exactly
cancelling the order of `Sigma'` in `(8)`.  Therefore `eta` has neither a
finite pole nor a finite zero, including above every critical value of
`Sigma`.

### 1.3 The complete infinity invoice

Put

```text
delta=gcd(A,B),
M=(NB+A)/delta,                  gamma=gcd(D,M).          (11)
```

At each branch at infinity of `(5)`,

```text
ord(c)=-B/delta,          ord(e)=-A/delta,
ord(omega)=(AB-A-B)/delta-1.                             (12)
```

The map `(7)` has ramification index `D/gamma` there, while
`ord(b)=-M/gamma`.  Pulling back a differential contributes the different
`D/gamma-1`, and division by `Sigma'(b)` contributes
`(D-1)M/gamma`.  Consequently every place above infinity has exact order

```text
ord(eta)=K/(delta gamma)-1,                              (13)

K=DAB-A+B((D-1)N-D).                                    (14)
```

For `N>=2`, this integer is nonnegative.  Indeed

```text
K=A(DB-1)+B((D-1)N-D)>0,                                (15)
```

because `D,N>=2`.  Write `A=delta a`, `B=delta q`, so
`M=Nq+a`.  Modulo `gamma`, both `D` and `Nq+a` vanish, and

```text
K/delta =D delta aq-a+q((D-1)N-D)=0 mod gamma.          (16)
```

Thus `delta gamma|K`; positivity makes `K/(delta gamma)>=1`, proving the
claim from `(13)`.

The infinity ledger is the whole canonical divisor, not merely a pole
bound.  Riemann--Hurwitz gives

```text
2g_X-2
 =D(AB-A-B-delta)+(NB+A)(D-1)+delta(D-gamma)
 =K-delta gamma,                                         (16a)

div(eta)=(K/(delta gamma)-1) sum_(x at infinity) x.      (16b)
```

There are exactly `delta gamma` infinity points.  Equality of the two lines
also independently verifies that no finite zero was omitted.
Here the three Riemann--Hurwitz contributions are respectively the degree-`D`
pullback of the base canonical degree `AB-A-B-delta`, the finite ramification
`(NB+A)(D-1)`, and the infinity ramification `delta(D-gamma)`.

When `N=1`, formula `(14)` becomes

```text
K=DAB-A-B.                                               (17)
```

It is positive except for the unique tuple `(D,A,B)=(2,1,1)`, and the same
divisibility argument applies in every positive case.  In the exceptional
case, write the leading terms as

```text
f(c)=a c+...,       g(e)=d e+...,       Sigma(b)=s b^2+... .
```

There are two unramified infinity places with `b/c -> +/-rho`, where
`s rho^2=-a/d`.  Formula `(8)` has simple-pole residues

```text
+1/(2ds rho),                     -1/(2ds rho).         (18)
```

Both are nonzero, while every exact differential has zero residue.  Thus
the exceptional rational-fibre case is also nonexact.

Outside `(18)`, the differential `eta` is therefore nonzero and holomorphic
on the smooth projective normalization of the generic fibre.  In characteristic zero,
the differential of a rational function has a pole at every pole of that
function.  A rational function without poles on a projective curve is
constant.  Hence a nonzero holomorphic differential cannot be exact, in
contradiction with `(9)`.

## 2. The cheapest elliptic controls

Take

```text
N=D=2,                  Sigma=b(b+1),
Q=c+e.                                                   (19)
```

The generic fibre is the plane cubic

```text
c^2(w-c)=b(b+1),                                        (20)
```

and `eta=-dc/(2b+1)` is its regular differential.  Formula `(13)` gives
order zero at infinity.  At `w=1`, the projective cubic is smooth; thus
`(20)` is already an elliptic rational-mate obstruction in total coordinate
degree one.

For the A13 row

```text
N=D=3,                  Sigma=b(b^2+1),
Q=c^2+e^2,                                               (21)
```

one has `(delta,M,gamma,K)=(2,4,1,16)`.  There are two infinity
places, and `(13)` gives order seven at each.  There are no finite zeros, so

```text
div(eta)=7P_+ + 7P_-,                 genus=8.           (22)
```

This is a first-kind obstruction with no residue signal: a residue-only
search would miss it completely.

## 3. Adding a linear b-channel still fails in the smallest target

The preceding obstruction is not merely an artefact of omitting `b`.
Keep `N=2` and let `Sigma` be any squarefree quadratic.  For

```text
Q=lambda b+mu c+nu e,                 nu!=0.            (23)
```

there is again no rational constant-bracket mate.

On `Q=w`, eliminate `e` to obtain

```text
F_w(b,c)=c^2(w-lambda b-mu c)-nu Sigma(b)=0.            (24)
```

If `(lambda,mu)!=(0,0)`, this is a plane cubic.  Its Hamiltonian field satisfies

```text
D_Q(c)=F_b,                       D_Q(b)=-F_c,           (25)
```

so the forced time form is the cubic residue

```text
eta=dc/F_b=-db/F_c.                                      (26)
```

The generic projective cubic is smooth.  At infinity, `[1:0:0]` is smooth
because `F_Z=-nu*lc(Sigma)`.  If `lambda!=0`, the only other infinity point
has `F_B!=0`; if `lambda=0`, it is absent.  Affine singular points with
`c!=0` must satisfy

```text
lambda c^2+nu Sigma'(b)=0,
mu c^3-2nu Sigma(b)=0,                                  (27)
```

a finite system, and hence occur for only finitely many `w`.  Points with
`c=0` are smooth by squarefreeness.  Thus `(26)` is a nonzero holomorphic
differential on a smooth genus-one generic fibre and cannot equal `dP`.

The concrete coordinate `Q=b+c+e` is therefore already rationally blocked.
Its singular-value polynomial for `Sigma=b(b+1)` is

```text
16w^4+16w^3-8w^2+24w-11,                               (28)
```

so `w=1` is an exact smooth control.

If `lambda=mu=0`, `(24)` is instead a smooth conic.  At its two infinity
points `b/c -> +/-rho`, where `lc(Sigma)rho^2=w/nu`, the form `(26)` has
opposite nonzero simple-pole residues.  It is again nonexact.  Thus every
affine-linear coordinate with `nu!=0`, including `Q=e`, is rationally
blocked.

## 4. Sharp boundaries and assumption challenges

1. Both transverse channels in `(3)` are load-bearing.  For `Q=c`,
   `P=b/c^N` satisfies `{P,Q}=1` rationally.
2. The degree condition on `Sigma` is sharp.  If `Sigma=b`, then `Y=A2_(c,e)`
   with `{c,e}=-1`, and `P=-c,Q=c+e` is a polynomial Darboux pair.
3. The affine-linear corollary needs `nu!=0`: for `Q=b`, the rational mate
   `P=1/c` satisfies `{P,Q}=1` on the exponent-two target.
4. The exponent boundary is sharp.  At `N=0`, the relation is
   `e=Sigma(b)` and `{b,c}=1`; hence `P=b,Q=c+e` is a polynomial Darboux
   pair.  If instead `Sigma=a in C*` is constant, put
   `H={b,c+e}=c^N-Nc^(N-1)e`; then `{H,c+e}=0`, so
   `P=b/H,Q=c+e` is a rational Darboux pair.  Squarefreeness is the smooth
   target/type boundary; no singular extension is claimed.  Characteristic
   zero is load-bearing: in characteristic `p`, `Sigma=b^p-b`,
   `Q=c+c^p+e`, `P=-e` is a polynomial Darboux pair.
5. Holomorphicity, not genus alone, is load-bearing.  A high-genus fibre
   with a meromorphic time form may still have a different exactness ledger.

## 5. Preservation and loss ledger

```text
source       generic Q-fibre in the smooth Danielewski surface
target       smooth projective normalization of its plane/fibre-product model
map          restriction of the Hamiltonian mate equation
preserved    {P,Q}=kappa becomes dP=kappa eta
obstruction  eta is nonzero holomorphic, or logarithmic with nonzero residues
lost         affine arm labels and polynomial degree data
sidecar      all finite critical-value and infinity valuations of eta
cheap fail   Q=c (one channel) has rational mate b/c^N
cheap core   N=D=2, Q=c+e gives an elliptic time form
```

The theorem is a full-field strengthening of the separated polynomial
hostiles in THM-3591 and a Danielewski analogue of THM-3552.  It does not
exclude a genuinely mixed coordinate whose forced time form acquires
cancellable poles.

## 6. Exact verification contract

The companion checks:

- the Poisson/time-form contraction identities;
- the local critical-point cancellation model;
- all `1<=N<=12`, `2<=D<=12`, `1<=A,B<=12` infinity orders,
  divisibility gates, and the unique residue boundary;
- the elliptic rows `(19)--(20)` and affine rows `(23)--(28)`;
- the A13 order-seven/genus-eight invoice `(21)--(22)`;
- the one-channel, degree-one, `nu=0`, exponent-zero, constant-`Sigma`, and
  characteristic-`p` hostiles.

The proof is universal; finite rows are exact controls, not extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_danielewski_separated_timeform_rational_nonentry_thm3595.py
python3 -O 04-computation/jc2_danielewski_separated_timeform_rational_nonentry_thm3595.py
```
