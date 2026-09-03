# Exceptional quartic: unrestricted normal motion pays, admissible motion does not

**FINITE-EXACT + INDEPENDENTLY AUDITED representative-level slice,
2026-09-03.**  This scratch note does not promote a theorem and does not
identify its map with the full THM-4067 graph-family transgression.  It
records a positive-minor computation, an independent direct-matrix replay,
and the exact typing boundary that remains.

## Outcome

Let

```text
T=K[b,c,e]/(c^2e-b(b+4)),
N=K[x],
gamma_0:T -> N,
gamma_0(H)=H(B(x,Q),C(x,Q),E(x,Q)),
S=gamma_0(T).
```

Vary the Russell compiler at fixed `x` by `q=Q(x)+s` and define the chosen
target-representative normal sidecar

```text
gamma_1(H)=d/ds H(B(x,Q+s),C(x,Q+s),E(x,Q+s))|_(s=0).
```

Put `d=d/dx` and

```text
E_x=N/dS,
U_quad=(dS+gamma_1(T))/dS subset E_x.
```

THM-4381 gives `S^sn=S+K r`.  Differentiation injects its one-dimensional
defect into `E_x`:

```text
S^sn/S -> E_x,                 [r] |-> [r'].
```

Indeed, if `r'` belonged to `dS`, then `r-s` would be constant for some
`s in S`; constants already lie in `S`, contradicting `r notin S`.

The scratch verifier and independent audit prove

```text
boxed: U_quad=E_x.
```

Consequently `[r']` is zero after quotienting by this representative-level
first normal sidecar.  In this precisely stated quadratic compiler slice,
the THM-4381 class is **killed**, not surviving.

## Exact finite certificate

THM-3703's monic Apéry basis gives one monic element of `S` in every
semigroup degree.  Differentiating those elements shows that `E_x` has the
89-dimensional monic complement

```text
{x^(g-1): g is one of the 89 THM-3703 gaps}.
```

The three generator sidecars are

```text
gamma_1(B)=3x^2D(D+2),
gamma_1(C)=2x^3(D+1),
gamma_1(E)=2(D+1),                 D=1+x^2Q.
```

The verifier reduces the sidecars of all actual target monomials
`B^i C^j E^k` of restriction degree at most 240 against the differentiated
monic basis.  At

```text
p=421,                   alpha=126 mod 421,
```

the exceptional quartic vanishes, every encountered denominator is a unit,
and `p` exceeds every differentiated degree.  The rank ladder is

```text
cutoff : columns : rank : rank after adjoining r'
    30 :       3 :    3 :  4
    60 :      12 :   11 : 12
    90 :      28 :   24 : 25
   120 :      54 :   41 : 42
   150 :      93 :   64 : 65
   180 :     146 :   89 : 89
   210 :     217 :   89 : 89
   240 :     307 :   89 : 89
```

Thus 89 columns already present by cutoff 180 give a nonzero minor after a
good specialization.  If the corresponding characteristic-zero minor were
zero, every good reduction would be zero.  The reduction is nonzero, so the
exact rank is at least 89; `dim_K E_x=89` makes it exactly 89.  This is the
sound positive-minor direction.  No modular rank failure is used as a
characteristic-zero obstruction.

The selected-column hash is

```text
e2d88c2341b620cbac9248d682c14509bcd94d115d70a68592d3942d1835d0bd.
```

The first sidecar order already sees the retained missing line.  At
`x=(-1,0,1)`, the generator values are

```text
gamma_1(B): ( 0,0, 0),
gamma_1(C): ( 2,0,-2),
gamma_1(E): (-2,4,-2).
```

For the retained functional `Lambda=(5,-18,13)/18`, their responses are

```text
0, -8/9, -6.
```

Thus `gamma_1(C)` is already locally transverse to the derivative plane of
`S`.  The degree-180 full-rank calculation additionally pays every shifted
gap coordinate of `E_x` and proves equality in that quotient.  These 89
coordinates are not thereby identified with THM-4067 graph periods or node
coordinates.

The independent audit imports only the audited THM-3703 restriction-algebra
certificate, not this scout or the THM-4034 conductor.  It avoids the scout's
gap normal form: in the ordinary coefficient space through degree 172 it
builds `dS` directly from derivatives of the monic Apéry basis and obtains

```text
ambient dimension =173,
rank(dS)           = 84,
rank(dS+gamma_1(T))=173.
```

Thus the quotient rank is again `173-84=89`.  The two implementations select
the same 89 target monomials, with the same hash displayed above.

The audit also proves genuine non-descent rather than relying only on a
warning.  If `gamma_1` factored through `S`, it would induce a derivation
`S -> K[x]`.  Since `Frac(S)=K(x)`, that derivation would extend to
`v(x)d/dx`, forcing

```text
gamma_1(B) C' - gamma_1(C) B' =0.
```

The exact left side is nonzero; its good reduction has degree 42 and
coefficient hash

```text
5ff697cdc557fc1db4c7e8c06420cb7d37091f350122f30aff0099989e08e2e5.
```

Consequently some target representatives with zero restriction have nonzero
first normal sidecar.  This is non-descent to `S`; no claim about a chosen
quotient or graph-family gauge follows from it.

## Why this is not yet the full THM-4067 answer

The computed source, target, and map are exact, but they differ from the
hypotheses of THM-4067 in three load-bearing ways.

1. `gamma_1` is the fixed-`x` derivative for the source-normal deformation
   `q=Q+s`.  THM-4067 writes a fixed-edge-coordinate first-output derivative
   for a supplied moving graph family.  No theorem currently identifies
   these operators modulo exact edge derivatives for all 87 conductor
   fibres.
2. `gamma_1` acts on actual target representatives in `T`; it does not
   descend to a map on their special-fibre restrictions `S`.  THM-3737
   explicitly warns about this.  Kernel representatives can have zero
   restriction and nonzero first normal sidecar.
3. A literal THM-4067 specialization still needs a `K[[s]]` graph family,
   moving endpoint sections for the node and triple fibres, an exact edge
   primitive module, and a proof that its special-fibre image `U` is exactly
   `U_quad`.

Therefore the strongest present statement is:

```text
PROVED in the representative-level quadratic normal-sidecar slice:
    [r'] is killed.

OPEN as an identification with the full THM-4067 transgression:
    U_quad =? U_graph after coordinate change and specialization.
```

The missing coordinate is no longer a rank or conductor calculation.  It is
the first-order deformation/endpoint comparison map between the
`q=Q+s` compiler and THM-4067's graph complex.

## Hostile comparison: admissible kernel motion does not pay the class

The unrestricted source `T` is too large for the stagewise lifting question.
THM-4039 already identifies the solution-preserving kernel of `L_0` as

```text
{(P,Q)=(C'U,E'U): U in K[x], C'U in S, E'U in S}.
```

An exact bridge scout now extracts these `18` channels rather than merely
counting them.  Write `c=Lr`, where `L=x(x^2-1)`.  At an ordinary node,
distinct full target tangents have nonproportional `(C',E')` projections
unless the projection of the Russell surface

```text
c^2e=b(b+4) -> (c,e)
```

ramifies; its ramification divisor is `B=-2`.  At the inherited good fibre
`(p,alpha)=(137,44)`, both polynomial degrees survive and

```text
gcd(h_172,B+2)=1 mod 137.
```

This is the sound positive-resultant direction, so no one of the `86`
geometric nodes lies over `B=-2`.  The two kernel equations consequently
force `U` to vanish on every conductor branch, hence `U=rH`.  Node-local
membership is then automatic; at the retained triple it is exactly

```text
Lambda(r'H)=0,                 Lambda(E'r'H)=0.          (*)
```

The two independent equations in `(*)` act only on `H mod L`.  Their common
line has the especially small monic representative

```text
H_0=x^2-(4/27)x-4/3,             U_0=rH_0.
```

Therefore the potential module and one explicit free channel basis are

```text
M={U:C'U,E'U in S}=K U_0+cK[x],

M=K[Z]U_0 direct_sum direct_sum_(j=0)^16 K[Z](cx^j),
v_U=(C'U,E'U).
```

The leading degrees are `177,178,...,194`, one in every residue class modulo
`18`.  Freeness follows from those distinct residues.  For generation, use
`Z=LW` with `W` monic of degree `15`: then `ZU_0=c(WH_0)`, and
`1,x,...,x^16,WH_0` is a triangular `K[Z]`-basis of `K[x]`.  The exact
restriction-pair hash is
`72d394a366977b69e3d02cdce4380b14e8d289fab93da2bb9d9af0442a65dec0`.
This extraction is finite-exact bridge data, not a new theorem claim.

This smaller source has the opposite retained answer.  Since `P,Q` have
common retained values while

```text
C'_ret=(3,3,3),             E'_ret=(-9,4,9),
```

one first obtains `U_ret=0`.  Put `v=U'_ret`.  The derivative rows of both
`P` and `Q` lie in the retained tangent plane, so

```text
Lambda(v)=0,                  Lambda(E'_ret v)=0.
```

These two planes meet in exactly the constant line; hence
`v=lambda(1,1,1)`.  For any actual target representatives of `P,Q`, their
retained normal values are uniquely determined by these derivative rows:
`B'=gamma_1(B)=0` at the triple and the rows of `C',E'` are independent.
Therefore

```text
gamma_1(P)_ret=lambda gamma_1(C)_ret,
gamma_1(Q)_ret=lambda gamma_1(E)_ret,
```

independently of representative.  The fixed-`x` first-output term is then

```text
tau_ret=C' gamma_1(Q)-E' gamma_1(P)
       =12 lambda (1,1,1),
Lambda(tau)=0.
```

By contrast `Lambda(r')` is nonzero, and `Lambda` annihilates `dS` at the
retained triple.  Thus the admissible kernel-transgression span cannot hit
`[r']`.  This is precisely THM-4039 equations (15), (24)--(29), re-read in
the seminormal language; it needs no new canonical theorem.

The sharp lesson is therefore a domain boundary:

```text
unrestricted target-representative motion:  image in E_x has rank 89;
solution-preserving L_0-kernel motion:        retained [r'] gate survives.
```

Neither line is yet THM-4067's graph-family `U`.

## Relation to THM-4046

The sidecar parameter `s` is the even stable parameter in the quadratic fold
`q=Q+t^2`.  This makes the computation directly relevant to the THM-3737/
4046 compiler.  It does **not** remove THM-4046's `J_8` obstruction.  Paying
the seminormal derivative class by unrestricted target-representative motion
says nothing about the source of the order-eight scalar `kappa`, because
those directions need not preserve the preceding stable equations.  The
retained order-eight identity contains the additional closed-two-form and
stagewise compatibility data, and the admissible kernel comparison above has
zero `Lambda` response.

The unexpectedly strong signal is that the raw first normal sidecar spans
all 89 dimensions of `N/dS`, not just the single seminormal line.  The
THM-4039 hostile comparison explains why this does not tune the lift:
unrestricted representative motion includes directions that are not
solution-preserving.  A future graph-family comparison must first construct
a descended base derivative and then identify its image in the exact
THM-4067 quotient; the raw 89-column minor cannot substitute for that step.

## Reproduction

```bash
python3 -B 04-computation/jc2_exceptional_quartic_normal_sidecar_transgression_scout_transgression_niche_s615.py
python3 -B -O 04-computation/jc2_exceptional_quartic_normal_sidecar_transgression_scout_transgression_niche_s615.py
python3 -B 04-computation/jc2_exceptional_quartic_normal_sidecar_transgression_independent_audit_s615.py
python3 -B -O 04-computation/jc2_exceptional_quartic_normal_sidecar_transgression_independent_audit_s615.py
PYTHONHASHSEED=17 python3 -B 04-computation/jc2_exceptional_quartic_normal_sidecar_transgression_independent_audit_s615.py
python3 -B 04-computation/jc2_exceptional_quartic_kernel_first_output_obstruction_thm4388_drafter_s615.py
python3 -B -O 04-computation/jc2_exceptional_quartic_kernel_first_output_obstruction_thm4388_drafter_s615.py
```

Normal, optimized, and fixed-hash-seed independent runs byte-match the
stored audit transcript.  Its raw LF hashes are

```text
primary scout:      8055613271c4100c2cb3e8b4e76528b3eb8963fc5f9e5ebb9d7b94ad62cb0f3d
independent script: 8668b446bc937492ca2c88c71af48bcd0f832af28424fa15ba2461e55dd56a07
independent output: cfbd0d956d0bd46e14786ee6c753c3e28847f4b48f8dcd777bcd9215a88718f2
kernel script:      2bb70db0ab6dc9d119bf97993f0d686069eeb1bc23858e94978194cef3e56376
kernel output:      4bb7bacca1b42f6ce04496dadeca8a3eb45cafe264fb0a4fa7b84cfdda58d05a
```

The computation is deliberately stored as a session scout rather than a
canonical theorem companion.  It proves no graph-family identification,
coherent all-order lift, convergence, algebraization, Keller map, or
consequence for `JC(2)` or `DC(2)`, which remain open.
