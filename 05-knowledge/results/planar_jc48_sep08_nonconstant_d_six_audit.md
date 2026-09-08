# Independent audit: nonconstant-D finite-six / infinity-double

**Status: INDEPENDENT ANALYTIC + SOURCE + EXACT-REPLAY AUDIT PASS.**
No mathematical correction is required. This audit accepts the complete
nonconstant-`D` polynomial-mate exclusion in the frozen
[primary proof](planar_jc48_sep08_nonconstant_d_six.md). The stronger
rational exclusion is accepted only after its polynomial gate has imposed
`c=nu=0`. Constant-`D` rational mates are retained.

The audit read the entire primary and source, independently reconstructed
the original-coordinate field Jacobians and the quadratic trace, checked
the generic residue implication in an ordinary local coordinate, and
replayed the source normally and with `-O`. Both replays match all 399
bytes of the frozen 64-gate output. A separate inline computation without
producer imports passed ten original-coordinate, conjugation, residue,
asymptotic and hostile identities. The source/output were not edited.

## 1. Accepted universe and original global rows

The entry is exactly `F=H^2+L` with `H in L2`, `L in L1` on the fixed
DG surface, leading polynomial `N=u^6`, `u=x-p`, and nonconstant `F|D`.
It does not reduce arbitrary Keller pairs to this entry. The finite
position `p` is arbitrary; the change to `u` is not used as a surface
automorphism.

I independently checked completeness of the pre-active row equations.
The three coefficient rows of the quadratic section in the original
`x=u+p` are

    Q,  P-2Qx^2,  N-Px^2+Qx^4,

all of degree at most four. If `deg Q>=3`, the second condition makes
the highest part of `P` equal to `2Qx^2`; the last then retains the
uncancellable term `-Qx^4` of degree at least seven, since `deg N=6`.
Thus `deg Q<=2`, `deg P<=4`. Cancelling the remaining degrees six and
five gives exactly

    Q=(P4-1)u^2+(P3-2pP4+4p)u+Q0.

The linear section rows `R,M-Rx^2` similarly give the displayed `R`
coefficients. These conditions are also sufficient: after the two high
cancellations, every original row has its required degree. There is no
hidden condition on `P0,P1` or `M0,M1` before the active gate.

The full active rows in the primary are therefore complete. Direct
substitution into the original second chart gives slopes `-A` and
`-lambda` for `H|D` and `L|D`. A nonzero `A` leaves a quadratic term in
`F|D` that no affine `L|D` can cancel. If `A=0`, its remaining slope is
`-lambda`. This verifies the exact equivalence

    F|D nonconstant iff (A,lambda)!=(0,0).

## 2. All local entry branches and the polynomial-ring gate

The original infinity has multiplicity two times a unit,
`r^2(1-pr)^6`. The primary retains the actual factor `r^2` in the
volume. A normal unit or M-unit is covered by the proved local suppliers.
In the other case the full generic tangent quartic has nonzero constant
term and four simple nonzero roots. The repeated-root eliminant is
independent of the generic fibre coefficient and has constant
`-4a^2`, so it has only finitely many exceptional levels. Every actual
branch has unweighted form order `-1`, improved to `+1` at infinity.
This exhausts the possible infinity branches; it does not import a
finite logarithmic obstruction without its weight.

At the finite sixfold point, the M-unit branch has total primitive
capacity at most two, including its balanced Morse case. All other
boundary points are regular. A generic transverse point of nonconstant
`F|D` makes a primitive have local degree three, exceeding that capacity
on its own compact component. The argument does not require the generic
fibre to be connected. If `M` vanishes but the normal derivative is a
unit, every relative form is regular instead. The remaining shared
first-jet gate then forces exactly the active values and derivatives
stated in the primary. These are proved dependencies, and the retained
constant-`D` alternatives are not silently used here.

The polynomial-weight argument in Section 3 is valid for arbitrary
finite polynomial support. Under `w=u^2t`, an original monomial
`u^i t^j` becomes `u^(i-2j)w^j`. If the least exponent `ell` of a
nonzero polynomial mate is negative, every monomial in its leading
coefficient has `j>=1`; hence that coefficient is divisible by `w`.
If `F=f(w)+O(u)` with `f'!=0`, the unique leading bracket term is

    -ell f'(w) g_ell(w) u^(ell+1).

Terms involving `F_u` start one order later. If `ell<-1` this is an
uncancellable negative power; if `ell=-1` its coefficient cannot be one
because it is divisible by `w`. For `ell>=0` the bracket has positive
`u` order; the `ell=0` leading derivative vanishes, leaving order at
least two. This proves the all-degree gate `c=nu=0`.

There is no use of the separately pending general weight-gate bundle:
the complete needed proof is written locally. The stronger second gate
is correctly marked as unused. Polynomiality is essential; the rational
hostile at the end of the primary has nonconstant `f(w)` and is not
subject to this argument.

## 3. Independent actual field and quadratic-trace reconstruction

I recomputed in the original `(u,t)` coordinates, without substituting
the producer's later field formulas. Put `v=u+u^3t`. After the first gate,

    H=v^2+bv+d+A u(v-2p),
    L=lambda u(v-2p)+mu v+e,
    J(u,v)=u^3,
    J(H,v)=A(v-2p)u^3.

Thus `omega=u^-3 du wedge dv`. For `A!=0`, the genuine field inverse
`u=(h-h0)/[A(v-2p)]`, `h0=v^2+bv+d`, gives

    omega=A^2(v-2p)^2/(h-h0)^3 dh wedge dv.

These denominators are allowed in a rational field argument; no point
on them is assumed regular. The generic quadratic in `h` is precisely

    h^2+delta h-delta h0+mu v+e-f=0,
    delta=lambda/A.

Its discriminant has the irreducible simple factor linear in the
transcendental `f`, with slope four, and is not a square in `C(f)(v)`.
This remains valid when `lambda=mu=0`, where the extension can be a
constant quadratic extension. The proof correctly avoids asserting
geometric generic integrality in this step.

Here is an independent trace derivation using conjugate roots rather
than the source's companion matrix. Set `z=h-h0`,
`B=2h0+delta`, `C=h0^2+mu v+e-f`. If the roots are `z1,z2`, then
`z1+z2=-B`, `z1z2=C`, and

    sum_i 1/[z_i^3(2z_i+B)]
      = (1/z1^3-1/z2^3)/(z1-z2)
      = -(z1^2+z1z2+z2^2)/(z1^3z2^3)
      = 1/C^2-B^2/C^3.

This is the **full** trace; the normalized trace is half of it. I also
verified it by putting the other root equal to `-B-z`, clearing
 denominators, and taking the polynomial remainder modulo
`z^2+Bz+C`; the remainder is identically zero.

With the primary's convention `dF wedge eta=omega`, the sign and
carrier are consequently

    Tr eta=-A^2(v-2p)^2[B^2/(P-f)^3-1/(P-f)^2]dv,
    P=h0^2+mu v+e.

A rational mate gives an exact trace relative to `C(f)`. Only necessity
is used, so neither a trace-zero converse nor a normal-versus-full trace
ambiguity enters the proof.

## 4. Independent residue calculation and complete parameter split

For generic `f`, the roots of `P(v)-f` are simple and avoid `P'` and
`v=2p`. Their points are generic in the `v` line, since `P` is a
nonconstant quartic. Define `D_P=(1/P')d/dv` and

    L1=(v-2p)^2(2h0+delta)^2/P',
    L2=(v-2p)^2/P'.

The residue in the `P(v)` parameter is
`D_P^2 L1/2-D_P L2`. I independently recovered this in an ordinary local
coordinate. If the Taylor data at the root are

    P-f=p1 x+(p2/2)x^2+(p3/6)x^3+...,
    n1=n10+n11 x+(n12/2)x^2+...,
    n2=n20+n21 x+...,

then the residue of `n1 dv/(P-f)^3-n2 dv/(P-f)^2` is

    [n12/2-3a n11+(6a^2-3b0)n10]/p1^3
       -(n21-2a n20)/p1^2,
    a=p2/(2p1), b0=p3/(6p1).

Literal differentiation of `D_P^2(n1/P')/2-D_P(n2/P')` gives the same
expression. This checks the factor two and both pole signs through a
path independent of the source's Taylor control.

The resulting rational function is independent of `f`. Its vanishing at
a generic inverse root forces an identity in `C(v)`. Integrating there
is legitimate and gives

    L1'-2(v-2p)^2=kappa P'.

I checked its asymptotics directly by clearing denominators and comparing
polynomial degrees. Its left side has degree difference two and leading
coefficient one; the right side is either cubic, if `kappa!=0`, or zero.
Thus the `A!=0` branch is impossible for all `p,b,d,mu,lambda`, including
`lambda=0` and all lower-coefficient degenerations.

If `A=0`, nonconstant `F|D` gives `lambda!=0`. Direct differentiation in
original coordinates gives `J(F,v)=lambda(v-2p)u^3`, paying the actual
rational field and

    eta=lambda^2(v-2p)^2 dv/(f-P)^3.

Its residue condition is `D_P^2 L2=0`, hence `L2=alpha P+beta`.
The nonzero rational function `L2` has degree difference `-1` and tends
to zero at infinity. An affine function of a quartic cannot do this
unless it is zero. This closes the complete complementary field case.
No value of `b,d,mu,p` or branch with constant quadratic discriminant was
lost. The case `(A,lambda)=(0,0)` is exactly the excluded constant-`D`
entry and remains covered by its separate proved theorems.

## 5. Scope, hostile and dependency acceptance

The independent original-coordinate calculation also reproduced

    H=u^2(1+u^2t)^2+(1+u^2t),
    L=u(1+u^2t),
    G=1/[2u(1+u^2t)]-H,
    J(H^2+L,G)=1.

This is the actual constant-`D`, nonconstant-`L` rational hostile. It
confirms that neither the first weight gate nor the theorem's polynomial
conclusion can be applied to arbitrary rational coefficients. Its
source poles are not erased by a field change.

I checked the current status of both named constant-`D` suppliers:
[constant-D shared six](planar_jc48_sep08_constant_d_shared_six.md) and
[translated constant-D M-unit six](planar_jc48_sep08_const_d_translation.md)
are `PROVED` and independently audited on disk. The finite-double
[leading residue](planar_jc48_sep08_leading_exactness.md) is also proved.
Accordingly the stated all-location binary `(6,2)` polynomial consequence
has the required scoped inputs once the audited primary is promoted.
This audit does not claim a reduction of arbitrary planar Keller maps
or a blanket rational exclusion in the constant-`D` boundary.

## 6. Frozen source, replays and pins

Reproduction commands are

```sh
python3 -B 04-computation/planar_jc48_sep08_nonconstant_d_six.py
python3 -B -O 04-computation/planar_jc48_sep08_nonconstant_d_six.py
```

Both independent executions pass **64 always-active gates** and agree
byte-for-byte with
[the frozen output](planar_jc48_sep08_nonconstant_d_six.out).
The source was read in full; it contains no producer imports, sampled
parameter claims, removable assertions, or dependence on the unused
second polynomial weight gate. Its finite negative-row monomials are
correctly controls for the separate unbounded proof. The trace matrix
is applied only after the actual field coefficients have been checked.

| Artifact audited | Bytes | SHA-256 |
| --- | ---: | --- |
| Primary proof before status promotion | 14,392 | `6ebc3102a69b46bf4820b17e4525768c38841ebab815863a1cc613ff671a10e4` |
| Source | 7,363 | `946cae045bfd1fd3b4ce8c054cf444196c2a2e99d00d0265ef7e070684bdaa3f` |
| Frozen output and both audit replays | 399 | `16d2315be218c959d98169b153b6f5c2395f3fb12cbf4f8804c45496278693c6` |

The semantic record is
`5821fc3ad3231a3db28a6c9cc3652432ba765959a9c219b76354a2cc578082a8`.
The additional ten independent symbolic controls used original-coordinate
Jacobians, conjugate-root polynomial reduction, ordinary local Taylor
residues and degree comparison; they did not import or execute producer
functions. The primary, source and output remain unchanged by this audit.
