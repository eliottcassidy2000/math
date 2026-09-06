# Independent audit of sharp finite depth recognition

**Status: ANALYTIC AUDIT PASS + FINITE-EXACT INDEPENDENT RECONSTRUCTION.**
The audited claims are the complete capped depth description, the full
truncated annihilator, and the sharp recognition cutoff in
[the depth report](planar_jc_long_20260906_depth.md), including its
all-height extension. No mathematical correction was required. A suggested
rewrite of the upper-bound inequality distinguishes its continuous linear
bound from a stepped expression involving ceilings.

The proof inputs were read directly: THM-4308's actual depth generators,
[THM-4369, Section 4](../../01-canon/theorems/THM-4369-source-packet-pascal-circuit-kernel-and-boundary-basis.md)
and its integral principal image, and
[THM-3973, Section 5](../../01-canon/theorems/THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion.md)
and its entire graded pieces. The principal image and all-height graded
formula are inherited; the finite support recognition theorem is the new
consumer being audited.

## Analytic audit

For height two, a literal generator `x^a u^b p^c y^e` has intercept
`ell=2c+3e-a`, initial row `b+c+2e`, and final row `b+2c+3e`. Under
`a+b<=d`, the final row is at most `ell+d`. Its signed diagonal is a
monomial times `(1-z)^(c+e)`, with `c+e>=ceil(ell/3)`. This checks both the
factor and the finite degree cap. The complete Bernstein base at depth
zero, followed by multiplication by `u`, really does give all coefficients
over the integers. The exceptional intercepts `ell=1` and `ell<=0` are
handled correctly, including the zero module at `ell=1,d=0`.

The projection rank follows from the initial unit triangular block of
multiplication by `(1-z)^rho`. The consecutive simplex bank has exactly
that codimension, annihilates the actual source generators, and has a unit
minor. Thus it is the full annihilator, including its integral assertion;
this is stronger than comparing rational dimensions alone.

For recognition, the crucial input is the fixed polynomial's *zero later
coefficients*. Dividing its diagonal polynomial by `(1-z)^rho` gives the
reciprocal-series criterion stated in the report. When its degree is below
`rho-1`, the proposed square coefficient matrix is exactly the displayed
binomial Toeplitz matrix. The determinant product is nonzero for
`k<=rho`, `B>=0`; the Desnanot–Jacobi recurrence has the correct adjacent
minors, and its `B=0` boundary is triangular with determinant one. The
remaining full-degree case simply tests the complete diagonal. This proves
recognition at the claimed diagonal row, rather than merely a rank heuristic.

For height two, the resulting diagonal row is at most
`floor(4T/3)+d+1`. For general height `h>=2`, the same argument applies to
the *specified capped filtration* of `B_h`, using its inherited graded
pieces. The upper bound can be written transparently as

```
n_ell <= T+d+1+ell*(1-1/h-1/(h+1))
       <= h^2*T/(h+1)+d+1.
```

Taking the integer part gives the claimed
`N_h(T,d)=floor(h^2*T/(h+1))+d+1`. This argument does not require monotonicity
of the individual ceiling expression. Here `h` indexes completion surfaces,
not the number of coordinates in the Jacobian conjecture.

The sharp witness `t^T` is valid. Its diagonal polynomial is the constant
one, which has no required root at `z=1`. Truncating the reciprocal of
`(1-z)^rho` produces an *actual source image*, and its first mismatch is at
exactly `N_h(T,d)`, with a nonzero integer coefficient in characteristic zero.
This verifies the necessity of the extra support information and the exact
endpoint, including `T=1`, `d=0`, and the `T=0` exception.

The characteristic-three hostile is also an actual source construction:

```
p^3-3y^2=t^3-3x^4t^5-2x^6t^6.
```

Modulo three this matches `t^3` through row five, although `t^3` has no
required double root on its diagonal. It first differs at row six. Thus the
characteristic-zero hypothesis of the sharpened recognition cutoff cannot
simply be removed. The integral diagonal image itself is unaffected.

## Independent exact path

The [audit source](../../04-computation/planar_jc_long_20260906_memory_depth_audit.py)
imports no repository mathematics or producer helpers. Its
[output](planar_jc_long_20260906_memory_depth_audit.out) reports **828,684
always-active gates**. It reconstructs literal monomial spans and compares
them with Hasse-derivative kernels at `u=-1`. For fixed support it computes
the dimension of the intersection of the actual projected source span with
the subspace having those zero later rows. Equality with the full Hasse
kernel proves agreement of the two tested membership subspaces. This route
uses neither the producer's reciprocal-series matrices nor its determinant
formula to decide recognition.

The complete declared controls are:

- 530 literal source-image universes: `h=2..6`, `d=0..3`, `ell=-d..24`,
  omitting negative ambient dimensions.
- 1,980 fixed-support diagonal universes: `h=2..6`, `T=1..6`, `d=0..3`,
  with every nonempty intercept from `-d` through `h*T`.
- 120 sharp actual-source controls, checking success through `N_h-1`,
  failure at `N_h`, and nonmembership of `t^T` in the full source image.
- 504 complete simplex-bank universes at height two, with every source
  generator checked against every bank row.
- 140 exact determinant controls using a separate symbolic determinant
  algorithm, plus the explicit characteristic-three source polynomial.

These controls independently support the universal proof; they are not a
finite substitute for its quantifiers. The audit also accepts the distinction
between a projected depth image and its restriction to a bracket-compatible
family, and between finite recognition at a supplied support bound and any
unproved bound on the support of a hypothetical Keller pair.

Reproduce with

```
python3 -B 04-computation/planar_jc_long_20260906_memory_depth_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_depth_audit.py
```

The source and output hashes below freeze this independent path. The full
mathematical scope is the depth report's stated characteristic-zero
filtrations. This audit makes no JC(2), entry, degree-bound, or polynomial
termination claim.

Normal and optimized audit streams byte-match. SHA256 pins are

```
audited producer 15a885b1c0aa5bdd8b956414ac1acf79bc9bb392c53d483733971a2b8e688984
audit source c9a479683792b1278a87ba8f68a54b34171ec1555f876cef429d089ffe4ffb7d
audit output d582116fe1db0259c0121cbc4606984216fa6013396156625aee8d7309b13e29
semantic 640659ee96f28dc4cfcc754787e060dfdfb710f30b54c6fddda465b2ce5bdd7c
```
