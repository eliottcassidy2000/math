# Complete doubled-return signs: a proved quadratic family and the carry gauge

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED for
the family in Section 2; FINITE-EXACT for the fourteen named controls.**
The requested starting support `(-4,1,6)` is a positive control, not a
counterexample: its complete second return is strictly negative at both
zeros of its first return after the specified monomial normalization.
This sign statement holds throughout the primitive family
`(-4,g-4,2g-4)`, odd `g>=5`. General higher-channel constant-sign separation
remains **OPEN**. This does not newly close width-four coprimality.
The generated [lower-carry endpoint-fifteen family](trinomial_width15_empty_core_returns_sep06.md)
now extends two-rung detection beyond the prior endpoint-eight strip; its
independent trace/norm certificate is a separate unbounded result.

## 1. Inheritance and the exact starting test

The closest incoming mechanism is
[THM-4436, complete factorial-row simple negative roots and trinomial phase wall](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
It proves individual complete rows have simple negative scalar roots, and
explicitly leaves relations between different masses open. The
[width-eight report](overnight_20260906_moments_width8.md) already proves
two-rung coprimality for every support considered by the family below.
Its symbolic resultants and shifted positive coefficients are inherited
operations, not new discoveries here.

The corrected near miss is discarding channels or their index shifts.
The canonical hostile `(-13,1,8)` retains both doubling carries. The
[two-channel note](trinomial_two_channel_empty_core_returns_sep06.md),
canonically routed to
[THM-4432](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md),
also records an adjacent-pair determinant equal to zero for `(-4,3,10)`
and positive for `(-4,5,14)`. Those refute pair compression, not the
complete-row root-sign statement: both complete quadratic rows pass here.
The least-used sidecar is the position of the Euclidean remainder's zero
relative to the full first-row roots.

Targeted searches on `origin/main` for signed remainder, constant-sign
remainder, root sign, the named supports, and the return/tunability routes
recovered these warnings and mechanisms. No prior uniform root-sign
theorem was found. No external-priority or literature-novelty claim is made.

For `f(u)=alpha*u^-4+beta*u+gamma*u^6`, put
`tau=alpha*gamma/beta^2` and `X=alpha*beta^4`. Direct complete-channel
enumeration gives

```text
CT(f^5)=X*P(tau),       P=5+30tau+10tau²,
CT(f^10)=X²*Q(tau),
Q=45+840tau+3150tau²+2520tau³+210tau⁴.
```

The first channels are `(1,4,0),(2,2,1),(3,0,2)`; the second are
`(2+j,8-2j,j)`, `j=0,...,4`. Thus no generated or carry term is omitted.
Exact division yields

```text
Q=(21tau²+189tau-525/2)*P + 7770tau+2715/2.
rho=-181/1036,
P(rho)=34305/536648>0.
```

The remainder zero `rho` lies right of the first quadratic's vertex
`-3/2`, so it lies strictly right of both roots. Its slope is positive;
therefore `Q` is negative at both roots of `P`. The proposed cheap
refutation would require `P(rho)<0`; it has the opposite sign here.
The sign is for the displayed real scalar quotient. The actual complex
moment includes the nonzero factor `X²`.

The active board is: individual negative-root geometry; cross-mass
remainders; canonical carry gauges; resultant signs; actual Laurent
realizability. The starting positive signal supports the second lane;
the examples below separate it from the first and fifth.

## 2. An unbounded uniform negative-sign certificate

Let `g>=5` be odd, and take nonzero complex coefficients on

```text
f(u)=alpha*u^-4+beta*u^(g-4)+gamma*u^(2g-4).
```

The support is primitive exactly because `g` is odd:
`gcd(4,g-4,2g-4)=gcd(4,g)=1`. Its complete first return has mass `g`
and three channels; all support returns have masses divisible by `g`.
Write

```text
tau=alpha*gamma/beta²,       X=alpha^(g-4)*beta⁴,
P_g(tau)=12tau²+12(g-2)tau+(g-2)(g-3),
Q_g(tau)=sum_(j=0)^4 (2g)_(8-j) tau^j/[(8-2j)! j!].
```

Here `(v)_k=v(v-1)...(v-k+1)`. The full identities are

```text
CT(f^g)=X*g(g-1)*P_g(tau)/24,
CT(f^(2g))=X²*Q_g(tau).
```

Their channels are `(g-4+j,4-2j,j)`, `0<=j<=2`, and
`(2g-8+j,8-2j,j)`, `0<=j<=4`. Both carries vanish.

**Claim.** At each of the two scalar zeros of `P_g`, `Q_g<0`.
The inequality is strict; there is no equality case. These scalar zeros
are attainable on the coefficient torus, for example by taking
`alpha=beta=1` and `gamma=tau`. Consequently tuning either first zero
gives first nonzero moment exactly `2g`, as already implied by the older
coprimality theorem. The added predicate is its uniform normalized sign.

**Proof.** Division over `Q(g)` gives `Q_g=C_g+D_g*tau mod P_g`, where

```text
C_g=g(g-3)(g-2)(g-1)(2g-3)(2g-1)
       *(503g²-1723g+1470)/30240,
D_g=g(g-2)(g-1)(2g-3)(2g-1)(3g-5)(11g-18)/180.
```

In particular `D_g>0` for `g>=5`. The remainder zero is

```text
rho_g=-(g-3)(503g²-1723g+1470)/[168(3g-5)(11g-18)].
```

The discriminant of `P_g` is `48(g-2)(2g-3)>0`. Its two roots are
negative, since their sum and product have the required signs. Its
vertex is `v_g=-(g-2)/2`. Define

```text
H(g)=5141g⁴-31762g³+69537g²-62820g+18900,
J(g)=2269g³-11468g²+19233g-10710.
```

Exact substitution gives

```text
P_g(rho_g)=(g-3)(5g-7)H(g)/[2352(3g-5)²(11g-18)²],
rho_g-v_g=J(g)/[168(3g-5)(11g-18)].
```

For `s>=0`,

```text
H(s+5)=5141s⁴+71058s³+364257s²+820900s+686100,
J(s+5)=2269s³+22567s²+74728s+82380.
```

Both are strictly positive. An upward quadratic is positive to the
right of its vertex exactly beyond its larger root. Hence `rho_g`
lies right of both first roots. Since `D_g>0`, the remainder, and thus
`Q_g`, is strictly negative at each. This proves the claim.

The displayed polynomial argument actually holds for every **real**
`g>=5` in the falling-factorial extension. For even integral `g`, the
displayed mass-`g` and mass-`2g` identities are still valid, but the
support is nonprimitive and `g` need not be its first return. The
first-return conclusion here is only for odd integral `g`.

## 3. A carry changes canonical sign without changing noncancellation

An independent audit by the certificate lane supplied the actual
three-channel control `(-15,1,9)`. Its first mass is eight,
`(A,B)=(2,3)`, and its canonical third residue is `z0=1`. The rows are

```text
P=56+560tau+56tau²,
Qcanonical=16+10920tau+400400tau²+1681680tau³
             +720720tau⁴+8008tau⁵,
Qcanonical mod P=-47087024-466126752tau.
```

The remainder zero is `rho=-2942939/29132922`. With
`p=P/56=1+10tau+tau²`, exact arithmetic gives

```text
p(rho)=23910838225/848727144258084>0,   rho>-5.
```

It again lies to the right of both first roots, but this remainder has
negative slope. Thus `Qcanonical` is **positive** at both roots.
This refutes canonical negativity across supports. It does not refute
one nonzero sign across the roots within a support.

For the exact gauge, let the first canonical anchor be `v`, the affine
step `d=(B-A,-B,A)`, and `tau=coefficients^d`. The lower carry is
`epsilon_z=floor(2z0/A)`. The second canonical anchor is

```text
v2=2v-epsilon_z*d,
CT(f^(2g))/coefficients^(2v)=tau^(-epsilon_z)*Qcanonical(tau).
```

Here `epsilon_z=1`; multiplying the positive canonical value by
`tau^-1<0` restores a negative first-anchor-normalized value. The same
sign flip occurs in the old two-carry control `(-13,1,8)`.
Any general negative-sign candidate must retain this factor. The upper
carry still changes the polynomial's degree and coefficients, even though
it does not change the left anchor.

## 4. The failed implication and strongest survivor

Individual simple negative roots do not imply constant root sign, even
when both polynomials have positive coefficients, are coprime, and their
degrees are in ratio two. A logical hostile is

```text
P=(tau+1)(tau+3),
Q=(2tau+1)(tau+2)(tau+4)(tau+5).
```

All roots are distinct and negative, and the polynomials are coprime.
Nevertheless

```text
Q mod P=-23-11tau,
Q(-3)=10,      Q(-1)=-12,
P(-23/11)=-120/121<0.
```

This pair is not claimed to come from one Laurent support. Its role is
to refute the first unsupported implication from THM-4436's individual
root geometry to a cross-mass sign. It also demonstrates that constant
nonzero root sign is stronger than coprimality.

The strongest current survivor is: all fourteen named actual controls
below have a negative second value after division by the square of the
first anchor, and Section 2 proves that predicate for one unbounded
family. The general candidate

```text
P(tau)=0  =>  tau^(-epsilon_z)*Qcanonical(tau)<0
```

is **OPEN**, as is its weaker general constant-sign version. No actual
support with opposite second-row signs at first-row roots was found in
this bounded probe; that is not a proof of either assertion.

For a quadratic first row `P=p0+p1*tau+p2*tau²`, `p2>0`, and linear
remainder `R=A*tau+B` with `A!=0`, the exact next test is

```text
N=p2*B²-p1*A*B+p0*A²=A²*P(-B/A).
```

With distinct real first roots, `N>0` is equivalent to a common nonzero
remainder sign there; `N<0` gives opposite signs; `N=0` gives a common
root. Coprimality only requires `N!=0`. A nonzero constant remainder is
an immediate positive control; the zero remainder fails noncancellation.
The concrete next symbolic family is the carry control extended to
`(-15,2g-15,3g-15)`, `g>=8`, `gcd(g,15)=1`. Its first row is quadratic
and its lower carry is one. Proving a sign for its `N` is a bounded
symbolic obligation with unbounded parameter, not a larger height census.

The map from full channels to `P,Q` preserves every cancellation after
the nonzero monomial factors are restored. Passing to real-rootedness
forgets relative root position. The remainder retains precisely that
sidecar; passing further to its nonzero resultant forgets the stronger
same-sign predicate. The original empty-hexagon paper motivates this
certificate-and-sidecar research move at the method level; no orientation
or geometric theorem is being transferred to Laurent moments.

## 5. Exact reproduction, controls and scope

The [source](../../04-computation/trinomial_root_sign_empty_core_returns_sep06.py)
and [frozen output](trinomial_root_sign_empty_core_returns_sep06.out) use
integer/rational arithmetic and SymPy. No repository producer or theorem
is imported. Thirteen named primitive supports have at least three first
channels:

```text
(-4,1,6), (-4,3,10), (-4,5,14), (-5,1,7),
(-6,1,8), (-6,5,16), (-7,1,9), (-7,2,11),
(-8,1,10), (-8,3,14), (-9,1,11), (-10,1,12), (-15,1,9).
```

The fourteenth is the old two-channel/two-carry control `(-13,1,8)`.
Original charge-equation enumeration reconstructs both complete rows;
independent repeated Laurent multiplication checks every coefficient and
the emptiness of all earlier moments. Exact Sturm root isolation and
rational interval Horner evaluation certify each remainder's sign at
every first root. A midpoint sign alone is never used. The closest-root
control `(-10,1,12)` requires refinement to `10^-16`.

The symbolic check reconstructs the entire family remainder over `Q(g)`,
both root-location identities, and all positive shift coefficients.
The logical hostile and the first-anchor carry transformation are also
checked. All **330 explicit gates** pass normally and with optimization;
the two outputs are byte-identical. These controls are a named list, not
a support census and not a claimed minimal counterexample search.

```sh
python3 -B 04-computation/trinomial_root_sign_empty_core_returns_sep06.py
python3 -B -O 04-computation/trinomial_root_sign_empty_core_returns_sep06.py
```

Semantic digest:
`13a831940fc814386bb016642ea73640b544cef983871042d6dadea4fbc04b9e`.
Raw SHA-256 values:

```text
source c0af333a03d57719bcccd4872993f508c2cdd9fcead7fe39392956a915a70258
output 1aeef00bf711d86fcfdec5610e3f482b4dd79f423c75304ce3d3b89657577660
```

The [companion certificate audit](trinomial_trace_sign_empty_core_certificate_sep06.md)
independently reconstructed the symbolic
remainder by the recurrence
`tau²=-(g-2)tau-(g-2)(g-3)/12`, rather than by Euclidean division. It
recovered both root-location identities and the positive shifts, and
audited the odd-integer first-return restriction: **PASS**. Its equivalent
trace/norm check gives negative sum and positive product of the two
remainder values, hence the same strict sign conclusion.
Its final-text audit also accepted the carry control, the distinction
between constant sign and coprimality, and the next-family quantifiers.

The work strengthens a sign predicate on an already closed smaller-endpoint
four family. It neither proves unrestricted two-rung coprimality nor
improves the existing width-eight return theorem, and it supplies no
general Gaussian-moment or LRC conclusion.
