# An open coefficient rectangle has negative response at every original phase

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**

For the fixed two-anchor coefficient pencil, the full carried response is
strictly negative at every positive original first-row root whenever

    83<=x=e3<=85, 34<=y=e4<=36, 0<=z=e5<=5.             (1)

This is a coefficient-prism theorem, even without beta-root or interlacer
assumptions. Every nonnegative-root beta polynomial in the displayed x,y
rectangle automatically satisfies z<540/109<5. Consequently the entire
nonnegative-root rectangle, including its weak and discriminant boundaries,
is closed. Neither interlacer is needed for this stronger statement. It
contains a genuine open rectangle around (84,35), rather than only the
previously proved one-dimensional z fibre.

The remaining two-anchor model outside this rectangle is OPEN. Neighboring
model coefficients are not asserted to be factorial rows of an integer
Laurent support.

## 1. Inheritance, exact object and connection

The closest proved mechanism is the all-root fixed-fibre certificate in
`05-knowledge/results/long_frontier_sep06_residue_tail.md`, with its exact
admissible-domain companion `long_frontier_sep06_fibre_domain.md`. The
newly proved decoder in `continuing4_20260906_moments_packet.md` identifies
the full weak model with two ordinary degree-eight moment conditions and
shows that its fixed-(x,y) z fibres are compact intervals. A general fibre
can have a positive lower endpoint; no general z=0 admissibility is assumed.

The canonical hostile is the positive-response degree-seven surrogate from
that decoder report. The corrected near miss is promoting truncated residue
positivity to beta-root geometry, or dropping the original first equation
after elimination. The underused sidecar is the first critical point of
the nonnegative-root beta polynomial, before projecting geometry to a few
moment determinants.

The board is: exact degree-eight domain; compact z fibres; first critical
point; original phase brackets; and coefficient positivity on compact and
unbounded intervals. The source-to-target map retains the original x,y,z
pencil and all carries, then uses root geometry only to enclose its z fibre.
The target sign proof deliberately covers the larger prism(1). It loses
the precise admissible fibre endpoints but preserves every actual original
phase. The degree-seven and z=6 hostiles below check both losses.

Define, with x,y,z real,

    B(v)=v^5-13v^4+55v^3-xv^2+yv-z,
    C(v)=v^4-12v^3+45v^2-(2x/3)v+3y/7,
    D(v)=v^4-11v^3+36v^2-(5x/12)v+y/7.

If B has five nonnegative roots, their sum is13 and square sum59. In the
original variable t keep

    beta=t^-1(1+13t+55t^2+xt^3+yt^4+zt^5),
    C_raw=t^-1(1+12t+45t^2+(2x/3)t^3+(3y/7)t^4),
    D_raw=t^-1(1+11t+36t^2+(5x/12)t^3+(y/7)t^4),
    O=sum_j binom(14,2j+1)t^j, E=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star(beta^2+2t C_raw D_raw).         (2)

Star means coefficientwise multiplication. The relative factor2t and
both Laurent shifts remain fixed. Q has every exponent -1,...,8 and
q_(-1)=28. No replacement phase equation is used.

## 2. Beta geometry gives a uniform fifth-coefficient cap

Write F(v)=B(v)+z and assume B has five nonnegative roots, with
x in[83,85] and y in[34,36]. Its derivative has four nonnegative real roots
by Rolle's theorem, including multiplicities. Since F'(0)=y>0 and

    F'(8/25)
      <=5(8/25)^4-52(8/25)^3+165(8/25)^2-166(8/25)+36
       =-146524/78125<0,

the first critical point eta lies in(0,8/25). It lies between the first
two B roots. On that closed root interval B is nonnegative, also when the
two roots coincide. Thus z<=F(eta).

For 0<v<=8/25, the term v^5-13v^4 is strictly negative and 55v<=88/5.
Therefore

    F(v)<36v-(327/5)v^2 <= 540/109 <5.                 (3)

It follows that **0<=z<540/109**, including zero-root and repeated-root
boundaries. This bound uses beta geometry alone, not a determinant sign or
either interlacer. Through the degree-eight iff decoder, the same conclusion
is a consequence of its exact two-matrix domain restricted to this rectangle.
The first-critical-point argument is separate from the sign calculation.

## 3. Four intervals exhaust the original phases uniformly

From(2), the original first equation is

    P(-s)=2002 p(s),
    p(s)=z s^4-(12/7)y s^3+x s^2-10s+1/11.             (4)

For the full prism(1), the following endpoint signs are strict:

|Interval|Left sign|Right sign|
|---|---|---|
|(1/102,1/100)|+|-|
|(1/9,13/100)|-|+|
|(1,8/5)|+|-|

Also p(19/2)<0. Each endpoint expression is affine in x,y,z, so checking
the eight corners proves these signs on the entire closed prism. The
exact extrema are in the certificate and transcript. For example,

    min p(1/102)=10751/13618836>0,
    max p(8/5)=-41189/9625<0,
    max p(19/2)=-2058747/1232<0.

For z>0 the positive leading coefficient gives a fourth root in
(19/2,infinity). Four disjoint sign changes exhaust degree four, so all
four roots are positive and simple. For z=0 the degree is exactly three,
and the first three intervals exhaust it with simple roots. Thus no
uncontrolled phase remains, including the root escaping as z decreases
to zero. This argument does not assume that z=0 is an admissible endpoint
of every actual model fibre.

## 4. The full response is negative on every retained branch

On the first interval, retain the uneliminated polynomial T=sQ(-s).
Every raw q_j is a polynomial with nonnegative coefficients in x,y,z.
For each positive monomial of T evaluate x,y,z,s at85,36,5,1/100; for
each negative monomial use83,34,0,1/102. This gives the uniform bound

    T <= -728550005046322718853208807
          /1704830652000000000000000 < -427 <0.         (5)

For the other three intervals, eliminate z using precisely(4):

    z=(12/7)y/s-x/s^2+10/s^3-1/(11s^4).

After forming the complete response first, define h by

    sQ(-s)=-(14/11)h(s;x,y)                            (6)

at every original zero. It is a degree-eight polynomial. For an explicit
compact specification, write h=-(11/14) sum_(i=0)^8 r_i s^(8-i), where

|i|r_i|
|---|---|
|0|-(26075790/7)y^2|
|1|(153780300/7)xy-(179344800/7)y^2|
|2|-16900975x^2+(647843760/7)xy-1329865290y|
|3|-53986980x^2+1122025905x-4282905900y|
|4|3467704710x+(5521932000/11)y-10070260200|
|5|-(3690469830/11)x-9902880y-30313505040|
|6|6175260x+3654364350|
|7|-1022439600/11|
|8|565082|

Set x=83+2X,y=34+2Y, with0<=X,Y<=1. For each finite interval[a,b] form

    (1+u)^8 h((a+b u)/(1+u);x,y),

and for the unbounded interval form h(19/2+u;x,y). Expand in powers u^k,
0<=k<=8, and express each coefficient in the tensor Bernstein basis of
degree(2,2) in X,Y. All243 rational coefficients are strictly positive:

|Phase interval|Number of coefficients|Minimum coefficient|
|---|---:|---|
|[1/9,13/100]|81|10396002758034286225411/24500000000000000|
|[1,8/5]|81|1992478168568/49|
|[19/2,infinity)|81|165789872820/49|

The complete lists are frozen in the JSON certificate. Their construction
uses only the binomial theorem and exact rational arithmetic. The Bernstein
basis is nonnegative and sums to one on the rectangle, so each u coefficient
is strictly positive. Hence each transformed polynomial is positive for
u>=0. For a finite right endpoint the leading coefficient supplies h(b)>0.
This proves h>0 on all three closed intervals, including the whole
unbounded branch. Equation(6) gives Q(-s)<0 at their original phases.

Together with(5) and the degree exhaustion of Section3, this proves the
coefficient-prism theorem(1). Combining it with(3) proves the entire
nonnegative-beta rectangle statement, without either interlacer.

## 5. Nonvacuity, exact-domain hostiles and boundaries

For every(x,y) in the rectangle, z=1 supplies a B with five positive roots.
Its signs at v=0,1/10,1,3,5,7 are respectively -,+,-,+,-,+; all four
coefficient corners satisfy these strict signs. Degree exhaustion gives
five simple positive roots throughout that two-dimensional section.
This is a nonvacuity statement about beta geometry, not an assertion that
both interlacers hold at every point of that section.

At the central point(84,35,1), the original rows in(2) are the literal first
and doubled factorial rows of support(-27,1,15). The source independently
enumerates all charge-zero compositions at masses14 and28, recovering all
five and ten channels. The known actual fixed-row result is inherited;
the new conclusion is the free x,y rectangle and all its admitted z values.

Two hostiles preserve the stopping scope:

- The degree-seven surrogate x=78071/1000,y=601/50,s=57/2 with its original
  z has positive Q. Both degree-eight ordinary Hankel determinants are
  strictly negative, reconstructed exactly here. It is rejected by the
  exact domain, not accepted as a counterexample to the model.
- At x=84,y=35,z=6, the original p has a root in
  (16693/2000,41733/5000). A separate negative-coefficient transform proves
  h<0 on that whole interval, hence positive response at its original root.
  Its D degree-eight determinant is -25668. Thus dropping the z cap is
  false even inside this x,y rectangle; the hostile lies outside beta
  geometry and the exact interlacer domain.

The universal sign proof includes coefficient boundaries x=83,85 and
y=34,36, z=0 and5, and all escaping phases. Actual beta discriminant
boundaries require no limiting sign argument because they already lie
strictly inside the coefficient cap. No claim is made about the full
free-x,y model outside the rectangle, or a general actual Laurent family.

## 6. Reproduction

The [standalone source](../../04-computation/continuing5_20260906_moments_rectangle.py),
[transcript](continuing5_20260906_moments_rectangle.out), and
[certificate](continuing5_20260906_moments_rectangle_certificate.json)
use only Python's standard library. From the source directory:

    python continuing5_20260906_moments_rectangle.py
    python -O continuing5_20260906_moments_rectangle.py

Both runs pass **438 always-active exact gates** and produce byte-identical
raw LF stdout. The source explicitly configures LF; the capture wrapper
checks for carriage returns and does not normalize output. The finite
certificate universe comprises the full symbolic original convolutions,
same-zero elimination,56 phase endpoint corners,243 positive coefficients,
the independent transform evaluations, domain and nonvacuity controls,
all15 literal actual channels, and both exact-domain hostiles. Numerical
probes located convenient intervals only; no numerical root or optimization
result is used by the proof or verifier.

- Source SHA256: `da6236b1828b05bf5d30fe09d97e95005db36fe1baaa390c1b422aba85f13a25`.
- Output SHA256: `4aab430109b707da57cd2654ee0138d1b141002a1d67d33164cbe66ba0bb4a28`.
- JSON SHA256: `1e86411634a8f62f98a512d741d4b8df8ba5a069fb262d92bafe64e40569a4f8`.

No repository, maintained file or Git state was edited by this producer.
The proof remains a candidate until independent audit.

Independent [proof and exact referee](continuing5_20260906_moments_rectangle_audit.md) passes.

Filed checkpoint provenance: the [raw-byte manifest](continuing5_20260906_manifest.json)
pins the final report, source and output. Reviewed candidate-report hashes
above refer to the pre-promotion bytes. Source-location defaults and local
links were made portable where necessary; all emitted outputs were replayed
as raw bytes. The independent audit supplies the stated promotion basis.
