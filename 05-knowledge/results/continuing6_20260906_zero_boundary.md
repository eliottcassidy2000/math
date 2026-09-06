# Every original phase on the zero-beta boundary has negative full response

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The entire zero-root boundary of the anchored model has
strictly negative full carried response at every original phase. A stronger
coefficient-packet theorem uses only two low-order C/B Gram constraints.
Thus the conclusion also holds on the C-only zero-beta boundary, with
repeated roots allowed. The second interlacer is not needed as a hypothesis,
but its original carrier remains in the response.

The remaining shared-root boundaries and the general anchored sign question
remain OPEN. This is not a claim that every coefficient pair is an actual
factorial row, or that the general Laurent noncancellation problem is solved.

[Independent proof and exact audit](continuing6_20260906_zero_boundary_audit.md) passes.

## 1. Inheritance and the retained coordinate

The closest proved mechanism is
[continuing5_20260906_simple_beta_boundary.md](continuing5_20260906_simple_beta_boundary.md), with
its independent audit. It excludes repeated beta roots in the full weak
two-interlacer model and gives exact zero-atom peeling. The incoming
[long_frontier2_sep06_boundary.md](long_frontier2_sep06_boundary.md) reduces fixed-phase extrema to zero or
shared-root boundaries after that exclusion. The larger coefficient-prism
theorem in [continuing5_20260906_moments_rectangle.md](continuing5_20260906_moments_rectangle.md) closes a neighborhood
of (84,35), but not the entire zero-root boundary.

The canonical hostile is the C-only repeated beta configuration (0,0,3,5,5).
It invalidates extending the full-model simplicity theorem to C alone; the
new sign theorem deliberately includes it. The corrected near miss is
dropping the original phase after eliminating a coefficient. The underused
sidecar is the ordinary three-by-three moment condition for C/B alongside
its shifted two-by-two condition, before retaining the full root geometry.

The live board is: fixed zero atoms; the original phase; the reduced
quadratic response; a necessary curved coefficient region; raw separate
convexity; and positive coefficient charts. The map solves for y at the
original zero with z=0. It preserves that zero, all relative factors, the
binomial parameter14, and the inverse carry. It forgets some beta geometry;
the resulting larger coefficient region is proved sufficient for the sign.

## 2. A stronger coefficient-packet theorem

For real x,y keep the prescribed carriers at z=0:

    B(v)=v(v^4-13v^3+55v^2-xv+y),
    C(v)=v^4-12v^3+45v^2-(2x/3)v+3y/7,
    D(v)=v^4-11v^3+36v^2-(5x/12)v+y/7,

    beta(t)=t^-1(1+13t+55t^2+xt^3+yt^4),
    C_raw(t)=t^-1(1+12t+45t^2+(2x/3)t^3+(3y/7)t^4),
    D_raw(t)=t^-1(1+11t+36t^2+(5x/12)t^3+(y/7)t^4),
    O(t)=sum_j binom(14,2j+1)t^j,
    E(t)=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star(beta^2+2t C_raw D_raw).

Star is coefficientwise multiplication. In particular q_(-1)=28; neither
the inverse shift nor the mixed response is removed.

**Theorem.** Suppose

    x>=75,     0<=y<=c(x):=(7/72)(x-75)(135-x).             (1)

Then every zero of P is strictly negative and simple, and

    Q(-s)<0 whenever s>0 and P(-s)=0.                      (2)

If y>0 there are three zeros; if y=0 there are two. No beta-root or
interlacer hypothesis is needed beyond the explicit coefficient packet (1).

**C-only zero-boundary corollary.** If B has five nonnegative roots with
the fixed sum13 and square sum59, has a zero root, and C weakly interlaces
B, then (2) holds at every original phase, including repeated-root cases.
In particular the entire zero-root boundary of the full C/D model is closed.

To prove the corollary, weak interlacing gives a positive probability
residue measure for C/B, including canceled/repeated beta roots as in
[continuing3_20260906_residue_floor.md](continuing3_20260906_residue_floor.md), Sections 1--2. Its moments begin

    mu0=1, mu1=1, mu2=3,
    mu3=x/3-16, mu4=16x/3-373-4y/7.

Nonnegative beta nodes make the shifted matrix PSD, so

    det(mu_(i+j+1))_(i,j=0)^1=(x-75)/3>=0.

The ordinary matrix is PSD, so

    det(mu_(i+j))_(i,j=0)^2
      =(x-75)(135-x)/9-8y/7>=0.                          (3)

Also y>=0 by beta geometry and z=0 by the zero root. These are exactly (1).
This is only the necessary direction from geometry to the packet: the
packet is not asserted to characterize interlacing or nonnegative roots.

The preserved zero-atom identities from the preceding theorem remain useful
for interpreting the full simple-zero boundary. The present proof needs
even less geometry and does not change O14 or E14 to lower-height carriers.

## 3. The packet supplies all original phase branches

The original phase polynomial is

    P(-s)=2002 p(s),
    p(s)=-(12/7)y s^3+x s^2-10s+1/11.                    (4)

From (1), 75<=x<=135 and

    0<=y<=175/2-(7/72)(x-105)^2<=175/2.

The rectangular bounds give

    p(1/110)>=81/13310>0,
    p(1/90)<=-7/1980<0,
    p(63/1000)<=-7207/2200000<0.                         (5)

For fixed s>0, substitute the upper bound c(x) for y. Completing the
square in x gives a lower bound whose quadratic coefficient is s^3/6 and
whose minimum occurs at x=105-3/s. In particular

    p(1/8)>=3/2816>0,
    p(9/16)>=3335/22528>0.                              (6)

For y>0, the negative leading coefficient gives p(s)<0 for large s.
There is therefore one root in each of the disjoint intervals

    I1=(1/110,1/90),
    I2=(63/1000,1/8),
    I3=(9/16,infinity).                                 (7)

The three sign changes exhaust the cubic degree, so these are all its
roots, each simple. If y=0, the degree is two with x>0; the first two
intervals exhaust it. No finite third root is asserted at y=0.

## 4. The exact restricted response and feasibility inequality

At an original zero, solve (4) exactly:

    y=Y(x,s)=7x/(12s)-35/(6s^2)+7/(132s^3).              (8)

Reconstructing the full carried response before substitution yields

    s Q(-s)=-F(x,s)/968,                                (9)

where the complete eleven-term polynomial is

    F= 5182210385 s^6 x^2+8439169200 s^5 x^2
       -235681249320 s^5 x-585418008010 s^4 x
       +2361406947900 s^4+38085709200 s^3 x
       +6069283306240 s^3-385825440 s^2 x
       -498310898015 s^2+8357151600 s-38651536.

Its x^2 coefficient is

    A(s)=264385 s^5(19601s+31920)>0.                     (10)

Thus the restricted response is strictly concave in x, so a mere check of
segment endpoints would not prove its negativity. The proof below retains
the location of its vertex through its discriminant and derivative.

Substituting (8) in y<=c(x), multiplying by the positive quantity
(72/7)s^3, gives the necessary inequality

    J(x,s):=s^3(x^2-210x+10125)+6x s^2-60s+6/11<=0.     (11)

This inequality is imposed only at the original phase. It is the additional
coordinate that rules out the potentially adverse portions of the quadratic.

## 5. Complete positive certificates on all three branches

For a univariate polynomial f of degree at most d on [a,b], the certificate
uses the exact polynomial

    (1+u)^d f((a+bu)/(1+u)).                             (12)

Strictly positive coefficients prove f>0 on the closed interval, including
the right endpoint by its leading coefficient. The complete rational
coefficient lists are retained in the JSON; positivity is not inferred from
samples. There are107 positive coefficients in total.

**First branch.** Use the raw T(x,y,s)=sQ(-s), before (8). Direct coefficient
extraction gives

    T_xx=397670 s^5(66-85s)>0,
    T_yy=(4011660/7)s^7(140-13s)>0

on I1. Separate convexity bounds its maximum on the containing rectangle
[75,135] x [0,175/2] by its four corner values. All four polynomials -T at
those corners have nine strictly positive coefficients under (12), d=8,
on [1/110,1/90]. This proves T<0 throughout I1, even without imposing the
phase equation inside that interval.

**Lower second branch, 63/1000<=s<=11/100.** Let Delta(s)=disc_x F.
The polynomial -Delta(s)/s^4 has seven strictly positive coefficients
under (12). Hence Delta<0 there. With A(s)>0 from (10), the whole quadratic
F(x,s) is strictly positive for every real x.

**Upper second branch, 11/100<=s<=1/8.** Three certificates prove

    J(105,s)>0,       F(105,s)>0,       F_x(105,s)<0.     (13)

Since J_x=s^2(2s(x-105)+6)>0 for x>=105, feasibility (11) forces x<105.
Because F_x is increasing in x, F_x(x,s)<0 for every x<=105. Therefore
F(x,s)>=F(105,s)>0 on the whole feasible upper-second branch.

**Third branch.** A four-coefficient certificate gives

    J(84,s)>0 for 9/16<=s<=3/5.

For x<=84 on that interval,

    J_x/s^2=2s(x-105)+6<=6-42(9/16)<0.

Thus J(x,s)>=J(84,s)>0, contradicting (11), so every feasible point there
has x>84. The polynomial F(84+X,9/16+u) has all21 coefficients strictly
positive for 0<=deg_X<=2 and 0<=deg_u<=6. Hence F>0 throughout this part
of the branch. Finally F(75+X,3/5+u) also has all21 coefficients positive,
covering every x>=75 and s>=3/5. This covers the entire unbounded branch.

For reproducibility, the univariate certificate sizes and minima are:

| Positive polynomial and interval | Count | Minimum coefficient |
|---|---:|---|
| -T(75,0,s), I1 closure | 9 | 340209247973/566899520 |
| -T(75,175/2,s), I1 closure | 9 | 8128988108758357/13718968384000 |
| -T(135,0,s), I1 closure | 9 | 6831521934461/14172488000 |
| -T(135,175/2,s), I1 closure | 9 | 260452122980509/548758735360 |
| -disc_x F/s^4, [63/1000,11/100] | 7 | 17987365736233620424877833708131/2500000000000000 |
| J(105,s), [11/100,1/8] | 4 | 40761/110000 |
| F(105,s), [11/100,1/8] | 7 | 32275785602081/262144 |
| -F_x(105,s), [11/100,1/8] | 7 | 260266674481272503/20000000000 |
| J(84,s), [9/16,3/5] | 4 | 2008239/45056 |

The two21-coefficient quadrant arrays complete the107 coefficients. Together
with the exhausted original phase intervals, these prove the theorem.

## 6. Scope controls and the remaining boundary problem

The full-model point (x,y,z)=(84,35,0) is nonvacuous: its remaining quartic
has four positive roots, and the C/B and D/B full ordinary moment matrices
have leading principal minors respectively

    (1,2,11,188,7125),     (1,3,26,705,30600).

The source independently verifies these facts. The C-only repeated-root
control (75,0,0) lies on (1), so its original phases are covered even though
D does not interlace it. No simplicity assumption was introduced for the
new theorem.

The packet itself cannot simply be discarded. At

    x=84, y=1050/11, z=0, s=1/6,

the original equation holds, but

    sQ(-s)=228261457669/313632>0,
    det(mu_(i+j))_(i,j=0)^2=-639/11<0.

This is a coefficient hostile, not an admissible C-interlacing model.
Conversely, at the admissible point (84,35,0) but the non-phase s=1/4,

    sQ(-s)=412444713751/32768>0,
    p(s)=335/176 !=0.

Thus the original-zero requirement is also essential. These two controls
separately test the geometry-packet and phase projections.

Combining the proved full-model simplicity theorem, the incoming boundary
extremum theorem and this zero-boundary closure gives a concrete reduction:
**if a nonnegative or positive response exists in the full model at an
original phase s, then a (possibly different) admissible shared-root shape
at that same s has a response with the same sign property.** The converse
is immediate because the shared-root shapes are in the full model. This
is an existence reduction, not a claim that every nonnegative point itself
lies on the boundary. The zero branch cannot host it. No sign on those remaining shared-root branches is
claimed here; neither an exterior resultant zero nor a lower-moment
surrogate is an admissible substitute.

## 7. Reproduction and freeze

The standalone source imports no producer or repository implementation.
It reconstructs P,Q by the actual Hadamard convolutions, the C/B moments by
formal division, the complete restricted response, feasibility inequality,
all107 coefficient entries, phase bounds, and positive/hostile controls.
All arithmetic is exact. The uniform theorem follows from the identities
and complete coefficient arrays, not from a finite sampling claim.

    python continuing6_20260906_zero_boundary.py
    python -O continuing6_20260906_zero_boundary.py

Both runs pass **168 always-active exact gates** and produce identical raw
LF stdout. The source configures newline handling explicitly; captures were
written as bytes without text normalization. The JSON destination is
adjacent outside the repository, or 05-knowledge/results after filing the
source in 04-computation.

- Source SHA256: `dd28188852e5fb955d19d886687f2cb791e25cbf1a6f1c3b2bf66158cc467e6b`.
- Output SHA256: `6fe1ddf5e2182828405968b34c8236bc8ccea8f5b046e843781dcf26362fef0c`.
- Certificate SHA256: `f6c169d4190f96465f93d195fe4517537a9c40e85370ae6a736e0e2886e9b12e`.

No previous frozen artifact, repository file, namespace or Git state was
changed. The independent proof and coefficient audit above passed; the manifest pins the filed artifacts.
