# Sharp stability when the two largest positive roots are macroscopic

**Status: PROVED ANALYTICALLY + FINITE-EXACT; independently audited analytically and computationally (root).**
Regional continuation of the frozen
[global bound](open_frontier_sep06_stability.md). A separate
[complementary proof](open_frontier_sep06_stability_complement.md) now
combines with this result to prove the sharp global constant. Both
proofs and their exact sources have passed mutual independent audit;
see the [complement audit](open_frontier_sep06_stability_complement_audit.md).
No theorem ID or priority claim.

## 1. Regional theorem and the scope of packing

Use the inherited normalization and quantities

    p_1=sum r_i=1, p_2=sum r_i^2=1, E=(1-p_4)/2>0,
    D=-[s^4]prod_i(1+r_i s)^2=(5-8p_3+3p_4)/6,
    J=D/E, c_*=(13-8sqrt(2))/3,
    a>=b>=0: the largest two positive entries,
    d^2=2-sqrt(2)(a+b).

**Theorem.** If b>=1/sqrt(3), then

    J-c_* > K_3 d^2,
    K_3=4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))].                 (1)

The constant is sharp as an infimum among these finite regional lists.
Both signs, zero entries and all integer multiplicities are retained.
No finite equality case exists. The region is an extra hypothesis;
this does not settle the best constant outside it.

The closest proved mechanism is the energy-adjusted signed envelope.
Its exact hostile relaxation permits a fractional number of equal roots
between two and three. The corrected near miss is treating a locally
positive Hessian as a global extremum theorem. The least-used sidecar is
the upper bound on the entire tail square mass in the present region.
The live board is: top two positive roots; energy; convex tail packing;
the changed first moment; three-root boundary geometry; actual dust lifts.

The map below replaces the tail by one positive root with the same square
mass and proves a one-sided inequality for a polynomial objective. It
preserves a,b,p_2 but changes p_1. It is **not** a transformation into
another admissible normalized cancellation list. The original first
moment makes the comparison strict and is restored exactly in the finite
sharpness family of Section 5.

## 2. Convexity gives a rigorous three-root objective comparison

Write

    u=sqrt(2), v=sqrt(3), h=1/u, z=1/v,
    A=2-u, B=u-1,
    gamma=3K_3/(4u)=(2-u)(3-v)/4,
    C_0=A+u gamma=(1+u+v-sqrt(6))/2,
    t=a+b, C=C_0-gamma t=A+(3K_3/8)d^2.

The target in (1), with a non-strict sign, is exactly

    F_actual=B-(C-A)-p_3+C p_4
            =(3E/4)(J-c_*-K_3 d^2) >=0.                  (2)

Put c^2=sum_(i>=3)r_i^2=1-a^2-b^2. The regional assumption gives
0<=c<=z<=b and t>=2z. Consequently

    0<C<=C_3:=C_0-2gamma z=(3-v)/2 <3v/8.

For 0<x<=c^2<=1/3, the function H_C(x)=x^(3/2)-Cx^2 obeys

    H_C''(x)=3/(4sqrt(x))-2C >=7v/4-3>0.                  (3)

Strict convexity and H_C(0)=0 imply
sum H_C(x_i)<=H_C(sum x_i) for any finite nonnegative list: apply
H_C(x_i)<=x_i H_C(sum x_i)/(sum x_i) and sum. A negative root loses
strictly when its cubic term is replaced by that of its absolute value.
It follows that

    p_3-Cp_4
     <= a^3-Ca^4+b^3-Cb^4+c^3-Cc^4.                     (4)

There is necessarily a negative tail entry: a+b>=2/sqrt(3)>1 whereas
p_1=1. Thus (4) is strict for every actual list in the theorem. The
surrogate has exactly three nonnegative entries a>=b>=z>=c>=0 and
square norm one. Its first moment is not claimed to equal one.

## 3. Concavity leaves two explicit boundaries

For this surrogate, retain t=a+b and c. Since
ab=(t^2+c^2-1)/2, its moments are

    p_3=(3t-t^3-3tc^2)/2+c^3,
    p_4=(3c^4-2c^2t^2-2c^2-t^4+2t^2+1)/2.

Let F(t,c) be (2) with these moments. At fixed 0<=c<=z, the
eligible interval for the top-root sum is exactly

    z+sqrt(2/3-c^2) <= t <= sqrt(2(1-c^2)).                (5)

The left endpoint has b=z; the right endpoint has a=b.
Differentiation gives

    F_tt=3t-2C_0(3t^2+c^2-1)+2gamma t(5t^2+3c^2-3).      (6)

This is strictly negative throughout (5). Here is a rational certificate
rather than a numerical concavity check. We have

    0<gamma<3/16, C_0>27/32, 8/7<2/v<=t<=u<10/7.

The coefficient of c^2 in (6) is
-2C_0+6gamma t< -27/16+(9/8)(3/2)=0. Thus (6) is at most

    Q(t)=10gamma t^3-6C_0t^2+(3-6gamma)t+2C_0.

Its derivative satisfies

    Q'(t)<(3/8)(15t^2-27t+5)<0,

because the convex quadratic in parentheses is negative at both
8/7 and 10/7, with values -307/49 and -145/49. Finally

    Q(2/v)<35/(4v)-81/16<0,

the last sign being 140^2<81^2*3. For completeness the parameter
bounds follow from 24/17<u<17/12 and 19/11<v<26/15:
gamma<35/187<3/16, and
C_0=1-(u-1)(v-1)/2>61/72>27/32. Every radical comparison is a
comparison of positive rational squares.

Therefore F is concave in t on (5), and its minimum lies at an
endpoint. It remains to check the two boundary configurations.

## 4. Positive boundary factorizations

On the boundary b=z, let z<=a<=u/v and
c=sqrt(2/3-a^2)<=z. Since C<=C_3, f_C(s)=s-Cs^2 is increasing on
[0,z]: its derivative is at least 2-v>0. Replacing the third contribution
by c^2 f_C(z) bounds F below by the signed tail envelope at b=z:

    F_env(a,z)=gamma(1-a)(a-z) P(a),
    P(a)=a^3-(1+u)a^2+(2/3)a+2/3.                         (7)

On this interval,
P'(a)<=8/3-2(1+u)/v<0, while
P(u/v)=(2u/3)(2/v-1)>0. Hence (7) is nonnegative.

On the other boundary put a=b=x, z<=x<=h, and c=sqrt(1-2x^2)<=z.
The envelope and its exact missing contribution are

    F_env(x,x)=2gamma(1-x)(x-z)(x-h),
    F(x,x,c)=F_env(x,x)+c^2(x-c)[1-C(x+c)].               (8)

The identities x-c=3(x-z)(x+z)/(x+c) and
h-x=c^2/[2(h+x)] turn this into

    F(x,x,c)=c^2(x-z) {
      3(x+z)/(x+c) [1-C(x+c)]-gamma(1-x)/(h+x)}.           (9)

Its braces are positive. Indeed x+c<=2z (differentiate
(1+s)/sqrt(2+s^2) for s=c/x in [0,1]), so
1-C(x+c)>=1-2zC_3=2-v>1/4. Also
(x+z)/(x+c)>=1, (1-x)/(h+x)<=1, and gamma<1/4.
The braces are therefore greater than 3/4-1/4=1/2.
Equation (9) is valid at both endpoints without division by its vanishing
prefactors, and is nonnegative throughout.

Concavity and these two boundary checks prove F(t,c)>=0. Combining
this with the strict inequality in (4) proves (1). The surrogate
objective vanishes only at the three equal roots or at the two equal
roots with c=0; neither is itself an actual p_1=1 list.

## 5. An actual finite family proves regional sharpness

For any integer n>=3, take the three raw positive roots

    x_n=1+1/n, x_n=1+1/n, y_n=1-2/n,
    L=n^2, Q_n=3+6/n^2,
    q_n=(9-Q_n)/[L(3+sqrt(Q_n+(9-Q_n)/L))].

Append exactly L negative roots -q_n, and divide every root by
S_n=3-Lq_n>0. The tuning quadratic

    9-Q_n-6Lq_n+L(L-1)q_n^2=0

is exact, so the normalized finite list satisfies p_1=p_2=1 and has
positive energy. Its two largest positive roots are a=b=x_n/S_n.
Because q_n<2/L,

    S_n^2=Q_n+Lq_n^2<3+10/n^2<=3(1+1/n)^2.

Thus a,b>1/sqrt(3) for every n>=3: this is an actual family inside
the region, unlike the equal-three-plus-dust family that approaches the
boundary from below. As n tends to infinity, S_n tends to sqrt(3), the
three positive roots tend to 1/sqrt(3), and the dust square mass tends
to zero. Its signed first moment is retained by the exact normalization.
The inherited moment formula then gives

    (J-c_*)/d^2 -> K_3.

This proves sharpness of (1) as a regional infimum without a finite
equality claim. Outside the region the tail square mass and coefficient
C need not remain in the convexity range, and a single packed root may
exceed the second-largest root. No global packing statement is inferred.

## 6. Reproduction and audit manifest

[Source](../../04-computation/open_frontier_sep06_stability_packing.py)
uses only the standard library, with exact quadratic and biquadratic
arithmetic. It checks the complete second-derivative polynomial identity,
both boundary factorizations, the diagonal correction after reduction by
c^2=1-2x^2, and the radical signs used on the full intervals.

The declared actual-row universe consists of the 24 prefixes
(1,x,y), where x is in {1/2,3/4,1} and y is in
{-1,-3/4,-1/2,-1/4,0,1/4,1/2,3/4}. Each receives its unique fourth
entry making e_2=0 and is normalized by its first moment. The output
prints all 17 exclusions; seven rows have positive energy and lie in
the regional domain. Literal product coefficients check their strict
target, signed tail and zero-padding invariance. The source additionally
checks the actual sharpness construction at n=3,4,8,20 (including an
independent literal-factor path at n=3,4), the two surrogate equality
shapes, and a strict rational three-root control.

Reproduce from the repository root:

```bash
python3 -B 04-computation/open_frontier_sep06_stability_packing.py
python3 -B -O 04-computation/open_frontier_sep06_stability_packing.py
```

All 72 explicit gates pass; the normal and optimized outputs are
byte-identical. The [frozen output](open_frontier_sep06_stability_packing.out)
has the following manifest:

```text
source SHA256 a31fc847573b9fbc2799863057870c0be90c7645bf2c98b76da0c516c202c09b
output SHA256 98bd68e9ffb7023bc838764768d6f56a3b840b21651381540f33ca9148ab654c
semantic SHA256 b3f72cc3c4a34f0c97fd3a41e29d2e71ad346b59f3e3a240a28cf3cab3a512f8
```

Root independently audited Sections 1–5 and reported **PASS**: the strict
negative-tail comparison, convexity range, exact fixed-c interval,
concavity certificate, both boundary factorizations and the actual
regional sharpness family were checked. Root also read the complete
frozen source and independently reproduced all 72 gates in normal and
optimized execution, obtaining byte-identical saved output and the two
declared file hashes: **PASS**. No numerical root search or finite-degree
extrapolation is a proof input.
