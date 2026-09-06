# An explicit fourth-coefficient floor and an effective Laurent phase tail

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent referee](continuing2_20260906_effective_anchor_audit.md)
accepts the complete proof without repair; 124 independent exact gates pass. This sharpens the
incoming [qualitative two-anchor theorem](open_frontier_sep06_quadratic_anchor.md)
within exactly its stated two-anchor, two-interlacer model. The incoming
proof supplies a qualitative positive floor and a qualitative uniform tail.
Here the new consequences are the explicit, nonoptimal bounds

    e4 > 1/100,
    s >= 118163898523 and P(-s)=0  ==>  Q(-s)<0.

The remaining finite-phase sign question and general actual Laurent
noncancellation remain OPEN. No new actual endpoint-family closure follows.

## 1. Inheritance and the retained predicate

The closest proved mechanism is the incoming compact weak-interlacing
argument: its hypothetical two-zero boundary is accepted by C and rejected
by D. The canonical hostile is the linearly anchored model with escaping
positive response. The corrected near miss is a finite search mistaken for
a uniform tail. The underused sidecar here is the exact sign at the fourth
and fifth ordered roots, before replacing roots by coefficient compactness.

The map retains all five nonnegative roots, both exact anchors, both
contiguous rows, the same original first zero, and the entire carried
response. It destroys no phase normalization. The new operation is a
quantitative perturbation of the rejected boundary, followed by cancellation
of the two largest carried coefficients. Its cheapest hostile is still the
exact repeated-root tuple (0,0,3,5,5). This is a model theorem; an arbitrary
real-rooted carrier without both interlacers does not inherit it.

## 2. A uniform explicit lower bound

Let 0<=a<=b<=c<=d<=e with sum 13 and square sum 59; equivalently e1=13,
e2=55. Write epsilon=e4 and f=e5 for their elementary symmetric functions.
Let

    B(v)=v^5-13v^4+55v^3-e3 v^2+epsilon v-f,
    C(v)=v^4-12v^3+45v^2-(2/3)e3 v+(3/7)epsilon,
    D(v)=v^4-11v^3+36v^2-(5/12)e3 v+(1/7)epsilon.

Assume C and D both weakly interlace B. In particular C(d)<=0 and D(e)>=0,
including repeated roots. Suppose epsilon<=1/100, toward a contradiction.

First c>2/3: Cauchy bounds d+e<=sqrt(118)<11, so a+b+c>2. Also d>1:
otherwise a+b+c+d<=4 and e>=9 contradict square sum 59. The nonnegative
term bcde<=epsilon therefore implies b<=(3/2)epsilon. Consequently

    t=a+b <= 3epsilon <= 3/100,
    f=a(bcde) <= (3/2)epsilon^2.

Eliminating B at the fourth root gives the exact identity

    C(d)=d^2(d-5)^2/3 - (5/21)epsilon + 2f/(3d).

The last term is nonnegative. Since C(d)<=0 and d>1, we obtain
|d-5|<1/10. Put delta=d-5. Then c<=13-2d<16/5, hence 5-c>9/5.
The remaining three roots have sum 13-t and pair sum
55-13t+t^2-ab. Substituting d=5+delta yields

    (5-c)(5-e)=-3t+2delta+t^2+delta t+delta^2-ab.

Using ab<=t^2/4, the absolute value is at most

    9/100 + 1/5 + 9/10000 + 3/1000 + 1/100 + 9/40000
      =2433/8000 <9/25.

Thus |e-5|<1/5. At the largest root, the second exact elimination is

    D(e)=e^2(7e^2-67e+157)/12 -(23/84)epsilon+5f/(12e).

Writing eta=e-5 gives 7e^2-67e+157=-3+3eta+7eta^2<=-53/25.
Since e>24/5 and f<=(3/2)(1/100)^2, it follows that

    D(e) <= -2544/625 +1/76800 <0.

This contradicts D(e)>=0. Hence e4>1/100 throughout the whole weak class.
The floor is not asserted optimal. The argument includes zero roots and
repeated roots, so no compactness or limiting genericity is needed.

## 3. Retain the leading cancellation to obtain an explicit tail

Use exactly the incoming original coefficient normalization. Define

    G_B=(1,13,55,e3,e4,e5),
    G_C=(1,12,45,(2/3)e3,(3/7)e4),
    G_D=(1,11,36,(5/12)e3,(1/7)e4),
    beta=t^-1 G_B, C_raw=t^-1 G_C, D_raw=t^-1 G_D,
    O=sum_j binom(14,2j+1)t^j, E=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star (beta^2+2t C_raw D_raw).

Here star is coefficientwise multiplication, and tuple notation means the
ordinary polynomial with those coefficients. Write Q(t)=sum_(j=-1)^8 q_j t^j.
All q_j are nonnegative. The original first equation is

    P(-s)=182-20020s+2002e3 s^2-3432e4 s^3+2002e5 s^4.

At any positive root s, including a multiple root, this gives

    e5 s=(12/7)e4-e3/s+10/s^2-1/(11s^3)
         <=(12/7)e4+10/s^2.

The leading carried coefficients, with the crossing term retained, are

    q8=binom(28,18)e5^2,
    q7=binom(28,16)(2e4 e5+(6/49)e4^2).

Since (12/7)binom(28,18)-2binom(28,16)=-38346750<0,

    q8 s-q7 <= -(6/49)binom(28,16)e4^2
                +10binom(28,18)e5/s^2.

For the remaining coefficients, all shifts preserve coefficient mass at 1.
AM-GM gives G_B(1)=product(1+a_i)<=(18/5)^5 and e5<=(13/5)^5.
Coefficientwise G_C(1),G_D(1)<=G_B(1). The multiplier coefficient at j is
binom(28,2j+2)<=binom(28,14). Consequently

    sum q_j <= 3binom(28,14)(18/5)^10.

For s>=1 every lower term j<=6, including the negative-power carry, has
absolute contribution at most q_j/s after division by s^7. Combining with
the strict fourth-coefficient floor gives, at every original positive root,

    Q(-s)/s^7 < -delta0 + H/s,
    delta0=(6/49)binom(28,16)/10000 =2607579/7000,
    H=3binom(28,14)(18/5)^10+10binom(28,18)(13/5)^5
      =17194291313831660508/390625.

The integer 118163898523 is strictly greater than H/delta0. This proves
the displayed effective tail. It does not replace the finite-phase interval
by a verified census. Neither the threshold nor the floor is optimized.

## 4. Controls, independent review and reproduction

The finite universe is the formal root-elimination and carry identities,
the displayed rational inequalities, one strict rational model with roots
(1,3,9,22,30)/5, and the weak boundary (0,0,3,5,5). Exact alternating signs
verify both strict interlacings in the positive control. At the boundary,
C=v(v-2)(v-5)^2 still interlaces B, while D(5)=-25/4 rejects it. These
controls are not a census of the compact model; the inequalities above
supply all universal quantifiers.

The [standalone producer](../../04-computation/continuing2_20260906_effective_anchor.py)
imports no repository implementation; its [output](continuing2_20260906_effective_anchor.out)
is frozen with the independent audit. It checks
the eliminations as formal polynomial identities and reconstructs all
carried coefficients from the original Laurent shifts and convolution.
The producer passes 78 always-active exact gates. Root replayed normal,
optimized and frozen output agreement for both the producer and the
independent standard-library referee. The latter reconstructs the E/O
multiplier independently before compiling every original carry coefficient.
The universal inequalities are proved above; the finite controls check their
algebra and failure boundaries.

    python -B 04-computation/continuing2_20260906_effective_anchor.py
    python -B -O 04-computation/continuing2_20260906_effective_anchor.py

The next target is a finite-phase sign argument on 0<s<118163898523,
or a stronger boundary estimate that interacts with the exact first-root
equation. The positive-algebraic amplitude repairs on individual actual
rows do not supply uniform coefficients over this model class.
