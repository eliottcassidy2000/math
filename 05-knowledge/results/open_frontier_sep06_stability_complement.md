# The complementary region closes the sharp global stability constant

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
Together with the independently constructed
[regional packing theorem](open_frontier_sep06_stability_packing.md),
this proves the exact best constant in the finite signed-root stability
problem. The [independent complementary audit](open_frontier_sep06_stability_complement_audit.md)
and [third-reader global audit](open_frontier_sep06_stability_global_audit.md)
both passed; root separately audited the full regional proof and source.
The consolidated canonical theorem is
[THM-4454 / sharp global signed-root duplication stability](../../01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md).

## 1. Exact target and inheritance

For any finite real list with p1=p2=1, set

    E=(1-p4)/2>0, D=(5-8p3+3p4)/6, J=D/E,
    c_*=(13-8sqrt(2))/3,
    a>=b>=0: the largest two positive roots, padded with zeros,
    d^2=2-sqrt(2)(a+b).

The claimed sharp global theorem is

    J-c_* > K3 d^2,
    K3=4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))].                 (1)

Every finite eligible list satisfies the strict inequality. The constant
is the infimum of the quotient (J-c_*)/d^2, approached by actual finite
lists from the regional theorem. There is no finite equality case.
Also d^2>0: equality in a+b<=sqrt(2) would force exactly the two
positive roots1/sqrt(2), with every other root zero, contradicting p1=1.
This solves this quantitative signed-root problem; it does not prove
actual binomial two-rung separation or LRC(14).

The closest proved mechanism is the
[energy-adjusted envelope](open_frontier_sep06_stability.md), whose constant
K0 is optimal for its unrestricted fractional-tail relaxation. The incoming
regional result retains the whole tail square mass when b>=1/sqrt(3).
The complementary observation is that the old envelope already suffices
at K3 when b<=1/sqrt(3). The formal obstruction to the old method lies in
the other region, where the new packing argument repairs it.

The live concepts are energy, signed tails, the second-positive-root
threshold, continuous envelope, integer packing and actual finite dust.
The connection keeps the objective and normalized target identical and
partitions the domain at the same threshold. No transformed surrogate
is asserted to preserve p1. The largest equality shapes in the relaxed
domain have first moment greater than one; signed feasibility makes the
actual inequality strict.

## 2. Signed secant on the complementary domain

Write

    u=sqrt(2), v=sqrt(3), r=1/v, h=1/u,
    A=2-u, B=u-1, alpha=3K3/8, gamma=u alpha,
    C0=A+2alpha=(1+u+v-uv)/2,
    C=C0-gamma(a+b)=A+alpha d^2.

The target gap has the exact form

    F_actual=B-alpha d^2-p3+C p4
            =(3E/4)(J-c_*-K3 d^2).                     (2)

Suppose b<=r. Since a>=b, C<=C0-2gamma b. The function
b(C0-2gamma b) is increasing on [0,r], because
C0-4gamma r>0. Also C>=A>0. Consequently

    1-2Cb >= 1-2r(C0-2gamma r)=2-v>0.

The scalar f_C(x)=x-Cx^2 is strictly increasing on [0,b]. Every
positive root after a is at most b; every nonpositive root contributes
at most zero to p3-Cp4. Thus the signed secant gives

    F_actual >= F(a,b),
    F=B-alpha d^2-a^3-(1-a^2)b
      +(C0-gamma(a+b))[a^4+(1-a^2)b^2].                 (3)

All signs and arbitrary finite root counts are retained. This is a
one-sided objective bound, not a change of admissible root list.

## 3. Monotonicity and three boundaries

At fixed a, differentiate (3):

    F_b=(1-a^2)T,
    T=gamma(1+a^2-2ab-3b^2)+2C0 b-1.

On a>=b>=0, a^2+b^2<=1,

    T<=U(b)=gamma(2-6b^2)+2C0 b-1,
    U-T=gamma[(1-a^2-b^2)+2b(a-b)]>=0.

Here U is strictly increasing on [0,r], since
U'(b)>=2C0-12gamma r>0, while U(r)=2C0/v-1<0.
The latter sign is equivalent to uv>1+u, whose squared comparison
reduces to9>8. All auxiliary radical inequalities are also checked
exactly by the standalone biquadratic-field verifier.

For a<1, F therefore decreases strictly with b. It suffices to put

    b=min(a,r,sqrt(1-a^2)).                              (4)

**First boundary: 0<=a<=r, b=a.** Exact factorization yields

    F(a,a)=2gamma(1-a)(a-r)(a-h)>=0.                     (5)

It is positive for a<r and vanishes only at a=r in this interval.

**Second boundary: r<=a<=u/v, b=r.** Here

    F(a,r)=gamma(1-a)(a-r)H(a),
    H(a)=a^3-(1+u)a^2+2a/3+2/3.                        (6)

On this interval H'(a)<=8/3-2(1+u)/v<0, and

    H(u/v)=2u(2v-3)/9>0.

Thus (6) is nonnegative, with equality only at a=r.

**Third boundary: u/v<=a<=1, b=sqrt(1-a^2).** At a=1,
F(1,0)=0. Otherwise write t=a+b, E2=a^2b^2 and
g2=B-a^3-b^3+A(a^4+b^4). The exact two-atom identity is

    g2/[E2(2-ut)]=u(At+1)/(t+1)^2.

Here 1<t<=u/v+1/v<u, so both denominator factors are positive.
Moreover the associated quotient satisfies

    4g2/[3E2(2-ut)]
      >=4u(3-u)/[3(1+u)^2]
       =4(13u-18)/3 >1/2>K3.                           (7)

The first bound uses t>=1 and t<=u separately, which is sufficient
on this boundary. The strict lower comparison in (7) follows from
104sqrt(2)>147, checked by21632>21609. The final comparison follows,
for example, from v/(1+v)<2/3 and u>1, giving K3<4/9<1/2.
Equation (2) for the surrogate therefore gives F(a,b)>0 when a<1.

These cases prove F>=0 over the entire complementary domain. Its only
zeros are (a,b)=(r,r) and (1,0). The latter forces E=0 and is excluded.
In the former, a+b>1 forces a negative tail entry, so the secant in (3)
is strict. Thus every actual eligible list with b<=r satisfies (1).

## 4. Global conclusion, equality and actual sharpness

The [packing theorem](open_frontier_sep06_stability_packing.md) supplies
the same strict inequality for b>=r. The domains overlap at the exact
same threshold and exhaust all eligible lists. Its proof has a strict
signed-tail comparison because a+b>1 in that region. Consequently (1)
holds globally, with no finite equality case.

That theorem also constructs, for every integer n>=3, raw positive roots
(1+1/n,1+1/n,1-2/n), exactly n^2 negative dust roots with a displayed
algebraic tuning, and a positive normalization making p1=p2=1. The two
largest normalized roots exceed r; all three tend to r, while dust square
mass tends to zero. The limiting quotient is exactly K3, with positive
limiting distance and energy. Hence no larger global constant is possible.
The sharp bound is an infimum over actual finite lists, not over a
fractional-multiplicity relaxation or an infeasible positive-only list.

## 5. Source, audit and finite universe

[Source](../../04-computation/open_frontier_sep06_stability_complement.py)
uses only the standard library. It implements exact Q(sqrt2,sqrt3)
arithmetic and sign comparison, then checks formal polynomial identities
for the derivative, domain remainder, diagonal, fixed-r boundary and
independent two-atom moment elimination. Radical inequalities certify the
full intervals in the proof; no phase grid supplies a proof input.

The separate actual-row controls take the49 ordered rational prefixes
from {-2,-1,-1/2,0,1/2,1,2}^2, append the unique third root making e2=0
when the prefix sum is nonzero, and normalize by the full first moment.
Every exclusion is printed:19 zero-prefix or zero-energy cases. The30
eligible actual rows include16 complementary and14 packing cases.
Literal product/square coefficients independently check the strict
global target. These are controls, not an exhaustion argument.

Reproduce:

```bash
python3 -B 04-computation/open_frontier_sep06_stability_complement.py
python3 -B -O 04-computation/open_frontier_sep06_stability_complement.py
```

All90 explicit gates pass, with normal and optimized outputs byte-identical.
The [frozen output](open_frontier_sep06_stability_complement.out) is pinned by:

```text
source SHA256 bc1c612f840d5e9f1e8135fa17086f3244442d13962e8e37b633080c77c2d5a1
output SHA256 c43d2e8da2947154415299d155766777a3384754221d4526dff92b00ed3dda9d
semantic SHA256 81ddb22d7f5948976e3f7e11912bc534c6fab16249c6f3a84a8cdfee415bbf3c
```

The independent analytic/source audit and normal/optimized/frozen replay
passed with exactly these hashes. A third reader separately checked both
regions, global strictness, actual sharpness, and both90/72-gate bundles.
