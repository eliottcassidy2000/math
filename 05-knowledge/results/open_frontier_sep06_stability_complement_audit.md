# Independent audit of the sharp global stability completion

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT — PASS.**

Auditor: `three_ray_geometry`, independent of the complementary proof's
author, root. The auditor authored the separately checked
[regional packing theorem](open_frontier_sep06_stability_packing.md).
Root independently audited that theorem and its frozen exact source.
This note records the reciprocal audit of the
[complementary theorem](open_frontier_sep06_stability_complement.md),
which closes the global sharp constant

    K3=4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))].

The target is the finite real-list problem p1=p2=1, E=(1-p4)/2>0,
J=(5-8p3+3p4)/(6E), with a,b the largest two positive entries and
d^2=2-sqrt(2)(a+b). It asserts the strict inequality
J-c_*>K3 d^2, where c_*=(13-8sqrt(2))/3. All signs, zero entries
and integer multiplicities remain allowed. The constant is an infimum
over actual finite normalized lists; no finite equality is asserted.

## 1. Analytic checks

The target identity was independently reconstructed:

    (3E/4)(J-c_*-K3 d^2)
      =B-alpha d^2-p3+(A+alpha d^2)p4.

Here A=2-sqrt(2), B=sqrt(2)-1 and alpha=3K3/8. The identity
retains the energy factor and uses the actual normalized moments.
The distance is positive: Cauchy gives a+b<=sqrt(2), and equality
would force exactly two nonzero entries 1/sqrt(2), contradicting p1=1.

For b<=r=1/sqrt(3), the secant estimate is valid on every positive
tail entry because x-Cx^2 increases on [0,b]. For a nonpositive entry,
x^3-Cx^4<=0<=x^2(b-Cb^2); at a negative entry the comparison is
strict. The bound therefore preserves the signed finite-list scope.
It is an objective comparison and makes no false assertion that the
surrogate still has first moment one.

I differentiated the envelope, checked its displayed derivative factor,
and independently verified

    U-T=gamma[(1-a^2-b^2)+2b(a-b)]>=0.

The signs needed for U to increase to a negative endpoint are exact.
For example, with u=sqrt(2), v=sqrt(3),

    U'(r)=7-2u-(5-2u)v>0.

Both terms being compared are positive, and their squared difference
is 32u-42>0, certified by 2048>1764. The additional secant monotonicity,
negative U endpoint, boundary derivative and two-atom constant signs
were checked through the auditor's independent biquadratic arithmetic.

The upper endpoint min(a,r,sqrt(1-a^2)) exhausts the feasible b
interval. Its three cases are exactly a<=r, r<=a<=sqrt(2/3),
and sqrt(2/3)<=a<=1. I checked the diagonal factor and fixed-r
cubic factor; the latter agrees with the independently derived regional
boundary identity. I also reconstructed the two-atom moments and the
cleared quotient identity. Its coarse lower bound

    4(13sqrt(2)-18)/3>1/2>K3

uses only 1<=a+b<=sqrt(2), so it holds throughout the required branch.
The branch endpoints are handled separately before division.

The only zeros of the complementary envelope are (r,r) and (1,0).
The second forces E=0. At the first, a+b>1 and p1=1 force a negative
tail entry, making the secant strict. Thus the actual complementary
inequality is strict. The regional theorem covers b>=r with the same
target and constant. The overlap is harmless and the two regions exhaust
all eligible lists.

Finally, the actual regional sharpness sequence has three positive
macroscopic entries tending to r and exactly n^2 negative entries whose
square sum tends to zero. Its exact tuning retains p1=p2=1, positive
energy, and b>r. The limiting energy and distance are positive, and its
quotient tends to K3. This proves global optimality over actual lists,
without fractional multiplicities or an inadmissible limiting first moment.

## 2. Source and independent replay

I read the complete
[source](../../04-computation/open_frontier_sep06_stability_complement.py):
field arithmetic and exact signs, formal polynomial derivative and
boundary identities, the independent two-atom elimination, the complete
49-prefix rational control universe, every exclusion, normalization,
and the literal coefficient consumer. Gates use explicit exceptions and
remain active under optimized execution. No numerical sweep is used to
prove a continuous-domain sign.

I independently ran normal and optimized Python and compared both outputs
byte-for-byte with the saved output. All **90 gates PASS**. The 49 ordered
prefixes give 19 explicit exclusions and 30 eligible rows, with 16 in
the complementary region and 14 in the packing region. The file hashes
independently match the author's manifest:

```text
source SHA256 bc1c612f840d5e9f1e8135fa17086f3244442d13962e8e37b633080c77c2d5a1
output SHA256 c43d2e8da2947154415299d155766777a3384754221d4526dff92b00ed3dda9d
semantic SHA256 81ddb22d7f5948976e3f7e11912bc534c6fab16249c6f3a84a8cdfee415bbf3c
```

No mathematical correction or source defect was found. The result solves
the stated quantitative signed-root constant problem. It asserts no
transport to actual Laurent two-rung separation or to LRC(14) without
an additional map.
