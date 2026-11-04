# Exact ordered-cover reduction for an equal-profile cyclic port core

**PROVED:** an exact formula for every ordered-cover coordinate of the physical two-port cyclic comparison, and a cancellation of the separated-port-slot contribution. **NOT PROVED:** coordinatewise domination for every tournament profile or lifting the equal-profile theorem through every strongly connected ambient context. An explicit nonnegative-sequence obstruction shows why the unconditional scalar positivity theorem alone is insufficient for that extension.

Use the notation and physical query of the adjacent `covered_response/PROOF.md`: T=C3[A1,A2,A3], U=C3^op[A1,A2,A3], all three child profiles equal q, x->y, a=V1 and c=V2. Write T+ and U+ for these actual extensions. Put P=XYZ, R=1-P and

    E0=X+Y+Z+XY+YZ+ZX+3XYZ,
    K=X²Y²+X²Z²+Y²Z² +X²Y+X²Z+XY²+XZ²+Y²Z+YZ²
       -XYZ(X+Y+Z)-6XYZ.

Let L_q(X^rY^sZ^t)=q_r q_s q_t, with q_0=0 and zero extension beyond the child order.

## 1. All ordered-cover coordinates

For every h>=1,

    q_h(T+) - q_h(U+)
       = (h/3) L_q( K E0^(h-1) / R^(h+1) ).             (1)

The coefficient extraction is finite. The right side is always an integer for actual tournament profiles, despite its displayed division by three.

To prove (1), use the cyclic old-run matrix G=N/R from the companion proof. The generating series for one nonempty old-only path slot is E=1^t G 1=E0/R, symmetric in X,Y,Z and unchanged by quotient converse. Let F be the two-port word series, and let

    Fx=(1+1^t G a)(1+(1-a)^t G 1),
    Fy=(1+1^t G c)(1+(1-c)^t G 1).

These are single-port slot series, allowing the port alone. An ordered h-cover either has both ports in one slot, giving h F E^(h-1), or has ports in distinct slots, giving h(h-1) Fx Fy E^(h-2). The latter term is absent for h=1. All remaining slots are nonempty old paths. Reading each slot in order and then the slots in order supplies the child-run ordering, so the same L_q converts these series into exact covers of physical names.

Write the numerator of F-F^op as D over R^3 and the numerator of FxFy-(FxFy)^op as D2 over R^4. The companion identity is sum_(pi in S3) D(pi X)=2RK. Also

    sum_(pi in S3) D2(pi X)=0.                            (2)

Indeed, exchanging X and Y reverses the cyclic quotient while exchanging the two single-port masks. Their product is unchanged by exchanging the ports, so D2 is antisymmetric under that variable exchange. Pair the six permutations in this involution. The exact polynomial checker also verifies (2).

Because L_q is symmetric and E0,R are symmetric, average the ordered-cover contrast over six variable permutations. The distinct-slot term vanishes by (2); the common-slot term becomes precisely (1). This is an all-order counting and polynomial argument, not a finite reconstruction inference.

## 2. What would suffice for arbitrary ambient lifting

Suppose for one actual common tournament profile the right side of (1) is nonpositive for every h. Its h=1 value is strictly negative by the companion positive-square theorem. Therefore the entire ordered-cover profile of U+ dominates that of T+, with strict Hamilton coordinate. In every fixed rooted ambient quotient R0, substituting these port-extended modules gives a positive actual query difference: the substitution coefficients are nonnegative and the one-run coefficient is H(R0)>0. The ports use the original module's exterior incidences.

Thus (1) is a precise sufficient structural reduction of this ambient problem to explicit finite profile inequalities. It is not a proof that those inequalities always hold.

## 3. A genuine obstruction to unrestricted sequence positivity

For the artificial finite nonnegative sequence

    q=(0,9,0,1,1,1,0),

one has

    L_q(K E0/R^3)=219,
    (2/3)L_q(K E0/R^3)=146>0.

Hence the h=2 coordinate in (1) has the opposite sign from its h=1 coordinate. This sequence is not claimed to be a tournament profile. Its zero q_2 with positive higher entries already excludes that interpretation. It nevertheless refutes an attempted proof that applies to every finite nonnegative sequence, as the scalar positive-square theorem does. No exact response collision or actual ambient cancellation follows from this artificial example.

## 4. Finite controls and scope

There are fifteen distinct ordered-cover profiles among all 1,099 labelled tournaments of orders one through five. The exact formula has the domination sign for every h up to the full port-core order for each of those fifteen profiles. Literal subset-cover calculations for TT1, TT2, TT3 and C3 children independently match all coordinates. For TT2 children, q_h(U+)-q_h(T+) at h=0,...,8 is

    (0,36,168,288,192,0,0,0,0).

For C3 children it is

    (0,2052,23220,105300,254880,354240,272160,90720,0,0,0,0).

These are exact finite checks. No statement about all tournament profiles or arbitrary strongly connected ambient quotients is inferred from them. Full named response coefficients are preserved in the companion packet; this reduction neither drops cubic/quartic information nor replaces named response equality by scalar equality.

Run `python verify.py > receipt.json` and `python -O verify.py > receipt.optimized.json`. The sibling `covered_response` and `response` packets must be present. The verifier uses only the standard library and explicit exceptions.
