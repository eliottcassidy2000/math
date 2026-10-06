# Controlled unbounded Poisson witnesses for Collatz source floors

2026-10-05. **PROVED:** two explicit weighted contraction estimates, a
Poisson dual theorem for controlled unbounded functions, its residual bill,
and finite path extraction in the normalized weight. **FINITE-EXACT:** the
declared row and boundary controls. **OPEN:** a positive independently proved
dual for every source. Enlarging the valid witness space does not prove
all-source positivity or a new convergence basin.

[Script](../../04-computation/experiments/collatz_weighted_poisson_growth_20261005.py)
and [output](collatz_weighted_poisson_growth_20261005.out).

## 1. Inheritance and the precise extension

The closest mechanism is the [bounded Poisson source dual,
P1-P2](collatz_poisson_source_dual_20261005.md): a globally proved adjoint
inequality with positive ROOT evaluation supplies a source floor. Its
canonical hostile is the inverse ray with weights1/2 and harmonic function
phi(j)=2^j: the omitted boundary stays1, giving the false readout1<=0 if
unbounded functions are admitted without a tail condition.

The extension below replaces absolute boundedness by boundedness relative to
an explicitly proved growing weight. The least-used coordinate is ordinary
source height *inside the inverse operator*, rather than along one forward
orbit. The anchor remains independent source positivity; the niche is a
weighted operator space; the wildcard is an exact integer-height weight
that avoids irrational evaluation in a certificate checker.

Board: **selected source / inverse sibling phase / growth envelope / weighted
contraction / signed residual / ROOT boundary**. The META-PATTERNS card
**Find the hidden second coordinate in a nearly true theorem** applies:
an unbounded function is admissible only together with a controlled weighted
tail. The card **Controlled forgetting and unlabeled quotients require a
sidecar** applies to the normalization: the weight is retained at each base.

No novelty or literature-priority claim is made for weighted operator duality.
The work here is the explicit Collatz arithmetic and certified constants.

## 2. Operator and shifted height

Use the exact base space and killed-root operator from the companion:

    B={b>0 odd: v2(3b+1) in {1,2}},     S(x)=4x+1,
    U(b)=S^k G(b),
    (Af)(b)=2^(-k-1)f(G(b)) for b!=1,     (Af)(1)=0,
    P=A*,     g=(I-A)^(-1)delta_1.

Every odd source has a unique form S^j b. The inherited unweighted bound is
||P||_infinity=6/7. The nonnegative g has g(1)=1 and has positive mass exactly
on the actual ROOT component. These statements do not assume Collatz.

At a base c, an inverse child at sibling depth k is

    y=S^k c,     b=(2^a y-1)/3,
    a=2 if y=1 mod3,     a=1 if y=2 mod3.

Depths with y=0 mod3 and the child b=1 are omitted. Each other child is a
base with exact parent(c,k). Since S^k c=c+k mod3, exactly one depth class
modulo3 is forbidden, before the extra root-loop deletion. The identity

    3b+1=(2^a/3)(4^k(3c+1)-1)                         (1)

retains the ordinary-height coordinate through this infinite inverse sum.

## 3. A polynomially growing envelope

Put

    V_p(b)=((3b+1)/4)^(1/12),       V_p(1)=1.

**WG1 — PROVED.** On every base,

    P V_p <= kappa_p V_p,       kappa_p=3280/3367<1.    (2)

By (1), V_p(b)/V_p(c)<(4/3)^(1/12)4^(k/12). Set
t=2^(-5/6). Any two retained phases in each period have weighted sum at most
(1+t)/(1-t^3). Thus

    (P V_p)(c)/V_p(c)
       <= (1/2)(4/3)^(1/12)(1+t)/(1-t^3).

Two integer checks certify rational upper bounds:

    3*41^12 >= 4*40^12,       9^6 >= 2^19.

They give (4/3)^(1/12)<=41/40 and t<=9/16. Substitution yields exactly
3280/3367. Root-loop deletion only decreases the sum. This constant is a
convenient certified bound, not an asserted optimum.

The exponent is not unrestricted. At exponent1/4, base25 already has a
finite partial normalized row sum greater than1. This refutes one-step
contraction in that particular weight, not all possible norms or iterates.
At exponent1/2, the infinite row diverges: on admissible phases the individual
weighted terms tend to sqrt(2^a/3)/2>0.

## 4. An exact rational height envelope

Let bitlength be the ordinary integer binary length and put

    V_l(b)=(25+bitlength(3b+1))/28,       V_l(1)=1.

**WG2 — PROVED.** On every base,

    P V_l <= kappa_l V_l,       kappa_l=641/686<1.     (3)

Equation (1) and 2^a/3<2^(a-1) imply

    bitlength(3b+1)<=bitlength(3c+1)+2k+1.

The extra cost over the three possible phase patterns is exactly

    sum_(admissible k) 2^(-k-1)(2k+1)
      in {95/49,106/49,93/49}.

Root deletion can only reduce it. The constant part costs at most6/7,
and 25+bitlength(3c+1)>=28. Consequently

    P V_l/V_l <= 6/7 + 106/(49*28) = 641/686.

This envelope is unbounded but has exact rational values at every source.
The integer formula helps a symbolic or interval checker retain its growth
budget without requiring algebraic-number arithmetic.

## 5. Weighted duality and residual payment

For either envelope V, write kappa for its certified contraction. Define

    ||f||_(1,V)=sum_b V(b)|f(b)|,
    ||phi||_(infinity,V)=sup_b |phi(b)|/V(b).

**WG3 — PROVED.** A is a contraction of norm at most kappa on the weighted
ell1 space, and P is a contraction on the displayed weighted sup space.
For the nonnegative Green vector,

    sum_b V(b)g(b) <= 1/(1-kappa).                    (4)

Indeed the column estimate (PV)(c)<=kappa V(c) is exactly the weighted ell1
bound. Summing the Neumann series against the unit root value V(1)=1 proves
(4). Its vector agrees with the inherited unweighted g by uniqueness there.

If q and phi have finite weighted sup norms and satisfy the **global**
coordinatewise inequality

    (I-P)phi <= q,

then

    <q,g> >= phi(1).                                 (5)

All pairings are absolutely convergent by (4), so moving P across the pairing
is legitimate. There is also a unique exact Poisson solution in this weighted
space. Its uniqueness, like the bounded-space version, does not force a
positive root value.

For an approximate inequality with a proved eta>=0,

    (I-P)phi <= q+eta V,

the quantitative conclusion is

    <q,g> >= phi(1)-eta/(1-kappa).                    (6)

The two explicit bills are 3367/87 for V_p and686/45 for V_l. A refined
measurable residual can instead be paired directly with g if that pairing
has its own independent upper bound. A numerical grid is not a proof of the
global inequality or the residual bill.

For q=delta_b, a positive right side gives g(b)>0 and an exact source floor.
For n=S^j b, the fixed-price floor is multiplied by2^-j; [Poisson source
dual P3](collatz_poisson_source_dual_20261005.md#5-a-fixed-price-floor-gives-a-logarithmic-counter-deadline)
then gives a finite ROOT deadline, a mixture-weight floor, and the existing
signed localized measurement. The growth envelope itself supplies none of
these positive signs at an arbitrary target.

## 6. Strict extension of the valid witness class

For target3 the actual base edge has G(3)=1,k=1. The finite packet

    psi(3)=1,     psi(1)=1/4,     psi=0 elsewhere

satisfies (I-P)psi=delta_3, including all boundary equations. Now set

    phi=psi-V_l/8.

This phi is unbounded below, so it is outside the preceding bounded theorem.
Nevertheless (I-P)V_l>0 by (3), giving

    (I-P)phi<=delta_3,       phi(1)=1/8.

The weighted theorem rigorously certifies g(3)>=1/8. The actual value is1/4;
this example proves admissibility of a larger witness class, not a new
convergence result. The same construction applies to any independently
checked path packet, and does not pretend to provide an unknown packet.

More generally, if |phi|<=M V, q=delta_b and phi(1)=delta>0, then at a
positive non-target vertex the normalized value phi(c)/V(c) is dominated
by a nonnegative weighted average of child values whose total mass is at
most kappa. Some child therefore has normalized value strictly greater than
a phi(c)/V(c), for any fixed1<a<1/kappa. Choosing such children must hit b
before a^T delta>M. For a computable witness and envelope, strict interval
comparisons find a suitable child by enumeration; no exact maximum is needed.
The proof still yields a finite actual connecting path, with a different
admissible analytic representation of the witness.

The ray hostile explains the boundary. With g_j=2^-j and Pphi(j)=phi(j+1)/2,
V_j=2^j has weighted contraction1, not less than1, and g_j V_j=1 never tends
to zero. The claimed dual readout1<=0 fails. Taking V_j=(3/2)^j instead gives
weighted contraction3/4 and a vanishing boundary. Strict contraction is a
sufficient condition used here; other justified tail conditions are not
excluded.

## 7. Connections and the next precise obligation

| Source | Map and preserved predicate | Lost information and required sidecar | Decisive test |
|---|---|---|---|
| Inverse Collatz edges | Shifted height (1) bounds their weighted column sums | Residue phase, root deletion, and height cannot be discarded | All three forbidden-phase patterns; root1 separately |
| Bounded Poisson witness | Replace absolute norm by phi/V; retain valid adjoint pairing | Absolute boundedness is lost; weighted Green integrability replaces it | Critical ray with nonzero boundary |
| Guarded run-block family | Its known source floor can serve as a proved patch datum | Membership and decreasing source rank remain necessary | New23 family, with residual27 still not certified by that grammar |
| Port elimination | Apply the [patch calculus](collatz_poisson_patch_20261005.md) in a weighted norm | Exact forcing and boundary residual must survive elimination | A zero-target completion or disconnected component stays unpaid |

The next obligation is a concrete source-indexed analytic ansatz satisfying
(5) or a paid version of (6) with positive root value on sources outside the
proved families. The new envelopes make unbounded ansatzes legal; they do
not supply the required global inequality. This is more specific than asking
for an arbitrary positive measure, and it retains a falsifiable residual.

## 8. Reproduction and finite scope

```text
python -B 04-computation/experiments/collatz_weighted_poisson_growth_20261005.py
python -B -O 04-computation/experiments/collatz_weighted_poisson_growth_20261005.py
```

There are16,037 exact checks:192 actual bases below512, sibling cuts24 with
analytic infinite-tail enclosures, integer-certified radical intervals,
all three phase costs, a positive unbounded witness, the exponent1/4 and1/2
boundaries, and24 critical/subcritical ray controls. Python arbitrary-size
integers and Fractions are used; no floating comparisons or ROOT census
enter the row certificate. The universal bounds have the proofs above;
finite controls are not promoted to all-source positivity.
