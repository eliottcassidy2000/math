# The effect cube and the full rooted return observer

**PROVED, scoped; exact matrix certificates and independent actual-source tests.**
Same-agent proof review, not an independent-agent audit.

## JC1. Domain and six coordinates

Start with an actual finite strong tournament, two distinct retained ports,
and a specified root. Apply the inherited two-birth NR seal, retaining the
root and any one of the three available second-port choices. Thereafter use
the actual C3, C5 and NR macros of
[RC2](../rooted_carry_closure_20261103/ROOTED_RETURN.md) and
[RC4](../rooted_carry_closure_20261103/SINGULAR_RETURN.md).
The domain is this sealed-source language, not every arbitrary balanced
profile. Graphs, physical names and native path programs remain sidecars.

At the current root r put

    v = (h,s,e,u,p,d)
      = (H(T), starts_r(1), ends_r(1), starts_r(2), P_T(2), P_(T-r)(2)).

P(2) counts ordered two-bin spanning path covers, allowing empty bins.
In starts_r(2) the first bin is nonempty and starts at r. The balanced
profile is (h,s,s,h). The full golden triple is

    (h,X,Y) = (h,(u-e)/2,(u+e)/2).

There are exact row-action laws v'=v A_j. The matrices are:

```
A3 = [ 61  9  0     0     0     0
       -2  0  0     0     0     0
        0  0  4     0     0     0
        0  0  0  -104  -544  -112
        0  0  0    81   468    86
        0  0  0   100   488   108 ]

A5 = [421 33  0     0     0     0
       -2  0  0     0     0     0
        0  0 16     0     0     0
        0  0  0  -944 -6688  -976
        0  0  0   597  4956   614
        0  0  0   940  6488   972 ]

AN = [  3  1  0     0     0     0
        0  0  0     0     0     0
        0  0  1     0     0     0
        0  0  0    -1    -8     0
        0  0  0     2     7     2
        0  0  0     0     6     0 ].
```

Proof certificate: [effects.py](effects.py) constructs the inherited exact
46-state endpoint-cover operators. Under the polarity-zero births used
here, the number of unnamed terminal endpoints cannot decrease. Discarding
states with more than two therefore cannot affect these queries. The three
NR seal images span a 32-dimensional invariant subspace. On it the four
initial queries h,s,e,u close in dimension six, with the displayed exact
intertwining equations. This is a finite identity proving every word and
every source in the stated domain, not extrapolation from bounded orders.

The six coordinates are linearly minimal over Q. Six actual strong
order-seven sources in [controls.json](controls.json) have determinant
-82944=-2^10*3^4. Moreover, writing u3=u(C3), u5=u(C5),

    U3=u3+104u=81p+100d,
    U5=u5+944u=597p+940d,
    p=(940U3-100U5)/16440,
    d=(81U5-597U3)/16440.

Thus every natural coordinate is a linear combination of future outputs.
Actual source rank gives the matching lower bound. The same six witnesses
have rank six over F7. This is a minimum for a **linear rooted count/golden
observer**, not for arbitrary encodings or the full eighteen-polynomial
two-port attachment language.

## JC2. Exactly eight invariant linear types

All three matrices preserve the blocks (h,s), (e), (u,p,d). Their generated
unital associative algebra is exactly

    M2(Q) direct-sum Q direct-sum M3(Q).

[algebra.py](algebra.py) gives fourteen actual matrix words and rational
coefficients for every block matrix unit. The upper bound fourteen is
visible from the blocks; the constructive units prove equality. Those
fourteen words also have rank fourteen modulo seven.

Every invariant subspace of the six-dimensional row module is consequently
a sum of some of these three simple summands. Indeed the three central
idempotents split any invariant subspace, and the full matrix units make
each nonzero block intersection the whole block. Hence there are exactly
eight invariant subspaces, and exactly eight quotient types, with dimensions
0,1,2,3,3,4,5,6. Their inclusion diagram is a Boolean cube.

This is a cube defined by **observable effects of compositions**. It refines
the earlier [operation-language cube](../two_selector_boolean_20261102/CARRIER.md):
the eight subsets of {C3,C5,NR} are a different cube. For our four initial
queries, the empty operation set has dimension four and every nonempty set
already closes in dimension six. Eight vertices need not mean eight distinct
dimensions. No typed-lambda semantics or octonion multiplication is asserted.

## JC3. Integral closure removes the apparent odd-prime obstructions

Let L be the Z-span of all future query columns A_w q, for
q in {h,s,e,u}. Exact Hermite closure gives

    length <=1: index 40,
    length <=2: L = <h,s,e,u,p,2d>, index 2,
    length <=3: unchanged.

[lattice.py](lattice.py) certifies integer invariance of the final lattice,
so it is the all-word closure. In particular, over **every odd modulus**,
queries through length two recover all six natural residues. The odd factors
3,5,137 in the convenient two-probe denominator 16440 are coordinate-chart
artifacts, removed by mixed queries. The only formal integral defect is 2.
On the formal Z^6 module, reduction modulo an even m leaves the possible
hidden shift d -> d+m/2. We do not assert that every such shift is realized
by actual tournaments, whose coordinate lattice has additional congruences.

## JC4. A literal quarter-turn, not an equality of unrelated algebras

The full matrix orders modulo seven are ord(A3)=12 and ord(A5)=48.
For the three-dimensional cover block B3,

    char(B3) = (x-3)(x^2+2) in F7[x].

The projector onto the scalar line is 2(B3^2+2I); its complement is the
quadratic plane. On that plane gamma^2=-2=5 and gamma^7=-gamma.
Thus J=B3^3 satisfies J^2=-I. The saved two-by-three matrix L satisfies

    L J = [0 1; -1 0] L

and has rank two. It explicitly intertwines an ordinary signed quarter-turn
with a real return action on the cover observer. This transfers the C4
action from a plane in the earlier
[signed lambda/octonion lift](../session143_magic_square_cyclic_bridges_20261103/lambda_octonions/PROOF.md).
It preserves that group action only: multiplication, magic-square equations,
integral sources and typing rules are not transported.

The quadratic eigenvalues in the C5 cover block are the negatives of those
in its pair block: their polynomials are x^2+x+3 and x^2-x+3 respectively.
The field is shared, but the sign and physical coordinate map remain data.
This is why a pair-phase repair need not repair the full golden phase.
