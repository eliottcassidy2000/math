# A finite central boundary and an arithmetic obstruction to resolving its points

2026-10-26. **PROVED** by the arguments below, with **FINITE-EXACT** controls
in `check.py`. The Gaussian norm classification is classical; its elementary
proof is included. This is an arithmetic refinement and correction of
[THM-4773](../../../../01-canon/theorems/THM-4773-fixed-twist-rational-lifts-on-the-hourglass-boundary-central-conic-outer-genus-five.md),
not a new claim about the generic genus-five boundary or integer hourglasses.

## Inheritance and the failed implication

The inherited [fixed-twist theorem](../../fixed_twists_conductors_reverse_potential_20261025/twists/FIXED_TWIST_LIFTS.md)
correctly classifies the **open** central line by three opposite-pair
squareclasses and a conic. Its complete-central-line "none or infinitely
many" sentence overlooks the three excluded triple points. Indeed its own
example `d=(1,1,1,1,-1,1,-1)` has rational lifts at `(0:0:1)`.
The closest proved mechanism is its pointwise squareclass test; the hostile
is a rational node with conjugate branches. The corrected near miss is
extending an open-stratum dichotomy to its closure. The least-used sidecar
is the coefficient of the transverse central root, invisible on that closure.

Anchor: classify fixed-twist central points, including the missing stratum.
Niche: arithmetic of exceptional conics. Wildcard: whether changing only the
invisible transverse coefficient changes rational resolution lifts.
The four live concepts are the central fibre, its component normalization,
the normal surface, and the exceptional conic. The computation below keeps
all four separate; the wildcard gives an exact change invisible to the first.

## C1. The exact trichotomy on the finite normal surface

Let `k` be any infinite field of characteristic different from two and
`d_i in k*`. Use the seven branch forms

    L=(E,E+x,E-x,E+z,E-z,E+x+z,E-x-z),   s_i^2=d_i L_i.

The equations define the finite normal root surface `X_d`; points are
projective root coordinates. Put

    q1=-d2/d1, q2=-d4/d3, q3=-d6/d5.

Square means a nonzero square in `k`. On the open central line
`E=0, x z (x+z)!=0`, rational lifts exist exactly when all three `qj`
are squares and the smooth conic

    A^2/d1 + B^2/d3 = C^2/d5                         (C1)

has a `k`-point. If so, there are infinitely many central points. If not,
the only possible central points lie above the three triples, with tests

| Base point | Necessary and sufficient square tests | Number of lifts |
|---|---|---:|
| `(0:0:1)`, roots `012` vanish | `q2,q3,d5/d3` | 8 |
| `(0:1:0)`, roots `034` vanish | `q1,q3,d5/d1` | 8 |
| `(0:1:-1)`, roots `056` vanish | `q1,q2,-d3/d1` | 8 |

Thus the complete central divisor has **zero, exactly eight, or infinitely
many `k`-points**. The finite nonempty case occurs iff exactly one `qj`
is nonsquare and the triple opposite that missing pair passes its table test.
This includes every number field, including fields where `-1` is square.

**Proof.** At a generic central point each opposite root pair is nonzero,
so its ratio squares to `qj`. If the three ratios are available, choose their
signs and put `(A,B,C)=(s1,s3,s5)`; the remaining equation is exactly (C1).
It is smooth since its three diagonal coefficients are nonzero. A rational
point parametrizes the conic by lines, giving infinitely many points over
an infinite field; exclude its finitely many points with `ABC=0` to recover
the open central line. Conversely these formulas construct every root point.
At `x=0`, normalize `s3=1`: the three other nonzero coordinate squares are
`q2,d5/d3,q3*d5/d3`. This proves the first test and the count `2^3=8`.
The other rows are identical coordinate calculations. No point has two of
`x,z,x+z` zero. If two distinct table rows pass, all `qj` are squares and
one of those points supplies a point of (C1), hence the infinite case.
If all `qj` are squares and even one table row passes, the same argument
applies. Fewer than two square `qj` allow no table row. This proves the
classification and both directions of the finite-case criterion. QED.

In the minimal witness `d=(1,1,1,1,-1,1,-1)`, exactly `q1=-1` fails
over `Q`. The eight points are

    [0:0:0:1:eps4:eps5:eps6],   eps4,eps5,eps6 in {+1,-1}.  (C2)

They lie above `(0:0:1)`. There are no other rational central points.
The strongest survivor of the old dichotomy is the open-central statement;
the first invalid implication was that a rational point of the singular
closed divisor must lift to its component normalization.

## C2. A coefficient invisible on the central fibre controls resolution lifts

Over `Q`, keep the one-parameter family

    d_delta=(2/delta,1,1,1,-1,1,-1),   delta in Q*.

Its entire central scheme is independent of `delta`: setting `s0=0`
removes the only changed coefficient. C1 gives the same eight rational
points (C2) for every `delta`. All are ordinary double points of `X_delta`.
On the minimal resolution `pi: Xtilde_delta -> X_delta`, the exceptional
curve above **each** point (C2) is the conic

    Q_delta: a^2+b^2=delta c^2 in P^2.                (C3)

Consequently the inverse image of the central divisor has rational points
iff `Q_delta(Q)` is nonempty; in that case it has infinitely many, all on
the eight exceptional conics. Otherwise it has none. In particular:

* `delta=2`: `[a:b:c]=[1:1:1]` splits each exceptional conic;
* `delta=3`: no exceptional conic has a rational point;
* both finite normal surfaces still have exactly the same eight rational
  central points and the same rational central scheme.

**Local proof.** Work at any (C2) in chart `s3=1`, with
`c=s0,a=s1,b=s2`. The opposite-pair equation is exactly

    a^2+b^2=2c^2/d0=delta c^2.

The other equations become

    s4^2=1-delta c^2,
    s5^2=1+a^2-delta c^2/2,
    s6^2=1+a^2-3delta c^2/2.

At every chosen sign triple their derivatives in `s4,s5,s6` are nonzero.
Thus the surface is etale-locally the displayed nondegenerate quadratic
cone; its completed local ring is
`Q[[a,b,c]]/(a^2+b^2-delta c^2)`. Blowing up its vertex resolves the node
and its exceptional curve is its projectivized tangent cone, (C3).
The same formulas prove this for all eight sign triples, without assuming
that a geometrically rational exceptional curve is rational over `Q`.

All rational points of the inverse-image divisor map to one of (C2).
Above each node the fibre is (C3), proving the resolution assertion.
Equivalently the central divisor itself has local branches `a=+i b` and
`a=-i b`; its normalization has no rational point above any of the eight
nodes. Its eight geometric conic components are paired by conjugation,
and the strict transforms supply no rational point. The rational points
created when (C3) splits lie on exceptional curves, not on a hidden
rational component of the original central divisor. QED.

This statement distinguishes three operations: `X_delta` is already the
normalization of the complementary quotient product; normalizing the
divisor inside it separates its branches; resolving `X_delta` inserts the
tangent conic. They have different rational-point behavior.

## C3. Exact number-theoretic decision and maps

The exceptional conic is soluble precisely when

    delta > 0 and v_p(delta) is even for every prime p=3 mod 4.  (C4)

This is the classical rational Gaussian-norm criterion, not a new prime
classification. A point of (C3) has `c!=0` over `Q`, and maps to
`alpha=(a+ib)/c in Q(i)*` with norm `delta`. Conversely a norm witness
`alpha=u+iv` gives `[u:v:1]`. This preserves the actual twist, not merely
its geometric conic class.

Here is a self-contained algebraic proof of the criterion. The Gaussian
integers are Euclidean for `N(a+ib)=a^2+b^2`: round both coordinates of
a quotient to obtain remainder norm at most one half of the divisor norm.
Hence they are a UFD. For an odd prime `p=3 mod4`, `-1` is nonsquare modulo
`p` by Fermat's little theorem, so divisibility of `A^2+B^2` by `p` forces
both `A,B` divisible by `p`. Removing squares shows its valuation is even.
Clearing rational denominators preserves this parity; positivity is immediate.
For `p=1 mod4`, Wilson's theorem gives a square root of `-1 mod p`:
`((p-1)/2)!^2=-1 mod p`. Wilson follows by pairing each non-self-inverse
unit with its inverse. If `p` were Gaussian-prime it would divide one
factor of `(r+i)(r-i)`, impossible because the imaginary coefficient is
`+/-1`. UFD therefore factors `p` nontrivially; the factors' integer norms
multiply to `p^2` and each exceeds one, so both norms are `p`. Thus `p`
is a sum of two integer squares. Also `2=N(1+i)` and each even power of
a `3 mod4` prime is a norm of an integer. Multiplication and division of
these norms builds every positive rational satisfying (C4). QED.

For `delta=3` there is also a direct hostile proof: a primitive integer
solution of `a^2+b^2=3c^2` has `3|a,b` by reduction modulo three, hence
`3|c`, contradicting primitivity. The positive control `delta=2` and this
negative control have identical central scheme and rational node set.

For an explicit parametrization, from `u^2+v^2=delta` multiply the norm
witness by `(1+it)/(1-it)` to obtain

    a(t)=[u(1-t^2)-2vt]/(1+t^2),
    b(t)=[v(1-t^2)+2ut]/(1+t^2), c(t)=1.             (C5)

This parametrizes the conic with the omitted limiting value at `t=infinity`.
An actual tangent-direction point is therefore available for every soluble
twist, not just a numerical solubility predicate.

There is a precise local F1 consequence. A soluble `Q_delta` is
`P^1_Q`, hence has the split-torus decomposition `G_m` plus two rational
points. An anisotropic `Q_delta` has no rational point and cannot have a
finite decomposition into nonempty split `Q`-tori: every such torus has a
rational identity point. Thus (C4) also decides split torifiability of
these individual exceptional conics. This is not a torification claim
for the whole resolved surface, nor for its finite central fibre.

Primary check: Keith Conrad, [The Gaussian Integers](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/Zinotes.pdf),
Theorems 9.6--9.10, read 2026-10-26, gives the same classical norm argument.
The proof above and the map to the hourglass exceptional fibre are explicit.

## Controls, connection contract, and stopping boundary

`check.py` uses exact rational arithmetic and imports no producer code.
It exhausts all 4,096 outer-twist tuples from `{-2,-1,1,2}^6`; checks the
three special fibres directly; for every tuple lacking a square pair class
checks all primitive central base points of height at most five; verifies
all eight displayed nodes and the local equations for representative deltas;
and compares the norm criterion with direct two-square enumeration for
positive integers through 500 and a rational denominator grid. It also
checks the `delta=3` obstruction by primitive residues modulo nine and the
positive parametrization. Explicit exception gates survive `python -O`.

Source -> target: a rational node of the fixed-twist root surface -> its
exceptional conic and Gaussian norm torsor. Map: projectivized tangent cone,
then `[a:b:c] -> (a+ib)/c`. Preserved: the coefficient `delta` and rational
lifting under the resolution. Lost by the central scheme: that transverse
coefficient. Required sidecar: the second-order local equation. Decisive
test: `delta=2` versus `delta=3`.

The concept-board return is exact: the constant central fibre does not
determine the resolution's rational fibre; its divisor normalization also
does not model the resolution. The relevant meta-pattern is **classify the
target over the exceptional stratum before closing by density**. No new
meta-pattern card is proposed from one thread. No claim is made about
interior rational points, integral points, all outer-line twists, broader
F1 frameworks, or a torification of the whole resolved surface. A natural
continuation is to keep these tangent
conics as arithmetic sidecars in any resolved norm-seven fibre-product atlas.

Reproduction from repository root:

```sh
python3 -B 05-knowledge/results/normseven_arithmetic_f1_20261026/boundary/check.py
python3 -B -O 05-knowledge/results/normseven_arithmetic_f1_20261026/boundary/check.py
```

The JSON outputs must agree with `normal.json` and `optimized.json`.
They are byte-identical, SHA-256
`b9c212deb975c73d5ea2a97a8242f76e6db90dafeba81ec81c4e402b39854cf4`.
The receipt records 184,940 explicit gates: 144 finite-eight twist cases,
3,888 empty cases forced by a missing pair class, and 64 all-square-pair
cases deliberately left to the conic criterion rather than bounded search.
