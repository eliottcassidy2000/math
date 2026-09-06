# The complete mixed (2,2,1) Smith law and its intrinsic dyadic bit

**Status: PROVED in the stated scope; primary exact verification PASS;
independently proof-audited and verified by the root agent.** This classifies every Smith factor for
three arbitrary primitive projective directions with multiplicities `(2,2,1)`.
There is exactly one shape where the metric needs an extra bit. The result
has no theorem ID and makes no external novelty or priority claim.

## Inheritance, hostile first, and the concept board

The closest proved mechanism is
[THM-4443 / arbitrary-jet precision and dyadic unit boundary](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md):
reciprocal Hermite coefficients give the largest Smith factor. Its explicit
`(2,2,1)` hostile supplies two equal-metric observers with different precision,
but does not classify all primes, shapes, depths, or intermediate factors.
The new
[primitive-direction transport](creative_20260906_smith_bridge.md)
provides the lawful projective change of coefficients and local jets.

The canonical hostile is the labelled pair with doubled nodes `0,8` and
simple node `4` versus `-4`, whose full dyadic lists are
`(0,0,4,7,9)` and `(0,0,4,8,8)`. The corrected near miss is MISTAKE-547:
uniform two-jet metric precision does not survive deleting one derivative.
The least-used sidecar is the **reference dependence of a projective
cross-ratio before its active precision scale**.

The five live concepts are the remaining derivative multiplicities; the
three-by-three residual observer; its first determinantal ideal; reciprocal
numerator cancellation; and projective bracket residues. The first hostile
reproduced the inherited unequal lists. A second hostile found that the
obvious bit can change under a lawful projective chart at shallow depth,
although the Smith partition does not. This led to the intrinsic threshold
proved below, rather than treating an arbitrary fourth direction as canonical.

The higher-jet Deuring and p31 intermediate-ideal results were recovered.
The [p31 packet](overnight13_20260906_jets_p31_intermediate.md) and
[adjacent-factor continuation](overnight14_20260906_jets_p31_adjacent_even.md)
show why determinant mass or one largest factor does not generally determine
the full observer or full kernel. Here the residual rank is exactly three,
so its first ideal, determinant, and largest factor do determine every factor.
Targeted current-canon/results searches found the hostile examples and general
inverse formula, but no complete `(2,2,1)` classification. The p31 partition
is not being claimed as a consequence of the present smaller observer.

## 1. The observer and complete classification

Let `v_1,v_2,v_3` be primitive integer vectors with pairwise nonzero brackets
`Delta_ij=det(v_i,v_j)`. The first two directions are doubled and the third
is simple. Choose integer tangents `w_i` with `det(v_i,w_i)=1` and observe
homogeneous binary forms of degree four by

```text
E(F) = (F(v_1), [T]F(v_1+T w_1),
        F(v_2), [T]F(v_2+T w_2), F(v_3)).
```

This is an integral square map of rank five. Its complete Smith partition
is determined prime by prime as follows. Fix a prime `p` and put

```text
A=v_p(Delta_12), B=v_p(Delta_13), C=v_p(Delta_23),
e=min(A,B,C),    f=max(A,B,C).
```

The minimum depth occurs at least twice. The following lists are in increasing
order and give all five `p`-Smith exponents.

**Case I: the doubled pair is not uniquely closest (`A=e`).** Put
`delta=min(2e,e+[p=2])`. Then

```text
(0,0, delta, 4e-delta, 2e+2f).                      (1)
```

This includes every equilateral case, including three distinct dyadic
projective residue classes with `e=f=0`.

**Case II: the doubled pair is uniquely closest (`A=f>e`) and `p` is odd.**
Put `delta=min(f,2e)`. Then

```text
(0,0, delta, f+3e-delta, 3f+e).                     (2)
```

**Case III: the doubled pair is uniquely closest at `p=2`.** Write `f=e+d`,
where `d>=1`. If `d>=2`, put `delta=min(2e,e+d+1)`; then

```text
(0,0, delta, 4e+d+1-delta, 4e+3d-1).                (3)
```

If `d=1`, the two shallow boundaries are

```text
e=0: (0,0,0,2,2),
e=1: (0,0,2,5,5).                                 (4)
```

For `d=1,e>=2`, choose any primitive reference direction `w` satisfying
`det(v_i,w)` odd for all three `i`, and form the half cross ratio

```text
tau_w = Delta_12 det(v_3,w)
        / [2 Delta_13 det(v_2,w)] in Z_2^*,
epsilon = [tau_w = 1 modulo 4].                    (5)
```

Such a reference exists; a unimodular tangent to `v_1` is one choice. The bit
is independent of the reference, of primitive representative units, of the
order of the two doubled points, and of projective `GL_2(Z_2)` changes.
The complete exponent list is

```text
(0,0, e+2, 3e+2-epsilon, 4e+epsilon).               (6)
```

Therefore the complete partition is metric-only at every odd prime, and at
two outside precisely the shape in `(5)`. Within that shape each of the two
bits occurs for every `e>=2`, and their full lists differ. One additional
binary coordinate is both sufficient and necessary for the partition in
this active regime. These claims concern the weighted, multiplicity-labelled
metric; forgetting which direction is simple is not permitted.

The global integer invariant factors are obtained by multiplying the prime
powers in each of the five positions. Only primes dividing a pair bracket
can occur. In every case the last exponent is the sharp extra observation
precision for recovering all five coefficients modulo `p^N`, for every
`N>=1`.

## 2. Reduction to three elementary invariants

First work at affine nodes `0,a,b`, with `0,a` doubled and `b` simple, on
polynomials of degree below five. Clearing the value and first derivative
at zero gives two unit factors and the residual matrix

```text
R = [ a^2    a^3     a^4 ]
    [ 2a     3a^2    4a^3]
    [ b^2    b^3     b^4 ].                          (7)
```

The gcd of all residual entries is exactly

```text
gcd(a^2,2a,b^2).
```

Indeed these three entries are present and divide every other entry after
taking their common gcd. Write `A=v_p(a)`, `B=v_p(b)`, `C=v_p(b-a)`. Thus
the smallest residual exponent is

```text
rho=min(2A,A+v_p(2),2B).                            (8)
```

The confluent determinant gives total exponent

```text
T=4A+2B+2C.                                        (9)
```

Specializing THM-4443's exact inverse formula to multiplicities `(2,2,1)`
gives the sharp largest exponent

```text
L=max(2A+B, 2A+C, 2B+2C,
      3A+2B-v_p(a+2b),
      3A+2C-v_p(3a-2b)),                            (10)
```

where a zero numerator contributes minus infinity. For example the two
reciprocal order-one coefficients have respective numerators `a+2b` and
`3a-2b` over denominators `a^3 b^2` and `a^3(a-b)^2`, up to sign. The
simple direction contributes `1/[b^2(b-a)^2]`.

There are exactly three residual factors. Their first exponent is `rho`,
their last is `L`, and their sum is `T`. Hence the **entire** list is

```text
(0,0,rho,T-rho-L,L).                                (11)
```

No assumption that a determinant alone determines Smith factors enters this
step. It uses the exact residual rank and the separately attained largest
inverse denominator. Consecutive-factor ordering is inherited from the
actual Smith module; it can also be read directly from `(1)--(6)`.

## 3. Exhaustive valuation analysis

If `A=e`, at least one of `B,C` is `e`, and the other is `f>=e`.
The simple-point candidate in `(10)` is `2e+2f`. Every other candidate
is at most this value: both linear numerators have valuation at least `e`,
and the order-zero candidates are at most `3e+f<=2e+2f`.
Thus `L=2e+2f`. Formula `(8)` becomes `rho=min(2e,e+v_p(2))`, proving `(1)`.
This argument includes the equilateral case.

Now suppose `A=f>e`, so `B=C=e`. At an odd prime, both `a+2b` and
`3a-2b` have valuation exactly `e`, since the `2b` term has strictly
smaller valuation than the other term. Thus `L=3f+e` and
`rho=min(f,2e)`, giving `(2)`. This remains true at `p=3`.

At two write

```text
a=2^(e+d)u,       b=2^e v,       u,v odd.
```

If `d>=2`, both linear numerators have valuation `e+1`. Their candidates
are `4e+3d-1`, dominating both `3e+2d` and `4e`; `(8)` gives
`rho=min(2e,e+d+1)`. This proves `(3)`.

If `d=1`, the two numerator valuations are

```text
e+1+v_2(u+v),       e+1+v_2(3u-v).
```

Their normalized sum is `4u`. Because `u,v` are odd,

```text
min(v_2(u+v),v_2(3u-v))
 =1  if u/v=1 mod4,
 =2  if u/v=3 mod4.                                 (12)
```

The second line is exact: both terms are divisible by four and their sum
has valuation two, so they cannot both have larger valuation. Individual
terms can vanish or cancel more deeply; the common minimum remains fixed.
Write `epsilon=[u/v=1 mod4]`. Then `(10)` becomes

```text
L=max(3e+2,4e+epsilon),
rho=min(2e,e+2),        T=8e+4.                     (13)
```

At `e=0,1`, the maximum masks the bit, giving `(4)`. At `e>=2`,
`L=4e+epsilon` and `rho=e+2`, giving `(6)`.
This proves all affine branches, including signed nodes and arbitrary unit
lifts. The analysis is valuation-exact, not a finite-depth extrapolation.

## 4. Why the extra bit is an intrinsic projective coordinate

For any two unit-separated references `w,w'`, the Plucker determinant
identity gives

```text
det(v_3,w)det(v_2,w')-det(v_2,w)det(v_3,w')
 = +/- det(v_2,v_3)det(w,w').                        (14)
```

All four reference brackets are units. Their cross-ratio quotient is therefore
`1 modulo 2^e`, because `v_2(Delta_23)=e`. The half cross ratios in `(5)`
are themselves units, so

```text
tau_w = tau_(w') modulo 2^e.                        (15)
```

At `e>=2` this proves reference independence modulo four. Scaling any
primitive local representative by a unit cancels exactly in `(5)`, and the
four determinant factors make projective `GL_2(Z_2)` covariance exact.
Swapping the doubled endpoints sends `tau` to `tau/(2tau-1)`, which is
congruent to `tau` modulo four because `tau` is odd. Thus epsilon is intrinsic
to the declared weighted projective configuration.

In an affine chart with the reference at infinity, translate the first doubled
node to zero. Then `tau=a/(2b)=u/v`, so `(5)` is exactly the bit used in `(12)`.
The new coordinate is not an arbitrary affine numerator disguised as a
projective invariant; equation `(15)` proves the required chart independence.

The depth condition cannot simply be erased. At affine points `0,4,2`,
with the first two doubled, `tau=1`. The integral determinant-one projective
change `x -> x/(1+x)` gives points `0,4/5,2/3`, where the reference at the
new infinity gives `tau=3/5=3 mod4`. The common depth is `e=1`, and the
full Smith list is `(0,0,2,5,5)` in both charts. The first failed implication
would be to treat the raw modulo-four coordinate as intrinsic before it is
observable. The repaired theorem uses it only for `e>=2`, exactly where
`(13)` ceases to mask it.

## 5. Projective transfer and completeness of the cases

The projective observer has complete local jet multiplicities `(2,2,1)` on
homogeneous degree four. The coefficient and target transformations in the
[primitive-direction transport theorem](creative_20260906_smith_bridge.md)
apply: substitution by `GL_2(Z_p)` is unimodular on source forms, unit changes
of representatives and local parameters give invertible triangular jet blocks,
and an integer tangent shift adds `4k` times the value to the derivative.
Thus the full `p`-Smith module is independent of those choices.

For three projective directions, a common affine chart with unit denominators
exists unless all three residue classes of `P^1(F_2)` occur. Indeed choose
the chart pole in an unoccupied residue class; there are `p+1` available
classes and at most three occupied classes. In the exceptional three-class
case every pair bracket is a unit and the full observer is unimodular.
This follows either by unit-separated interpolation or by the weighted
determinant formula

```text
|det E|=|Delta_12|^4 |Delta_13|^2 |Delta_23|^2.       (16)
```

For completeness, `(16)` follows on an affine chart by ordinary confluent
Vandermonde. Passing from `(x_i,1)` to a primitive representative with second
coordinate `b_i` contributes a value/derivative block determinant
`-b_i^6` at each doubled direction and `b_3^4` at the simple direction.
These powers cancel the corresponding denominators of the bracket products.
An integer shear can make all second coordinates nonzero, as in the previous
transport proof, so directions at infinity are included without assuming
unit denominators for this rational identity.

In a common integral chart, normalize the three directions to `(x_i,1)` and
translate `x_1` to zero. The affine proof applies over `Z_p`; it used only
integral matrix operations and exact valuations, so it extends directly
from integers, or by sufficiently close integer approximation above the
known determinant valuation. Bracket valuations become difference valuations.
In the sole unit-sensitive case, Section 4 identifies the coordinate
independently of the chosen chart. The exceptional three-class observer
has all exponents zero, already included in `(1)` with `e=f=0`.
This completes the classification for all primitive integer directions.

## 6. Actual finite-precision distinction and stopping boundary

Take doubled nodes `0,2^(e+1)` and simple node `+2^e` or `-2^e`, for any
`e>=2`. They have identical weighted dyadic metrics, and their bits are
respectively one and zero. Let `K_epsilon(N)` denote the kernel of the
corresponding full observer modulo `2^N`. Because all five invariant factors
are now classified, Smith coordinates give the exact consequence

```text
|K_0(N)| / |K_1(N)|
 =2  if 3e+2 <= N <= 4e,
 =1  otherwise,                         N>=1.       (17)
```

Indeed `log_2 |K_epsilon(N)|` is the sum of the five capped exponents
`min(N,alpha_i)`. Only the last two exponents differ, their sum is fixed,
and direct comparison gives `(17)`. There are exactly `e-1` affected integer
precisions, an unbounded interval of distinguishability produced by one bit.
Unlike a claim from an isolated pair of unclassified p31 factors, `(17)` is
about the entire kernel because the remaining factors are known.

The result advances an explicit previously unclassified mixed observer. It
does not infer a general higher-jet projective partition, a p31 classification,
or an LRC consequence. A next operation is to replace the simple node by a
third doubled node: that returns the already proved metric-only two-jet
three-node theorem. Replacing it by multiplicity three is a different open
classification; the two numerator cancellation budget in `(12)` no longer
controls every reciprocal jet or intermediate ideal.

## 7. Reproduction

Run from the repository root:

```text
python3 -B 04-computation/second_20260906_smith.py
python3 -B -O 04-computation/second_20260906_smith.py
```

The [source](../../04-computation/second_20260906_smith.py) and
[matching output](second_20260906_smith.out) freeze the exact universe:
bounded signed affine nodes, every labelled closest-pair placement, independent
outer and inner depths, arbitrary selected unit lifts, projective transforms,
reference changes, a dyadic three-residue-class control, both inherited twins,
and the new shallow coordinate hostile. The source compares complete rational
DVR Smith lists with the closed formulas and a separate reciprocal precision
calculation. In small cases it reconstructs every determinantal divisor from
all minors and also inverts the literal matrix. It checks the actual full
kernel-size ratio, not just the changed factor pair.

Normal and optimized outputs agree byte for byte. There are 8,732 complete
Smith/metric rows, 648 projective/reference symmetry rows, 1,506 independent
all-minor determinants, and 215,537 exact gates. The source uses no floating
point or Python assertion statements. The analytic proof supplies the
universal quantifiers. The root's independent
[SymPy verifier](../../04-computation/second_20260906_root_audit.py) and
[matching output](second_20260906_root_audit.out) reconstruct 1,059 literal
matrices, compare 4,236 prime partitions, and separately check kernel windows
through depth 15, without importing this producer.

```text
source_sha256:
e8a5f1e057e4d3c1261d0a4750a9e3291ba8bf8de6b472d41454a72074880b20
output_sha256:
7e703fe507560ee14d6a928fc6158e0bffe1a9735ba83339cb38c9fc5d3ba321
semantic_sha256:
bdc00c74a3232590f3be4031adac63126a662f1ff427bce8446eda77429a16a7
```
