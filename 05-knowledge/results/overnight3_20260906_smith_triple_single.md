# Four ternary nodes need a unit residue beyond the distance tree

**Status: PROVED by a finite symbolic minor certificate + VERIFIED-EXACT +
[INDEPENDENT MATHEMATICAL AND EXACT AUDIT PASS](overnight3_20260906_smith_triple_single_audit.md).** The theorem below covers one
complete unit family at a fixed nested depth pattern. It refutes a general
four-node metric-only Smith law. A separate exact finite reduction proves
that the displayed integer counterexample has minimal diameter. No theorem
ID or external novelty claim is made.

**Concurrent connection, fetched at `d2f64b809`:**
[THM-4435 / four-node metric blindness and universal Hermite precision](../../01-canon/theorems/THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md)
independently proves a **dyadic** four-node metric hostile at all outer
scales. The present **ternary** nested unit-family law, intrinsic residue,
and minimal integer diameter are complementary results. No priority claim
is made for the shared general observation that four-node metrics can fail.

## Inheritance and the first decisive test

The closest proved mechanism is the full weighted-minor certificate in
[the ternary double-pair law](overnight2_20260906_smith_double_pair.md),
which proves metric sufficiency for occupancy `(2,2)` at arbitrary depth.
The next shape `(3,1)` has a nested inner triple. Its metric is fully
parameterized by an outer depth `e`, an inner common depth `d`, and an
extra depth `z>=0` of the closest inner pair: the six pair valuations are
`e,e,e,e+d,e+d,e+d+z`. At `z=0` the three normalized inner residues are
all distinct; at `z>0` there is a unique closest pair.

The now-proved
[THM-4429 / arbitrary-three-node-two-jet-smith-form-and-metric-precision](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md)
shows that every three-node two-jet Smith list depends only on its metric.
It is a legitimate proved dependency here, including in the minimality
reduction. It does **not** imply a four-node result.

The canonical hostile is
[THM-4419 / twojet-prime-wall-precision-and-dyadic-triple-smith-law](../../01-canon/theorems/THM-4419-twojet-prime-wall-precision-and-dyadic-triple-smith-law.md)'s
failure of an unchanged prefix after adjoining a node. The corrected near
miss is
[THM-4010 / confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall](../../01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md)'s
separation of a determinant from all intermediate Smith invariants. The
least-used sidecar is the *leading residue* of a weighted minor at a tie
between two equally valued terms.

The concept board is: nested metric depth; unit directions; divided
derivative rows in characteristic three; competing intermediate ideals;
and finite-precision kernels. Varying unit lifts in the first small depth
boxes found a decisive hostile at `(e,d,z)=(2,1,1)`. The rest of this report
proves the entire unit family at that depth, rather than extrapolating an
all-depth formula from the search.

## 1. The exact family and minimal hostile

Let `a,b` be arbitrary integers prime to three, and observe values and
first Hasse derivatives of polynomials of degree below eight at

```text
x(a,b)=(0,9,27a,81b).
```

The full labeled pairwise valuation matrix, in this order, is always

```text
      0  1  2  3
0:    -  2  3  4
1:    2  -  2  2
2:    3  2  -  3
3:    4  2  3  - .
```

Set `kappa=1` if `a=2 mod 3`, and `kappa=0` if `a=1 mod 3`.
The complete three-adic Smith exponent list is

```text
(0,0,2,6,7,12+kappa,15-kappa,22).                   (1)
```

Thus both `b` and all higher unit digits of `a` disappear, but the first
unit residue of `a` is essential. In particular:

```text
(0,9,27,81): (0,0,2,6,7,12,15,22),
(0,9,54,81): (0,0,2,6,7,13,14,22).                  (2)
```

Their entire labeled metric matrices agree, not just the multiset of
distances or the determinant valuation. Both determinant valuations are
64. Their largest exponents are both 22, so their sharp uniform coefficient
precision losses also agree.

THM-4435's independently incoming Hermite formula explains why the largest
exponent survives here. At either closest-pair node, the sum of the three
distance valuations is `S=9` and its unique closest depth is four, so the
local largest-factor contribution is `2S+4=22`. The inner singleton has
`S=8`, closest depth three, hence contribution at most 19; the outer
singleton has `S=6`, closest depth two, hence at most 14. Thus the preserved
worst precision is structural while an intermediate determinantal ideal
still detects `chi`. This use of the incoming theorem is explanatory; the
complete-minor proof of (1) below is independent of it.

Nevertheless their finite-precision kernels differ. Modulo `3^13`, the
kernel orders are `3^53` and `3^54`; modulo `3^14` they are `3^55` and
`3^56`. This follows from the exact Smith-coordinate formula
`|ker J mod 3^N|=3^(sum_i min(N,alpha_i))`. Thus even keeping both the
determinant and the sharp precision loss would miss a real observer
distinction.

Every corresponding three-node subobserver has the same spectrum in the
two examples, by THM-4429. Those are observers on their own degree-below-six
modules; their agreement cannot reconstruct the four-node integral lattice
on degree below eight.

## 2. The smallest residue that repairs this family

The metric tree uniquely identifies an outer singleton, a closest inner
pair, and the remaining inner singleton. If `y` is either point of the
closest pair, define

```text
chi = (x_inner_single-y)/(3*(x_outer-y)) mod 3.      (3)
```

The fraction is a three-adic unit, and (3) is independent of the choice
between the two closest points. In the displayed coordinates it is `a`;
choosing `y=81b` instead of zero gives
`(a-3b)/(1-9b)=a mod 3`. It is also unchanged by affine changes of
coordinate. Thus `chi` is an intrinsic residue once the roles in this
metric tree are identified. It is lost when only valuations are retained.

For this family, the repaired invariant is exactly the distance tree plus
`chi in F_3^*`. The spectrum takes its exceptional branch precisely at
`chi=-1`. This is a sufficient and necessary sidecar for the two possible
spectra here, not a proposed complete invariant for every nested cluster.

## 3. Full weighted-minor proof

Clear the value and derivative at zero by integral row operations and set
`A=3a`, `C=9b`. The residual six-by-six matrix has columns of degrees
`2,...,7`, and rows

```text
0: V(1), 1: D(1), 2: V(A), 3: D(A), 4: V(C), 5: D(C),
V(u)_q=u^q,                  D(u)_q=q*u^(q-1).
```

For row set `R` and degree set `Q`, let `w=sum(Q)-number_of_derivative_rows`.
Its actual minor is exactly `3^(2w) P(A,C)`, where `P` is the corresponding
residual determinant. A polynomial term `c*A^i*C^j` consequently has
valuation at least

```text
2w+i+2j+v_3(c).                                    (4)
```

There are 923 nonempty residual minors. Their complete symbolic expansions
contain 5,587 nonzero coefficients. The following finite polynomial lemma
gives the minimum weights and every normalized residue that survives modulo
three. All minors not listed in a row have zero normalized residue there.
Degrees in this table are actual monomial degrees, rather than array indices.

| Rank | Coefficient floor L | Rows | Degrees | Actual minor / 3^L mod 3 |
|---|---:|---|---|---|
| 1 | 2 | `1` | `2` | `2` |
| 2 | 8 | `0,1` | `2,3` | `1` |
| 2 | 8 | `1,3` | `2,3` | `a` |
| 3 | 15 | `0,1,3` | `2,3,4` | `2a` |
| 4 | 27 | `0,1,3,5` | `2,3,4,5` | `a^2 b(1+a)` |
| 5 | 42 | `0,1,2,3,5` | `2,3,4,5,6` | `2a^6 b` |
| 6 | 64 | all six | `2,3,4,5,6,7` | `a^8 b^4` |

The accompanying self-contained certificate verifies the lemma exactly.
It constructs every polynomial twice, first by cached Laplace expansion
with integer polynomial arithmetic, then by a separate literal permutation
expansion. In the second route, a permutation term is simply its signed
product of derivative degrees times `A^i C^j`; it does not use the first
route's polynomial matrix or determinant recursion. Every coefficient
dictionary must agree.

For each coefficient of every rank-`k` minor, the certificate checks
`2w+i+2j+v_3(c)>=L_k`. Terms at equality are reduced modulo three, giving
**exactly** the displayed residue table. These are finite integer polynomial
identities and divisibilities, so they prove the assertion for all units
`a,b`, including unbounded integer lifts and units in `Z_(3)`. No bounded
numerical substitution is used to establish the universal implication.

The residue table proves that the determinantal valuations at ranks
`1,2,3,5,6` are exactly `2,8,15,42,64`. At rank four, it proves valuation
27 when `a=1 mod 3`, and a lower bound of 28 when `a=2 mod 3`. To attain
that latter bound, take rows `0,1,2,3`, columns `2,3,4,5`. Its exact minor is

```text
3^24 * A^4(A-1)^4,
```

of valuation 28 for every unit `a`. Hence the residual determinantal
valuations are precisely

```text
(2,8,15,27+kappa,42,64).
```

Their successive differences, preceded by the two cleared unit factors,
prove (1). The factors are already ordered, and every minimum has an actual
attaining minor. In particular no inference from the determinant alone or
from a continued Smith prefix appears in the proof.

## 4. Why the nested residue becomes visible

The decisive rank-four residual minor factors, up to sign, as

```text
2AC(A-1)(C-1)(A-C)*(10AC-5A-5C+3).                  (5)
```

The last factor becomes

```text
3*(90ab-5a-15b+1),
```

whose quotient by three is `1+a mod 3`. The terms `3` and `-5A` in (5)
have the same valuation at the chosen nested scale and can cancel. The
metric records their valuations, but cannot decide their residue sum.

This is consistent with the derivative obstruction that motivated the
case: the second divided difference of first-derivative rows approaches
`D^[2](P')=3D^[3]P`, which loses a direction in characteristic three.
Formula (5), together with the complete competing-minor table, identifies
the actual surviving residue in this concrete observer. A generic rational
rank calculation would discard it.

It would also be wrong to track only (5). With `b=1`, its last normalized
factor is `85a-14`. At `a=2,8,17` its valuations are respectively one, two,
and three. In fact its valuation can be made arbitrarily large by choosing
`a` modulo higher powers of three, because 85 is a unit. Yet the full
rank-four ideal has valuation exactly 28 throughout this exceptional
branch. The other Hermite minor caps the change. This is why both residue
cancellation and every competing minimum must be retained.

## 5. Minimal diameter among integer four-node ternary counterexamples

For a finite set of integers write its diameter as `max(x)-min(x)`.
The minimum possible value of the larger of the two diameters in a pair
of four-node integer configurations with the same three-adic metric tree
but different Smith spectra is **81**.

The upper bound is (2). For the lower bound, translate the least node to
zero; translation is unimodular on the polynomial coefficient module.

If the common depth is zero, separate the residue classes by integral
CRT. The monic products of squared node factors in different classes are
coprime over `Z_(3)`, so the observer is integrally equivalent to the direct
sum of its residue-class observers on their respective degree boxes. Every
class has at most three nodes. The proved metric-only laws through three
nodes therefore determine this direct sum completely from the distance
tree. No counterexample can have common depth zero.

If the common depth is positive and the diameter is below 81, all translated
nodes are multiples of three, so the complete universe is

```text
(0,x,y,z),  {x,y,z} subset {3,6,...,78},
binomial(26,3)=2600 configurations.
```

An exhaustive exact check computes each full Smith list and the canonical
pairwise valuation matrix, minimizing its six upper-diagonal entries over
all 24 vertex permutations. There are 11 metric trees, and no tree has two
different spectra. This is a complete finite reduction, not an estimate of
how often a counterexample occurs. Combined with the CRT argument it proves
the lower bound 81. It is specific to integer nodes, four nodes, prime
three, and this diameter criterion.

## 6. Validity boundary and synthesis

The precise failed implication is:

```text
same entire pairwise valuation tree
  -> same degree-below-eight two-jet Smith spectrum.
```

It remains true through three nodes, and for the separately proved ternary
`(2,2)` family, but fails for the nested `(3,1)` tree in (2). The strongest
repair proved here keeps exactly the first residue `chi` in addition to
that tree. The unresolved all-depth `(3,1)` law must allow residue data at
ties between weighted-minor terms.

The source is the full observer. The first map replaces its minors by
weighted valuations; that map preserves lower bounds but loses leading
residue sums. The sidecar restores the critical normalized residue, and a
competing minor supplies the upper bound. The exact kernel comparison above
checks an actual consequence of the reconstructed lattice, not merely a
minor statistic. This is a concrete instance where the metric quotient is
too small even though both the determinant and worst precision survive.

## 7. Reproduction and evidence manifest

Run from the repository root:

```text
python3 -u -c 'import runpy; runpy.run_path("04-computation/overnight3_20260906_smith_triple_single.py", run_name="__main__")'
```

The [source](../../04-computation/overnight3_20260906_smith_triple_single.py)
and [preserved output](overnight3_20260906_smith_triple_single.out) include:

- all 923 symbolic minors, 5,587 coefficients, two independent polynomial
  determinant algorithms, the complete normalized residue table, and the
  exceptional-branch ceiling witness;
- all 2,916 unit pairs with `-40<=a,b<=40` and neither divisible by three;
- 144 independent rational DVR controls with signed lifts, translation,
  and unit scaling;
- every minor of the original eight-by-eight observer in four controls,
  including both branches of (2);
- the full 2,600-configuration minimal-diameter check with canonical metric
  matrices, plus deeper critical-minor cancellation controls.

Modular elimination uses precision one greater than the known determinant
valuation, so all eight invariants are visible. Normal and optimized Python
runs give identical output and 1,282,215 explicit gates. No floating-point
constants or optimization-disabled assertions occur. The symbolic proof
does not import a repository theorem as a numerical oracle.

SHA-256 values over LF bytes:

```text
source:
f639fc43d51e36fcbedcbe0350a2fb238e19a9c9a2ec172532350e07d06aba11
output:
26878c358ef657ab6feb4faef207ba5e85ded531d422cc7f0c09f03520c7c4b7
complete symbolic-minor semantic digest:
172cc189edeae9ad4742909b7dacf23b3499c22928b2d38bc331a4e492326610
unit-pair semantic digest:
2757cd27c4c540e6b9fb2b1bfdba5059ac7dd7318b2824287922efdf2b757f3f
```
