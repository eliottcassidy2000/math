# The exact roles of 12, 31, 41, 49 and 5

Session242, 2026-11-10. **PROVED elementary identities and denominator bound;
FINITE-EXACT uniqueness certificate; OPEN tournament transfer.** No literature
novelty is claimed. This note consolidates classical arithmetic and the
repository's earlier constructions, and separates three incompatible notions
of descent.

## 1. Inheritance and provenance

Closest proved mechanism: [AP square-sum classification](../../ap_square_pythagorean_20261019.md), §§1–2 and 8, and
[congruent square-gap classes](../../congruent_square_gap_classes_20261019.md), §§2–6.
The former explicitly identifies `31,41,49` as half of the pair-sum roots
`62,82,98` of the minimum positive AP square-sum clique `(482,3362,6242)`.
The latter gives the class-5 triangle and the elliptic point below.

The value `12` is recovered algebraically as the least common denominator
of the rational square roots of a gap-5 AP. We do not claim to have found an
earlier text explicitly assigning it that role: the inherited sources
explicitly display denominators 6 and 24 in the triangle and graph-gap maps.
The denominator-12 assertion below now gives the exact missing relation.

Canonical hostile: doubling the class-5 point gives a triangle outside the
positive graph-label chamber. Corrected near miss: an AP square-sum clique is
not a Euclidean triangle of those graph labels. Least-used sidecar: the actual
elliptic point, including sign and multiplier, rather than merely its area
squareclass. Relevant tree source is THM-3756,
`01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md`;
branch conventions must be retained.

Concept board: square AP / denominator; triangle / area squareclass;
elliptic point / signed multiple; Berggren node / primitive height;
tournament / Hamiltonian-path count and rooted response. The calculations
below explicitly identify the losses between the first four; no map to the
fifth is asserted.

## 2. Lossless arithmetic packet

The fundamental equality is

`41² - 31² = 49² - 41² = 720 = 5·12²`.

Consequently the three rational squares

`(31/12)², (41/12)², (49/12)²`

form an AP of common difference 5. Their roots are not themselves an AP.
The rational right triangle obtained by sum and difference of the outer roots
has sides

`a=(49-31)/12=3/2`, `b=(49+31)/12=20/3`, `c=2·41/12=41/6`.

It has area 5. Multiplying the sides by 6 gives primitive `(9,40,41)`,
area `180=5·6²`. The positive graph-label construction gives

`2(C²-4AB), 2C², 2(C²+4AB) = (482,3362,6242)`.

Its gap is `2880=5·24²`, and its pair-sum roots are `(62,82,98)`.
Thus 12, 6 and 24 are three different normalization scales in one exact
commuting arithmetic construction. Neither 31 nor 41 nor 49 is being
asserted to be the area; 49 is a square numerator, not evidence that the
area-squareclass equals 1.

## 3. The denominator 12 is forced and sharp

**Proposition.** Suppose `N` is a positive squarefree integer and
`(p/d)²,(q/d)²,(r/d)²` is a strictly increasing AP of positive rational
squares with difference `N`, where `p,q,r,d` are positive integers. Then

`(4 if N is odd, otherwise 2) · (3 if 3 does not divide N, otherwise 1)`

divides `d`. In particular `12|d` when `N=5`. No reduced-fraction assumption
is needed.

**Proof.** The squared integers have common step `D=N d²`. The only square
residues mod16 are `0,1,4,9`. Checking these four residues gives the complete
list of `(first residue, step)` for a three-square AP:

`(0,0),(1,0),(1,8),(4,0),(9,0),(9,8)`.

Hence `8|D`. For odd squarefree N this implies `v₂(d)>=2`; for even
squarefree N it implies `v₂(d)>=1`. Modulo3, a three-square AP must have step
zero, since the square residues are just 0 and 1. If `3∤N`, then `3|d`.
Combining the coprime divisibilities proves the assertion. ∎

At `N=5`, the displayed denominator12 example proves sharpness. Its positive
increasing numerator triple is unique at that denominator: factor

`(q-p)(q+p)=720`.

For each divisor `u<=sqrt(720)`, set `v=720/u`, require equal parity, and test
`p=(v-u)/2>0`, `q=(v+u)/2`, and `q²+720` a square. The only survivor is
`(p,q,r)=(31,41,49)`. This finite divisor certificate exhausts the universe;
it does not impose an arbitrary numerator cutoff. `verify.py` executes it,
and independently verifies the absence at denominators 1 through 11.

This is a minimal-denominator statement, not a classification of all
class-5 witnesses, nor a solution of the congruent-number problem.

## 4. Two natural elliptic maps differ by signed doubling

This is the useful map distinction hidden in the original numbers.
Let `0<p<q<r` and `q²-p²=r²-q²=N d²`. Set

`a=(r-p)/d`, `b=(r+p)/d`, `c=2q/d`.

The triangle map used by the inherited notes gives

`P=(N(a+c)/b, 2N²(a+c)/b²)` on `E_N: y²=x³-N²x`.

The direct square-AP map instead gives

`Q=(q²/d², pqr/d³)`.

Then **`Q=-2P`**, not `Q=P`.

**Proof.** The membership of Q follows by multiplying
`(q²/d²-N)(q²/d²)(q²/d²+N)`.
Since `N=ab/2`, write `x(P)=a(a+c)/2` and `y(P)=a x(P)`.
Using `b²=c²-a²`, the tangent slope is `(2a+c)/2`. The usual tangent formulas
therefore simplify to

`x(2P)=c²/4`, `y(2P)=c(a²-b²)/8`.

Substitution of a,b,c gives `x(2P)=q²/d²` and
`y(2P)=-pqr/d³`. ∎

For our packet,

`P=(25/4,75/8)` and `Q=(1681/144,62279/1728)=-2P`.

The inherited class-5 double has the negative y-coordinate, as required.
Applying the triangle map inversely to Q gives a *different* right triangle,
not `(3/2,20/3,41/6)`. The negative-first-label example in the inherited note
therefore does not contradict positivity of the original clique: it concerns
the different triangle obtained after doubling. Retaining the map and the
signed multiplier resolves the apparent conflict.

The proof is algebraic for every such rational square AP. The companion
checks the general doubled-coordinate identities independently for all 331
primitive Euclidean-parameter triangles with parameter `2<=m<=40`.

### 4.1 The two positivity chambers must also be distinguished

For the direct AP point Q=(X,Y), take square roots of X-N, X and X+N.
The rational square-sum graph labels reconstructed from these pair-sum roots
are

`((X-2N)/2, X/2, (X+2N)/2)`.

Thus this **direct AP chart** is positive exactly when `X>2N`. At our Q,
`X=1681/144>10`; multiplying these three labels by `576=24²` gives precisely
`(482,3362,6242)`. Q is admissible in this chart.

In contrast, first apply the inverse **triangle chart** at a point (X,Y),
and then its quadratic clique construction. Its positivity condition is

`(X²+N²)² > 8NX(X²-N²)`.

This fails at our Q, as in the inherited doubled-triangle example. These are
different maps and different inequalities. The signed-doubling identity does
not identify their positivity chambers.

### 4.2 Exact halving obstruction without citing a general descent theorem

For any rational R=(x,y) with y nonzero on E_N, duplication gives

`x(2R) = ((x²+N²)/(2y))²`,
`x(2R)-N = ((x²-2Nx-N²)/(2y))²`,
`x(2R)+N = ((x²+2Nx-N²)/(2y))²`.

These identities follow by expansion using `y²=x(x²-N²)`. Therefore every
finite rational double with nonzero y-coordinate has all three of its shifted
x-coordinates rational squares. This is the needed elementary direction of
the classical full two-descent criterion; no converse theorem is invoked.
For a positive rational square AP, §4 already supplies an explicit rational
half, namely `-P` for Q.

At `P=(25/4,75/8)`, however, `x(P)-5=5/4` is not a rational square. Hence P
has no rational half, and the same is true of -P. The chosen inverse step
`Q -> -P` leaves the direct square-AP locus.

One can strengthen this to **Q is in 2E_5(Q) but not in 4E_5(Q)**, checking
all possible halves rather than only the chosen one. If `Q=4R`, then
`2R=-P+T` for some rational two-torsion point
`T in {O,(0,0),(5,0),(-5,0)}`. The x-coordinates of these four translates
are, respectively,

`25/4, -4, 45, -5/9`.

The first fails the shifted-square test above; the other three already fail
the square test on x itself. Thus none is a rational double. The translate
formulas follow directly from the chord law: translation by (0,0) sends
x to `-N²/x`, by (N,0) to `N(x+N)/(x-N)`, and by (-N,0) to
`-N(x-N)/(x+N)`. This is a concrete finite halving obstruction, not a
proof of finite halving depth for every elliptic point.

## 5. Why the obvious descent does not transfer the predicate

On the inherited repeated outer Berggren branch the forward chain is

`(3,4,5) -> (5,12,13) -> (7,24,25) -> (9,40,41)`.

Canonical parent descent reverses this chain and strictly decreases the
hypotenuse. But its area squareclasses in forward order are

`6,30,21,5`.

In particular the first descent from class5 leaves class5. It also leaves
the positive graph-label chamber: `25²=625<4·7·24=672`, whereas
`41²=1681>4·9·40=1440`. Therefore this canonical descent preserves the
Pythagorean equation but neither the required congruent class nor the graph
positivity predicate. No conclusion about all possible descent procedures is
being asserted.

Conversely elliptic multiplication preserves class5 but is not the same
well-founded height descent. The inherited exact valuation formula yields
`v₂(x(2P))=-4`, `v₂(x(4P))=-6`, and subsequently a drop by 2 on every doubling.
This is denominator growth and a non-torsion certificate, not a decreasing
positive-integer rank. Moving towards the real point at infinity may improve
the leg-ratio chamber while making rational height large. Real closeness,
valuation and numerator/denominator height must not be interchanged.

Fermat's square-area descent, cited in the inherited congruent-gap note §8,
works because every assumed square-area solution produces a *smaller solution
of the same forbidden type*. Its key exportable feature is predicate
preservation, not the presence of squares or the word descent.

## 6. A precise Chan–Pak transfer obligation, not a proof

Let A denote the target set in the particular Chan–Pak spectrum formulation
being studied, with its genuine known exceptions removed. Let B be a proved
base set. A sufficient induction mechanism would consist of, for every
`h in A\B`, finitely many smaller targets `h_i in A`, auxiliary conditions
on their actual tournament witnesses, and a proved constructor producing a
tournament with Hamiltonian-path count exactly h. All auxiliary conditions
must be supplied by the induction hypothesis or established independently.

The distinction from Fermat descent matters: a hypothetical missing h has
*no tournament to delete*. One must descend the target and construct from
smaller witnesses, or encode a hypothetical obstruction that actually
exists. Deleting vertices from already existing tournaments alone cannot
supply the missing target's predecessor object. Moreover a scalar witness
`H(T)=h_i` may fail the rooted-response condition required by a constructor;
that response must be part of the induction state, not silently inferred.

Connection contract:

- Source: rational square AP with denominator and point multiplier retained.
- Target: a prospective witness-carrying recursive tournament constructor.
- Transferable mechanism: well-founded rank plus preservation of the exact
  recursive predicate, with an explicitly proved lift from children.
- Not transferred: squareclass, elliptic group law, graph square-sum edges,
  or rational denominator to tournament Hamiltonian-path count.
- Missing map: a construction whose exact path count and required rooted
  response satisfy the proposed recursive equations.
- Cheapest hostile: the class5 Berggren parent already proves that retaining
  the broad ambient equation does not retain the desired subclass. In the
  tournament lane, test whether two witnesses of the same H have different
  required rooted responses before using scalar induction.

The denominator bound offers a practical analogy only: build congruence
filters from the actual recurrence, retain its extra state, then prove a
strict rank decrease on every admitted inverse branch. The two modular tests
here are necessary and happen to be attained by a witness at N=5; necessity
alone supplies no global coverage theorem. No Chan–Pak assertion follows
from the arithmetic packet.

## 7. Reproduction and status boundary

From repository root:

```
python3 05-knowledge/results/session242_congruent_descent_20261110/arithmetic/verify.py
python3 -O 05-knowledge/results/session242_congruent_descent_20261110/arithmetic/verify.py
```

The outputs are byte-identical and stored in `verify.json`. Checks use explicit
exceptions, not removable assertions. Residue universes, divisor universe,
small-denominator negative controls and the finite triangle parameter universe
are specified above and in source. General claims follow from the algebraic
proofs, not from the 331-triangle test.

Priority boundaries: classical AP/right-triangle/elliptic correspondences are
not claimed new. The inherited graph positivity theorem is credited above.
The denominator sharpening and signed-map reconciliation are derived here
as research synthesis, with no firstness claim. Chan–Pak completeness and
any arithmetic-to-tournament witness map remain **OPEN**.
