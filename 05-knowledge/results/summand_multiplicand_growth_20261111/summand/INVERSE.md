# Deleting doubles has an exact binary inverse

Date: 2026-11-11. **PROVED:** the observation inverse, its coefficient bound,
exact row size, and sequential reconstruction. **FINITE-EXACT:** small
integer-matrix and impulse controls in`check.py`. The parallel multiplicand
lane supplies the square-deletion analogue; this note owns the additive side.

## The recovered object and the new observer

The May3 strict summand graph has `a→n` iff`1≤a<n` and`n≠2a`.
The September repair
`05-knowledge/results/arithmetic_braids_20260917_summand.md`, §1,
proves that its complement inside the ascending relation is exactly the
doubling forest. This does not identify two-parent constructibility with
ordinary reachability. The current note studies a **linear incidence
observer** on that precisely specified labelled shadow.

Assign an integer weight f(a) to every positive vertex and put f(0)=0.
Observe the new vertex's own weight plus all weights on its strict incoming
arcs:

    g(n)=sum_(1≤a≤n) f(a) − [2|n] f(n/2).                (1)

Identity is included, making the finite observation matrix unit triangular.
No pair is counted twice here: this is a weighted vertex-incidence sum,
not the number of unordered summand pairs. The latter needs its companion
sidecar and is a different observable. Over any commutative coefficient
ring, the same identities hold for ring-valued weights.

The closest old mechanism is the deleted swap-fixed diagonal in THM-2422,
`01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`.
The July reflection
`07-reflections/the-summand-alternation-is-an-involution-parity-law-kps-S131.md`
suggested promoting that fixed locus to a new object. Here it becomes the
exact correction tree of an inverse, rather than only an odd/even count.

## 1. One-vertex reconstruction follows binary ancestry

Set g(0)=0 and Δg(n)=g(n)−g(n−1). Subtract adjacent instances of(1):

    f(n)=Δg(n)+(-1)^n f(floor(n/2)).                     (2)

For even n, the newly deleted diagonal is f(n/2); for odd n, the previous
deleted diagonal is restored. This proves(2) directly, including n=1.
It is an exact online algorithm: the next observation needs only its
predecessor and the already reconstructed half-index weight.

Iterating along the binary ancestors gives

    f(n)=sum_(j=0)^floor(log2 n)
       (-1)^popcount(n mod2^j) Δg(floor(n/2^j)).          (3)

The empty remainder at j=0 has sign+1. Each next sign records the parity
of the current integer, so their product is precisely the parity of the
removed low binary digits. Formula(3) uses O(log n) observations to recover
one coordinate without reconstructing all smaller coordinates.

For comparison, retaining equal-parent arcs gives the weak observer
`g(n)=sum_(a≤n)f(a)`, whose inverse is simply Δg(n). Thus the entire extra
memory in(2) is caused by the deleted doubling forest.

## 2. The ordinary-series inverse is a Thue–Morse product

Let F(z)=sum_(n≥1)f(n)z^n and similarly G(z). Equation(1) is exactly

    G(z)=F(z)/(1−z)−F(z²).

Consequently, as formal power series,

    F(z)=sum_(j≥0) P_j(z) G(z^(2^j)),
    P_j(z)=product_(r=0)^j (1−z^(2^r)).                  (4)

For each coefficient only finitely many j contribute, since F and G have
zero constant coefficient. Iterating
`F(z)=(1−z)G(z)+(1−z)F(z²)` proves(4): after terms j=0 through J, the
remainder has no terms below degree 2^(J+1). Unique binary expansion gives

    P_j(z)=sum_(k=0)^(2^(j+1)−1) (-1)^popcount(k) z^k.    (5)

This is an actual finite binary-sign operator at each scale, not a
numerical resemblance to the Thue–Morse sequence.

## 3. Every inverse matrix entry is −1,0,or1

Let K(a,n) be the coefficient of g(a) in f(n). Equations(4)–(5) give

    K(a,n)=sum_{j≥0: 0≤n−2^j a<2^(j+1)}
                    (-1)^popcount(n−2^j a).             (6)

**Coefficient theorem.** For all positive a,n,

    K(a,n)∈{−1,0,1}.                                    (7)

*Proof.* The scale-j contribution is supported on the half-open interval
`2^j a≤n<2^j(a+2)`. If a≥2, these intervals are disjoint: the next interval
starts at2^(j+1)a, at least the preceding upper endpoint. There is at most
one contribution.

For a=1 there can be two contributions, at adjacent scales j,j+1. Write
`n=2^(j+1)+r`, where0≤r<2^j in the overlap. The two remainders in(6) are
`2^j+r` and r, with opposite binary signs. They cancel exactly. No three
scales overlap. This proves(7). QED.

The small source a=1 is the entire cancellation boundary; it corresponds
to the first missing diagonal arc1→2. In particular K(1,1)=1. For n≥2,
K(1,n)=0 when the binary expansion begins10, and equals
`(-1)^(popcount(n)−1)` when it begins11.

The multiplicand lane's analogous correction is indexed by repeated square
roots and binary prime-exponent signs. On prime-power inputs the exponent
map`a↦p^a` identifies that incidence matrix with(1), so its squarefree-source
cancellation at a single prime must agree with(7). The full multiplicative
matrix has additional independent prime coordinates and can have entries±2;
this additive proof does not erase those coordinates.

### A synchronized-scale tensor, rather than separate coordinate inverses

The precise common factor is

    κ_j(a,e)=(-1)^popcount(e−2^j a)
              if 0≤e−2^j a<2^(j+1), and 0 otherwise.

Then the additive inverse is `K_+(a,n)=sum_j κ_j(a,n)`. The proved
square-deletion formula in [the multiplicand packet](../multiplicand/PROOF.md)
has the exact prime-exponent form, for d,n≥2,

    K_×(d,n)=sum_(j≥0) product_(p prime) κ_j(v_p(d),v_p(n)).

Only primes dividing dn matter, since κ_j(0,0)=1. The allowance a=0 in
this coefficient kernel records an absent prime; it does not add a zero
vertex to the positive summand graph. This identity follows by multiplying
the binary-sign factors of each Euler factor at the same scale j. It is
finite because at least one exponent v_p(d) is positive.

The order of operations is essential: first tensor the prime coordinates
at one common scale, then sum the scales. For d=6,n=36, the two exponent
coordinates both have source 1 and target 2. Each separate additive
inverse is `K_+(1,2)=−1+1=0`, whereas

    K_×(6,36)=(−1)(−1)+(+1)(+1)=2.

Thus taking a product of completed one-coordinate inverses destroys the
synchronized scale. This is the exact extra sidecar that separates the
double-deletion bound 1 from the square-deletion bound 2. On a single
prime-power ray there is only one nontrivial coordinate, recovering the
additive result without loss.

## 4. The exact observation budget at every new node

Let R(n)=sum_(a≤n)|K(a,n)|, which by(7) is the number of nonzero entries in
inverse row n. Then

    R(1)=R(2)=1,
    R(n)=2+R(floor(n/2)) for n≥3.                        (8)

Indeed in(2), the new columns n,n−1 are disjoint from all half-index
columns when n≥3. At n=2 the g(1) terms cancel, leaving f(2)=g(2).
Writing k=floor(log2 n), equation(8) yields

    R(n)=2k−1 if n≥2 begins10 in binary,
    R(n)=2k+1 if n≥3 begins11 in binary.                 (9)

Thus the frontier alternation becomes a precise two-block row-complexity
law at every binary scale. The inverse has no coefficient amplification,
but its exact number of contributing observations grows logarithmically.
An observation error bounded by ε therefore changes reconstructed f(n)
by at most R(n)ε. This bound is sharp for a single row by choosing each
observation error with the sign of its nonzero coefficient.

The inverse also gives a complete positivity test: an integer observation
sequence is the image of nonnegative weights exactly when the recursively
reconstructed values in(2) are all nonnegative. This is an exact predicate
with its full input retained, not a positivity conclusion from a mean or
from the parity of g alone.

## Reproduction and transfer contract

`python3 -B check.py` forms finite observation matrices directly from the
strict-arc definition, inverts them by integer triangular substitution,
and compares every entry with(6). It also checks(2) on signed weights,
all impulse columns through128, the cancellation boundary, and(9).

Source→target: labelled strict summand shadow plus identity→linear observer.
Preserved: all weighted vertex incidences. Lost: co-parent pairings and any
two-parent generation requirement. Sidecar for those questions: the paired
operation fibres, as in[CLOSURE.md](CLOSURE.md).

Source→target: observation sequence→binary-ancestry inverse. Preserved:
every weight exactly. Lost: nothing, with f(0)=g(0)=0 fixed. Cheapest
hostile: n=2, where failing to cancel the source1 terms produces the wrong
answer despite both scale contributions being individually nonzero.
