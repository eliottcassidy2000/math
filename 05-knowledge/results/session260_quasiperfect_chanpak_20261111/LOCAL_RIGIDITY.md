# A complete one-edge barrier beside every consecutive excess-two source

**PROVED, using the new crossing-rank formula and inherited positive path-cover substitution.**
This is a statement about a named graph and every possible single-edge reversal,
not about all graphs with the same Hamilton count.

Let T(a,b)=C3(TT_a,TT_b,1), with A→B→z→A, and use one-based positions
A_i and B_j within their increasing transitive blocks. The arithmetic lane
proves H(T(a,b))=(2^a−1)(2^b−1)+2. Its boundary crossing reversal has
increment2−2^(a−1)−2^(b−1).

## 1. Complete singleton-edge formulas

Set f(0)=0 and f(t)=2^(t−1) for t≥1. Set

    g(0,t)=f(t),
    g(l,t)=2^(l−1)(2^t−1), l≥1.

Reversing z→A_i changes H by

    Δ_A(i)=2f(i−1)−2g(a−i,b).                  (1)

Reversing B_j→z changes H by

    Δ_B(j)=2f(b−j)−2g(j−1,a).                  (2)

Proof. Delete z. The remaining graph is the transitive chain A,B. Its
endpoint capacity from position x to y>x is c(y−x), with c(1)=1 and
c(d)=2^(d−2) for d≥2: choose an even-sized subset of the intermediate
vertices to form an even-vertex path, and the complementary tournament has
one Hamilton path. The sum of capacities over distances1 through t is f(t);
over distancesl+1 through l+t it is g(l,t).

The newborn z beats exactly A. Moving A_i from the outgoing to the incoming
shore removes all capacities from A_i to B and adds those from earlier A
vertices to A_i. Moving B_j in the opposite direction adds capacities from
B_j to later B and removes those from A to B_j. The exact endpoint-cut law
H=h+2q gives (1) and (2). These arguments retain every longer even path.

## 2. Internal module flips cannot decrease H

For a tournament R on m vertices, the number of ordered spanning r-path
covers with all r paths nonempty is at least r!S(m,r): for each ordered
partition of the vertices into r nonempty sets, Rédei gives at least one
Hamilton path in every induced set. TT_m attains equality. The general
substitution formula has nonnegative skeleton coefficients. Therefore
replacing a transitive module by any tournament of the same size cannot
lower H. If the new module is nontransitive, the one-run contribution
strictly increases, since its H is at least3 rather than1.

Thus reversing an adjacent arc inside A or B gives an isomorphic transitive
module and changes nothing; reversing a nonadjacent internal arc creates a
triangle and strictly increases H. This covers every internal arc.

## 3. All crossing arcs are classified

The [crossing-rank proof](bridges/CROSSING_FLIP.md) derives the exact response
for every A_i→B_j. For a,b≥2, its conclusion is stronger than the needed
bound: the boundary A_a→B_1 is the only crossing arc whose reversal decreases
H. Every other crossing reversal strictly increases H. The boundary decrement
is 2^(a−1)+2^(b−1)−2.

## 4. Sharp minimum for consecutive block sizes

**Theorem.** For every a≥2, every single-edge reversal of T(a,a+1) either
preserves or increases H, or decreases it by at least

    d_a=3·2^(a−1)−2.

Equality occurs both at the crossing A_a→B_1 and at B_2→z.

Proof. Internal and crossing arcs are settled above. In (1), f(i−1) increases
with i and g(a−i,a+1) decreases, so Δ_A(i) strictly increases. Its largest
value is Δ_A(a)=−3·2^(a−1), still negative and below −d_a. In (2), Δ_B(j)
strictly decreases with j. At j1 it is0; at j2 it is
2^(a−1)−2(2^a−1)=−d_a; later values are more negative. This exhausts all
arcs and proves sharpness. The j1 zero move simply changes the block sizes
(a,a+1) to(a+1,a), preserving the symmetric product. QED.

The formal excess-two source n_a=2^a(2^(a+1)−3) has nearby graph count
H(T(a,a+1))=n_a+3. Its desired Chowla neighbor n_a+1 lies exactly two below.
Because d_a≥4, **no single-edge reversal of this graph realizes n_a+1 at
any a≥2**. This upgrades the finite a2..5 observation to an all-order result.
Yet the arithmetic lane explicitly realizes n_a+1 for every a≥3 by changing
the graph. At a2 the target21 is genuinely impossible. Local rigidity and
spectrum membership have now been separated for the entire family, not
only for isolated numerical hostiles.

The theorem does not obstruct two-edge changes, extra vertices, source
replacement, or a different realization of n_a+3. The cheap next constructive
probe is the complete two-edge neighborhood, retaining the exact changed
arcs rather than assuming a second numerical decrement can be composed.
