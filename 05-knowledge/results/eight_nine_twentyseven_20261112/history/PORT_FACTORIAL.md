# Selective vertices shift the excursion factorial filtration

Date: 2026-11-12. **PROVED:** the all-order divisibility theorem and its
precision cutoff. **FINITE-EXACT:** literal path-layer controls in check.py.
This gives a universal divisor, not an exact context ideal or positive
target-coverage theorem.

**Concurrent proof credit.** The incoming
[Session295 construction atlas](../../session295_construction_atlas_20261113/atlas.md#new-selective-port-filtration),
read at `origin/main` commit `9138a603e2`, independently proves exactly the
same all-r filtration by exceptional-excursion contraction and free block
permutations. Neither result strengthens the other. This packet also
records the general (k,s) top-layer formula and the indexed-reflection
provenance; Session295 explicitly exhibits its sharp mod27 cutoff case.

## 1. The object and statement

Let T be a tournament, and let C be a nonempty induced subtournament.
An exterior vertex is **uniform on C** when it beats every vertex of C
or loses to every vertex of C. Suppose at most r exterior vertices are
not uniform. Such vertices are called selective; their individual arcs
to C are retained, not replaced by a scalar port count.

Let h_k(T;C) count Hamilton paths of T whose C vertices occur in exactly
k nonempty consecutive runs. Write x_+=max(x,0). Then for every k≥1,

    ((k−r)_+)! ((k−r−1)_+)! divides h_k(T;C).          (1)

For r=0 this is the familiar k!(k−1)! excursion factor. For r=1 it is
(k−1)!(k−2)! when k≥2. The proof explains why one selective vertex costs
at most one whole block of factorial divisibility, rather than erasing
all the permutation freedom at once.

## 2. Contract the exceptional excursions, not their endpoint data

Split a Hamilton path uniquely into its k C-runs, k−1 nonempty internal
exterior excursions, and a possibly empty exterior prefix and suffix.
Call an exterior excursion exceptional if it contains a selective vertex.
Let e be the number of exceptional internal excursions, and let a be the
number of exceptional prefix/suffix excursions, so

    e+a≤r.                                           (2)

Every ordinary internal excursion starts at a uniform vertex beaten by
C and ends at a uniform vertex beating C. Its endpoint requirements
therefore accept every C vertex.

Absorb each exceptional internal excursion into its two adjacent C-runs.
Consecutive absorptions merge transitively. This gives

    L=k−e

compound blocks, each beginning and ending in C, separated by L−1 wholly
uniform internal excursions. Every compound block retains its complete
actual path, including its exceptional exterior vertices and all internal
seams. No nonuniform incidence is forgotten.

If the exterior prefix is exceptional, anchor the first compound block;
if the suffix is exceptional, anchor the last. Let b be the number of
distinct anchored blocks. Then b≤a. The other L−b compound blocks can
be permuted arbitrarily: all new external seams are through uniform
vertices. Independently, all L−1 ordinary internal excursions can be
permuted arbitrarily. Leave the prefix, suffix and anchored blocks fixed.

These two actions commute. Their joint action is free because all blocks
and excursions have nonempty, pairwise disjoint actual vertex sets, and
the alternating decomposition recovers their orders from the resulting
Hamilton path. A combined nonidentity permutation cannot return the same
literal vertex word.

The orbit therefore has size

    (L−b)! (L−1)!.                                   (3)

The contraction data are recovered uniquely from each path, so these
orbits partition the k-run layer. Equations (2)–(3) give

    L−b≥k−e−a≥k−r,
    L−1=k−e−1≥k−r−1.

Factorial divisibility proves (1). If L=1, both boundary anchors may hit
the same compound block. There is still one orbit position and b≤a;
the requested divisor in (1) is then 1, so this boundary is harmless.
Equivalently, one can dispatch k≤r+1 first, where (1) is already trivial.
QED.

This action is on literal path blocks, not an automorphism action on C.
It can force prime factors when C is asymmetric.

## 3. Exact finite prime-power cutoffs

For a prime p and s≥1, define

    J(p,s)=min{j≥1 : v_p(j!)+v_p((j−1)!)≥s}.          (4)

Every layer with k≥r+J(p,s) vanishes modulo p^s. Hence

    H(T) ≡ sum_(1≤k<r+J(p,s)) h_k(T;C) (mod p^s).    (5)

This holds for every exterior size and orientation satisfying the stated
selective-vertex bound. Only nonempty C-runs are counted. The surviving
layers still require their actual endpoint constraints; (5) does not
identify them with unmarked ordered cover counts.

For p=3 the initial values are

| precision | J(3,s) | last potentially nonzero k, r=0 | r=1 |
|---|---:|---:|---:|
| 3 | 3 | 2 | 3 |
| 9 | 4 | 3 | 4 |
| 27 | 6 | 5 | 6 |
| 81 | 7 | 6 | 7 |

For an eight-vertex C, the fully split k=8 layer is automatically divisible
by 81 when r≤1, and by 27 when r≤2. At least three selective vertices are
necessary before this particular layer can contribute nonzero modulo 27.
This is a statement about the k=8 layer, not about the whole count: low-run
layers can already break divisibility with one selective vertex.

## 4. A literal top-layer family

Let C be a transitive k-vertex core, k≥2. Add k−2 uniform pairs
x_i,y_i, with C→x_i and y_i→C, and exterior total order
`x_1,y_1,...,x_(k−2),y_(k−2)`. Add one selective vertex z, beating a
chosen s-element subset of C and losing to the other k−s vertices, with
1≤s≤k−1. Orient its remaining exterior arcs arbitrarily.

In the k-run layer all C-runs are singletons. There are k−1 internal
excursions. The exterior has only 2k−3 vertices, forcing z to be a
singleton excursion and every ordinary excursion to be one x_i,y_j pair;
prefix and suffix are empty. The allowed pair relation i≤j has only the
identity perfect matching, since summing i≤sigma(i) forces equality.

The exceptional singleton may occupy any of k−1 internal slots. Its
adjacent ordered C vertices have (k−s)s choices. The remaining C vertices
and the ordinary excursions each have (k−2)! orders. Consequently

    h_k=(k−1)s(k−s)((k−2)!)^2.                       (6)

This independently exposes the shifted factorial in (1). For k=5,s=1,
the layer has 576 paths and v_3=2; a universal extra factor of 3 at that
layer would be false. This control concerns a declared path layer, not
the full Hamilton count.

At k=6,s=1, the same formula gives

    h_6=5*5! *4!=14400,    v_3(h_6)=2.

This is exactly Session295's independently derived sharp one-port model.
Thus the mod27 cutoff k>=7 for r=1 cannot be lowered: a physical six-run
layer can remain nonzero modulo27. The general formula proves this
specialization; the local receipt's bounded family controls stop at k=5.

## 5. Relation to the asymmetric eight-vertex core

The inherited S293 core has rows

    (56,117,129,22,100,76,141,59)

and ordered path-cover counts

    (513,6102,32508,99432,187200,216000,141120,40320).

Their gcd is 9, while the exact uniform-context ideal has generator 27
after multiplying the k-th entry by (k−1)!. Its proof is
`05-knowledge/results/session293_cycle_returns_20261113/control/MODULE_IDEAL.md`.
Those facts are inherited here. The new result is what remains when
uniformity is broken at a bounded number of actual exterior vertices.

The full cover gcd does not price internal excursion permutations, while
an automorphism census does not see them at all. The compound-block proof
retains precisely which nonuniform excursions prevent those permutations.
It supplies a uniform precision compiler for selective contexts; it does
not assert that their low-run coefficients attain an arbitrary residue
or that a prescribed +2 repair exists.

Source→target: a Hamilton path with selective exterior incidences→its
compound blocks and uniform excursions. Preserved: every actual vertex,
every directed adjacency, the number of C-runs, and the whole nonuniform
interface. Lost only by taking the orbit size: which particular block
ordering is used. The needed sidecar is the labelled ordered content of
each exceptional compound block and any anchored endpoint block.
