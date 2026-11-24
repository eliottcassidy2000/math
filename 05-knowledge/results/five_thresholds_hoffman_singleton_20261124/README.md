# Five at a threshold: converse pairs, Heegner reduction, Hoffman–Singleton

2026-11-24. Extends the [five-bit Heegner selector](../chanpak_heegner_selector_20261122/HEEGNER.md),
[Fermat and solid constraints](../chanpak_fermat_solids_seams_20261122/README.md),
[MRB paired pressure](../chanpak_collatz_mrb_fibers_20261116/MRB_PRESSURE.md),
and the [native Moore spine](../moore57_graph_spine_tiles_20261124/).
The main new construction is an axis-symmetric 56-site tiling of an
exact blow-up of Hoffman–Singleton. The native 50-site version cannot
be axis-symmetric. The matching arithmetic mechanism is reduction
at 5, not a numerical equality of unrelated constants.

## 1. The tournament peak is global, and the quotient explains it

Let t_n be the number of tournament isomorphism classes, s_n the
self-converse count, and p_n=(t_n−s_n)/2 the number of distinct
converse pairs. Then s_n/p_n has its unique maximum **4 at n=5**
wherever defined (n≥4). The corresponding proportion s_n/t_n=2/3
also peaks there for n≥4; it is 1 at n≤3, where the pair ratio
itself is undefined.

Let m=floor(n/2), ε=n−2m. A self-converse tournament has an
anti-involution (Sylow 2 in its automorphism/anti-automorphism
group; a tournament automorphism group has odd order). Each
anti-involution has at most one fixed vertex. Its conjugacy
class within the tournament's extended group is unique by Sylow,
and since conjugation by an anti-automorphism can be replaced by
conjugation by the product with the marked anti-involution, there
is exactly one marked class over each self-converse class.

Normalize each of the m pairs by orienting its distinguished
edge from a_i to b_i. Choose independently a tournament on the
a_i, a simple graph whose bit records a_i→b_j iff a_j→b_i,
and ε center bits for each i. The b_i layer and opposite center
bits are forced. The only residual action is simultaneous
permutation of the m pairs. This gives a folded recurrence:
for an odd cycle partition λ of m, put r(λ)=Σ_i(λ_i−1)/2+
Σ_(i<j)gcd(λ_i,λ_j), and z_λ=∏_a a^(m_a)m_a!. Then

    s_(2m+ε) = Σ_(λ odd) 2^(2r(λ)+εℓ(λ)) / z_λ.

Indeed an even permutation cycle reverses an unordered pair
under a power and cannot preserve the quotient tournament;
otherwise r(λ) counts pair orbits shared by both bit layers.
In particular s_n≤2^(m(m−1)+εm). On the other side t_n≥2^C(n,2)/n!,
so s_n/t_n≤n!/2^floor(n²/4). At n=8 this is 40320/65536<2/3
and the bound strictly decreases thereafter: the multiplier is
(n+1)/2^floor((n+1)/2)<1. For n=4,5,6,7 the (s_n,p_n) pairs are
(2,1),(8,2),(12,22),(88,184). At n=5 the folded data itself is
two index pairs plus a center; the center provides exactly the
extra bits at the peak. The small counts agree with
[McKay's data](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html).

## 2. The analytic five and two different parity effects

For radius r, V_n(r)=π^(n/2) r^n / Γ(n/2+1). Strict log-concavity
in the continuous dimension follows from positivity of trigamma.
Consequently dimension 5 is the unique integer maximizer iff

    15/16 < r < 16/(5π).

At the endpoints the adjacent dimensions tie. At unit radius
the parity recurrence V_(n+2)/V_n=2π/(n+2) instead gives separate
odd/even winners **5 and 6**; the real maximum is not exactly 5.

In [Burns's cube interpretation](https://marvinrayburns.com/ThegeometryV11.pdf),
f(n)=n^(1/n) is the side length of the n-cube of volume n, and the
MRB constant is the convergent grouped sum Σ_(m≥1)
(f(2m)−f(2m−1)). The f(n) themselves peak at 3, not 5. But
**the largest negative correction among the MRB pairs is 5→6**.
Here f''(x)/f(x)=[(1−log x)^2+x(2 log x−3)]/x^4>0 for x≥5.
Thus adjacent drops decrease from 5 onward. The only earlier
negative pair is 3→4 and rational power bounds give
f(3)−f(4)<0.0281<0.0317<f(5)−f(6).
If all adjacent steps are allowed instead of the MRB pairing,
the largest drop is 4→5. These are comparable parity-sensitive
thresholds, not a relation between numerical values of the
MRB constant, volumes, or tournament counts. In particular the
volume peak moves when the radius moves; the discrete ratio does not.

## 3. Five is the actual residue characteristic

The inherited lists label the same five rows: d=4p−1=8b+3 and
b=(p−1)/2. The norm in the class-number-one ring of discriminant
−d is N(x+yω)=x²+xy+py², where ω=(1+√(−d))/2.
The least rational prime splitting in this ring is **p**.
A split or ramified prime ℓ<p would have a principal prime above
it and hence be the norm of a nonreal algebraic integer; the
displayed norm is at least p on y≠0 (complete the square, treating
|y|=1 separately). Conversely p has norm witness ω and splits
because its discriminant is 1 mod p. This proves a first-split
threshold, stronger than merely reading the five sets as labels.
See [Granville–Mollin on Rabinowitsch](https://citeseerx.ist.psu.edu/document?doi=5a6dfc882e0d54fc4d145823370a595171cfce28&repid=rep1&type=pdf).

| p | d | b | −d mod 5 | ring mod 5 | CM j mod 5 |
|---:|---:|---:|---:|---|---:|
| 3 | 11 | 1 | 4 | F5 × F5 | 2 |
| 5 | 19 | 2 | 1 | F5 × F5 | 4 |
| 11 | 43 | 5 | 2 | F25 | 0 |
| 17 | 67 | 8 | 3 | F25 | 0 |
| 41 | 163 | 20 | 2 | F25 | 0 |

In the inert rows conjugation is the Frobenius x↦x^5.
Reduction modulo 5 partitions the five CM j-values
−32³,−96³,−960³,−5280³,−640320³ exactly as the 5-row of the
repo's Fermat sign basis. The last three are the supersingular
rows; the first two are ordinary. The j calculation is on
the algebraic values, not rounded Ramanujan exponentials.
For potential good reduction these labels agree with
[the Deuring criterion as used by Elkies–Ono–Yang](https://arxiv.org/abs/math/0512350).

Now take the exact characteristic-5 curve Γ: X^6+Y^6+Z^6=0
that [Shimada uses to build Hoffman–Singleton](https://arxiv.org/pdf/1405.4643),
and its degree-four quotient

    Γ → E: (X:Y:Z) ↦ (X²:Y²:Z²),
    E: U³+V³+W³=0.

The quotient group consists of projective coordinate sign changes.
The Fermat cubic has j=0, by its order-three elliptic automorphism
fixing a point at infinity. It is therefore the same geometric
j-class as the reductions of all three inert Heegner rows. This
does not assert that their chosen elliptic models have the same
twist over F5, or that any CM lift is canonical.

There is an exact 36/126 transfer. Over F5, cubing is bijective,
so E has 6 rational points and trace 0. Over F25 it has 36 points
by the trace recurrence. Nine points have a zero coordinate.
The 126 points of Γ(F25) attain the upper count
4·(36−9)+2·9=126, so every nonbranch rational point has four rational
lifts and every branch rational point has two. This is a quotient
of the actual curve in the graph construction, not just a
match of the two numbers. Shimada takes components of the conic
incidence geometry to obtain the 50 vertices; E's 36 points
are not the graph vertices. Moreover, Γ has genus 10 and is maximal
over F25, so its 20 Frobenius eigenvalues are all −5. The elliptic
quotient has eigenvalues −5,−5 there. [Tate's isogeny criterion](https://pazuki.perso.math.cnrs.fr/index_fichiers/Tate66.pdf)
therefore gives Jac(Γ) isogenous over F25 to E^10. The Heegner
relationship sees the entire isogeny type, not merely one point count.

## 4. Platonic five, the Fermat five, and the star step

Spherical regular Schläfli pairs satisfy 1/p+1/q>1/2. For
integers p,q≥3 this gives the five ordered pairs, and the
icosahedron/dodecahedron pair is the only pair with a 5.
The rotational group A5 identifies with PSL2(F5), acting on the
six projective points (or six axes through opposite icosahedral
vertices). Its order-five element x↦x+1 fixes ∞ and cycles the
other five points. On that affine line, the scalar x↦2x swaps
the two unoriented step colors {±1} and {±2}; its determinant
class is nonsquare, hence this operation is in PGL2(F5) but not
PSL2(F5). It is the **color swap**, not a Platonic rotation.
Shimada's conic stabilizer is PGL2(F5), with its index-two A5,
so these exact group names are in the characteristic-5 construction
as well (Theorem 1.14 and the preceding paragraph).

Fermat number F_1=5 supplies the constructible regular pentagon
and its golden-ratio subfield. The same cyclotomic polynomial
has Φ5=(x−1)^4 in characteristic 5, and the real trace
t=ζ5+ζ5^(-1) obeys t²+t−1=0, with double root 2 mod 5.
On the other hand **every Fermat number except F_1 is nonsquare
mod 5**: F_0≡3 and F_i≡2 for i≥2. This is an all-index statement
and does not presume exactly five Fermat primes. For the five
known prime examples see the repo's THM-4861 for their
Hamiltonian-count realizations; do not infer that a 5-fold
tournament automorphism certifies a prime count.

## 5. Hoffman–Singleton rigidity and a symmetric 56-site encoding

See [Hafner's pentagon/pentagram chart](https://www2.math.ethz.ch/EMIS/journals/JACO/Volume18_1/np3862514r41w675.fulltext.pdf).
Permit arbitrary simple intra-fiber graphs on both halves of the
biaffine chart over F_q; assume no same-side cross-fiber edge and
use precisely the cross matching y=x+hi. **Even with those arbitrary
fiber graphs, a nontrivial diameter-two Moore graph forces q=5.**
If both halves have at least two fibers, neither side's fiber-index
set can be thinned: two vertices in distinct same-side fibers can
only meet through the other side, and varying their vertical difference
requires every opposite index i (respectively h).

With full indices fix δ≠0. The pair P(h,x),Q(i,x+hi+δ) is a
nonedge. Its unique common neighbor must use either the internal
edge (x,x+δ) in P_h or the internal edge (x+hi,x+hi+δ) in Q_i.
The two indicators sum to 1. For two h's the arguments differ
by (h−h')i, which traverses the field, forcing each indicator to
be constant in the vertical coordinate. So all P fibers have the
same symmetric connections and all Q fibers their complement.
Both colors are triangle-free. The elementary R(3,3)=6 argument
gives q≤5; regularity forces equal fiber degrees, so q is odd,
and q=3 has no half-size symmetric nonzero step set. At q=5 the
two symmetric sets, up to scaling, are precisely {±1},{±2}.
They give the Moore graph: different same-side fibers have one
cross common neighbor, a same-fiber nonedge has its unique
common neighbor on C5, and cross nonedges were already accounted
for. Thus changing intra-fiber arrows or thinning the full
matching chart cannot make it a degree-57 construction.

### Native symmetry obstruction, quantitative version

For any diameter-two Moore graph on v vertices of degree k>1,
a fixed-point-free involution cannot preserve all edges. Otherwise
every paired vertex must be adjacent to its partner: a nonadjacent
pair would have a fixed common neighbor. But an edge between two
different paired edges would make a 4-cycle together with its
image, so k≤1, a contradiction.

More quantitatively for any pairing ρ, let B be the number of
present edges whose images under ρ are absent, and let m be the
number of adjacent pairs. Each of the v/2−m nonadjacent pairs
forces a broken edge on its unique length-two path. A broken
edge can serve at most two such pairs, one from each endpoint.
On the 2m sites in adjacent pairs every between-pair edge is broken.
The degree sum gives ≥2mk−vk/2 internal edges there, of which
at most m are pairing edges. Therefore

    B ≥ max((v/2−m)/2, (2k−1)m−vk/2)
      ≥ v(k−1)/(2(4k−1)).

The edit distance to the nearest ρ-invariant simple graph is B,
one edit per mismatched two-edge orbit. In native backward-chord
representations of a Moore graph on an even Hamilton spine,
including spine edges in the decoded graph if needed, reflection
is fixed-point-free. So a native axis-symmetric tiling is
impossible; any attempt to symmetrize the graph needs at least
6 edits for Hoffman–Singleton and 401 for a hypothetical degree-57
graph. This is a marked-spine statement, not a claim that an
isomorphism class has no other self-converse representation.

### Six twin copies are necessary and sufficient in this encoding

For an involution let F be its fixed vertices, A the moved sites
paired with an adjacent site, and B' all other moved sites.
A has no fixed neighbors and no edges between different pairs
in A; each B' has exactly one fixed neighbor. Count A–B' edges
to get |A|≤|B'|. If d_F denotes the induced-fixed degree sum,
then |B'|=k|F|−d_F and v≤(2k+1)|F|−2d_F. Thus any nonidentity
involution of a hypothetical degree-57 graph fixes ≥30 sites
(use even parity; the right-hand bound without d_F already suffices).
In Hoffman–Singleton, |F|=2 or 4 would force the induced fixed
graph K2 or K1,3: it has diameter≤2 by the unique-common-neighbor
property, girth≥5 and odd degrees (moved neighbors occur in pairs).
Those degree sums give 26 or 48 on the right, strictly below 50.
No involution fixes fewer than six.

The chart has the involution P(h,x)↦P(−h,−x), Q(i,y)↦Q(i,−y)
with exactly six fixed sites, inducing K1,5. Replace each by
a pair of independent twins and swap the twins. This makes the
extended involution free on 56 sites. Pick one representative
of every reversal pair; the graph on the 28 representatives has
max degree ≤7, since it is induced from the original graph.
Its complement is Hamiltonian by Dirac and we can take a
complement path ending at a twin representative. Follow that
path by its reflected reverse on the other 28 sites. The central
twin pair is a nonedge. Start with the transitive tournament on
that spine and reverse the graph edges of the blown-up graph.
The spine is still Hamiltonian, and reflection reverses every
tournament arc. Merging false twins in the backward graph
recovers exactly Hoffman–Singleton, since distinct original
vertices cannot have equal neighborhoods. The replay gives an
explicit labelled spine. The 56 bound is optimal **for this
fixed-site duplication scheme**; not a lower bound over arbitrary
self-converse encodings. A count of Hamiltonian paths of the
56-site output is not claimed.

### Why the other 126 is not an unsigned identification

The native backward-chord count along a Hamilton spine is
175−49=126. The same integer is the cycle-space dimension
175−50+1. But the signed action on that space is not the
permutation action on the 126 Hermitian points: for the
explicit six-fixed-point involution, the graph has five
fixed edges with individually fixed endpoints, and six
fixed edges with swapped endpoints (A has 12 sites).
The traces on oriented edge chains and vertex chains are
therefore −1 and 6, respectively, giving trace −6 on H1
by the graph chain sequence (H0 trace 1). A permutation action
on curve points has nonnegative trace at every element.
The equality of cardinalities cannot produce this equivariance
without a sign or a different module.

Replay: run python3 replay.py in this directory. It checks the
converse ledger by Burnside, the root comparisons by rational
powers, the finite-field fibers, all Hoffman–Singleton neighbor
conditions, and every arc of the marked self-converse tiling.
