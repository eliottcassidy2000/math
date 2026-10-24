# Fifteen Richelot kernel choices: incidence, outer duality, and weighted orbits

**Status: PROVED elementary finite geometry; FINITE-EXACT complete certificates;
CITED algebro-geometric interpretation.** Session 52, 2026-10-24.

This packet recovers a genuine fifteen-element geometric carrier, rather than
identifying unrelated lists of length fifteen. Its concrete addition is a
fully explicit syntheme-graph/Kneser-graph map and a replay of the published
Bolza kernel action. The latter turns fifteen kernel choices into four
weighted isogeny classes with weights **6,4,4,1**. The weight-six class, not
the singleton, has the cited `j=8000` elliptic-square codomain.

## 1. Inheritance and scope

The main carrier is already proved in
[cardinality_bridges_20261018m.md, §3](inputs/cardinality_bridges_20261018m.md):
perfect matchings of six labels correspond to the fifteen Lagrangians of a
binary symplectic four-space. That note also records the 45 point-line flags,
the non-equivariance witness, the marked-A5 bridge, and the product/gluing
split 9+6. These are recovered results, not discoveries of this session.

The six-pentad outer automorphism is likewise inherited from
[information_observers_collatz_20261014a.md, §3](inputs/information_observers_collatz_20261014a.md)
and classical work. [Howard--Millson--Snowden--Vakil, §1.4](https://arxiv.org/html/0710.5916v1)
defines the syntheme/pentad construction and its outer action. Our new
standalone checker independently generates all relevant finite objects;
missing historical checker files in this sparse checkout are not claimed as
rerun dependencies.

The geometric anchor is
[THM-4012-weighted-leading-face-good-elliptic-factor-observer, §5](inputs/THM-4012-weighted-leading-face-good-elliptic-factor-observer.md).
It independently derives the Bolza elliptic quotients and `j=8000`, and keeps
the higher-weight stable-specialization hypothesis explicit. This packet
neither strengthens its conditional claims nor addresses the Jacobian
conjecture.

Inheritance checklist:
- Closest mechanism: the explicit even-subset symplectic carrier.
- Hostile: equal cardinalities with different fixed-point characters.
- Corrected near miss: conflating kernels, their automorphism orbits, and
  distinct quotient isomorphism classes.
- Least-used sidecar: the actual six pentads and the published Bolza action.

Concept board: counted object; incidence; relabelling action; outer transport;
weighted quotient. The finite graph below is an incidence graph on kernels
at one source, **not** an isogeny graph between abelian surfaces.

## 2. Fifteen points and fifteen Lagrangians

Let `V=F2^4` with

    B(x,y)=x0*y1+x1*y0+x2*y3+x3*y2.

A Lagrangian is a two-dimensional totally isotropic subspace. Every nonzero
`v` has a three-dimensional orthogonal complement. Its six elements outside
`<v>` fall into three pairs `{w,w+v}`, so exactly three Lagrangians contain
`v`. Each Lagrangian has three nonzero vectors. Counting flags gives

    15 * 3 = number_of_Lagrangians * 3,

hence fifteen Lagrangians and 45 flags. Equivalently there are `15*6` ordered
isotropic independent pairs and six ordered bases per plane. This also
proves completeness, rather than merely exhibiting fifteen examples.

For a six-label set `Omega`, use the equivalent model

    V = {even subsets of Omega} / <Omega>,
    B([S],[T]) = |S intersection T| mod 2.

The even-subset space has dimension five, with radical `<Omega>`; quotienting
gives the four-dimensional symplectic space. Each nonzero class has a unique
two-element representative, called a duad. Distinct duads are orthogonal
exactly when disjoint. A perfect matching `M={e1,e2,e3}` gives the plane

    L_M = {0,[e1],[e2],[e3]},

because the three pair classes sum to zero. Conversely the three nonzero
vectors of an isotropic plane have disjoint duad representatives and form a
matching. This is a natural `S6`-equivariant bijection. The matching count
`6!/(2^3*3!)=15` is therefore the same carrier count, not a coincidence.

The checker constructs an explicit symplectic basis of the quotient, saves
every duad's four-bit vector, and checks both addition and pairing for every
pair of even subsets. It enumerates planes independently of matchings and
checks equivariance under all 720 relabellings.

## 3. The syntheme graph and its explicit outer duality

Let `G` have the fifteen perfect matchings as vertices. Distinct matchings
are adjacent when they share a duad.

**Proposition.** `G` is strongly regular with parameters `(15,6,1,3)`.

**Proof.** A matching has three edges; each lies in two other matchings, so
its degree is six. If `M,N` share edge `e`, the unique third matching
containing `e` is their only common neighbor. A matching containing distinct
edges from `M\{e}` and `N\{e}` would need those edges disjoint, but the two
different matchings of the remaining four labels have no disjoint cross
pair. If `M,N` share no edge, their union is a six-cycle. Each edge of `M`
is disjoint from exactly one edge of `N`. These three disjoint pairs each
have a unique completing edge, giving exactly three common neighbors.
No other common neighbor is possible. QED.

Thus its adjacency matrix obeys `A^2=3I-2A+3J`. The spectrum is
`6^1,1^9,(-3)^5`: positivity of the nonadjacent common-neighbor count gives
connectedness; on the all-ones complement the polynomial is
`(A-I)(A+3I)=0`; dimensions and trace determine the multiplicities. These
are combinatorial spectra, not superspecial-isogeny spectral assertions.

A pentad is a partition of the fifteen duads into five matchings.

**Lemma.** There are six pentads, every matching lies in two, and any two
edge-disjoint matchings lie in a unique pentad.

**Proof.** The complement in `K6` of the six-cycle `M union N` is a
triangular prism: two triangles and three cross edges. Its perfect matchings
are the three cross edges together, or one cross edge with the forced edge
in each triangle. The all-cross matching cannot belong to a factorization
of the prism, since the leftover two triangles have no perfect matching.
The other three matchings partition its edges uniquely. This gives the
unique pentad through `M,N`.

A fixed matching has eight edge-disjoint alternatives (all fourteen others
minus its six neighbors), and a pentad through it uses four. Hence it lies
in two pentads. Double counting gives `15*2/5=6` pentads. Two pentads cannot
share two matchings, by the uniqueness just proved. QED.

Let `P` be the **set of six pentads**, not the original label set. Define

    F(M) = {the two pentads containing M}.

The lemma makes `F` an injective map between two fifteen-element sets, hence
a bijection. Two matchings sharing a duad cannot occur in one pentad, so
their pentad pairs are disjoint. Two edge-disjoint matchings occur in exactly
one common pentad. Therefore

    M adjacent N in G  iff  F(M) disjoint F(N).

This is an explicit graph isomorphism `G -> KG(P,2)`. The receipt saves the
six pentads and all fifteen images. Its construction is canonical relative
to the original six-label set, but its target labels are pentads.

Relabelling `Omega` induces `Phi:S6 -> Sym(P)`. A permutation fixing all
pentads fixes every matching (each is recovered from its pair of pentads),
then every duad (intersect two matchings containing it), then every original
label (intersect two duads). Thus `Phi` is faithful; both groups have size
720, so it is an isomorphism. A transposition induces cycle type `2+2+2`
on pentads. After choosing an ordering of `P`, `Phi` is an **outer**
automorphism of `S6`, since inner automorphisms preserve cycle type.

The finite checker verifies the complete action and composition along every
adjacent-transposition generator. The graph map is equivariant for the
induced pentad action, not for the ordinary six-label action on duads.

## 4. The ordinary equivariance obstruction is structural

The stabilizer of a matching in the natural `S6` action is

    (S2)^3 semidirect S3, of order 8*6=48.

It can swap within pairs and permute pairs, and is transitive on all six
labels. A duad stabilizer is `S2 x S4`, also of order 48, but its natural
point orbits have sizes two and four. The subgroups cannot be conjugate
inside this `S6`, because conjugation preserves their point-orbit sizes.
A bijection of transitive `S6`-sets would require conjugate stabilizers.
Therefore no ordinary `S6`-equivariant duad/syntheme bijection exists.

An independent hostile gives the same obstruction: a transposition fixes
seven duads (`the swapped pair` plus six pairs among the four fixed labels)
but only three matchings (those containing the swapped pair). The checker
saves fixed counts for every conjugacy type and directly generates the
48-element wreath stabilizer.

This is a stronger diagnosis than “no canonical bijection was found.” It
rules out **every** bijection respecting the specified ordinary action,
while the preceding section constructs an actual relation-preserving
bijection after transporting the action.

## 5. What Richelot geometry contributes, and what it does not

**CITED interpretation.** Lagrangian two-torsion kernels define polarized
`(2,2)`-isogenies; see
[Florit--Smith, Definitions 2.1--2.4](https://eprint.iacr.org/2021/012.pdf).
Their graph uses kernel-orbit weights, not a simple-neighbor count, and
Lemma 3.2 gives

    |RA(A)| * w(phi_dual) = |RA(A')| * w(phi).

This is an automorphism-weight balance law for isogeny edges under those
hypotheses, not symmetry of the syntheme adjacency matrix.

[Florit--Smith's Atlas, §§1,2.1,3.1](https://arxiv.org/html/2101.00917v2)
interprets the fifteen matchings as quadratic splittings and gives weighted
out-degree fifteen in the superspecial graph over `F_(p^2)`. Over a general
base field some kernels may not be rational. The Atlas assumes
characteristic outside `{2,3,5}`.

[Atlas §4.14, Table 12](https://arxiv.org/html/2101.00917v2#S4.SS14)
identifies the Bolza normal form `y^2=x(x^4+1)`, its reduced `S4` action,
and the following kernel orbits:

| Kernel orbit | Weight | Stabilizer order | Generic codomain description |
|---|---:|---:|---|
| 1,4,8,9,12,13 | 6 | 4 | elliptic square `E^2`, `j(E)=8000` |
| 3,6,10,11 | 4 | 6 | Type IV |
| 2,7,14,15 | 4 | 6 | Type IV |
| 5 | 1 | 24 | loop |

At `p=7`, both weight-four classes become loops, and the product has
`j=8000=1728` in the residue field. At `p=13,29`, one Type-IV codomain
becomes Type V, depending on the chosen root of unity. Characteristics
`2,3,5` are excluded throughout this cited Atlas specialization. In
particular, this abstract group computation cannot determine the full
geometric automorphism group at `p=5`.

The finite replay validates the three published permutations, their
24-element closure, these four orbits, and faithful four-point actions.
It does not independently calculate curves or isogenies. Four edge classes
need not mean four distinct codomains. No pairing of these labels with
primes is supplied.

The useful bridge to THM-4012 is consequently concrete but modest: the
same Bolza Jacobian already present in its `j=8000` analysis has this cited
fifteen-kernel port set, with nonuniform orbit sizes. Its product neighbor
is reached by **six** kernel choices in one edge class. Neither the fifteen
nor the four is a new formula for a CM constant or a stable-model theorem.

## 6. Transport ledger and stopping boundary

| Source | Map / target | Preserved | Lost or extra sidecar |
|---|---|---|---|
| six-label matching | nonzero vectors of its Lagrangian | incidence, isotropy, ordinary relabelling | curve and field are not part of the finite carrier |
| matching | pair of pentads containing it | adjacency as disjointness | ordinary label action becomes pentad action; retain `Phi` |
| Bolza kernel | reduced-automorphism orbit | cited isogeny edge class | individual kernel; retain orbit member and weights |
| edge class | codomain isomorphism class | endpoint | isogeny and kernel data; target alone need not recover an edge |

There is **no** constructed map from the Ogg prime set to these kernels,
no privileged cyclic ordering of fifteen kernels, and no distinguished
reduction modulo fifteen. In particular, a count of fifteen, a modulus of
fifteen, and a fifteen-element prime set are different mathematical types.
No relation to the Marvin Ray Burns constant, no MRB defect formula, and no
identification of Ogg primes by residues modulo fifteen follows here.

The bounded lane stops at a complete finite carrier, an explicit twisted
relation-preserving map, and an authenticated published geometric quotient.
Further arithmetic would require an actual proposed map and preserved
predicate, not another equality of cardinalities.

## 7. Reproduction and certificate coverage

Run from any directory:

    python3 check.py --receipt receipt.json
    python3 -O check.py --receipt receipt-optimized.json
    cmp receipt.json receipt-optimized.json
    python3 audit.py

`check.py` uses only the standard library and explicit runtime checks, never
removable Python assertions. The primary receipt records 28,604 gates:
all fifteen matchings and Lagrangians, all 32 even subsets and their ordered
pairs, 45 flags, every graph vertex-pair, all six pentads, all 720
relabellings, the 3,600 generator-composition checks, two nonconjugate
48-element stabilizers, and the published 24-element Bolza action.

All input files have validated SHA-256 digests in `source-manifest.json`.
`audit.py` compares normal and optimized output and repeats both from a
real relocated copy containing only declared inputs; it also confirms that
a tampered input is rejected in optimized mode. The audit receipt records
actual outcomes. It does not claim independent algebro-geometric
verification of the cited source.
