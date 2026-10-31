# Reflection overlap, Hamiltonian residues, and the black degree

**Status: PROVED statements below, using the cited odd-cycle collection formula;
FINITE-EXACT checks are stated separately. No literature-priority claim.**
Session 130, 2026-10-31. Source checkout commit
`5431278adcd6269b5eca7756acc7a842938cd8cf`.

## 0. Inheritance and notation

The inherited objects are directed Hamiltonian paths of a finite tournament
`T`, counted as ordered vertex words, and directed simple odd cycles of length
at least three, counted modulo cyclic rotation **not** modulo reversal and
**not** by vertex support. Write `H(T)` for the former number.

The conflict graph `Omega(T)` has the latter cycles as vertices; two are
adjacent exactly when their supports meet. The inherited OCF is

`H(T)=I(Omega(T),2)=sum_k 2^k alpha_k`,                         (0.1)

where `alpha_k` counts vertex-disjoint `k`-sets of odd cycles. This is
THM-002-ocf.md, attributed to Grinberg--Stanley and restated by Irving--Omar.
The version-specific primary citation is Irving--Omar,
[*Revisiting The Redei-Berge Symmetric Functions via Matrix Algebra*, v3,
Corollary 20](https://arxiv.org/html/2412.10572v3). Numbering/title differ in
v1, where the same specialization is Corollary 19. No change to the
mathematical attribution follows from that version difference.

Fix an involutory anti-automorphism `rho` of `T`: `rho^2=id` and
`u->v iff rho(v)->rho(u)`. Define its involution on odd cycles by

`rho_*(C)=reverse(rho(C))`.

Let

* `q` count fixed cycles of `rho_*`;
* `i` count nonfixed two-element orbits `{C,rho_*C}` whose supports intersect;
* `j` count such orbits whose supports are disjoint.

These are counts of orbits, not counts of individual members.

An anti-automorphism fixes at most one vertex: two fixed vertices would
require their mutual arc to point both ways. Thus an anti-involution has no
fixed vertex at even order and one at odd order. Every invariant odd support
contains that fixed vertex. In particular any two fixed odd cycles intersect.

Closest inherited mechanisms:

* THM-791-H-companion-laws-to-the-transitivity-flow.md already proves the
  even-order self-converse stratum `H=1 (mod 4)` by odd-cycle pairing.
* THM-643-blue-black-line-parity-structure.md and THM-4731-tiling-line-
  metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-
  tribonacci-partner.md supply odd fibres, fixed blue fibres, and black twins.
* The 2026-10-20 odd-blue-spectrum/hollow-cubes packet supplies the precise
  centralizer quotient for blue path counts. That denominator must not be lost.
* Session 124's cut-response and parity-induction packet supplies an exact,
  elementary weighted-cut identity for vertex insertion; it supersedes using
  a false equality between Type-II and orphan totals.

## 1. The third binary digit is an intersecting-mirror-cycle charge

**Theorem 1.** For every finite tournament with an involutory
anti-automorphism,

`H(T) = 1 + 2q + 4i (mod 8)`.                                 (1.1)

**Proof.** We have `alpha_1=q+2i+2j`. The involution induced by `rho_*` on
independent two-element sets fixes exactly the disjoint mirror pairs:
a fixed pair either exchanges its two cycles, giving one of the `j` pairs,
or fixes them separately; the latter is impossible because fixed odd cycles
intersect. Nonfixed independent pairs occur in pairs, so `alpha_2=j (mod 2)`.
Truncating OCF gives

`H = 1+2alpha_1+4alpha_2 = 1+2q+4i+8j = 1+2q+4i (mod 8)`.

All terms from `alpha_3` onward are divisible by eight. QED.

For even order, `q=0`, hence

`H(T)=1+4i (mod 8)`.                                         (1.2)

This refines, but does not replace, inherited THM-791. It says that the next
bit records *overlap* of an odd cycle with its reflected cycle. Disjoint
mirror pairs contribute to two successive OCF coefficients and cancel
modulo eight. Along an even-order blue line with endpoints `T,U`,

`H(U)-H(T)=4(i(U)-i(T)) (mod 8)`.                             (1.3)

There is no all-order mod-eight constancy along blue lines: at order four
the transitive node and its all-backward-tile partner have counts `1,5` and
charges `0,1`. This is the smallest nontrivial even-order failure.

## 2. Central insertion: cycle charge equals a mirror-cut charge

Suppose `|T|` is odd, and let `z` be the fixed vertex of `rho`. Let `r_z`
count nonfixed mirror pairs of odd cycles containing `z` (both members then
contain `z`). Put `B=T-z`, with its induced fixed-point-free anti-involution.

**Theorem 2.**

`H(T)-H(B)=2q+4r_z (mod 8)`.                                 (2.1)

**Proof.** A nonfixed mirror pair avoiding `z` is exactly a mirror pair in
`B`, with unchanged overlap. A pair containing `z` always intersects.
Therefore `i(T)=i(B)+r_z`. Apply (1.1) to `T` and (1.2) to `B`. QED.

This also follows directly from the inherited OCF deletion identity
`H(T)-H(B)=2 sum_(C containing z) H(T-V(C))`: a fixed odd cycle has an
even-order self-converse complement whose `H` is `1 (mod 4)`, while a
nonfixed pair has equal odd complement counts and contributes `4 (mod 8)`.

Now define the canonical pair weights on `B`. Let `W_uv` count ordered
pairs of nonempty directed paths partitioning `V(B)`, the first ending at
`u` and the second starting at `v`. Put `b_uv=(W_uv+W_vu)/2`. Session 124
proves that `b_uv` is a positive integer and, for the new center's set of
out-neighbors `S`,

`H(T)-H(B)=2 K(S)`,
`K(S)=sum_(u in S,v outside S,u->v) b_uv`.                    (2.2)

Because `rho` is an anti-automorphism of the extension, `S` contains exactly
one member of each `rho`-pair. Reversing both cover paths and applying
`rho` exchanges the order of the two paths and proves

`W_uv=W_(rho(v),rho(u))`, hence `b_uv=b_(rho(u),rho(v))`.

The out-cut arcs are preserved by `(u,v)->(rho(v),rho(u))`. Its fixed arcs
are exactly `u->rho(u)` with `u in S`; all other arcs occur in equal-weight
pairs. Consequently

`K(S)=sum_(u in S,u->rho(u)) b_(u,rho(u)) (mod 2)`.

Combining (2.1) and (2.2) proves the **cycle--cut bridge**

`q = sum_(u in S,u->rho(u)) b_(u,rho(u)) (mod 2)`.             (2.3)

Thus the parity of *all fixed odd cycles* in the centered extension is
detectable from the opposite-pair two-path-cover weights of the smaller
tournament. This is not a triangle-only assertion or an equality of insertion
and orphan counts. The full cut weight additionally satisfies

`K(S)=q+2r_z (mod 4)`.                                      (2.4)

## 3. Automorphism-normalized black degree

For a self-converse isomorphism class, let `a=|Aut(T)|`, `b` be its blue
tiling fibre, and `d` its degree in the fundamental black multigraph
(the quotient of black twin lines; a loop contributes two).
The inherited fibre/twin laws give

`H/a=b+2d`.                                                (3.1)

Hence Theorem 1 gives the all-order refinement

`2ad+ab=1+2q+4i (mod 8)`,
`ad=(1-ab)/2+q+2i (mod 4)`.                                (3.2)

The half-integer-looking term is integral because both `a` and `b` are odd.
Modulo two this becomes

`d=(1-ab)/2+q (mod 2)`.                                   (3.3)

For even order the `q` term vanishes. The full formula (3.2) retains the
automorphism factor; silently replacing `a` by one is invalid.

To compute `b` from one chosen reflection, the inherited exact formula is
`b=F_rho/|C_Aut(T)(rho)|`, where `F_rho` counts Hamiltonian paths fixed by
`P->rho(reverse(P))`. It is not generally `b=F_rho`. All anti-involutions
are Aut-conjugate (Sylow conjugacy in the signed automorphism group), so
`q,i,j` are independent of the chosen anti-involution up to that conjugacy.
This is a recovered sidecar, not a new claim.

## 4. Order-five table

Codes below are the minimum integer upper-triangle adjacency mask over all
vertex relabellings; low bit is `(0,1)`, then `(0,2)`, etc. A code is not an
explorer tile mask. Converse classes `2,11` merge, as do `4,9`.

| Code | H | a | fibre | blue b | q | i | fundamental-black d |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 | 1 | 0 | 0 | 0 |
| 2,11 | 3 | 3 | 1 per class | 0 | — | — | 1 merged |
| 4,9 | 5 | 1 | 5 per class | 0 | — | — | 5 merged |
| 5 | 3 | 3 | 1 | 1 | 1 | 0 | 0 |
| 8 | 9 | 1 | 9 | 1 | 2 | 1 | 4 |
| 10 | 9 | 1 | 9 | 3 | 2 | 1 | 3 |
| 12 | 11 | 1 | 11 | 3 | 3 | 1 | 4 |
| 40 | 15 | 3 | 5 | 1 | 3 | 2 | 2 |
| 41 | 13 | 1 | 13 | 3 | 2 | 2 | 5 |
| 76 | 15 | 5 | 3 | 3 | 3 | 2 | 0 |

The row `76` is the regular cyclic tournament: all three tilings are blue,
but it is not a leaf. There are eight blue lines and twenty-four retained
black lines, paired into twelve fundamental black edges. The latter graph
has one loop each at codes `12,41`. `normal.json` records every edge.

At `n<=5`, two odd cycles cannot be disjoint, so `j=alpha_2=0` and the
congruence strengthens to the exact equality

`H(T)=1+2q+4i`                                              (4.1)

for a self-converse tournament with a chosen reflection. It still does not
determine `b` or `d`.

## 5. The precise failure boundary: equivariant OCF is not blue incidence

The tournaments with canonical codes `8` and `10` both have

`H=9, a=1, score=(1,1,2,3,3), c3=3, c5=1, q=2, i=1`.

In fact their entire **cycle-length-decorated conflict graphs with the
reflection action** are isomorphic: both are `K4`, with one fixed triangle,
one fixed five-cycle, and the two remaining triangles exchanged. But their
blue masses are `1` and `3`, and their fundamental-black degrees are `4`
and `3` respectively. The unique anti-involutions, in vertex order `0..4`,
are `(4,3,2,1,0)` and `(0,4,3,2,1)`.

The first tournament has the unique fixed Hamiltonian path `(4,3,2,1,0)`.
The second has three:

`(2,1,0,4,3)`, `(3,1,0,4,2)`, `(4,3,0,2,1)`.

The separate inheritance worker's packet independently verifies this
equivariant-conflict-graph obstruction. Our raw independent witness lists
the adjacency matrices, all nine paths for each tournament, every odd cycle,
its reflection, and the fixed paths in `equivariant_ocf_hostile.json`.

**Conclusion (finite-exact counterexample to a universal factorization):**
even adding scores and automorphism order to the equivariant OCF conflict
graph does not recover the blue fibre. This quotient retains pairwise cycle
intersection as a Boolean predicate but discards the vertex-incidence data
needed for fixed Hamiltonian paths. The symmetry-resolved OCF charge controls
Hamiltonian residues; it does not identify the colored fibre or retained line
pairing. At orders below five, exhaustive class enumeration gives no such
self-converse collision; the witness is smallest by order in that domain.

## 6. Verification and scope

The stdlib checker independently computes Hamiltonian counts by subset DP,
literal odd cycles by subset/permutation enumeration, pair covers by cutting
permutations, and the order-five class/line graph by exhaustive relabellings.
It checks all blue tilings through order eight, not a random sample:

| n | blue tilings | H=1 mod8 | H=3 mod8 | H=5 mod8 | H=7 mod8 |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | 1 | 0 | 0 |
| 4 | 4 | 1 | 0 | 3 | 0 |
| 5 | 16 | 5 | 4 | 3 | 4 |
| 6 | 64 | 28 | 0 | 36 | 0 |
| 7 | 512 | 124 | 110 | 132 | 146 |
| 8 | 4096 | 1852 | 0 | 2244 | 0 |

All 4,694 presentations pass Theorem 1; all 530 odd-order presentations
through seven additionally pass central deletion, full weighted-cut equality,
and (2.3). Order-five DP agrees with direct Hamiltonian-path permutations;
all twelve automorphism/fibre products, eight blue edges, twelve black-twin
edges and loop-aware degrees are checked directly. Normal and `python -O`
outputs agree byte for byte. No assertion statement carries correctness.

These are congruence and quotient-refinement theorems. They do not classify
all black connections, prove metagraph connectivity at arbitrary order,
identify Hamiltonian-path counts from scores, or establish a global forbidden
`7/21` claim. The latter is deliberately not used as a dependency.

The original spectrum goal has an explicit stopping boundary here: a
mod-eight obstruction alone cannot rule out `H=7`, because the pure-blue
regular cyclic five-tournament has `H=15=7 (mod 8)`. Nor can it rule out
`H=21`, because the strongly connected self-converse four-tournament has
`H=5=21 (mod 8)`. The new cycle overlap and central-cut coordinates are
structural refinements, not forbidden-spectrum closure.

Reproduce from any directory:

```sh
python3 check_hamiltonian_bridges.py --max-blue-n 8 > normal.json
python3 -O check_hamiltonian_bridges.py --max-blue-n 8 > optimized.json
cmp normal.json optimized.json
```
