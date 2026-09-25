# Non-pair tournament decoders: intrinsic blocks, switching, and triangle ports

Status: **PROVED** for the elementary propositions below; **FINITE-EXACT** for the explicitly bounded controls. No universal Collatz root generator is claimed. Date: 2026-09-25.

## 1. Inheritance and target

The source object is the tournament

\[
F_Q(H)=Q[H,H,H,1],\qquad |H|=q\geq3\text{ odd},\quad H\text{ regular},\quad |F_Q(H)|=3q+1.
\]

Its vertices are the three actual copies of the seed and one singleton. The observable is the original directed arc between each distinct pair; there are no ties. No orientation gauge is silently identified. The target is a reversible structural description, and, where stated, directed reachability to a marked vertex. Arithmetic root certifiability additionally requires legal labelled Collatz edges and their certificates.

Closest proved mechanisms and boundaries:

- [decoder_halving_20260925.md](decoder_halving_20260925.md) supplies exact pair halving for the prescribed family `E(2^k q)=H_q[TT_(2^k)]`, and the hostile absence of two-vertex modules in `F_Q(H)`.
- [decoder_pair_repair_20260925.md](decoder_pair_repair_20260925.md) measures the individual arc reversals needed to repair that obstruction. The present switching result concerns the narrower operation of reversing a cut.
- [collatz_tournament_core_20260925.md](collatz_tournament_core_20260925.md) separates a four-block quotient from an arithmetic transition and retains the marked root.
- [THM-1960 — tournaments compose from regular seeds; spectral substitution law](../../01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md) is the relevant modular sidecar. Strong connectivity and modular primality are different properties.
- [THM-1415 — switching is the canonical star quotient](../../01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md) supplies the switching viewpoint. The invariant used here is derived directly below.
- [THM-3121 — path-cover walk-content substitution kernel](../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md) already proves the exact path-cover composition law. Its proved result supersedes the old open-kernel language in [THM-1975](../../01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md).

The live concept board is: intrinsic modules; cut gauges; strong-set contraction; marked endpoints; path-cover profiles; legal arithmetic proof edges. Each decoder below states what survives and which attachment data it needs.

## 2. The three seed blocks are intrinsically recoverable

**D1 — strong-block theorem (PROVED).** In every `F_Q(H)` above, the three original copies of `H` are strong modules. They are precisely the modules of cardinality `q`. Therefore the unmarked tournament determines their unordered collection and the remaining singleton, without a chosen matching.

Here a module is a set seen uniformly by every outside vertex. A strong module overlaps no other module: two intersecting modules must be nested.

First, a regular tournament of order `q>=3` is strongly connected. Otherwise a source strongly connected component `S` of size `s<q` dominates its complement. The average full outdegree in `S` is

\[
(s-1)/2+(q-s)>(q-1)/2,
\]

contradicting regularity. Now suppose a module `M` overlaps a seed block `H_i` without either containing the other. Choose `v` in `M\H_i`. For every `z` in `H_i\M`, the module condition forces the relation of `z` to every element of `M intersect H_i` to equal its relation to `v`. Since `v` is outside the block, its relation to all of `H_i` is uniform. Hence all arcs between the two nonempty parts `M intersect H_i` and `H_i\M` point the same way. This contradicts strong connectivity of `H_i`.

A module of size `q` must meet some seed block. Strongness forces nesting; equal cardinality then forces equality with that block. This proves the claimed uniqueness.

The decoder computes `q=(N-1)/3`, finds the three size-`q` modules, and returns their induced tournaments, their embeddings, and the marked four-vertex quotient. These data reconstruct every original arc exactly. The three blocks may be permuted by automorphisms; the construction does not select one of them equivariantly. A general recognizer for *identical* copies must additionally verify the three induced tournaments are isomorphic; the experiment tests the stated substitution family directly.

The scalar map `N -> (N-1)/3` is a **guarded reverse construction arrow**, with `N=3q+1` and `N=4 mod 6`. It is not a forward Collatz step, and does not supply `N/2`. A supplied ordinary Collatz certificate for `q>=3` begins with the actual edge `q -> 3q+1=N`; truncating that edge certifies `N`. Conversely a certificate for `N` can be prefixed by that edge to certify `q`. The structural decoder generates neither certificate by itself, and proves no compatibility `F_Q(H_q)=E(3q+1)`.

This reduction is genuinely useful before pair repair: it reveals the seed and the four-core with zero arc changes. Its stopping boundary is also exact. The regular cyclic tournament on five vertices has no proper nontrivial module, so module-only recursion can stop there. This says nothing against reductions that retain nonuniform ports.

## 3. Cut switching cannot manufacture a pair module

**D2 — switching obstruction (PROVED).** No tournament cut-switching-equivalent to any `F_Q(H)` above has a two-vertex module.

Write `S_uv=+1` for `u->v`, and `-1` for the reverse. A cut switch is

\[
S'_{uv}=d_u S_{uv}d_v,\qquad d_u\in\{\pm1\}.
\]

The pair `{u,v}` is a module after some switch exactly when the products `S_uw S_vw`, for all outside vertices `w`, are constant. Indeed switching multiplies every one of these products by the same `d_u d_v`. In the original tournament this is precisely the condition that the two outside rows are either identical or opposite.

If `u,v` lie in one seed, regularity excludes identical outside rows *within the seed*: the internal arc between them would make their full degrees differ by one. Thus some other seed vertex gives product `-1`. Every vertex outside that seed gives product `+1` by the substitution definition. If the pair crosses seed blocks, or contains the singleton, take `u` in a seed and `v` outside it. Vertex `u` has both an internal win and an internal loss, while `v` sees the seed uniformly, so the products again take both signs. Every case fails the criterion.

The missing information is not a choice of cut gauge. More general individual arc reversals remain available and are measured by the inherited pair-repair tables.

## 4. A cyclic four-core has two different switching types

**D3 — four-core gauge invariant (PROVED).** For a signed four-vertex tournament define

\[
P=S_{01}S_{23}-S_{02}S_{13}+S_{03}S_{12}.
\]

Then `|P|` is either `1` or `3` and is invariant under cut switching. Every summand is multiplied by `d_0 d_1 d_2 d_3`. It is also invariant under relabelling, since the Pfaffian changes only by the permutation sign.

Fix a marked root `r`. The unique cut gauge with `d_r=1` and `d_v=S_rv` makes the root a source. The other three vertices form either a transitive triple or a directed triangle. With the root labelled `0`, the normalized expression is `S_23-S_13+S_12`, whose absolute value is `3` exactly for a directed triangle. Consequently a four-core is cut-switchable to the transitive tournament exactly when `|P|=1`.

Among all 64 labelled four-cores:

| Sorted outdegree sequence | Count | `|P|` |
|---|---:|---:|
| `(0,1,2,3)` | 24 | 1 |
| `(1,1,2,2)` | 24 | 1 |
| `(0,2,2,2)` | 8 | 3 |
| `(1,1,1,3)` | 8 | 3 |

Thus a strongly connected cyclic four-core can become transitive under a cut switch. A source or sink attached to a directed triangle cannot. The mere presence of a cycle does not distinguish the two cases.

The three gauge bits and the marked root must be retained to restore the original object. Switching does not preserve reachability to that root or Hamiltonian path counts; switching the root cut turns a source into a sink. It therefore cannot silently stand in for an arithmetic proof edge.

## 5. Directed triangle contraction with attachment ports

**D4 — exact reachability quotient (PROVED).** In any tournament, contract a directed triangle `C` to a single vertex. For each remaining pair of quotient blocks, include an arc in a direction whenever any original arc has that direction. Both directions may occur. The result is a semicomplete directed graph: every distinct pair has at least one arc, and a pair may have a digon.

For every original ordered pair `(u,v)`, reachability from `u` to `v` is equivalent to reachability between their quotient blocks. Projection proves one implication. For the other, choose an original witness for each quotient arc; successive entry and exit ports within `C` are connected by a directed triangle path of length at most two. Singleton blocks require no internal movement. Concatenation gives an original walk and hence reachability. If the designated root belongs to `C`, its actual vertex is retained as the terminal port.

For complete reconstruction, retain the three crossing bits for each outside vertex, the cyclic internal triangle, the outside tournament, and the vertex embeddings. A mixed crossing pattern has six possibilities: merely recording its digon loses which triangle vertices supply which direction. Retaining these bits restores every arc, beyond the reachability guarantee of the coarse quotient.

Every cyclic four-core has a directed triangle. Contracting a marked triangle gives:

- a two-vertex digon for every strongly connected four-core;
- a single directed arc for a source or sink over a directed triangle.

The complete labelled controls contain 48 marked triangles of the first type and 16 of the second. This is an exact `4 -> 2` reduction of a cyclic core when the category permits digons. A two-vertex tournament cannot preserve both directions of reachability in a digon; choosing one direction necessarily loses the other.

Triangle contraction reduces order by two, so triangle contractions alone preserve parity. In particular, they cannot implement `10 -> 5` halving. They preserve reachability, not Hamiltonian paths, and they do not automatically turn structural arcs into legal arithmetic steps. To transport a Collatz root certificate, retained witness arcs and internal port paths must carry and verify their actual arithmetic guards. A directed triangle itself precludes a strictly decreasing integer rank on all its directed arcs; a proof DAG can select and certify an appropriate route, but vertex-count decrease of the compressed graph does not prove integer-orbit termination.

## 6. Hamiltonian structure survives through a richer inherited profile

This section applies the already proved [THM-3121 path-cover walk-content kernel](../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md), rather than reopening its old near miss. Let `pc_H(c)` count unordered spanning directed path covers with `c` components, and set `F_H(c)=c! pc_H(c)`. Let `W_Q(c)` count directed words in the loop-free quotient `Q` having content vector `c`.

\[
\operatorname{ham}(Q[H_i])
=\sum_{\mathbf c}W_Q(\mathbf c)\prod_i F_{H_i}(c_i).
\]

A Hamiltonian path splits uniquely into maximal runs within blocks. Those runs form an unordered spanning path cover inside each block; ordering its components at the occurrences of the block label contributes `c_i!`. The quotient word, covers, and assignments reconstruct the path. Counts alone do not constitute an invertible decoder.

If the root is the singleton block `r`, replacing `W_Q` by the number of quotient words **ending at `r`** counts Hamiltonian paths ending at that root exactly. A root inside a larger block requires an endpoint-refined path-cover profile.

For `H=C3`, `F_H(1),F_H(2),F_H(3)=(3,6,6)`. The experiment compares this macro formula with independent subset dynamic programming on all 64 ten-vertex substitutions, both for all Hamiltonian paths and those ending at the singleton. The latter count can be zero even when the former is 3159. Thus the marked endpoint is essential to any root-directed use.

The odd-cycle formula of [THM-002 — OCF](../../01-canon/theorems/THM-002-ocf.md) gives an independent structural explanation of odd Hamiltonian counts, and specializes on four vertices to `ham(Q)=1+2 c3(Q)`. The primary source is [Irving–Omar, Corollary 20](https://arxiv.org/html/2412.10572v3), which attributes the odd-cycle form to Grinberg–Stanley. Its positivity guarantees some Hamiltonian path; it does not specify the root endpoint or an arithmetic route.

## 7. Exact controls and remaining obligation

Reproduce from the repository root:

```text
python 04-computation/experiments/duck_tournament_20260925.py
python -O 04-computation/experiments/duck_tournament_20260925.py
```

The script uses explicit checks rather than optimization-disabled assertions. Its full output is [duck_tournament_20260925.out](duck_tournament_20260925.out).

The bounded universe comprises all 64 labelled four-cores; their substitutions by cyclic regular seeds of orders 3 and 5; all subsets of each order-10 object to check strongness; all 512 cut gauges modulo global sign for every order-10 object; all marked directed triangles in four-cores; and all Hamiltonian paths via subset dynamic programming on each order-10 object. Every decoder tests reconstruction or its exact stated predicate. The regular prime order-five seed and the order-parity obstruction are hostile controls.

An independent read-only agent audit checked D1–D4, the arithmetic scope, and a full `python -O` replay: **PASS**. Normal and optimized runs give identical frozen output; all local links in this note resolve.

The constructive survivors are: recover the three seeds intrinsically; retain a finite cut gauge where the four-core permits it; contract strongly connected triangles with explicit attachment ports; and transport Hamiltonian data using the proved full profile. Their sharp common limit is that structural compression alone supplies neither a legal Collatz transition nor a well-founded arithmetic rank. The next useful object is a port-labelled proof DAG whose edges are checked arithmetic certificates, with these exact decoders used to organize its data.
