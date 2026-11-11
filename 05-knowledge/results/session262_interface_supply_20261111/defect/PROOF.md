# Reusable defect reservoirs with a protected domination interface

**Session262, 2026-11-11. PROVED:** ordinal reservoir stabilization, fibrewise amplification, a simultaneous Boolean coefficient bank, and finite initial objects. **CONDITIONAL:** their use by a future legal correction operation. **OPEN:** Chan–Pak, arbitrary exact target synthesis, and a strong-tournament version of stabilization. No historical novelty claim is made.

## 1. Inheritance and the turn in viewpoint

Take the universal inequality and exact nonnegative fibres in
[Session261 EXTENSION4](../../session261_creative_chanpak_20261111/bridges/EXTENSION4.md) as inputs, not objects to re-audit. For a marked tournament `(T,z)`, put `U=T−z`, `A=N⁺(z)`, `B=N⁻(z)`, and let `G` have edges `ab` when `a∈A,b∈B,a→b`. Write

    a=H(U), d=d(G), Δ=H(T)−H(U), σ=Δ−d+1 ≥ 0.

The closest proved mechanism is the positive marked-cycle expansion

    Δ = 2 Σ_W p(W) H(U−W),
    σ_W = 2p(W)H(U−W) − ε(W)2^κ(W) ≥ 0,
    σ = Σ_W σ_W.

Here W ranges over nonempty even vertex sets, p counts A-to-B Hamilton paths on W, and ε indicates that every induced G-component has even order. The least-used sidecar is **the complementary factor H(U−W)**. Rather than making more paths on W or enlarging its component factorial, we enlarge the complementary factor while leaving W and every old edge untouched.

The hostile boundary is the empty cut graph: d=1 and there is no marked cycle through z. The corrected near miss is physical normalization from THM-5107: dividing a realizable count need not give a realizable source. Our construction never divides a tournament count. It adds actual vertices. The exchange-sidecar warning in `multiplicand_saturation_chanpak_20261111/exchange/DUAL_DEFECT.md` remains: equal scalar defects do not identify cycle fibres with path-cut incidence components.

Concept board: (i) protected cut graph, (ii) complementary OCF reservoir, (iii) ordered product of actual paths, (iv) Boolean marked derivatives, (v) target anchors versus selectable coefficients. Anchor: reusable positive source capacity. Niche: independent marked reservoirs, comparable to a tensor-product state space. Wildcard: a finite Boolean cube whose mixed derivatives are physical deletion experiments, even though its coefficients are not themselves arbitrary tournament counts.

## 2. The actual stabilization operation

Let L and R be arbitrary finite tournaments, possibly empty. Form

    T' = L → T → R,

with all edges from an earlier displayed block to a later one, retaining every internal arc and the old mark z. Put `h=H(L)H(R)` and `I=V(L)⊔V(R)`.

### Theorem A: source-preserving amplification (PROVED)

1. The new cut graph is exactly `G' = G ⊔ |I|K₁`.
2. `D ↦ D∪I` is a bijection between dominating sets of G and G'. In particular d'=d.
3. Every original marked vertex, arc, subset, and label is retained literally. After stabilization the same operation remains available, without finding a new port.
4. `H(T')=hH(T)`, `H(T'−z)=ha`, `Δ'=hΔ`, and

       σ' = hσ + (h−1)(d−1).

5. For each old even W, the new individual fibre is

       σ'_W = hσ_W + (h−1)ε(W)2^κ(W).

   Every new W meeting I has p'(W)=0 and ε'(W)=0.

**Proof.** All L vertices lie in B', all R vertices in A'. No A'→B' graph edge touches either exterior block, by the ordinal directions. Thus they are isolated graph vertices, which every dominating set must contain. A tournament path cannot leave a later block for an earlier block; a spanning path consequently concatenates one spanning path per block. Restricting a path to its blocks is the inverse map, proving the two Hamilton-count products at the level of actual paths. An A'-to-B' path cannot use L (it cannot start there or return there), nor R (it cannot leave R to end in B'). The old p(W), ε,κ are unchanged, while `H((T'−z)−W)=hH(U−W)`. Substitution proves the remaining formulas. ∎

**Why this is stronger than a scalar budget.** The growth can be addressed by the *same exact W* as before. No redistribution theorem, selection of a new source graph, or search for a new marked vertex is required. Every old nonempty admissible W acquires `(h−1)2^κ(W)` additional states on top of h copies of its old slack. Thus a source whose original inequality is sharp is not a dead end: every shell with h=H(L)H(R)>1 makes it a reservoir. A nonempty transitive shell has h=1 and adds no capacity.

**Guard.** The graph is preserved up to forced isolated vertices, not as a literally unchanged vertex set. The tournament becomes decomposable when both core and shell are nonempty. Neither strong connectivity nor a preassigned global numeric target is preserved. Other marked cut graphs involving exterior vertices are not asserted unchanged.

## 3. Uniform supply with finite initial objects

Let C be the directed triangle on labelled vertices c₀,c₁,c₂, with cyclic arcs. Its three paths are `(c₀,c₁,c₂)`, `(c₁,c₂,c₀)`, `(c₂,c₀,c₁)`. Let `R_k` be the ordinal sum of k copies of C, and set `S_k(T)=T→R_k`. The empty shell R₀ has H=1; each step costs exactly three vertices and retains the old z.

### Theorem B: unbounded reusable capacity from every nontrivial source (PROVED)

For every finite `(T,z)` with `d(G)>1` and every integer q≥0, choose any k satisfying

    3^k ≥ 1 + q/(d−1).

Then `S_k(T)` has the same source domination count and at least q units of nonnegative quantitative defect. In fact

    (a,d,Δ,σ) ↦ (3^k a,d,3^k Δ,3^k σ+(3^k−1)(d−1)).

The construction is uniform, uses `|T|+3k` vertices, and can be iterated from its own output. A sidecar containing the block labels and the original cut data suffices to recover the input interface by restriction. The number q is a guaranteed capacity, **not** a promise that an arbitrary q of those states can be deleted by an arc operation.

The hypothesis is exactly `E(G)≠∅`: an edge gives an admissible two-vertex set and d≥3; an edgeless graph has the unique dominating set V. In the latter case Δ=0: any nonempty even A-to-B path would have some A→B transition, hence a graph edge. Stabilization cannot generate a z-cycle from no z-cycle. This is a genuine initial-case guard, not a finite-search omission.

**Small initial objects.**

- Zero source: one marked vertex T={z}. `(a,d,Δ,σ)=(1,1,0,0)`. Any triangle shell leaves marked capacity zero, despite creating other Hamilton paths.
- First positive source: T=C₃, z=c₀. Then U is a directed edge, G=K₂, and `(a,d,Δ,σ)=(1,3,2,0)`.
- One recharge: `C₃→C₃`, mark in the first block. `(a,d,Δ,σ)=(3,3,6,4)`.
- Two recharges: three ordinal triangles. `(a,d,Δ,σ)=(9,3,18,16)`.
- In general the first positive source gives `(3^k,3,2·3^k,2(3^k−1))` and H=3^(k+1).

These are finite initialization *rules*, not a claim that finitely many seeds cover all source counts d. In particular the operation never changes d, so it cannot fix a missing domination-count source.

## 4. Explicit protected states rather than an unexplained excess integer

Use the inherited odd-cycle formula as a finite weighted set: each disjoint odd-cycle packing has a binary label for each cycle, hence weight 2^number-of-cycles. Let `N_W(T,z)` be its states whose unique z-cycle uses W∪{z}. Thus

    |N_W|=2p(W)H(U−W).

The source fibre `F_W` is the set of binary labels on the κ(W) components of G[W], when W is nonempty and admissible; otherwise F_W is empty. Its cardinality is ε2^κ. There is an explicit injection `j_W:F_W→N_W`:

1. Choose the lexicographically least A-to-B spanning path in each component, as supplied by the inherited endpoint lemma.
2. Enumerate the pairs `(cycle bit, permutation of components)` lexicographically.
3. Assign the lexicographically ordered binary component words to the first 2^κ pairs, possible because `2κ!≥2^κ`.
4. Concatenate the component paths, close through z, attach that cycle bit, and use the empty complementary packing.

Labels recover the chosen component order, so the map is injective. It is deliberately a map from the **positive even-component source**, not an unproved bijection from individual dominating sets. Their equal total cardinalities are the inherited positive formula.

The triangle shell has an explicit three-state OCF alphabet: `0=empty packing`, `+=triangle with first bit`, `−=triangle with second bit`. Ordinal blocks cannot share a directed cycle. Therefore

    N_W(S_k(T),z) = N_W(T,z) × {0,+,−}^k

as an actual labelled packing-set equality. Protect `j_W(F_W)×{0}^k`. Its complement decomposes disjointly into

    (N_W\j_W(F_W)) × {0,+,−}^k
    ⊔ j_W(F_W) × ({0,+,−}^k\{0}^k).

This proves the fibrewise slack formula by a literal disjoint union, and identifies the reserve address `(W, old state, shell word)`. The all-zero word preserves the chosen original source injection. Recharging appends a new ternary coordinate; it does not invalidate an old address. This is source **reusability under repeated amplification**, not a claim that spending states preserves the source.

The path factorization of Theorem A is independently literal. No unconstructed natural bijection from OCF states to Hamilton paths is required or asserted here.

## 5. Simultaneous Boolean source bank

Mark r_i=c₀ in every triangle shell. Deleting r_i turns that shell into its directed edge, with one Hamilton path. For J⊆[k], let T_J retain all shell vertices except r_i for i∉J; always retain the entire old core. Then

    H(T_J)=H(T)·3^|J|,
    H(T_J−z)=a·3^|J|,
    Δ_z(T_J)=Δ·3^|J|.

The unique multilinear polynomial agreeing with these Δ-values on the Boolean cube is

    F(x₁,…,x_k)=Δ ∏_i(1+2x_i).

### Theorem C: mixed coefficient bank (PROVED)

For every S⊆[k],

    Σ_{J⊆S} (−1)^(|S|−|J|) Δ_z(T_J) = 2^|S| Δ.

Here marks outside S are deleted in every term. The coefficient has a positive object interpretation: the original z-containing packing state, a nonempty labelled triangle in each shell indexed by S, and empty packing in the other shells. For a fixed original admissible W, the protected-source contribution to that coefficient has size `2^|S| 2^κ(W)`. Distinct S give disjoint supports in the full reservoir, so simultaneous demand is not double counted.

Every r_i remains marked in the undeleted final object. Its own cut graph is K₂ plus forced isolates, so d=3 at every such port, independently of k and the old core. Thus one actual tournament supplies arbitrarily many locally identical initialized domination ports, with a retained global interaction law rather than an unsupported independence assertion.

**Important distinction.** Mixed coefficients are positive packing subsets and alternating combinations of physical deletion experiments. They are not automatically counts of standalone tournaments, nor do they authorize subtraction of Hamilton paths. The whole physical shell has count 3^k; the support-exact coefficient 2^|S| is even when S is nonempty. Treating it as a tournament count would violate parity immediately.

## 6. What a Chan–Pak constructor can and cannot take from this

The inherited offset carriers retain `(D,U)` and require output-side data to recurse; an H-count alone does not close that interface. Our operation instead closes `(a,d,Δ,σ)` under a genuine graph operation, with all original labels and a recoverable source. This is a different interface, not a substitute proof for their missing `(D,U)` supply theorem.

A precise conditional transfer statement is available. Suppose a future legal operation accepts a protected-source interface and a finite list of distinct reserve addresses, changes the tournament according to a proved output-count formula, and preserves the core labels and a stabilization-compatible port. Then every bounded address-capacity requirement on a source with d>1 can be fulfilled before invoking it, with O(log q) additional vertices. If it requires several disjoint classes of addresses, give each demand a disjoint shell-word range; if it requires specific Boolean support S, use Theorem C and add shell indices as needed.

This implication is **CONDITIONAL on the operation's compatibility**. The current theorem supplies no operation that realizes H(T)−2m, no arbitrary residue control, and no initial source for every prescribed d. The graph output H=3^kH(T) is a specified target anchor, and the reserve count is exact, but an arbitrary desired target N is not anchored by selecting k unless its arithmetic and interface constraints already match.

A useful next open problem is therefore narrower than “find more slack”: construct a local operation that reads a reserved `(W,bit,shell-support)` address, alters a known count by a controlled amount, and returns the *same protected core* with a recharge port. The theorem here discharges the supply side of that proposed constructor under its explicit compatibility condition. It also separates three obstructions that should not be conflated: no nontrivial source d, insufficient reserve capacity, and no legal way to spend an available reserve. Only the middle one is closed.

## 7. Three random historical sources, with actual transfer contracts

Sampling seed **26211031**; pool is `.md` paths under `01-canon/` or `05-knowledge/results/`, in `git ls-tree -r --name-only HEAD` order. Head, pool size and exact paths are in `receipt.json`. All three were read; this is a creative connection pass, not a re-verification of them.

1. `05-knowledge/results/session17_prime_shapes_20261020/F5_PATH_PAIR_ICOSAHEDRON.md`. Source: a marked path product and its sign-equivariant finite-field labels. Target: the marked ternary shell product. Map: send each shell packing to 0,+,− and retain the actual block label. Preserved: product coordinates and a distinguished zero state. Lost: internal arc-incidence geometry under a bare word; restored by the labelled block decoder. The marking-hostile there prompted keeping the explicit alphabet and refusing to call arbitrary word permutations tournament automorphisms. Cheapest test: restrictions of concatenated Hamilton paths recover each block path uniquely.
2. `05-knowledge/results/session40_marked_ports_square_sum_20261023/audit-third/THIRD_AUDIT.md`. Source: ordered color-word collisions versus owner-labelled reconstruction. Target: reserve support S versus its actual ternary choices and old W. Map: shell word↦support. Preserved: Boolean grade; lost: 2^|S| independent cycle bits, and all old source ownership. Sidecar: `(W,old state,word)`. Cheapest hostile: one shell has two nonzero states with identical support. This becomes the exact coefficient 2, rather than an accidental double count.
3. `01-canon/theorems/THM-5107-physical-graft-content-excludes-cofactor-seven.md`. Source: physical replacement only when a normalized gadget is realizable. Target: distinction between a Boolean coefficient and a physical shell. Map: deletion cube→multilinear coefficient, preserving an exact integer but not standalone tournament realizability. Sidecar: the entire cube of physical tournaments. Cheapest hostile: coefficient2 cannot be a tournament H-count although its triangle shell is physical with H=3. This prevents interpreting abundant even reserve mass as a subtraction gadget.

The board changes coherently: the field-product source suggests explicit reusable labels; the owner-word source forces support and sign to stay separate; physical normalization blocks the false final step from selected coefficients to arbitrary H-values.

## 8. Finite-exact controls and reproduction

Run `python3 05-knowledge/results/session262_interface_supply_20261111/defect/probe.py` from the repository root. The new checker uses its own subset Hamilton-path DP and literal dominating-set enumeration. Universe: all labelled rooted tournaments on 1..4 vertices, each with k=0,1,2 shells, and every shell Boolean subset. It tests the actual new graph operation, cut isolation, source-count preservation, scaling, slack, and every mixed coefficient. It does not recensus the inherited domination inequality. Positive controls are the triangle initial objects; the hostile control is an empty cut with nontrivial shells. The receipt records exact counts and source sampling.

No external paper is invoked beyond the inherited OCF mechanism. The all-order conclusions are Theorems A–C and their explicit maps, not extrapolations from these controls.
