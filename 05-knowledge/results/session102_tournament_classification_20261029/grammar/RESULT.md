# Signed modular cycle indices and the single converse branch

**Status:** PROVED elementary consequences, relative to the **CITED** classical canonical strong-module decomposition. **FINITE-EXACT:** independently implemented direct-permutation/attachment checks in `check.py`, with `results.json`. No historical novelty claim. No complete classification of prime tournaments or tournament metagraphs is claimed.

## Inheritance, target and information ledger

Verified against source commit `3b6e3ef82382da222fc59aa02914eea72d5a105b`, not just historical scratch files:

- `01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md`, including its correction banner.
- `05-knowledge/results/tournament_composition_classification_20261026/CLASSIFICATION.md`, §2, and `fibre/RESULT.md`, F1–F3.
- `05-knowledge/results/cyclictriple_composition_20261026/README.md`, §1–2.

The inherited exact recursion retains the oriented child-isotype stabilizer. Its hostile is `C3[A,A,A]` versus `C3[A,A,A*]`, `A=C3→1`: automorphism orders 81 and 27 despite identical labelled cyclic triples. The original automorphism-invariance claim is REFUTED, not a dependency. The least-used sidecar here is the permutation cycle index, including its negative/converse coset. Classical modular decomposition is imported, not reproved. The results below are an explicit signed extension of the inherited recursion, not a replacement for an atom census.

Concept board: strong-module tree; child isotype/anti-isotype colours; signed permutation cycles; binary attachment orbits; converse-fixed roots. Source object: a tournament with its intrinsic complete oriented relation, no ties. The map to signed indices preserves permutation-cycle statistics and specified attachment orbit counts; it loses the full group action and tournament isotype. The full decorated tree remains the classifier. Global converse is not the free-arc complement of a Hamiltonian path presentation; the results here concern the former only.

## G1. Exact signed cycle-index recursion

For a nonempty finite tournament `T`, let `Gε(T)` be its vertex permutations preserving arcs when `ε=+1`, reversing every arc when `ε=−1`. Write `a(T)=|G+(T)|`. The negative set is empty or a coset of `G+(T)`. Define

    Z_T^ε(x_1,x_2,...) = (1/a(T)) Σ_{g∈Gε(T)} ∏_{cycles D of g} x_|D|.

Both indices of a leaf are `x_1`. For a canonical internal node `T=Q[T_1,...,T_k]`, define

    Hε = {σ∈Gε(Q): T_i is ε-isomorphic to T_σ(i), for every i}.

An ε-isomorphism reverses arcs precisely when ε is negative. Then

    a(T) = |H+| ∏_i a(T_i),

    Z_T^ε = 1/|H+| Σ_{σ∈Hε} ∏_{C∈cycles(σ)}
                         Z_{T_i}^{ε^|C|}(x_|C|, x_2|C|, x_3|C|,...),   (G1)

where any representative `i∈C` may be used. This holds for prime and linear nodes alike. At a linear node `H+` consists only of the identity; `H−` consists of the order-reversing child permutation if every opposite child is dual-isomorphic, and is empty otherwise. The root's linear children are its ordered strong components.

**Proof.** Every signed automorphism preserves the canonical strong modules and induces a permutation in `Hε`; conversely every such quotient map lifts by choosing signed isomorphisms along its child arrows. Fix a cycle `C=(i_0,...,i_(l−1))` of σ. The composite around C is in `G_{ε^l}(T_i0)`. All child automorphism orders on C equal a common `a`. For each specified composite, the first `l−1` signed isomorphisms may be chosen independently, giving `a^(l−1)` choices; the last is then uniquely forced. A composite vertex-cycle of length d becomes a lifted cycle of length ld. Therefore the total cycle polynomial for this quotient cycle is

    a^l Z_{T_i0}^{ε^l}(x_l,x_2l,...).

Different quotient cycles are independent. Multiply, sum over σ, and divide by `|H+|∏a(T_i)`. This proves (G1). The positive order formula is the inherited stabilizer formula. Choosing a different representative around a quotient cycle conjugates the composite by a signed isomorphism, preserving cycle lengths. ∎

**Boundary.** For arbitrary noncanonical substitution partitions this computes the partition-preserving signed maps, not necessarily all signed automorphisms. Canonical maximal proper strong modules are essential. Child group orders alone are insufficient even for the positive index.

## G2. Sign changes follow at most one recursive branch

Every positive tournament automorphism has only odd-length cycles. Every negative automorphism has cycles of lengths `2d` with d odd, together with at most one fixed point; there are no other cycles.

**Proof.** A positive automorphism with an even cycle has a power interchanging its antipodal vertices and would reverse that pair's arc. For a negative f, its square is positive, ruling out cycle lengths divisible by four. If an odd cycle of f had length l>1, `f^l` would reverse arcs while fixing at least two vertices of that cycle, impossible. The same pair argument rules out two fixed points. ∎

Consequently, in every negative term of (G1), each nontrivial quotient cycle is even and contributes a **positive** child index. Only a possible unique fixed child contributes a negative index. Thus a specified anti-automorphism has a single possible chain of negative states through the modular tree; the off-chain factors are positive. This is a concrete answer to how recursive patterns behave under sign change, rather than an assertion that all branches independently dualize.

At each level, all children off that chain contribute an even number of vertices. Hence the parity of the child on the chain equals the parity of the root. If the tournament has odd order, the chain terminates at its unique fixed leaf. If it has even order, the chain must terminate at a quotient node without a fixed child. The chain depends on the chosen anti-map and is **not** an extra canonical invariant of an unmarked tournament.

## G3. Exact marked one-vertex extension and converse counts

For `S⊆V(T)`, form `E(T,S)` by adding a distinguished apex v, with `v→u` iff `u∈S`. Isomorphisms must preserve the apex. The number of such extensions up to marked isomorphism, for fixed base isotype `[T]`, is

    A(T) = Z_T^+(2,2,2,...).                               (G3a)

The number of these marked isomorphism classes which are self-converse **by an anti-isomorphism fixing the apex** is

    R(T) = Z_T^−(0,2,0,2,0,2,...).                        (G3b)

In particular `R(T)=0` for every odd-order base, and for any non-self-converse base.

**Proof.** Positive maps act on subsets by `S↦g(S)`; a permutation fixes `2^(number of cycles)` subsets. Burnside proves (G3a). If the negative coset is empty, (G3b) is zero. Otherwise choose h negative. Converse acts on positive subset-orbits by `S↦V\h(S)`. Its square is the positive action of `h²`, so it is an involution on the orbit set. The number of its fixed orbits is

    (1/a(T)) Σ_{f∈G−(T)} #{S : S=V\f(S)}.

For completeness, a stable positive orbit O contributes exactly `a(T)` pairs `(f,S)` to this sum: for each S∈O, there are `a(T)/|O|` maps from the negative coset satisfying the equality. An unstable orbit contributes zero. Along a permutation cycle the subset indicator must alternate; there are zero choices for an odd cycle and two choices for an even cycle. This proves (G3b). On odd n every negative permutation has a fixed point by G2, so its contribution vanishes. ∎

These are **marked** counts. Removing the apex mark from arbitrary extensions can identify inequivalent attachment orbits; no general unpointed-growth formula is asserted.

### Hostile: same size and group order do not determine growth

Let `L_k` denote the transitive tournament of order k. The order-six tournaments

    U=C3[L_2,L_2,L_2],       V=C3→L_3

both have automorphism group of order three. Yet

    Z_U^+ = (x_1^6+2x_3²)/3,       A(U)=24,
    Z_V^+ = (x_1^6+2x_1³x_3)/3,    A(V)=32.

In U the nonidentity automorphisms move all six vertices in two cycles; in V they fix three vertices. Thus size plus automorphism order cannot replace the cycle-index sidecar. Also `R(U)=8` and `R(V)=0` (V is not self-converse). These claims are checked directly, not inferred from orders alone.

### Hostile: even both signed indices do not classify tournaments

The transitive four-tournament `L_4` and the strong tournament `C3[L_2,1,1]` are nonisomorphic (one is strong and the other is not), but both have

    Z^+=x_1^4,       Z^−=x_2²,       A=16,       R=4.

Thus the signed pair is an exact sufficient sidecar for these orbit counts, not a complete tournament signature. Retain the full child-coloured quotient tree for classification. This hostile is directly verified by the checker.

## G4. A special unpointed descent: odd self-converse tournaments

In any odd-order self-converse tournament W, the vertices fixed by some negative automorphism form exactly one `Aut(W)`-orbit. Hence, if `sc(n)` denotes the number of ordinary self-converse tournament isotypes of order n, then for m≥1

    sc(2m+1) = Σ_[T], |T|=2m R(T).                         (G4)

The sum may be restricted to self-converse even-order bases. This is a consequence of the signed grammar, not a historical claim of a new enumeration identity.

**Proof.** Tournament automorphism groups have odd order: an involution would transpose a pair and reverse its arc, and Cauchy's theorem excludes even group order. The signed extension has order `2a(W)`; its Sylow 2-subgroups have order two, are generated by negative involutions, and are conjugate. Any negative f has order `2d` for odd d because f² has odd order; `f^d` is a negative involution. Its fixed vertex equals that of f, since G2 gives a unique fixed vertex on odd order. Given two negative involutions, a conjugator in the signed extension can be made positive by multiplying it on the right by the source involution if necessary. Their unique fixed vertices therefore lie in the same positive orbit.

Every odd self-converse W has a negative involution and thus a fixed vertex. Deleting it gives an even self-converse base and a marked self-converse extension counted by G3b. Conversely every extension counted by G3b is odd and self-converse. The unique orbit of admissible roots shows that forgetting the root identifies no two inequivalent marked self-converse classes. Bases of nonisomorphic deletions cannot create multiplicity, because equivalent marked classes have isomorphic deletions. ∎

G4 does not mean every vertex can serve as the fixed root, or that every deletion of a self-converse tournament is self-converse.

## Exact controls and reproducibility

Run `python3 check.py` in this directory. No dependencies beyond Python's standard library. The checker:

1. exhausts every labelled tournament of orders 1–5 (1,099 inputs);
2. constructs all modules from the definition, then the canonical maximal proper strong modules;
3. compares recursive signed cycle polynomials with independent direct permutation enumeration;
4. compares G3a/G3b against explicit binary attachment orbits and direct complement-converse action;
5. checks nine named examples/composites through order eight, including the growth hostile and anti-cosets with genuine 6-cycles.

6. verifies the unique anti-fixed-vertex orbit on every odd-order self-converse input, and G4 against independent isotype enumeration: `sc(3)=2`, `sc(5)=8`.

Both ordinary and `python3 -O check.py` runs pass and produce byte-identical JSON; all validation gates use explicit exceptions, not optimization-sensitive assertions. The polynomial coefficients use exact rational arithmetic. The direct comparison on prime inputs largely checks the base signed enumeration; composite controls specifically test nontrivial substitution. This is a finite exact check, not a universal computational proof. G1–G4 have proofs above. Empty tournaments are outside the implementation's domain; the order-one case in G4 can be added separately as `sc(1)=1`.
