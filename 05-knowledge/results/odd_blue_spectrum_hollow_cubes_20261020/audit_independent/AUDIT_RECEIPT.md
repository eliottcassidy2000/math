# Independent audit: odd-order blue-degree pairing and hollow-frame cubes

Date: 2026-10-20. Status: **APPROVED / PROVED**, for the precise theorem below. The finite computations are independently reproducible checks, not the all-order proof. No repository or remote files were modified by this audit.

## Verdict and exact scope

For every odd `n=2m+1 >= 3`, self-converse tournament isomorphism classes admit a fixed-point-free involution preserving the blue-tiling count `B` and reversing the parity of the number `c3` of cyclic triples. More strongly, each hollow-frame stratum decomposes into free cubes of dimension `r`, the number of pair-orbits of its hollow-frame automorphism group. Every cube has constant `B` and exactly half its classes on each `c3`-parity shore.

Every pair-orbit has odd size, so `r ≡ m (mod 2)`. Consequently, for `n ≡ 1 (mod 4)` and `n>=5`, each `B`-multiplicity on either shore is even, and its multiplicity among all SC classes is divisible by four.

No gap remains in the descent from marked to unmarked classes, the canonical orbit selection, the centralizer denominator, or the freeness of the strengthened cube action. `n=1` is genuinely excluded.

## Inheritance and independent route

The audit read `math/AGENTS.md`, the bounded startup packet and relevant current navigation/guardrails, the research protocol, THM-643-blue-black-line-parity-structure, THM-644-gridsym-fiber-law-and-anti-redei, THM-647-anti-redei-theorem, and the current staged corrections to THM-643 and HYP-3809-merged-metagraph-bucket-conjecture-atlas. The old all-order spectrum pairing is already refuted at `n=8`; nothing here revives it at even orders. The pure-blue-implies-leaf assertion is already refuted at `n=5`.

The independent implementation imports no producer or inherited census module. It uses an explicit `m`-vertex tournament `A`, an independent symmetric off-diagonal binary matrix `K`, and an `m`-bit matching completion `b`, with full vertex-permutation checks as a second route. This is different code and a different bit convention from the producer.

Live concepts audited: anti-involution conjugacy; hollow-frame automorphisms; fixed paths versus fixed quotient orbits; canonical orbit choice; free cube action; actual antipodal graph versus auxiliary cube graph.

## Decisive proof obligations

### 1. Existence, conjugacy, and forgetting the marking

For `n>=2`, the automorphisms and anti-automorphisms of an SC tournament form a signed group `G±` of order `2|G|`, where `G=Aut(T)`. A tournament has no involutory automorphism, because such an automorphism exchanges the endpoints of some arc. Hence `|G|` is odd. Sylow 2-subgroups of `G±` have order two and give anti-involutions.

All anti-involutions are conjugate in `G±`, and in fact under `G`: if an anti-automorphism `a` conjugates `rho` to `rho'`, the automorphism `a rho` has the same conjugation effect. Thus forgetting `rho` gives a bijection from isomorphism classes of marked pairs `(T,rho)` to unmarked SC classes. This is the essential missing step in an informal marking-level proof.

An anti-involution fixes at most one vertex. At odd order it therefore fixes exactly one vertex `f` and exchanges `m` disjoint pairs.

### 2. The hollow frame is intrinsic and has odd automorphism group

Polarize each pair by `x_i -> f -> y_i`, with `rho(x_i)=y_i`. Delete only the matching arcs `{x_i,y_i}`. The resulting partial tournament `D` itself remembers the marking: `f` is its unique vertex with no missing incident edge, and the missing-edge matching recovers `rho`. Thus `Aut(D)=Aut(D,rho,f)`.

Every automorphism preserves `X={x_i}`, and restriction to `X` is faithful because its action on `Y` is forced by the matching. Therefore `H=Aut(D)` embeds in the tournament group `Aut(T[X])`, has odd order, and has only odd-sized pair-orbits.

The representation by `(A,K)` is exact: `A` gives arcs within `X`, arcs within `Y` are reversed copies, and `K_ij=[x_i -> y_j]=[x_j -> y_i]` is symmetric. These coordinates are unrestricted.

### 3. Completion quotient and centralizer denominator

Completions form `{0,1}^m`, where `b_i=[x_i -> y_i]`. Marked completion classes are precisely `{0,1}^m/H`. They are also unmarked SC classes in this hollow-frame stratum: an arbitrary tournament isomorphism can be adjusted by an automorphism of its target so that it carries the chosen anti-involution to the chosen anti-involution. The adjusted map then preserves the hollow frame. Distinct hollow-frame types cannot overlap on one unmarked class.

The completion stabilizer is exactly

    H_b = C_Aut(T_b)(rho).

Let `P_rho(T)` be the ordered Hamiltonian paths fixed by `rho` followed by positional reversal. A blue Hamiltonian-path orbit under `G` has positional reversal equal to some anti-involution. Conjugacy puts one representative in `P_rho`. Two representatives in `P_rho` differ by a member of `G` precisely when that member centralizes `rho`. The centralizer acts freely on ordered Hamiltonian paths. Therefore

    blue path-orbits ≅ P_rho(T) / C_G(rho),
    B(T) = |P_rho(T)| / |C_G(rho)|.

This is not `B=|P_rho|` unless the centralizer is trivial. The aggregate `H_anti` from THM-644 is the sum over all anti-involutions, not the count for one fixed anti-involution.

### 4. Why the symmetric path set survives all matching changes

In a `rho`-symmetric Hamiltonian path, matched vertices occupy positions `j` and `n+1-j`. At odd `n`, these positions cannot be adjacent. Consequently the path uses no deleted matching edge. Its validity depends only on `D`, so the entire set `P(D)` is independent of completion `b`.

The action of `H` on `P(D)` is free. The twisted reversal involution on all directed Hamiltonian paths and Rédei parity give odd positive `h(D)=|P(D)|`. Hence

    b0(D)=h(D)/|H| is an odd positive integer,
    B(T_b)=b0(D) * |H.b|.

Thus the blue-count dependence on matching orientation is entirely the completion-orbit size, not a change in the fixed-`rho` path set.

### 5. Canonical single-orbit switch

To avoid an arbitrary label choice, encode rooted hollow frames `(A,K;i)` and minimize the encoding over relabelings that put `i` in the first position. Two root codes agree exactly when the roots are in the same `H`-orbit. The least root code therefore identifies one intrinsic odd orbit `O`.

Set `b'=b+1_O` over `F2`. Since `1_O` is `H`-invariant, `H_b'=H_b`. The hollow frame and the intrinsic chosen orbit are unchanged, giving an involution. The quotient and conjugacy arguments above make it independent of representative and chosen anti-involution. The denominator and path numerator are both preserved, so `B` is preserved.

A choice based on the completed tournament rather than the hollow frame would not suffice: it could change after the switch. An arbitrary canonical *vertex* label also needs its full orbit, not a singled-out vertex that the frame does not distinguish.

### 6. Triangle parity, most directly

`rho` pairs all nonfixed unordered triples and preserves cyclicity. Its only fixed triples are `{x_i,f,y_i}`. Such a triple is cyclic precisely when `b_i=0`. Thus

    c3(T_b) ≡ m - sum_i b_i (mod 2).

Toggling any odd pair-orbit reverses `c3` parity. This also proves the switched unmarked class cannot be the original class.

### 7. Free cube action

For all pair-orbits `O_1,...,O_r`, the translations by their indicators commute, preserve `H_b`, and descend to the completion quotient. Define

    p([b])_j = sum_{i in O_j} b_i (mod 2).

This is invariant under `H`, while toggling `O_j` changes precisely its `j`th coordinate because `|O_j|` is odd. Therefore every nonidentity cube translation changes an invariant of the quotient class. The action is free, not merely an action by commuting involutions.

Each cube has constant blue count and equally many vertices on both `c3` shores. The fixed-frame total blue mass is `b0(D) 2^m`, with `b0(D) 2^(m-1)` on each shore.

## Boundaries and hostile controls

- `n=1`: no matching pair and no available orbit. The one SC class has `B=1,c3=0`; there is no pairing. The signed-group argument also should not be used at this degenerate order.
- `n=3`: the unique cube pairs the transitive tournament and the directed triangle. Their full automorphism orders are `1` and `3`, but their selected centralizers and blue counts are both `1`.
- `n=5`, cyclic `T5`: independently verified pair has source `H=15, |Aut|=5, fibre=3, B=3, c3=5` and partner `H=11, |Aut|=1, fibre=11, B=3, c3=4`. Thus purity, total fibre size, full automorphism group order, and ordinary Hamiltonian count already fail to be preserved at five vertices.
- `n=7`, `C3 -> 1 -> C3`: the nontrivial orbit switch gives `h_rho=3`, centralizer `3`, and `B=1` on both sides, but fibre sizes `1` and `45`. Full counts are `(H,|Aut|)=(9,9)` and `(135,3)`. This detects the incorrect formula `B=h_rho`.
- Auxiliary cube edges are not claimed to be antipodal blue edges or an automorphism of the actual blue metagraph. Equal shore spectra and the cube decomposition do not determine endpoint adjacency.
- No even-order theorem follows: a symmetric path there can use a matched edge as its central arc, and there is no unique fixed root providing the above polarization.

## Independent finite results

`independent_pairing_audit.py` tested every polarized completion at `n=3,5,7,9`: respectively `2,16,512,65536` raw completions, giving `2,8,88,2752` unmarked SC classes. It checks every skeleton, every completion, canonical choice consistency, exact symmetric path-set equality, stabilizer equality, triangle parity reversal, class-level involutivity, and the absence of fixed classes.

For all 98 representatives through `n=7`, it additionally enumerates *all full vertex permutations*. It directly identifies all automorphisms and anti-involutions, checks conjugacy, enumerates all directed Hamiltonian paths and their full adjacency-word orbits, and verifies the blue-count/centralizer formula. The order-nine class census does not claim this full-`S9` check for every class.

The parity-shore spectra are:

- `n=3`: `{1:1}`.
- `n=5`: `{1:2, 3:2}`.
- `n=7`: `{1:5, 3:6, 5:8, 7:16, 9:9}`.
- `n=9`: `{1:12, 3:14, 5:26, 7:36, 9:50, 11:26, 13:42, 15:66, 17:56, 19:120, 21:60, 23:120, 25:120, 27:88, 29:136, 31:136, 33:138, 35:64, 37:48, 39:18}`.

`independent_cube_audit.py` separately enumerates the hollow-frame classes and completion cubes, checks freeness by the invariant parity vector, checks constant stabilizers and `B`, and verifies frame mass formulas. Cube dimension counts are:

- `n=3`: one dimension-one cube; one hollow frame.
- `n=5`: two dimension-two cubes; two hollow frames.
- `n=7`: four dimension-one and ten dimension-three cubes; twelve hollow frames.
- `n=9`: sixteen dimension-two and 168 dimension-four cubes; 176 hollow frames.

These numbers and spectra agree exactly with the producer, without importing the producer's classification or path-count code.

## Additional exact correction: the odd-only 3-power cap also fails

The historical proposed bound `B <= 3^floor((n-2)/2)` was already false at even order eight. The order-nine data also rule out its odd-order restriction. The independent direct full-`S9` witness has adjacency row masks

    [356, 273, 274, 263, 169, 206, 158, 15, 240].

Bit `j` in row `i` means `i -> j`, with vertices `0,...,8`. Its directly enumerated values are

    H=2659, |Aut|=1, fibre=2659, B=39, c3=28,
    number of anti-involutions=1, chosen centralizer=1.

Thus `39>27` refutes the odd-order cap without relying on the class census. `independent_n9_cap_witness.py` reproduces this by all `9!` permutations.

## Reproduction and remaining structural obligations

From the workspace root:

    python3 session14-odd-blue-pairing-audit/independent_pairing_audit.py --max-m 4
    python3 session14-odd-blue-pairing-audit/independent_cube_audit.py
    python3 session14-odd-blue-pairing-audit/independent_n9_cap_witness.py

The first run takes roughly 90 seconds here; the others are short. JSON and stdout receipts are retained alongside the scripts.

The actual unresolved structural questions now separate cleanly: determine which hollow frames and orbit sizes realize which blue degrees; find the correct all-order blue-degree extremal bound after the cap refutation; and connect the auxiliary hollow-frame cubes to actual antipodal endpoint incidence. None of these is discharged by the class pairing alone.

### Added integer-triangle-law audit

The producer's final draft adds the exact refinement

    c3(T_b)=c3(T_0)+sum_i (2m-1-2d_D^+(x_i)) b_i.

This is **APPROVED**. Starting from `b_i=0`, the full score pair changes from `(d_i,2m-d_i)` to `(d_i+1,2m-d_i-1)`; applying the standard score formula for `c3` gives exactly the displayed coefficient. Different matching bits affect disjoint vertex-score pairs, so the changes add. The independent cube verifier now checks this integer identity for every completion class of every hollow-frame type through order nine. No new dependency or scope exception is needed.
