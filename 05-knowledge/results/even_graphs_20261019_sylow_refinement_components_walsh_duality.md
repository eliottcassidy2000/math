# even_graphs_20261019 — tournaments versus even graphs: the Sylow refinement, components, and the Walsh duality (mac-mini-2026-10-19-evengraphs)

Owner's prompt: explore the repo's past sequence work in the vein of k-tournaments and hypertournaments; find an exact bijection between tournaments and even graphs; merge the two attached blueprints (Erdős 173; countable colourings of C-linear graphs) loosely; maybe Minkowski's ?.

Sibling notes of 2026-10-19 (same machine): [two_three_axis](two_three_axis_20261019_numbers_fixed_points_relations.md) (which lists the Mersenne, Farey, duality and Apollonius notes).

Canon: [THM-4721](../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md), [HYP-9340](../hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md). Scripts and outputs: `04-computation/experiments/even_graphs_20261019/`.

## 0. Status table

| item | statement | status |
|---|---|---|
| P | parity lemma: `P` odd normal in `G = <P,h>`, `h` of even order: the F_2-sum of the `P`-translates of the antipodal graph of an even cycle of `h` is `G`-invariant and `h`-odd | PROVED; checked for 440 pairs (n <= 7) |
| R | restricted twisted Burnside: `#Gamma`-orbits on `Fix_T(P)` = `#Gamma`-orbits on `Gamma`-even members of `Fix_G(P)`, every odd `P`, every `Gamma <= N(P)` | PROVED; checked for `Gamma = C(P), N(P)`, all odd subgroup classes, n <= 7 |
| S | Sylow refinement: `tau_T(P) = tau_E(P)` for Sylow `p`-subgroups, `p` odd | PROVED |
| C | prime circulants: `t_p = e_p = (1/(p-1)) sum_{d odd | p-1} phi(d) 2^{(p-1)/2d}`; refinement per odd multiplier subgroup; an `AGL(1,p)`-equivariant bijection exists | PROVED; p <= 200 by formula, p <= 31 brute force, criterion vs networkx p <= 13 |
| G' | `Gamma`-refined identity for every `Gamma <= S_n`; canonical equivariant bijection iff `|Gamma|` odd | PROVED |
| K | coloured: even `2k`-colourings = `k`-arc-coloured tournaments (THM-4720), HYP-9339 (3) | PROVED; exact for 2, 4, 6 colours, n <= 4 |
| E, A | even iff all components even; A334335 = connected even graphs; `EULER(A334335) = INVERT(A051337)` | PROVED; census n <= 8; b-file to n = 128 |
| W | Walsh basis of invariants over `Z[1/2]`, two orthogonality relations, `det M = +-2^k`, odd parts of `prod |Aut|` agree | PROVED; n <= 5 matrices, n <= 8 products |
| Q | class identity `sum_T m_T = sum_{X even} m_X + sum_{X odd, N(P)-even} m_X` | PROVED; n <= 7 element form |
| HYP-9340 (a) | `tau_T(H) = tau_E(H)` for every odd-order `H` | OPEN; exhaustive n <= 8 (all odd subgroup classes) |
| sequences | odd graphs, labelled even graphs, connected odd graphs not in OEIS | FINITE-EXACT / lookups 2026-10-19 |

## 1. What was already in the repo

THM-4531 (2026-10-01, tbij lane) already proves Theorem G (the `U`-subspace form of Royle et al.) and, via Lemma N1 (an `S_n`-equivariant map that is bijective on types exists iff the stabiliser-containment bipartite graph has a perfect matching), that **no `S_n`-equivariant bijection tournaments -> even graphs exists for `5 <= n <= 15` and sixteen larger `n`** (deficits 1, 2, 5, 12, 44 for n = 5..9; the census here reproduces 1, 2, 5, 12). HYP-9172 conjectures this for all `n >= 5` and reduces it to additive representations over twist-rigid block sizes. THM-4524 E1 is the switching-class form. So "find an exact bijection" cannot mean an `S_n`-equivariant one; this session asked what it can mean, and what finer identities hold.

The relevant past sequence program is the S48 Burnside enumerator (`burnside_enum_v2.c`: A000568 to 200 terms, hypergraphs, inverse Euler transforms for connected variants) and THM-4720 (k-valued tournaments `T_k`, 3-hypertournaments). The even-graph family joins it: even graphs = A000568; connected even graphs = InvEuler(A000568) = A334335; odd graphs = A000088 − A000568; even `2k`-coloured complete graphs = `T_{2k}`.

## 2. The principle and the lemma

Everything rests on one computation (THM-4721 (0.1)): for a group `G` of permutations and `h in G`, the signed count `sum_{X in Fix_G(G)} sgn_X(h)` is `2^{c_2(G)}` if every orbital graph `X_O` is `h`-even and `0` otherwise, while `|Fix_T(G)|` is `2^{c_2(G)}` for `|G|` odd and `0` otherwise. The two agree whenever "`|G|` even ⇒ some `G`-invariant graph is `h`-odd". Lemma P proves this for `G = P<h>` with `P` odd normal: with `Y_C` the antipodal graph of an even cycle `C` of `h`, the symmetric difference `X = sum_{k in P} kY_C` is `G`-invariant and `r(h, X) = sum_k r(h, kY_C) = |P| r(h, Y_C) = 1 (mod 2)` by the cocycle identity and the bijection `k -> hkh^{-1}` of `P`. Theorem R (restricted principle) follows by Cauchy–Frobenius on both sides.

## 3. The Sylow refinement and prime circulants

For a Sylow `p`-subgroup `P` (`p` odd) two extra facts make `N(P)`-orbits on `Fix(P)` equal to isomorphism types with a `P`-symmetry, on both sides: Sylow (two conjugates of `P` inside `Aut(x)` are `Aut(x)`-conjugate, so each type meets `Fix(P)` in one `N(P)`-orbit) and Frattini (an odd graph in `Fix(P)` has an odd automorphism normalising `P`). Hence Theorem S: `tau_T(P) = tau_E(P) = (1/|N(P)|) sum_{h in N(P), |<P,h>| odd} 2^{c_2(<P,h>)}`.

`n = p`: `P = Z_p`, `N(P) = AGL(1,p)`, types = multiplier classes (Turner's theorem, re-derived). Elements `x -> ax + b` with `ord(a) = d` odd give `<P,h> = Z_p x| <a>` with `(p-1)/(2d)` pair orbits, so `t_p = e_p = (1/(p-1)) sum_{d | p-1, d odd} phi(d) 2^{(p-1)/(2d)}`: 1, 1, 2, 4, 6, 16, 30, 94, 586, 1096, 7286 for p = 3..37 (`circulant_prime.out`, all 45 primes to 200 EQUAL; brute force to 31). The brute-force side needs an evenness criterion for circulants: by Alspach, `Aut(circ(p,S)) = Z_p x| Stab(S)` for `S` not empty or full, and `x -> ax` with `ord(a)` even reverses `|S|/ord(a)` edge orbits, so `circ(p,S)` is even iff `2 ord(a)` divides `|S|` for every even-order `a` in `Stab(S)`; checked against networkx for p <= 13. Refined per odd multiplier subgroup `K` (apply Theorem R to `Z_p x| K`): e.g. p = 31: `|K| = 3`: 4/4, `|K| = 5`: 2/2, `|K| = 15`: 1/1. Inclusion–exclusion on the divisor lattice gives Hall's condition for THM-4531's stabiliser-containment graph, so an `AGL(1,p)`-equivariant bijection exists for every prime `p`; for p = 7 it is `QR_7 -> empty`, rigid circulant -> complement of `C_7`. No uniform rule found: the multiplier derivative `D_c(S) = {x : [x in S] != [cx in S]}` (automatically symmetric, with an even number of classes) induces a bijection for `p = 7, 11, 19, 23` with `c = 4` but fails at 5, 13, 17; every fixed `c in {+-2, +-3, +-4}` fails somewhere (`circulant_prime.out`). The obstruction is arithmetic: the `K_max`-invariant tournaments must go to the empty graph, which `D_c` does iff `c in K_max`.

## 4. Components, the Euler transform, and the sequences

Lemma E: a graph is even iff every component is even (an odd automorphism `g` with its antipodal-edge cycles is traced to a power `g^k` on a single component). Hence even graphs are multisets of connected even graphs and A334335 (von Broemssen 2020, 19 terms, "it appears that a(n) is the number of connected even graphs (verified for n <= 10)") is exactly that; Royle et al. proved the total count, the component lemma was the missing step. The repo's A000568 b-file (129 terms) gives A334335 to n = 128 (`connected_even_graphs_bfile.txt`). Also `EULER(A334335) = A000568 = INVERT(A051337)`: multisets of connected even graphs and sequences of strong tournaments are equinumerous, two different "connected" decompositions of the same numbers (no bijection between the pieces: 2 connected even graphs vs 1 strong tournament at n = 4).

Census (n <= 8, `census_n7.out`, `census_n8.out`; nauty geng/gentourng + networkx automorphisms): even = tournaments (2, 4, 12, 56, 456, 6880); odd graphs 1, 2, 7, 22, 100, 588, 5466 (Burnside formula to n = 60 in `odd_graphs_bfile.txt`); labelled even graphs 1, 4, 20, 396, 16432, 1312056, 200582824 (of `2^{C(n,2)}`); connected even 1, 2, 8, 43, 398, 6413; connected odd 1, 4, 13, 69, 455, 4704; rigid even graphs 0, 0, 0, 8, 152, 3696 vs rigid tournaments 1, 2, 7, 41, 399, 6556; `sum_{X even} (-1)^{|E|}` = 2, 2, 8, 12, 88, 176 = self-converse tournaments (THM-4531 G2 confirmed). OEIS (2026-10-19): odd graphs, labelled even graphs, connected odd graphs: no match.

## 5. The Walsh duality: the exact linear bijection

For even `X`, `Phi_X(T)` = (copies of `X` in `K_n` on which `T` induces the `+` orientation class) − (copies in the `−` class); well defined exactly because `X` is even (`Aut X` preserves the classes). `Phi_empty = 1`, `Phi_{P_3} = #directed 2-paths − #2-stars`. THM-4721 Theorem W: the `Phi_X` are a `Z[1/2]`-basis of the isomorphism-invariant functions on tournaments; `sum_T |orb T| Phi_X(T) Phi_{X'}(T) = 2^{C(n,2)} |orb X| delta` and `sum_X |Aut X| Phi_X(T) Phi_X(T') = 2^{C(n,2)} |Aut T| delta`; `det M = +-2^k` with `k = mN/2 + (1/2) log_2(prod_T |Aut T| / prod_X |Aut X|)` (n = 3, 4, 5: `-2^2, -2^8, 2^49`, `walsh_check.out`); so the products of the automorphism group orders over all tournaments and over all even graphs have the same odd part (`3^328 5^12 7^4` at n = 8). This is the precise sense in which the two sets are "dual bases"; it is canonical and it is not a set bijection.

## 6. HYP-9340 (a): the refined conjecture and its reduction

`tau_T(H) = tau_E(H)` holds for every odd-order subgroup class of `S_n`, n <= 8 (`subgroups_fix.out`, `subgroups_fix_n8.out`: enumerate the `H`-invariant labelled objects as unions of pair orbits, canonicalise with labelg, count types / even types; 10 s at n = 8). Proved for Sylow `H`. For general odd `H` the Sylow proof's two steps fail separately and cancel: tournaments may have several `Aut`-classes of `H`-subgroups (n = 6, the two-3-cycle-module tournament, `Aut = Z_3 x Z_3`), and odd graphs may be `N(H)`-even (`K_{3,3}`: its only odd automorphisms swap the sides and do not normalise a 3-cycle). Proposition Q (the multiplicity-weighted identity) is what Theorem R gives directly.

Reduction. Write `iota(M) = [M contains a conjugate of H]` and expand `[x has an H-symmetry] = sum_K c(K) [K <= Aut x]` over subgroups `K` generated by conjugates of `H`, `c(K) = sum_{K' <= K, iota(K')=1} mu(K', K)`. Cauchy–Frobenius on both sides turns the conjecture into

    sum over even-order L of 2^{c_2(L)} sum_{a in L^+} phi_a(L) = 0,   phi_a(L) = sum_{M in [<a>, L]} mu(M, L) iota(M),

with `L^+ = {a in L : chi_O(a) = 1 for every L-orbit O}` (a subgroup, the kernel of the orbital characters). The odd-order `L` cancel termwise; the even-order ones must cancel in total. The termwise per-permutation identity is false (THM-4721 section 6: n = 4, type (3,1), labelled 16 vs 10, cancelled by the three double transpositions).

## 7. Negative results

- No `S_n`-equivariant bijection for 5 <= n <= 15 (THM-4531; deficits reproduced).
- No equivariant map from `Gamma`-even graphs to tournaments for any `Gamma` of even order (the empty graph).
- `|Aut|` multisets differ even in odd parts (n = 7: 15, 45, 315 occur only on the graph side); only the products agree.
- The edge-count generating function of even graphs is not a signed Burnside sum over odd-order permutations (the factors `1 − y^{|O|}` of self-paired orbits do not vanish for `y != 1`), so no natural tournament statistic can match `|E(X)|`; at `y = −1` it is the self-converse count (THM-4531 G2).

## 8. The attached blueprints and Minkowski

Erdős 173 (Plasma AI blueprint): every non-equilateral triangle has a monochromatic copy in any 2-colouring of the plane; the equilateral exception is the strip colouring by the parity of `floor(y/h)`, a parity class preserved by the triangle's symmetry. Countable colourings of C-linear graphs (second blueprint): sign assignments and their classes, bipartization as the edge cost of a sign class. The shape shared with this session is "a parity class preserved by every symmetry": even graphs are the graphs whose automorphisms preserve the two orientation-parity classes, and the signed Burnside sums are the bookkeeping of that preservation. Nothing transfers as a theorem. Minkowski's ?: the Walsh characters used in section 5 are the Fourier basis of the dyadic group, and ? conjugates Farey to dyadic structure; no honest statement about tournaments came out of it.

## 9. Scripts

`census.py` (n <= 8; evenness, spectra, cycle-type counts, Hall deficit and violator, matrix M for n <= 5), `classes_identity.py` (element-class identity), `subgroups_fix.py` (all odd-order subgroup classes, n <= 8), `parity_lemma_check.py` (Lemma P and Theorem R), `circulant_prime.py` (t_p, e_p, refinement, rules), `coloured_even.py`, `walsh_check.py`, `odd_graphs_bfile.py`; b-files `connected_even_graphs_bfile.txt` (A334335 to 128), `odd_graphs_bfile.txt` (to 60). All runs are seconds except `circulant_prime.py` (minutes; the networkx check must skip the empty and complete graphs).

## 10. Next obligations

1. Prove HYP-9340 (a) for non-Sylow `H` (pair the two defects, or the lattice identity of section 6); the first open case is `Z_3` of type `(3,1,1,1)` at n = 6.
2. A uniform natural rule for the prime-circulant bijection, or a proof that none is first-order definable in the ring `Z/p`.
3. OEIS: submit odd graphs, labelled even graphs, connected odd graphs; extend A334335 (b-file to 128) with the component lemma as the proof of its comment; cross-reference THM-4720's k-tournament sequences.
4. Audit THM-4721 (pre-audit): Lemma P's cocycle bookkeeping, the Frattini step, Theorem W's `Z[1/2]` argument.
5. Does Theorem S extend to Hall-type subgroups (products of Sylow subgroups for several odd primes, when they exist in `S_n`)? The data for `Z_15 <= S_8` (2/2) suggests yes.
