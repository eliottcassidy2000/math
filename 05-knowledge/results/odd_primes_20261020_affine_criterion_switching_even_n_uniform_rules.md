# odd_primes_20261020 — the affine fixed-point criterion, the even-n switching half of HYP-9172, and uniform circulant rules (mac-mini-2026-10-20-oddprimes)

Owner's prompt: prove HYP-9340 beyond Sylow (first open case Z_3 at n = 6), find a uniform rule or definability obstruction for the prime-circulant bijection, then prove or disprove HYP-9172; think in odd primes and their powers.

Sibling notes (same machine): [even_graphs_20261019](even_graphs_20261019_sylow_refinement_components_walsh_duality.md) (which lists the earlier 2026-10-19 notes).

Canon: [THM-4724](../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md); [HYP-9341](../hypotheses/HYP-9341-what-remains-open-after-thm-4722.md) updated. Scripts: `04-computation/experiments/even_graphs_20261020/`.

## 0. What had happened since the previous close-out (read before planning)

Between the 2026-10-19 close-out and this session, other sessions pushed THM-4722 (HYP-9172 (a): tournament half for every n >= 5 and switching half for odd n, by the Sylow refinement), THM-4723 (Sylow-profile form of HYP-9340 and a defect reduction for the 3-cycle type), HYP-9341, and the 2026-10-20 packet: `uniform_equivariance_obstruction_20261020` (HYP-9172 (a), both halves, every n >= 5, by monotone odd-prime valuations along an equivariant type bijection plus THM-4721's product identity), `prime_local_walsh_bijection_20261020` (HYP-9340 (a) PROVED for every odd p-group H by a signed p-group orbit-localisation lemma), `c15_symmetry_counterexample_20261020` (HYP-9340 (a) REFUTED for mixed primes: C15 on ten vertices, 32 vs 31, the missing graph K_{5,5}), and a complete mixed-prime census through ten vertices. So the owner's three questions were already answered in substance: (1) the 3-cycle case is a 3-group case, hence proved; the general statement is false; (2) only finite negative searches existed for uniform rules; (3) HYP-9172 (a) is proved, with one route (the Sylow one) incomplete at even n for switching classes.

This session: (A) audited those proofs (no gap found) and recounted C15 and the Sylow profiles independently; (B) closed the even-n switching half by the Sylow route with a cleaner general principle; (C) proved the obstruction for the derivative family of circulant rules and found the power-sum detector.

## 1. The affine fixed-point criterion (Theorem A)

Tournaments modulo an invariant subspace `U` form an affine space `A = T/U` over `V/U`. For odd-order `P` and `h ∈ N(P)`, `Fix_A(P)` is a nonempty affine space over `W = (V/U)^P` (odd order, characteristic 2), `h` acts on it affinely, and the `<P,h>`-fixed classes are the `h`-fixed points there: `|W^h|` of them if one exists, else none. Three facts finish it: `W ≅ ((U-perp)^P)^*` as `<h>`-modules (nondegenerate pairing on `V^P`, `(U-perp)^P = (U^P)^perp`), `dim M^h = dim (M^*)^h` for a cyclic group on an `F_2`-space, and the fixed-point criterion `c ∈ (h-1)W` ⟺ all `h`-invariant functionals kill `c`, where the functionals are the `X ∈ (U-perp)^{<P,h>}` and `<X, hT_0 - T_0> = r(h, X)` (Royle et al.'s Lemma 4.3 telescoping). Hence

    |Fix_{T/U}(<P,h>)| = Σ_{X ∈ (U-perp)^{<P,h>}} (-1)^{r(h,X)}

for every `U`, every `n`. Checked for `U = 0` and `U = cuts`, all odd subgroup classes `P` and all `h ∈ N(P)`, `n = 5..8`: 1832 pairs, no mismatch; at `n = 6` the switching-class fixed-point counts were also brute-forced (`affine_criterion_check*.out`). At even `n` the identity is often nonzero for `h` of even order (36 of 132 pairs at n = 6), which is exactly why THM-4722's Euler parity lemma fails there while the identity survives.

Consequences: the restricted twisted Burnside principle R(U) and the Sylow refinement S(U) for every `U` at every `n`; THM-4721's parity lemma as a corollary (tournaments have no even-order automorphisms).

## 2. HYP-9172 (a), switching half, even n (Theorem H'')

Lemma 0 (elementary, as in the packet): a class stabiliser has no involution fixing a vertex. Lemma L: for even `n >= 8`, `n ≠ 10`, two primes `p > q > n/2` (Ramanujan: `π(x) - π(x/2) >= 2` for `x >= 11`) — no class stabiliser contains a `p`-cycle and a `q`-cycle (block argument ⇒ primitive on the joint orbit; Jordan ⇒ `A_m` unless `m <= q + 2`, and `A_m` has involutions with fixed points, which lift; `m = p = q + 2` ⇒ Burnside ⇒ affine or 2-transitive of odd degree, both impossible); for `n = 6, 10` the pairs (5-cycle, Sylow 3) and (7-cycle, Sylow 3) with a short case analysis. Then Hall fails for the stabiliser-containment graph exactly as in THM-4722: `|S| = tau_C(P_1) + tau_C(P_2)`, `|N(S)| <= tau_EE(P_1) + tau_EE(P_2) - 1` (the empty graph), equality of the `tau`'s by S(cuts). With THM-4722 (odd n, tournaments) this is the second complete proof of HYP-9172 (a).

## 3. Uniform rules (Theorem W, Proposition D)

- Forcing: for `p ≡ 3 (mod 4)` and `p ≡ 5 (mod 8)` the only even circulant invariant under `K_max(p)` (odd part of `Z_p^*`) is the empty graph; at `p = 41` (`e = 3`) this fails.
- Theorem W: `D_c` sends the `K_max`-invariant tournaments to `∅` if `c ∈ K_max(p)` and to the odd complete graph if `c ≡ -1 (mod K_max(p))`; for every `c ≠ 1` the latter happens at infinitely many primes (non-squares: `p ≡ 3 (mod 4)` with `(c/p) = -1`; `c = d^2`: `p ≡ 5 (mod 8)` with `(d/p) = -1`, which is automatic for `d = ±2 m^2` — the identity `-4 = (1+i)^4` — so `D_4` fails at every `p ≡ 5 (mod 8)`; the cases `c ∈ {m^8, 16 m^8}` need Chebotarev in the 2-power Kummer tower). Table for `2 <= c <= 12`, `p < 400` and the coset facts for `p <= 43`: `dc_family.out`.
- A first draft of this section claimed a Weil-bound obstruction for all local dilation rules; that was wrong (dilation atoms are proportional, so Weil's independence does not apply) and was replaced before commit. What survives for general dilation rules is the forcing analysis: `φ` must vanish on every pattern realised by Legendre-symbol vectors `(χ_p(c_i))_i`; patterns outside that set are unconstrained, and general non-existence is OPEN.
- Proposition D (power sums): `σ_k(S) = Σ_{s∈S} s^k` transforms by `a^k`; `S` is traceless (all `σ_k = 0` for `gcd(k, p-1) = 1`) iff its multiplier stabiliser is nontrivial, when the odd part of `p-1` is a prime power or `p ≡ 3 (mod 4)` with two odd primes (`F_p`-Fourier analysis on `Z_p^*` and a `±1`-array argument); exact for all `5 <= p <= 31` (`circulant_power_sums.out`); false at `p = 61` (1444 traceless sets, 360 rigid, meet-in-the-middle, `traceless_p61.out`) and `p = 211` (slice-mixed construction, `traceless_p211.out`). In the good regimes the least unit index with `σ_k ≠ 0` normalises every rigid orbit canonically and uniformly in `p`.

## 4. Audits and independent recounts

THM-4722 read in full: sound. Packet: Lemma 1 (monotone valuations summing to zero) is the right use of the product identity; Lemma 4's 3-cycle support argument is correct; the Sylow-action addendum's signed localisation is correct; C15 recounted with an independent pipeline (labelg types, dreadnaut generators): 2 = 2, 8 = 8, 32 vs 31 (`c15_independent.out`). Sylow-p profiles (full conjugacy class of a Sylow p-subgroup of Aut) agree between tournaments and even graphs for p = 3, 5, 7 and n <= 8 (`sylow_profile.out`; e.g. n = 8, p = 3: 264 + 32 + 16 on each side), which is the packet's theorem seen from the census.

## 5. Odd primes and their powers: the dictionary this week produced

| object | role of the odd prime p |
|---|---|
| Sylow p-subgroup P of S_n | multiplicity one (Sylow) + Frattini ⇒ types with P-symmetry = N(P)-orbits ⇒ tau_T(P) = tau_E(P) (THM-4721 S; every U: THM-4724) |
| p-unit entries of the Walsh matrix | signed p-group orbit localisation ⇒ Sylow-p profiles agree ⇒ tau_T(H) = tau_E(H) for every p-group H (packet) |
| v_p(|Stab|) along an equivariant type bijection | monotone and summing to zero ⇒ equal ⇒ the empty graph's preimage would carry full Sylow 3 and 5 ⇒ HYP-9172 (a) (packet) |
| p-cycles with p > n/2 | Jordan + Burnside ⇒ no class stabiliser carries two of them ⇒ Hall violator at even n (THM-4724) |
| mixed primes 3·5 | C15 at n = 10: the N(H)-even odd graph K_{5,5} is uncompensated ⇒ HYP-9340 false in general (packet) |
| prime circulants, p ≡ 3 (mod 4) | QR_p forced to the empty graph; D_c fails iff c is a non-residue; power-sum detector exact when p-1 has at most two odd primes (THM-4724) |
| p ≡ 1 (mod 4), 2^e || p-1 | forcing to ∅ holds for e <= 2 and fails at p = 41; the detector fails at p = 61 (THM-4724) |

## 6. Next obligations

1. A general non-existence theorem for uniform dilation rules, or a uniform rule using the power-sum normalisation (the residue of HYP-9341 (c)).
2. The exact boundary of Proposition D (p = 151, 331: p ≡ 3 (mod 4) with a square odd prime factor).
3. HYP-9172 (b) (branch integrality for n ≡ 0 (mod 4)) is untouched by all of this.
4. Audit THM-4724 (pre-audit): Theorem A step (3) (the duality `W ≅ ((U-perp)^P)^*`), Lemma L's lifting of involutions, Theorem W's Kummer cases.
