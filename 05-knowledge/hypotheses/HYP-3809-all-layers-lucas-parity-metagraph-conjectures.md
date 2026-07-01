---
id: HYP-3809
title: THE ALL-LAYERS LUCAS PARITY LAW of the merged metagraph (unifies HYP-1772 d=1 + klein-S75 d=m into ONE Lucas-graded law) + fiber freeness + the rigidity-mod-2 answer + a conjecture menu. For the fixed-base-path tiling cube Q_m (m=C(n-1,2)), merged node u (mass M_u = tiling count), and Hamming layer d (flip d tiles): 2*lambda_u^{(d)} + tau_u^{(d)} = C(m,d)*M_u (each of M_u tilings has C(m,d) distance-d neighbors), so tau_u^{(d)} ≡ C(m,d)*M_u (mod 2). Since M_u ≡ [u is SC] (HYP-1772, PROVED: F(C)=H(C)/|Aut(C)| odd because Redei H odd + |Aut| odd [no involutory tournament automorphism]), and C(m,d) mod 2 = Lucas(m,d) = [d is a binary submask of m], the CROSS-INCIDENCE tau_u^{(d)} is ODD <=> (u is SC) AND (d submask of m). VERIFIED exhaustive n=4,5,6 (0 violations). SPECIAL CASES: d=1 (HYP-1772 wiggly edge-balance, C(m,1)=m) and d=m (klein-S75 blue/black lines, C(m,m)=1 => ALL SC nodes have odd d=m-cross-degree => #SC EVEN by handshake). PARITY-ACTIVE LAYERS = the binary submasks of m (a Boolean subcube of size 2^{s(m)}, s=digit-sum): {1,2,3}(n=4), {2,4,6}(n=5), {2,8,10}(n=6). FIBER: which tilings share a node = one Aut(C)-orbit of Hamiltonian paths (free action, LEM-003), F(C)=H/|Aut| tilings. RIGIDITY (owner's hypothesis 'metagraph = the constraint'): TRUE mod 2 -- the odd/odd/even category rule + the all-layers Lucas law FORCE the metagraph's entire Z2 parity skeleton; but NOT over Z -- the integer masses are F(C)=H/|Aut| (tournament data), which parity alone underdetermines. Plus a menu of ~13 conjectures
status: MIXED (main law verified n<=6 + near-proof; conjecture menu). VERIFIED exhaustive n=4,5,6: the all-layers law tau_u^{(d)}≡C(m,d)M_u mod2 (0 violations); parity-active layers = binary submasks of m; d=1/d=m specialize to HYP-1772/S75; F(C)=H/|Aut| holds and is odd. NEAR-PROVED: the law follows from the counting identity 2*lambda + tau = C(m,d)*M_u (exact, each tiling has C(m,d) layer-d neighbors) + M_u parity = [SC] (HYP-1772, proved all n). So the all-layers law is PROVED for all n GIVEN the identity (which is definitional) -- the only computed part is confirming no off-by-parity. HONEST: a unifying synthesis (mostly proof-backed) + fiber structure (LEM-003) + a precise rigidity answer + conjectures; not a closed enumeration of the metagraph. n>=7 not recomputed here.
source: klein-2026-07-01-S76
depends_on:
  - HYP-3808   # S75: the d=m blue/black line accounting (tripartite) -- the d=m case
  - HYP-1772   # merged-tiling bucket parity (F odd, SC odd/NS even), PROVED -- the M_u parity input
related:
  - HYP-3808   # merged-metagraph blue/black (d=m); THREE-WAY convergence klein-S75 + mac-mini-S83 (black=even graph/blue=odd, self-loop-conj refuted at n=6) + opus-S18 (GRID REFLECTION = COMPLEMENT, which explains conjecture 8: grid-sym tiling fixed by sigma=complement => SC => blue only touches SC); this all-layers law is the Lucas-graded extension of that cluster
  - THM-002    # OCF/Redei (H odd) -- the fiber oddness
  - HYP-2689   # half-address / half tiling model (grid-sym subcube = blue)
  - THM-551    # half tiling / grid-sym structure
  - HYP-3803   # flip-rank (same tiling cube; n=6 transitions)
external: Lucas' theorem (C(m,d) mod 2 = [d submask m]); LEM-003 (free Aut action on Ham paths); Redei
results:
  - 04-computation/all_layers_parity_fiber_klein.py
  - 05-knowledge/results/all_layers_parity_fiber_klein.out
---

# HYP-3809 — the all-layers Lucas parity law + fiber freeness + rigidity + conjecture menu

## The unifying law
Fixed base path, tiling cube `Q_m` (`m = C(n-1,2)`). For a merged node `u` (mass `M_u` = # tilings) and a
Hamming layer `d` (flip exactly `d` tiles), each of `u`'s `M_u` tilings has `C(m,d)` distance-`d`
neighbors, so
> `2*lambda_u^{(d)} + tau_u^{(d)} = C(m,d) * M_u`  (`lambda` = self-loops inside `u`, `tau` = cross), hence
> **`tau_u^{(d)} ≡ C(m,d) * M_u (mod 2)`**.
With `M_u ≡ [u is SC] (mod 2)` (HYP-1772) and **Lucas** `C(m,d) mod 2 = [d is a binary submask of m]`:
> **`tau_u^{(d)}` is ODD  ⟺  (`u` is SC) AND (`d` ⊆ `m` in binary).**
VERIFIED exhaustive `n=4,5,6` (0 violations). This is ONE law with two known special cases:
- **`d=1`** = HYP-1772's wiggly edge-balance (`C(m,1)=m`: odd cross ⟺ SC ∧ `m` odd).
- **`d=m`** = klein-S75's blue/black lines (`C(m,m)=1`: ALL SC nodes have odd `d=m`-cross-degree ⟹ **`#SC`
  even** by the handshake).
**Parity-active layers** = the binary submasks of `m` = a Boolean subcube of size `2^{s(m)}` (`s` = binary
digit-sum): `{1,2,3}` (n=4, m=3), `{2,4,6}` (n=5, m=6), `{2,8,10}` (n=6, m=10). At non-submask `d`, every
node has even cross-incidence.

## The fiber: which tilings share a node
`F(C)` = # tilings mapping to class `C` = `H(C)/|Aut(C)|` (HYP-1772), and **`Aut(C)` acts FREELY on the
Hamiltonian paths** of `C` (LEM-003), so the tilings in one node are exactly **one `Aut(C)`-orbit of
Hamiltonian paths** (each Ham path = one way to place the base path). Two tilings share a node iff a
relabeling carries one base-path-realization to the other (up to `Aut`/complement). Because `|Aut(C)|` is
odd (an involution would reverse an edge — impossible in a tournament), `F(C)` is odd, giving the SC-odd /
NS-even bucket parity. The **half tiling model** (HYP-2689/THM-551) parametrizes the grid-symmetric (blue)
tilings by half-addresses (one bit per `sigma`-orbit); the blue subcube has dimension `(m+fix)/2`,
`fix = floor((n-1)/2)`.

## The rigidity answer (owner's hypothesis: "the metagraph = the constraint")
> **TRUE mod 2, FALSE over Z.** The odd/odd/even category rule (SC pure-blue odd, SC mixed odd, NS
> pure-black even) together with the all-layers Lucas law FORCE the metagraph's *entire Z/2 parity
> skeleton* — every cross-incidence parity, `#SC` even, the 2-adic valuations (`v_2(M_u) = [u not SC]`,
> never ≥2). But the *integer* masses are `F(C) = H(C)/|Aut(C)|`, genuine tournament data that parity alone
> underdetermines (many parity-valid mass-vectors are not realized). So the constraint IS the metagraph's
> parity content, and the residual (the sizes) is exactly the Redei/`|Aut|` arithmetic.

## Conjecture menu (multitudes, small and large)
1. **All-layers Lucas law** (above) for all `n` [proof-backed via the counting identity + HYP-1772].
2. **Parity-active-layer count** `= 2^{s(m)} - 1`, `s(m)` = binary digit-sum of `m = C(n-1,2)`.
3. **`#SC` even for all `n`** (`d=m` handshake; HYP-1772: `#SC = 2,2,8,12,88` for `n=3..7`, all even).
4. **Fiber freeness**: `Aut(C)` acts freely on Ham paths, `F(C)=H/|Aut|`, odd (LEM-003 + Redei).
5. **No high 2-adic mass**: `v_2(M_u) ∈ {0,1}` always (`0` iff SC); never `≥2`.
6. **Pure-blue never self-loops** at `d=m` (verified `n<=6`); conjecture all `n`.
7. **NS-node self-loop onset at `n=6`** (`d=m` black): the first NS class admitting a flip-all-tiles
   symmetry (`T(flip(t)) ≅ T` or `T^op`); tie the onset to `s(m)` / the first non-rigid NS class.
8. **Tripartite color rule**: blue (`d=m`, grid-sym) lines touch only SC nodes; black touch only `NS∪MX`.
   Conjecture: grid-symmetry of a tiling ⟹ its tournament is SC-related (why blue avoids NS).
9. **`#pure-blue` = # SC classes all of whose `F(C)` tilings are grid-symmetric**; conjecture a `2`-power
   or `sigma`-fixed-point formula.
10. **Self-loop count** `lambda_u^{(d)} = (C(m,d)M_u - tau_u^{(d)})/2` = a symmetry census (flip-`d`
    symmetries of `u`'s tournaments); conjecture `lambda_u^{(m)}` counts flip-all-tiles automorphisms.
11. **Blue-line count** `= 2^{(m+floor((n-1)/2))/2 - 1}` exact (grid-sym / half-tiling dimension).
12. **n=6 as a universal threshold**: flip-rank shape death (HYP-3803), NS self-loop onset (HYP-3808), and
    the `m=10=1010_2` Lucas pattern all break at `n=6` — conjecture a common cause in `s(m)` / the first
    `m` with `≥2` binary `1`s beyond the top (`m: 1,3,6,10 = 1,11,110,1010`; two 1s from `n=4` on).
13. **Recursion** `(#PB,#MX,#KB) = (1,1,1),(3,5,2),(2,10,22)`; `#KB = (A000568−SC)/2 = 0,1,2,22,184`.

## Net
The merged metagraph's parity structure across ALL Hamming layers is one Lucas-graded law
`tau_u^{(d)} ≡ C(m,d)[u SC] (mod 2)`, with HYP-1772 (`d=1`) and klein-S75 (`d=m`) as endpoints and the
parity-active layers = binary submasks of `m`. The fiber is a free `Aut`-orbit of Hamiltonian paths
(`F=H/|Aut|`, odd). The owner's "metagraph = the constraint" is exactly right **mod 2** (the parity
skeleton is forced) and needs the `Redei/|Aut|` arithmetic **over Z**. Plus a menu of testable conjectures,
several proof-backed.
