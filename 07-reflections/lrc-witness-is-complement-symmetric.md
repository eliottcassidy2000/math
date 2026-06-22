# The Lonely-Runner optimum is complement-symmetric — the SC locus again

*kind-pasteur-2026-06-22-S35. A reflection on why the LRC(14) extremal witness lives on the same complement-fixed locus that the tournament maximizers do, and what that says about both.*

## The observation

The project's central involution is **complement / reversal**: `T ↦ T^op`, which on residues is
`v ↦ N − v` (the `a = −1` element of `(ℤ/N)*`). The whole merged metagraph is organized around its
fixed locus — the **self-complementary (SC) spine**, the principal line from the transitive class to
the regular maximizer, where `H` is maximized (Paley, the SC backbone).

The Lonely Runner gap problem has the *same* involution at its optimum. For LRC(14) the binding pairs
at the optimal witness `τ*` are `(1,13), (5,9), (3,11)` (HYP-2571c) — **each pair sums to N = 14**.
That is exactly the complement pairing `v ↔ N − v`. The optimal witness is **fixed by `v ↦ −v`**:
the configuration `{v·τ* mod 1}` is symmetric about `1/2`, because the speeds come in complement
pairs whose phases are negatives of each other. The dilation that realizes the optimum's symmetry is
`a = −1` in `(ℤ/D)*` — the same `T ↦ T^op` reversal.

So both extremal problems concentrate on the **complement-fixed locus**:
- Tournaments: `H` is maximized on the SC classes (`s_i + s_{n-1-i} = n-1`, fixed by `T ↦ T^op`).
- Lonely Runner: `M(S) = max_τ min_v ‖vτ‖` is extremized at a `v ↦ N−v`-symmetric witness.

## Why it is forced, not coincidental

In both cases the quantity being extremized is invariant under the involution, and the involution acts
on a compact configuration space. A standard averaging/convexity argument then pushes the extremum onto
the fixed locus:

- For `H`: complementation permutes tournaments preserving `H` (it is `T ↦ T^op`), and the SC classes
  are the fixed points; the maximizer is among them (proved cases through n≤8, and the SC-maximizer
  dichotomy).
- For `M`: `min_v ‖vτ‖` is invariant under `τ ↦ −τ` (since `‖−x‖ = ‖x‖`) and under any relabeling by
  `(ℤ/D)*`; the worst (tightest) sets are the ones whose binding structure is self-paired under
  `v ↦ N−v`, so the saddle sits where the two halves of each complement pair balance at `±M`.

The complement symmetry is the **balance condition**: the extremum is where no direction can be
improved because every speed has a mirror twin pulling the other way. That is the same reason SC
tournaments maximize `H` — perfect balance of the score sequence about its center.

## What this session added (and subtracted)

The session set out to test a *bounded-denominator* shortcut (Swing 3): that every covering 13-set has a
lonely point `τ = a/D` with `D` bounded (~41), which would close LRC(14) with no analysis. **It is
false** (HYP-2866): the witness denominator is *unbounded* — a far speed `w = 84m` divisible by the
small "AP-core-good" denominators (41, 53, 55, 65, 67, …) forces `D` up that ladder without limit. The
arithmetic of the far speed, not the tightness of the set, drives `D`.

But the *analytic* floor is real and robust (HYP-2867): `meas(GOOD ∩ G_P)` stays at 6.75–11.16× the
target for every small-part `P`. So the witness route (the measure floor + finite-V sampling) is the
necessary path, and the complement-symmetry is its organizing principle: the tightest covering sets —
the ones that actually stress the floor — are the complement-balanced ones, just as the H-maximizers
are the SC ones.

## The pointer beyond

If the LRC extremal locus is the complement-fixed locus, then the *right* coordinate for the uniform
floor is the one adapted to `v ↦ N−v` — the even/odd decomposition of the speed configuration under the
involution, exactly as the tournament invariants split into complement-even and complement-odd parts
(the merged metagraph `G_n/ℤ_2`). The conjecture to chase: **the uniform witness floor is controlled by
the complement-even part of the cluster alone**, the complement-odd part contributing a mean-zero
fluctuation that the L2/Parseval cancellation (HYP-2862) already handles. That would make the floor a
statement purely about the SC-projection of the speed set — and reduce the LRC hard core to its
self-complementary shadow, the same shadow the tournament project has been mapping all along.

— Related: [[lrc14-thread]], HYP-2571c (binding pairs sum to N), HYP-2866 (bounded-D refuted),
HYP-2867 (floor survives), HYP-2862 (cap-floor L2 duality); the SC spine in
`geometric-alignment-of-merged-metagraph.md` and `merged-metagraph-invariants.md`.
