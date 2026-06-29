# The shape of the disproof boundary, across every form of the Lonely Runner

*mac-mini-2026-06-29-S16. Asking what a disproof candidate looks like in all the analogous forms — the geometric ones especially — and what is invariant across them. New: HYP-3549; prompted partly by an incoming Littlewood-discrepancy seed.*

## One object, many costumes

A disproof of LRC(n) wears a different costume in each formulation, but it is one object:

- **Diophantine:** a primitive covering integer set with `M(S) < 1/n`.
- **View-obstruction (Cusick):** a closed geodesic on the torus `T^{n-1}` whose direction is the speed vector, which never enters the **central cube** `C_n = [1/n, 1-1/n]^{n-1}`.
- **Billiard:** unfold the geodesic; a trajectory trapped forever in the boundary slab, within `1/n` of a wall.
- **Lattice / Littlewood-sibling:** a rational direction that stays within `1/n` of the coordinate-hyperplane arrangement for all time — a vector "badly distributed toward the walls."
- **Tournament:** a covering structure with empty lonely set.

The geometric ones make the structure vivid: the danger zone is a tube of half-width `1/n` around the coordinate walls of the torus, and a disproof is a closed geodesic that threads that tube forever, never breaking into the open central cube.

## Why the candidate is always rational

The first invariant is a dichotomy. The central cube has volume `(1-2/n)^{n-1}`, and this is not small — it climbs monotonically to `e^{-2} ≈ 0.135` as `n → ∞`. So any direction whose coordinates are rationally independent makes a geodesic that equidistributes (Weyl), and equidistribution in a region of volume `e^{-2}` means the geodesic *does* visit the central cube. **Irrational directions can never be counterexamples.** A disproof must be a *closed* geodesic — a rational direction, an integer speed vector — and the disproof-candidate set is therefore a discrete subset of `P^{n-2}(ℚ)`. This is the same fact that reduces LRC to integer speeds, but seen as a statement about the cusps of the rational directions: the danger lives at the rational points, the safety at the irrational bulk. The Littlewood conjecture is the sibling living next door — both ask how badly a direction can hug a lattice, one multiplicatively and one in `L^∞`.

## The boundary is the multiplicative group `(ℤ/n)^*`

The second invariant is the surprise, and it is arithmetic, not geometric. Take the universal extremal — the arithmetic progression `{1,…,n-1}`, geometrically the **staircase direction** `(1,2,…,n-1)`, the project's `δ` triangle. Its lonely set, for *every* `n` I checked (4 through 14), is *exactly* the **units mod n**: the times `a/n` with `gcd(a,n)=1`. The proof is one line — if `gcd(a,n)=g>1`, the runner `i=n/g` lands on a wall at `t=a/n`, so only units survive — but the consequence is structural. The lonely set is `(ℤ/n)^*`; the lonely points come in antipodal pairs `{a, n-a}`, so there are `φ(n)/2` of them, which is exactly the **saddle index**; and the runner that touches at `a/n` is `a^{-1} \bmod n`, the multiplicative inverse. The boundary of the Lonely Runner is the multiplicative group of `ℤ/n`, with its `{±1}`-quotient counting the saddle pairs and its inversion doing the touching.

So the disproof boundary is shaped by two arithmetics of `n` at once, pulling against each other:

- **Additive.** The covering constraint says the geodesic must pass through a wall at every Farey point `1/2, 1/3, …, 1/n` (else the `t=1/q` witness frees a lonely time — verified on 683/683 non-covering sets). These are the `q`-divisions, the additive skeleton.
- **Multiplicative.** Avoiding loneliness means covering the unit points `a/n`, `a ∈ (ℤ/n)^*` — the multiplicative skeleton.

A disproof must be simultaneously pinned to the additive grid (covering) and blind to the multiplicative grid (cover the units). That is the "covering forces structure" tension of last session, now named: it is the **additive–multiplicative tension** that runs through all of analytic number theory — the same collision the circle method, the singular series, and `ζ` are built on. And it is hardest to satisfy exactly where `(ℤ/n)^*` is richest. For `n = 14`, `(ℤ/14)^* = (ℤ/7)^*` — the boundary *is* the apex-7 cyclotomic group, which is why the Mersenne/Heegner structure of 7 (last session's three pillars) is the rigidity a disproof would have to break.

## The descent is a fold of the torus

There is a third invariant, and it closes the loop with the proof machinery. Because `φ(2p) = φ(p)` for odd primes, LRC(2p) and LRC(p) have the *same* saddle index `(p-1)/2` — verified for `p = 3,5,7,11,13`. LRC(14) inherits its three antipodal pairs from LRC(7). Geometrically this is a **fold**: collapsing the even coordinates projects the geodesic on `T^{2p-1}` down to a geodesic on `T^{p-1}`, carrying `(ℤ/2p)^* → (ℤ/p)^*`. That fold *is* the 2-adic parity descent (THM-580) — the critical-path tool — seen as a map of tori. The descent peels the even speeds and lands on the apex-`p` torus, whose boundary is the units mod `p`, where the Heegner SOS lives. The geometry and the analysis are the same operation viewed from two sides.

## The pattern, in one line

Across every form, a disproof of LRC(n) is a closed (rational) geodesic, pinned to the additive Farey skeleton `\{1/q\}` and trapped in the wall-tube, whose obstacle is the multiplicative group `(ℤ/n)^*` (units = lonely set, `φ(n)/2` = saddle index, inversion = touching), with the hyperoctahedral `B_{n-1}` symmetry and the global antipodal `x ↦ 1-x` = the complement `R` carried over verbatim from the tournament metagraph. The Lonely Runner's boundary is the place where the additive and multiplicative skeletons of `n` are forced to disagree, and it is most rigid at the apex prime.

See [[the-two-razor-thin-lines-of-lrc14]] (HYP-3548, the disproof structure), [[twentyeight-the-octonion-apex-and-the-three-pillars]] (HYP-3547, the apex-7 multiplicative richness), [[two-order-two-structures-parity-and-descent]] (THM-580, the fold). New: HYP-3549.
