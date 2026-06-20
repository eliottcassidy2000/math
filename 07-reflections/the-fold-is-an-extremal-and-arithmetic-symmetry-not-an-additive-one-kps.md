# The converse fold is an extremal and arithmetic symmetry, not an additive one

**Source:** kind-pasteur-2026-06-20, the human owner's half-tiling framework (mirror over `y=x`
= reverse all arcs; half-tiling sizes `0,1,2,4,6,9,12,16,20,25,30`; even/odd recursions).
**Reads alongside:** mac-mini's `the-complement-folded-triangle-…` (the *geometric* fold +
cube-root/LRC-coverage `3+3+1`) and codex's `half-tiling-parity-address-quotient` (the address
quotient + the warning that `H` is **not** cell-affine). This reflection takes the third facet:
what the fold does to **extremes, parities, and counts**, where it turns out to be decisive.

## The fold, in one line

By THM-280, reflecting a fixed-HP tiling over `y=x` is exactly the converse `T -> T^op`
followed by the relabel `phi(i)=n+1-i`; on tiles it is the **pure** coordinate involution
`rho(a,b)=(n+1-b,n+1-a)` (no GF(2) bit-flip — arc-reversal and order-reversal cancel). Its
fixed cells lie on the anti-diagonal `a+b=n+1` (the SC spine), numbering `d=floor((n-1)/2)`;
its orbits number `half(n)=floor((n-1)^2/4)`. The half-tiling is its fundamental domain
(mac-mini THM-549, codex THM-550).

codex's warning is correct and worth restating: **`H` is not cell-affine** — the half-tiling
recursion `A+B-C(+D-E-F+G)` is a cell-count identity, and you cannot read `H` off the
half-region additively (THM-442). So at first glance the fold looks useless for the invariants
we actually care about. The point of this note is that the fold governs three *other* kinds of
structure — and there it is sharp.

## Layer 1 — Invariants: the fold is a free 2x (already known)

`H(T)=H(T^op)` and `c3(T)=c3(T^op)` (reverse a Hamiltonian path / a 3-cycle and you get one of
`T^op`). So every complement-invariant statistic is constant on `rho`-orbits and can be computed
on one side of the fold: a **2x** saving, the SC spine handled once (mac-mini). This is the
*invariance* layer — true but shallow.

## Layer 2 — Extremes: the maximizer localizes to the fixed locus (a 2^{(m-d)/2}, not 2x, win)

The deeper fact is about *where the maximum lives*, not its value. Exhausting the
`2^{half(n)}` grid-symmetric tilings directly (one free bit per `rho`-orbit), the **global
H-maximum is attained inside the fixed locus at every tested `n=3..9`**: `H_max = 3,5,15,45,189,
661,3357`, the grid-sym subcube reaching the true maximum each time (HYP-2688, independently
re-verified at n=8->661, n=9->3357). This is NOT the 2x of Layer 1 — it says the search can be
restricted from `2^m` to `2^{half}`, a **`2^{(m-d)/2}`** reduction (n=9: 4096x; n=14: ~7e10x).

But the mechanism is subtler than "symmetrize." `H=H^op` only makes the maximizer *set*
`rho`-invariant; a size-2 `rho`-orbit (a tiling paired with its mirror) has **no** fixed point,
and indeed the maximizer set is MIXED — at n=5,6,7 the non-grid-sym maximizers OUTNUMBER the
grid-sym ones (4/4, 12/18, 3/6). So the **strong** form ("all maximizers are grid-sym") is
**refuted**; only the **weak** form ("the set *contains* a grid-sym point") holds, which is
exactly what powers the speedup. The weak form is consistent with the canon SC/regular-maximizer
theorem (the global maximizer is regular and self-converse, so it admits a `phi`-self-converse
*relabeling*) — but translating that abstract self-converse symmetry into a grid-symmetric
*tiling* under the FIXED base path `P0` is an honest open gap. That is the one lemma between
HYP-2688 and a theorem (see the new OPEN-Q). The lesson: the fold's fixed locus is where the
extremum *sits*, even though `H` is blind to the fold additively.

## Layer 3 — Parity: the fold's fixed VERTEX pins c3 mod 2 (THM-552, proved)

`rho` (= `phi`) is an anti-automorphism, so it permutes directed 3-cycles; non-invariant ones
pair off, so `c3(T) ≡ #{phi-invariant 3-cycles} (mod 2)`. A 3-set is `phi`-invariant iff it
contains the `phi`-FIXED VERTEX `(n+1)/2`, which exists **iff `n` is odd**. Hence
**`phi`-self-converse tournaments have `c3` forced EVEN at even `n`, and reaching both parities
at odd `n`** (THM-552, proved; exhaustive over all grid-sym tilings n<=8). This corrects a
tempting but wrong reading of the geometry: even `n` is NOT distinguished by "missing the central
fixed cell" (the apex `(n,1)` is a fixed cell at every `n`) — it is distinguished by the fold
having **no fixed vertex**. The odd/even shape split (square `k^2` vs pronic `k(k-1)`) is the
*cell-space* shadow of this *vertex-space* fact: the lone median vertex.

## Layer 4 — Counts: the Burnside census (and a canon correction)

`rho` acting on the `2^m` tilings has exactly `2^{half}` fixed points (the grid-sym tilings), so
by Burnside the number of `rho`-orbits of tilings is `(2^m+2^{half})/2 = 2,6,40,544,16640,
1050624` (n=3..8) — the **fixed-HP-frame analogue of `V_merged=(A000568+SC)/2`**: the same `Z_2`
Burnside construction, on the labeled base-path-fixed cube instead of the `S_n`-quotient (so it
is much larger; equal only at `n=3`). In the full labeled space the anti-fixed count is
`Fix_anti(phi)=2^{half(n)+floor(n/2)}`, the extra `floor(n/2)` being the `phi`-self-paired
vertex-pairs `{u,n+1-u}` — NOT the `n-1` base-path arcs one might first guess.

Chasing the SC-halving bijection here turned up a canon error: the identity
`|SC(n)|=A000568(n-1)` stated in `two-models-staircase-recursion.md` is **false** — the true
self-converse counts are `2,2,8,12,88,176` (matching `unlocking-gn-at-all-n.md`), agreeing with
`A000568(n-1)` only at `n=4,6` (MISTAKE-081). `2^{half}` is a labeled fixed-point count, never an
iso-class count; the half-tiling supplies the Burnside *input*, not the class-level bijection
(still open).

## The thesis, and the engineering payoff

Put the layers together. The converse fold is **invisible to additive structure** (codex: `H` is
not cell-affine) but it **controls**:
- the **argmax** of `H` (Layer 2, conjecturally — the maximizer sits on the fixed locus),
- the **parity** of `c3` (Layer 3, provably — via the fixed vertex),
- the **census** of tilings (Layer 4, provably — Burnside `(2^m+2^{half})/2`).

So "the fold does nothing for `H`" is exactly wrong: it does nothing for `H` as a *linear/additive*
function on cells, and almost everything for `H` as an *extremal* object and for the *arithmetic*
of the configuration space. Symmetries that vanish from the value can still rule the extreme and
the count.

The clean engineering corollary: a self-converse tournament is determined by its `half(n)≈n^2/4`
orbit bits, so it compresses **2x** (vs the `C(n,2)` adjacency) with an exact round-trip, and a
uniform self-converse sampler is `O(half(n))` bits with no rejection — the `HalfTilingCodec`
(verified exhaustively n<=7). The fixed locus that holds the maximizer is also the natural
canonical fingerprint for the entire SC stratum of any antisymmetric pairwise-preference system.
