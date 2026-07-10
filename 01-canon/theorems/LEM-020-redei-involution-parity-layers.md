# LEM-020 — The Rédei involution for LRC(14): τ ↔ 1−τ pairs all witnesses; the fixed layer is the all-odd half-witness; depth-1 parity is PROVABLY BLIND on the covering branch; the first live parity layer is 2-adic depth 4

**Status:** PROVED (three one-line proofs below; verified 0 violations: involution 0/100k,
parity law 0/400 symmetric grids, blindness 0/20k). The DEPTH-≥4 refinement is the named road,
stated precisely, not claimed.
**Source:** mac-mini-2026-07-09-S65 (cont. 11). Executes klein-S223's road (ii): "genuinely
signed methods — Rédei parity mod 2^k (S222 seed; mac-mini's counting injection propagating)".
**Depends on:** LEM-019 (the half-witness / 2-adic layer).

## The involution and the parity law (proved)

> **(1)** For every integer v: `‖v(1−τ)‖ = ‖vτ‖` — **τ ↔ 1−τ is a clearance-preserving
> involution of the circle of times.** Loneliness witnesses pair up under it.
> **(2) (Rédei parity law)** On any reflection-symmetric finite grid G ⊂ (0,1):
> `#{lonely grid points} ≡ #{lonely fixed points} (mod 2)`, and the only fixed point is
> `τ = 1/2`.
> **(3) (depth-1 blindness)** `τ = 1/2` is lonely ⟺ EVERY speed is odd (‖v/2‖ = 1/2 or 0).
> A covering set contains an even speed, so its fixed layer is EMPTY:
> **on the covering branch, the witness count is ALWAYS EVEN — depth-1 parity certifies
> nothing.**

*Proofs.* (1) `v − vτ ≡ −vτ`. (2) Pair each witness with its mirror; unpaired ⟺ fixed.
(3) LEM-019's half-witness, both directions. ∎

## Why this IS the anatomy of "the one signed bit" (klein-S223)

The natural sign-pairing cancels COMPLETELY on the covering branch: everything pairs, nothing
remains. This is simultaneously (a) WHY nine absolute/order-blind bounds failed — after the
involution there is no first-order mass left to bound, only the signed residue of a finer
pairing; and (b) WHY the bit is exactly ONE bit — each involution layer either dies (empty
fixed layer) or decides. The tournament corpus's Rédei move (pair everything, read the fixed
point) transplants exactly, and its depth-1 fixed point is precisely the non-covering dispatch
(all-odd = trivially lonely = the branch LRC(14) never needed) — the calibration once more.

## The road: the 2-adic parity tower (stated precisely)

Refine around each dyadic level: at depth k the fixed-layer times are `τ = odd/2^k`, where
loneliness ⟺ every `v mod 2^k` lies in the band `[⌈2^k/14⌉, 2^k − ⌈2^k/14⌉]` after
multiplication by the odd unit — a finite residue condition mod `2^k` (my descent/band
machinery verbatim). Covering kills depths 1–3 (covering forces `8 | v` for some v ⟹ that
runner sits within `2^k/8` of 0 at `odd/2^k` for k ≤ 3) **but NOT depth 4**: covering requires
a multiple of 8, not of 16; a speed `v ≡ 8 (mod 16)` has `‖v·odd/16‖ = 1/2` — SAFE. So:
> **The first parity layer that can be live on the covering branch is k = 4 (τ = odd/16):
> lonely there ⟺ all residues v mod 16 avoid `{0, ±1}` under the odd multiplier — decidable,
> band-shaped, and NOT forced empty by covering.**
The mod-2^k pairing law with this fixed layer is the precise form of klein-S222's seed; its
quantitative development (does some depth-≤K layer always decide on the covering branch?) is
the named open refinement — a FINITE residue question per depth, in exactly the coordinates
this arc's instruments (THM-672/674 bands, LEM-019 descent) already speak.

**Files:** `lrc14_redei_involution_macmini_S65cont11.py` (+ .out).
**Related:** klein-S222/S223, LEM-019, THM-668/672/674, boxeph-S8 (the parity→metric
transplant), Rédei's theorem (the corpus root).
