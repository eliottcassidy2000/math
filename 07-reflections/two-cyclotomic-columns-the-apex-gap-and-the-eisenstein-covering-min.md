# Two cyclotomic columns: the apex gap (Q(cos 2π/p)) and the Eisenstein covering-min (Q(√−3))

*klein-2026-06-29-S23. Working the owner's reframe: the open conjecture binds off-cusp at n/Φ₆(n)=14/183, the Q(√−3) "existence column." The reframe slots cleanly against the recent apex-gap program and reveals two complementary cyclotomic columns.*

## The reframe's sharpest point: the inequality is trivial

`n/Φ₆(n) ≥ 1/n` unfolds to `n² ≥ n²−n+1`, i.e. `n ≥ 1`. So the analytic inequality the conjecture "reduces
to" carries **no content at all**. Everything is in the word "covering-min": the claim that `n/Φ₆(n)` is the
*tightest* off-cusp covering, that **no covering beats the construction**. This is the single most useful
thing the reframe does — it strips away the analysis and leaves a pure combinatorial-optimality statement.

## What Φ₆ is

`Φ₆(n) = n²−n+1` wears three hats, all the same hat:
- the **Eisenstein norm** `|n − ζ₆|²` (`ζ₆ = e^{iπ/3}`) — the field is `Q(√−3)`;
- the **point-count of `PG(2, n−1)`** (`= (n−1)² + (n−1) + 1`);
- hence the size of a **`(Φ₆(n), n, 1)` Singer difference set** (when `n−1` is a prime power).

For `n=14`: `Φ₆=183=3·61` (`61 ≡ 1 mod 6` splits in the Eisenstein integers, `3` ramifies), and `q=13` is
prime, so `PG(2,13)` exists — `14` speeds mod `183`, every nonzero difference once. The covering construction
IS the projective plane.

## Two columns

The recent program is really two complementary cyclotomic columns, and the open conjecture lives in the
first:

| | **covering / existence column** | **apex / measure column** |
|---|---|---|
| field | `Q(√−3)` (Eisenstein) | `Q(cos 2π/p)` |
| object | the Singer difference set `(Φ₆(n),n,1)` | the odd cycle `C_p` (the doublet) |
| atlas pole | the **SPREAD** pole (flat spectrum) | the **CONCENTRATED** pole (the gap) |
| value | `n/Φ₆(n)` (the covering-min) | `4 sin²(π/2p)` (the gap) |
| regime | OFF-cusp, positive measure, M tightest | the descent **attractor**; vanishes at the cusp |
| recent work | this (HYP-3705) | mac-mini-S41/HYP-3700 (edge isolation) |

The 2-adic descent renormalizes toward the **doublet attractor** (concentrated, apex). But the *conjecture's
binding* is not there — it is at the **difference set** (spread, Eisenstein), off the measure-0 cusp, where
the measure is positive and `M` is tightest. Two poles of the same multi-axis atlas (HYP-3610): the descent
end (doublet) and the covering end (difference set); two fields; two columns.

## Where the difficulty actually is

A `(v,k,1)` symmetric design — a projective plane / Steiner system `S(2,k,v)` — is the **optimal covering
design**: the covering number is attained exactly, no pair covered twice, no redundancy. So **"no covering
beats the construction" is TRUE at the combinatorial-design level.** The genuine open content is therefore
narrow and precise: the **bridge** from the LRC *continuous* covering floor `M` to the *discrete* design
covering number. Does the lonely-runner covering inherit the projective plane's design-optimality?

That is a cleaner open problem than any inequality. It says: the lonely runner at `n=14` is hard not because
of an analytic estimate (trivial) and not at the measure-0 cusp (existence handles that — the `φ(n)` units),
but because of a single combinatorial-optimality question about the projective plane `PG(2,13)`, sitting in
the Eisenstein field `Q(√−3)`, off the cusp, at the most-spread (difference-set) end of the core lattice.

## The picture

The mathematics keeps splitting into a concentrated/discrete/apex face (`Q(cos 2π/p)`, the doublet, the
descent attractor, the measure-vanishing isolated edge) and a spread/covering/Eisenstein face (`Q(√−3)`, the
Singer difference set, the off-cusp covering-min). The conjecture binds on the second. And on the second the
inequality is free; only the optimality of the projective plane remains — a finite, design-theoretic truth
whose only gap is its bridge to the continuous covering.

See HYP-3705 (this), HYP-3700 (mac-mini: the apex column's edge isolation), HYP-3610/3611 (the multi-axis
poles), HYP-3604 (the doublet/apex gap), HYP-3597 (off-cusp positive measure vs the cusp), THM-590.
