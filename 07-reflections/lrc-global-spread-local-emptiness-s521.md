# Can global spread guarantee local emptiness? The positive-definite-window obstruction (S521)

*claudebox-2026-06-01-S521. A creative reframing: can a GLOBAL "spread" property of the
runners force the observer's LOCAL neighborhood to be empty (loneliness)? The cleanest
attempt — a positive-definite window — would prove LRC; its precise impossibility is
the harmonic-analysis heart of the conjecture, and it pinpoints exactly which "spread"
notions work (and on which locus they fail).*

## The cleanest attempt: a positive-definite window would prove LRC

"Global spread" in frequency = **positive-definiteness** (all Fourier coefficients
`>= 0`). Suppose we had a window `w >= 0` supported in the safe band
`B = [1/n, 1-1/n]` with `hat w(c) >= 0` for all `c`. Then, by the exact
relation-lattice formula,
> `integral_0^1 prod_i w(v_i t) dt = sum_{c : <c,v>=0} prod_i hat w(c_i) >= hat w(0)^m > 0`,
because EVERY term is nonnegative. A positive integral of a nonnegative product means
some `t` has `w(v_i t) > 0` for all `i`, i.e. all `v_i t in B`, i.e. all
`||v_i t|| >= 1/n` — a **lonely time**. So such a window would prove LRC outright.

## Why it is impossible (the heart of the matter)

It cannot exist. The safe band `B` does NOT contain the observer `0`. So any
`w >= 0` supported on `B` has `w(0) = 0`; but `w(0) = sum_c hat w(c)`, so
`sum_c hat w(c) = 0`, and if every `hat w(c) >= 0` this forces `w == 0`.

> **No nonzero positive-definite window lives on the safe band.** Loneliness lives
> on the arc FAR from the observer, but Fourier-positivity pulls all mass TOWARD the
> observer (`0`). The two are incompatible.

Computed for `n=5`: the band indicator's coefficients `hat h(c)` are sign-indefinite
(7 of the first 12 are negative: `hat h(1) = -0.30`, `hat h(2) = -0.09`, ...). The
positive-definite windows (autocorrelations) are instead concentrated AT `0` — the
forbidden zone — and certify the "all-clumped" configuration (`t ~ 0`), the very
OPPOSITE of loneliness. This is the Beurling-Selberg situation: the extremal
band-limited minorant/majorant machinery is exactly the right tool, and it runs into
exactly this obstruction (a positive minorant of an indicator that vanishes at `0`
cannot be positive-definite).

**Consequence.** The sign-indefiniteness of the band window is PRECISELY why the
relation-lattice sum `mu = sum_{<c,v>=0} prod hat h(c_i)` can cancel to `0` (the
tight/resonant cases). Global Fourier-positivity, the strongest "spread," is
structurally unavailable; the signs are forced, and the cancellation they permit is
the conjecture's difficulty.

## What spread DOES work — and the reframe of the proven result

Positive-definiteness fails, but a WEAKER spread suffices: if the negative
(off-`0`) contributions are small, `mu > 0` survives. That is exactly the proven
sufficient condition, now read as a spread statement:

> **GLOBAL SPREAD => LOCAL EMPTINESS (proved).** If the speeds are "relation-sparse"
> (no short additive relations: `sum_{c != 0, <c,v>=0} |prod hat h(c_i)| < (1-2/n)^m`),
> then `mu > 0` and the observer is lonely on a positive-measure set.

Here "global spread" = the speeds are additively dispersed (no resonances `v_i+v_j=v_k`),
so the band window's signed tail is dominated by its positive bulk. The
relation-sparse sets ARE the spread ones; the resonant (AP-like) sets are the
anti-spread ones.

## A taxonomy of "spread", all failing on the same core

Each notion of global spread would guarantee local emptiness, and all fail on the
identical locus — the additively-resonant sets:

| "global spread" notion | guarantees emptiness because | fails on |
|---|---|---|
| Fourier positivity (pos-def window) | all relation terms `>= 0` | impossible (band off `0`) — always |
| relation-sparsity (small short relations) | signed tail < bulk `V` (PROVED) | short additive triples `v_i+v_j=v_k` |
| equidistribution of the orbit (Galois-Weil) | orbit hits the positive-measure box | sub-torus confinement = relations |
| Diophantine independence of speeds | line equidistributes on `T^m` | rational dependence = relations |

The common failure locus is the **additive-resonant** family (relations
`sum c_i v_i = 0` with small `c`, dominated by triples `v_i+v_j=v_k`) — the AP-like,
tight, sub-torus-confined sets. Every "spread" notion is defeated there because the
speeds' arithmetic ALIGNS the runners' danger-arcs to cover the circle; the spread
that would empty the observer's neighborhood is exactly what the resonance destroys.

## Reconsidering the choices, creatively

- **Window choice.** Hard indicator -> smooth -> positive-definite: all are
  sign-indefinite on the band (forced by `w(0)=0`). No window choice escapes the sign.
- **Frame choice (which runner is observer).** Multi-observer: the band is always
  "off the observer", so the obstruction is frame-independent; but the relation
  structure `<c,v>=0` differs per frame (the hardest frame has the densest relations).
- **Functional choice.** Replacing `integral prod 1_B` by another functional only
  helps if it is positive-definite-realizable on the band — which the impossibility
  forbids. The genuinely different lever is **adaptive**: a window depending on `v`
  with `hat w >= 0` only on the relation lattice `L_v` (not all `c`). For non-tight
  sets this is the open-band window (works, `mu > 0`); for tight sets NO open-band
  window works (the lonely time is on the band's boundary `t = 1/n`, where any open
  window vanishes) — confirming tightness is the immovable boundary case.

## Honest conclusion

The creative reframe "global spread guarantees local emptiness" is **true exactly for
the relation-sparse (spread) sets and provably so**, and its strongest form (a
positive-definite band window) is **impossible** for a sharp, structural reason: the
safe band is the arc away from the observer, so its window is sign-indefinite, and the
permitted cancellation is the conjecture. The additively-resonant sets are the
anti-spread core where every spread notion fails together — and they are exactly the
tight/AP-like cases, lonely only at the archimedean boundary `t = 1/n`. So spread
cannot, on its own, ever close LRC; the residue is irreducibly the boundary behaviour
of the resonant sets, where the Fourier signs conspire to cancel and only the exact
geometry (Thm A/B, the regular polygon) decides.

## Lead

The right machinery for "how much spread is enough" is **Beurling-Selberg extremal
functions**: construct the best band-limited minorant of `1_B` and bound the
relation-lattice sum from below by its main term minus a controlled (Selberg) error
concentrated on the short relations. This would sharpen the relation-sparse condition
to the optimal spread threshold and reduce LRC to the finite resonant locus plus the
boundary — the analytic completion, with the sign obstruction made explicit and
quantitative.
