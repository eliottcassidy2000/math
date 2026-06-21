# The LRC cycle space reinforces where the tournament's cancels

**Session:** mac-mini-2026-06-21-S9. The question was: how does the multi-dimensional Weyl error —
the one analytic gap left in the LRC(14) sector route — relate to tournaments? The answer is a clean
structural correspondence with one sharp surprise that redirects the whole program.

## The correspondence is real, and it is the cut/cycle seam

The carrier error is `corr(E) = measS7(E) − iid_k = Σ_{0≠n ∈ Λ(E)} K(n)`, a signed Fourier sum over
the **offset relation lattice** `Λ(E) = {n : Σ nᵢeᵢ = 0}`. That lattice is the LRC twin of the GF(2)
**cycle space** of `K_n` — both are the dependency kernel of the generators (offsets ↔ arcs). The
decorrelated `iid_k` is the single-particle **cut** side; `corr` is the **cycle** side. And the
stratification is by **support size**: support-2 relations form the 2-body cut layer (always ≥ 0,
small, the exact twin of THM-559's `c3 = 2-body line-graph Ising`), while support ≥ 3 is the
many-body cycle layer that carries the binding error. The leading cross-block term is a **support-3
additive triangle** `a + b = c` (e.g. `65 = 1 + 64`) — the LRC twin of a directed 3-cycle. So "cut
cheap, cycle dear" reappears, intact, as "support-2 cheap, support-≥3 dear." That much transfers
perfectly.

## The surprise: the OCF mechanism is the wrong sign

It is tempting to push further: the OCF keeps **odd** cycles and drops **even** ones, because a
directed cycle and its reverse carry opposite signs and cancel at even length (THM-560). The obvious
conjecture is that `corr` is likewise odd-support-dominated, an OCF on the offset lattice. It is not.
One dissociated set looked 98 % odd and nearly sold me on it — but the adversarial population kills
it: powers of two give even/odd = 2.0, Sidon sets 0.6, and a random population is split down the
middle. The reason is an exact, provable identity:

```
ĉ_T(−n) = conj(ĉ_T(n))  ⟹  K(−n) = conj(K(n))  ⟹  Re K(n) is EVEN under n ↦ −n.
```

A relation and its reverse **reinforce** — they add `2 Re K` — where a directed cycle and its reverse
**cancel**. The LRC cycle space is *conjugate-symmetric*; the tournament cycle space is
*reverse-antisymmetric*. Same lattice, opposite involution. So the OCF's defining move — drop the
even cycles by reverse-cancellation — has no LRC counterpart, and a reverse-pairing bound is exactly
the wrong tool for the signed Erdős–Turán estimate the carrier gap needs. The cancellation that
makes the absolute bound fail is real, but it lives in *cross-support sign mixing*, not in reverse
pairs. This is worth more than a positive analogy would have been: it tells everyone working the
gap (codex was about to start exactly this odd-support route) which door is bricked up.

## What does transfer, and the one lever

Two things survive. First, the cut/cycle = support-size seam: the program's existing support-based
covolume / relation-height machinery (THM-532/533, angle-F) is the right frame, and THM-559's
2-body Ising is literally the support-2 floor of `corr`. Second — the genuine new lever — **additive
energy is the multi-block splitting parameter.** Separating two blocks destroys their cross-block
additive triangles (a Schur triple `a+b=c` across blocks needs the far offset below `2w−2`), which
floors the cross-block additive energy and pins `|corr|` of a separated multi-block well below the
single coherent block's value: the consec block sits at 0.303, and *any* gap drops it to ≤ 0.093,
dissociation to 0.013. So the multi-block carrier gap — the last open regime — **reduces to the
single-block bound**, the one the program already targets. The triangle that is a 3-cycle is also a
Schur triple, and killing it is what separation does.

## The shape of the remaining wall

Everything now funnels to one place. With multi-block reduced to single-block, single-far closed,
and the dissociated tail handled, the entire route rests on consec maximizing the cover over the
bounded family — and through the cut/cycle lens that is the statement that the **support-3
Schur-triple terms, reinforced and phase-aligned mod 7, are largest at the arithmetic progression.**
It is the same mod-7 alignment that HYP-2704's slow/fast band trade kept pointing at, now wearing
the cycle-space costume. The wall has not fallen. But it is the same wall seen from a fourth side,
and each side has told us the same thing: the difficulty is the reinforced, arithmetic, many-body
sum — the cycle space doing what, in the tournament, the OCF does, but adding where the OCF cancels.
