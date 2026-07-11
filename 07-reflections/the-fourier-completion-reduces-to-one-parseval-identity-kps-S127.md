# The Fourier completion's analytic half is closed — one Parseval identity remains

*kind-pasteur-2026-07-11-S127. Owner: "work the measure-side closure, take the hardest available task, aim
to complete the formalization." A mapping pass (Explore agent) named the single biggest open sub-lemma; I
closed its analytic half and reduced the whole node to one clean algebraic identity. This note is the state
of that node.*

---

## The one open node on the measure route

A precise map of the LRC(14) measure-side surface found: the project is genuinely sorry-free, the entire gap
is a handful of named citation `Prop`s, and on the **Fourier / OffLine≤f(E3) route** exactly one link is
unproven — opus-S213's Stage **B.3**, the interval multiplicative pair-correlation bound

  `|C_w − b²/q| ≤ (q/4)·Σ_{h≠0} 1/(cdist h · cdist(w·h))  ≤  5q(log₂q+1)²/P(w).`

Everything around it is in Lean: opus's per-coefficient bound `‖B̂(h)‖ ≤ q/(2·cdist h)` (B.2), his
orthogonality `Σ_x e_q(hx) = q·1{q∣h}` (B.1), death-star's combinatorial `hyperbola_box_count` and dyadic
`harmonic_ratio_sum_mul_le`, my `offdiag_mcorr_sq_le`, opus's `E3_lt_choose`. B.3 is the sole missing joint,
and — the map's phrase — "proving B.3 lights up that whole chain."

## What B.3 actually decomposes into

B.3 is two independent pieces, and only one is analytic:

1. **The Parseval completion identity** (algebra): `C_w = b²/q + (1/q)·Σ_{h≠0} B̂(h)·conj(B̂(w·h))`.
   Fourier-invert `1_B`, expand the correlation `C_w = Σ_{x∈B} 1_B(wx)`, swap sums, collapse the inner
   character sum by orthogonality (B.1). Pure finite Fourier over `ℤ/q`.
2. **The analytic aggregation**: bound `‖Σ_{h≠0} B̂(h)conj(B̂(w·h))‖` using B.2 and the harmonic sum.

I closed **(2)**, kernel-pure, abstract in the coefficient function:

```lean
offDiag_bandSum_le      : ‖Σ_{h≠0} bc(h)·conj(bc(w·h))‖ ≤ (q²/4)·Σ_{h≠0} 1/(cdist h · cdist(w·h))
offDiag_bandSum_le_closed : … ≤ 5·q²·(log₂q+1)² / P            -- ∘ death-star's harmonic_ratio_sum_mul_le
```

The first is three moves — triangle inequality `‖Σ‖ ≤ Σ‖·‖`, multiplicativity
`‖bc(h)·conj(bc(w·h))‖ = ‖bc h‖·‖bc(w·h)‖`, and the termwise product of two B.2 bounds — with `w` a unit
keeping `w·h ≠ 0` so B.2 applies to both factors. The second composes death-star's `S·P ≤ 20(log₂q+1)²`;
the constants close exactly (`(q²/4)·20/P = 5q²/P`), and dividing the identity by `q` yields LEM-022's
target `5q(log₂q+1)²/P`.

## What this leaves — and why it is the honest boundary

With (2) done, **the entire Fourier-completion node reduces to the single identity (1)**. That identity is
not hard mathematics — it is textbook finite Fourier — but it is genuine *infrastructure* in Lean: a
character `e : ZMod q → ℂ`, its orthogonality over `ZMod q` (bridged from opus's `range q` / integer-`h`
form), Fourier inversion, the Fubini swap, and the unit reindexing `y = w·x`. Mathlib does not hand this
over ready-made for the concrete `exp` setup (which is exactly why opus proved B.1 by hand), so it is a
multi-session build, not a one-tactic call. I did not fake it, and I did not gamble a broken state on it
under budget — I banked the analytic half and scoped the rest precisely.

The abstraction of my lemmas in `bc` is deliberate: whoever proves identity (1) instantiates `bc = B̂` and
gets `|C_w − b²/q| ≤ 5q(log₂q+1)²/P` immediately, with no further analysis. The two halves compose at a
clean seam.

## The shape of progress at a genuine frontier

This is what advancing an open frontier looks like when you refuse to fake it: you don't close the node, you
*factor* it — separate the mechanical-but-infrastructural piece (the identity) from the genuinely analytic
piece (the aggregation), prove the one you can prove kernel-pure, and hand the other off with its exact
statement and derivation written down. The measure-side wall didn't fall this session. But the single
biggest attackable sub-lemma is now half its former size, and the remaining half is named to the theorem.

*Files: `LRCFourierAggregation.lean` (`offDiag_bandSum_le`, `offDiag_bandSum_le_closed`), building on opus's
`LRCFourierCompletionB` (B.2) and death-star's `LRCHyperbolaBox` (harmonic sum). The remaining identity is
opus-S213's Stage B.3 completion identity, derived above. Feeds `LRCMultCorrelation.offdiag_mcorr_sq_le`
(the per-cell `M` it assumes is `b²/q + 5qL²/P` once the identity lands).*
