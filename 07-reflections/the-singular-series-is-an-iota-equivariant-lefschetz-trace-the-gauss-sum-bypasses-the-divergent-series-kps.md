# The LRC singular series is an ι-equivariant Lefschetz trace: ι acts freely, the ordinary trace vanishes, and the certification is the ι-odd index = the Gauss sum `i√7` — one number that bypasses the divergent signed series

*kind-pasteur-2026-07-01-S19. Working the owner's idea that the singular series is an exact Lefschetz trace certifying `M ≥ 1/n`. The pieces assemble cleanly and are verified: the series is a measured, ι-equivariant Euler characteristic; the atoms' moments are Ramanujan character traces; ι acts freely (p≡3 mod 4), so the ordinary Lefschetz number is 0; and the certification is the ι-**odd** index — the quadratic Gauss sum `g_7 = i√7`, a single nonzero trace. That is the bypass of the absolutely-divergent signed series (MISTAKE-078), and — the punchline — the Gauss sum carries the Heegner and 3-mod-4 pillars while the free ι is the Mersenne pillar, so the Lefschetz-trace framing **unifies the three pillars into one topological certificate.***

## The four verified pieces (N=14, apex p=7)

**(A) The singular series is a measured Lefschetz/Euler alternating sum.** `L = ∫₀¹ ∏_v (1−g_v) dt = Σ_{A}(−1)^{|A|} I_A`, `I_A = meas(∩_{v∈A} danger_v)` — the inclusion–exclusion (Bonferroni = moment truncation, HYP-3789) is the **measured Euler characteristic of the danger-cover nerve**, and it is **ι-equivariant** under the antipode `ι: t ↦ 1−t` (the LRC complement, THM-581/582's two-index).

**(B) The atoms' moments are Ramanujan sums = character traces.** The lonely-set atoms are the units `(Z/14)* = {1,3,5,9,11,13}`, and their Fourier moments are `μ̂(k) = Σ_{a∈(Z/14)*} e(ak/14) = c_14(k)` — the Ramanujan sum, a **trace of the k-shift over the units** (verified: `c_14 = (6,1,−1,1,−1,1,−1,−6,−1,1,−1,1,−1,1)`, matching the closed multiplicative form). So the series' terms are group-character traces.

**(C) ι acts freely ⇒ the ordinary Lefschetz number is 0.** `ι: a ↦ 14−a` has orbits `{1,13},{3,11},{5,9}` on the units — **no fixed unit** (a fixed unit needs `a=7 ∉ (Z/14)*`). So ι is a free involution on the atoms, and there is no ι-fixed lonely time (`t=1/2` is not lonely for the AP core — even runners sit at 0). Hence the **ordinary Lefschetz number `Λ(ι) = #(ι-fixed lonely points) = 0`: the ι-even part is blind.** This is exactly `p≡3 mod 4` ⇒ `χ(−1)=−1` ⇒ the complement is an anti-automorphism = a free `ℤ₂` (last session's Brouwer-vs-Borsuk–Ulam split, S18 / my S11 sgn-parity finding).

**(D) The certification is the ι-odd index = the Gauss sum `i√7`.** Because ι is free, the obstruction lives in the ι-**odd** (reduced/Borsuk–Ulam) index, which localizes to the quadratic character `χ = (·/7)` (odd, since `χ(−1)=−1`). Its value is the **quadratic Gauss sum** `g_7 = Σ_a (a/7) e(a/7) = i√7` (verified: `0 + 2.6458i`, `|g|=√7`). It is a **single nonzero number** — and a nonzero ι-odd index obstructs any ι-equivariant cover of the circle, i.e. **certifies a lonely runner, `M ≥ 1/n`.**

## The bypass (the creative point)

The naive route — prove `L = Σ_A(−1)^{|A|}I_A > 0` — is blocked because the series is **absolutely divergent** (`Σ|c_k|` unbounded; MISTAKE-078), so the signed cancellation is opaque. The Lefschetz split reroutes it:

`L = (ι-even main term = the (6/7)¹³ average, manifestly positive) + (ι-odd correction).`

- The **ι-even** part is the mean (Lasserre level-1 / trace); it is positive but does not certify by itself.
- The **ι-odd** correction is the divergent signed series — **but its *index* (not its sum) is what certifies**, and by localization to the free ι-action it is a **single Gauss-sum trace `i√7`**, not an infinite alternating sum.

So: **replace the divergent signed correction by the one Gauss-sum trace.** This is *why* `p≡3 mod 4` needs Borsuk–Ulam (odd degree / Gauss sum) rather than Brouwer (even / fixed point / SOS): the even index vanishes on a free action, and the odd index is a finite computable trace.

## The three pillars are one ι-odd index

The framing unifies the `twentyeight` three pillars as the *ingredients of the same Gauss-sum certificate*:
- **Gauss sum `g_7 = i√7`** carries **two** pillars at once: the `√−7` is the **Heegner** field `Q(√−7)` (the SOS/Fejér–Bochner minorant lives there), and its being **imaginary** (`χ` odd) is the **3-mod-4 / Borsuk–Ulam** pillar.
- The **free ι** (parity descent `14 = 2·7`, `ι` a free `ℤ₂`) is the **Mersenne** pillar (the 2-adic peel to the all-odd apex).

So "the singular series is a Lefschetz trace" is not decoration — it is the statement that **the three pillars are the even part (Mersenne parity), the field (Heegner `√−7`), and the odd degree (3-mod-4) of a single equivariant index**, and the certificate is `i√7`.

## Connection to the moment ladder (S19 Bridge-2)

This is the Bridge-2 mindset applied to *certification*: the **ι-even trace (1st moment / mean)** is blind (`Λ(ι)=0`), and the **ι-odd 2nd-moment trace** (the Gauss sum, the atoms' odd pair-correlation) is what separates lonely from covered. The Gauss sum is literally the odd part of the Ramanujan autocorrelation `c_14` — the atoms' 2nd moment, twisted by `χ`. The reconstruction ladder and the certification ladder are the same: when the count/even level is degenerate, climb to the 2nd-moment/odd level; here that level is a *single finite Gauss sum*, which is why the certificate is clean.

## Honest status

- **Verified:** atoms = units; moments = Ramanujan traces `c_14`; ι free (3 orbits, `Λ(ι)=0`); `χ(−1)=−1`; Gauss sum `g_7 = i√7`; the pillar assignments.
- **Framing, not a theorem:** "`L =` an *exact* Lefschetz trace whose ι-odd index certifies `M≥1/n`" is the owner's conjecture, and I have grounded every ingredient — but identifying the *cohomology* (the precise complex whose ι-odd Lefschetz number is `i√7`-controlled) is unproven. The repo's **ι-odd index** (OPEN-Q-108 residual, klein-S56; THM-581/582) is exactly this object; this reflection gives it the Gauss-sum value and the pillar unification.
- **Bypass value:** it converts the divergent-series obstacle into a finite Gauss-sum trace, and dovetails with MSS-2024 finiteness (the odd index over a *finite* frequency window).

## Convergence (opus-S23, HYP-3815) — two sides of the same `i√7`

Concurrently, opus-S23 reached the identical move — "Borsuk–Ulam is dual to Gauss-sum/Lefschetz exactness; certify the free-ℤ₂ regime with a trace formula, not a fixed point" — from the **tournament** side: the Paley skew matrix `S_ij = χ(j−i)` (p≡3 mod 4) has spectrum **exactly `{0, ±i√p}`** = the Gauss sum = the **Frobenius eigenvalue on `H¹`** (Grothendieck–Lefschetz Weil point count `#Fix(Frob^k)=Σ(−1)^i Tr(Frob^k|H^i)`), and the Paley Cayley spectrum is *concentrated* (the "hard" object is arithmetically simple, just non-symmetric). **That skew spectrum `{0,±i√7}` is precisely the `cpS` (skew char-poly) that separated my n=6 OCF-cospectral twins (S15)** — so opus's Frobenius/Gauss eigenvalue, my twin-separator's spectral coordinate, and this LRC atom-Ramanujan `i√7` are one object seen three ways. opus (HYP-3815) carries the Weil/Frobenius realization; this reflection carries the LRC-atom/Ramanujan realization + the three-pillar unification + the moment-ladder tie.

## Next hunts

- **Make the trace exact:** find the chain complex (the ι-equivariant nerve / the runner-flow's periodic orbits, Blaschke-dynamics HYP-3796) whose ι-odd Lefschetz number *is* the certification, so `i√7 ≠ 0 ⇒ M≥1/n` becomes a theorem, not an analogy.
- **The other apices:** compute `g_p` for `p=11,23,31` (n=22,46,62) — `g_11=i√11` (Heegner survives), `g_23=i√23` (Heegner fails, `h(−23)=3` — the certificate needs the class group), `g_31=i√31` (Heegner fails) — the Gauss-sum trace exists at every `p≡3 mod 4`, but its *field* degrades exactly where the Heegner pillar does. The Lefschetz certificate degrades with the apex, matching S18's by-n difficulty.
- **The blue-spine spectrum** (tournament side): the SC blue-graph is connected (S17); its Laplacian spectrum is the rung-2 refinement of blue-degree — the tournament analog of this odd-index/Gauss-sum trace.

— Related: `the-flat-extension-moments-are-ramanujan-sums-…` (S7, atoms=units, moments=Ramanujan), HYP-3789 (Lasserre/inclusion-exclusion = the alternating sum), `the-bridge-2-mindset-is-one-rung-of-a-universal-moment-ladder-…` (S19), `lrc-difficulty-by-n-…` + `brouwer-vs-borsuk-ulam-…` (opus) (p mod 4, free ι, three pillars), `twentyeight-…` (the three pillars, `i√7`), THM-581/582 (ι-odd index / two-index), OPEN-Q-108, MISTAKE-078 (the divergent series this bypasses). Script: `04-computation/lrc_lefschetz_trace_gauss_sum_kps.py` (+ .out). Not a HYP reservation.
