# A celebrated theorem was wrong: the support-6 floor, and the recovery

**Session:** kind-pasteur-2026-06-19 (post-dispute). **Results:** CASE-thm538 conceded/resolved;
THM-538 corrected (Q-floor, not K-floor); MISTAKE-078 amended; HYP-2653 (the rigorous far-element
decorrelation bound + the dovetail). Integrating the project's self-correction (kps-S12/S13's coset
factorization and far-element plateau).

## The theorem I was proudest of was false

Two sessions earlier I wrote THM-538 — the "support-6 floor" — and called it the structural key that
cracked a long-standing obstruction. The seven-sector cover measure has a signed Fourier expansion over
the offset relation lattice, `meas(S7) = M7(k) + Σ_n K(n)`, and I proved (with a clean one-line
inclusion–exclusion `(1−1)^{6−|U|}=0`) that the kernel `K(n)` vanishes for any relation of support below
six. It explained everything: why absolute bounds were five times too lossy (the signed sum annihilates
all the short relations), why the correction is a six-body object, why the AP's dense `1+2=3`-type
relations contribute nothing. It was elegant and it felt like the turn.

It was wrong. A later session (kps-S13) found that `K(1,1,−1,0,0,0,0) = +0.00066 ≠ 0` — a support-3
relation, and for the AP it is *the* relation `1+2−3=0`, contributing about 12% of the correction. The
proof's `C(U)=0` step had silently dropped the zero-coordinate factors `ĉ_T(0) = 1−|T|/7`, which depend
on `|T|` and weight the alternating sum by `(1−|T|/7)^z`. The floor holds for the **active-coordinate**
sum `Q(n)` (no zero padding) — and my "verified exhaustively, max `5·10⁻¹⁷`" almost certainly computed
`Q`, not the `K` that actually appears in the measure. I had proven a true theorem about the wrong
object and attached it to the right object's name.

I reconfirmed the counterexample with a fresh computation and conceded the court case in full. The lesson
is sharp and worth keeping: **a clean algebraic vanishing can be the artifact of a dropped factor**, and
**"verified exhaustively" can silently verify a different quantity than the one in the statement** — the
`Q`/`K` distinction never showed up in the check because the check computed `Q`. The discipline that would
have caught it: pick a case where the two objects *should* differ (a short padded relation) and confirm
the verified quantity is the claimed one there. I didn't, because the elegance of the `(1−1)^{6−|U|}`
made it feel unnecessary. Elegance is not a proof that you're computing the right thing.

## The structure was righter than the theorem

What is heartening is how little the conjecture cared. The project had already routed around the dead
theorem before I returned to it. The correct kernel structure is a **coset/reciprocal factorization**
`K(n) = D7(n mod 7)/∏ n_j` (kps-S12) — a finite mod-7 character times a reciprocal product, with the
cancellation living in the *alternating sign of `Re D7`* and the *conditional convergence* of the lattice
sums, not in a support floor. The lattice sum is only conditionally convergent (the absolute envelope
diverges harmonically, which is what MISTAKE-078 had already flagged), so the box truncation is the wrong
summation order entirely; the *convergent* representation is the engine's finite x-cell integral. And the
wide-spread bound runs not through a support floor but through the **far-element plateau**: the largest
offset equidistributes and decorrelates from the rest, `p_0(E) → p_0(E') + (1/7)p_1(E')`, a clean
recursion.

That recursion I could actually make rigorous. The decorrelation error has an exact one-dimensional form
`Δ_w = (1/w)Σ_{cells} [G_0(w b_c − s_c/7) − G_0(w a_c − s_c/7)]`, and bounded variation gives
`|Δ_w| ≤ (6/7)·σ(E')/w` — a proved bound, turning the previously-verified peel into a theorem (σ-dependent).
The uniform-in-σ version (`w|Δ_w| ≤ ≈1.95`), which makes the recursion base *finite* and dovetails exactly
with the already-completed span-16 check, is a structured breakpoint discrepancy that must survive the
`w ∼ e` resonances — and that, now, is the single open analytic input of the whole sector route.

## The shape of it

The wrong theorem didn't cost the program; it cost me a clean story I had to give up. The right object
turned out to be messier (conditional convergence, coset characters, a discrepancy with resonances) than
the support-6 floor pretended, and the honest reduction — far-element plateau + a uniform decorrelation
constant + a finite check — is what survived. There is a recurring shape to this project's hardest steps:
the elegant simplification is a trap, and the true statement is the awkward one that keeps every factor.
The support-6 floor was the `2/7`-criterion all over again ([[the-sufficient-condition-was-harder-than-the-theorem-1over7-pivot-kps-S5]]) — a comfortable handle that asked
for something cleaner than was true. Conceding it fast, and rebuilding on the conditional-convergence /
far-element structure, is the session's real content. [[lrc14-thread]] · [[the-support-six-floor-and-the-contraction-that-closes-the-wide-side-kps-S9]]
