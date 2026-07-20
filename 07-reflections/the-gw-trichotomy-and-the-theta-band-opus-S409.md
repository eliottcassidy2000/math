# The GW trichotomy, its gate mechanism, and the θ-band where the twins separate

**Instance:** opus-2026-07-19-S409 (owner: work the next leads in conjunction). **HYP-8015.**
Script + frozen out: `lrc_conjunctions_gw_trichotomy_theta_opus_S409.py`. Two conjunctions,
both to verdict.

## A. The GW trichotomy law (λ-hole × (D,s)-rungs × gate, composed) — 16/16 exact

For GW_N = {1,…,N−2, N, 2(N−1)}, verified exactly on N = 5..20:

> **M(GW_N) = 1/(N+1)** (tight, slack 0) if **N ≡ 1 (mod 6)**;
> **M(GW_N) = 2/(2N+1)** via the escape pair **(3, 2(N−1))** if N odd, N ≢ 1 (mod 3);
> **M(GW_N) = 1/N** via the escape pair **(2, 2(N−1))** (s = 2N) if N even.

All three branches sit on slack-≤1 rungs, and all escapes run through the doubled ghost
2(N−1). **The gate mechanism, proved (two lines):** the D=2 escape needs a with
3a ≡ ±2 (mod 2N+1), which is solvable iff 3 ∤ 2N+1 iff **N ≢ 1 (mod 3)** — at gated N
the escape is arithmetically unavailable, not merely defeated. **The hole immunity,
verified 8/8 and provable:** the only residues that could kill an available escape are
±(N−1) mod (2N+1) — exactly the DELETED element and its mirror, i.e. the λ box-move hole
(S408 W2) — so whenever the escape exists it succeeds. The June gate (k ≡ 1 mod 6,
S405-F2), the box move, and the killer-absence are one mechanism: **the family frees the
±(N−1) classes by deleting its penultimate convergent; that same hole makes every
available escape unkillable; tightness survives only where mod-3 arithmetic denies the
escape a phase to live at.** Named residue for a full theorem: the upper halves (nothing
beats the stated M on each branch — ghost-packing-shaped, verified here) and the even
branch's parity story; one session, Lean-friendly throughout.

## B. The θ-band: the tight twins separate in the middle of the deformation

Exact M_θ at θ = p/8 (integer power-comparison, no floats): AP13 and GW are **value-tied
at θ ∈ {0, 1/8, 1/4} and again at {7/8, 1}, but SEPARATE throughout the middle band
θ ∈ {3/8, 1/2, 5/8, 3/4}**. Combined with S407 (c3 separates them at θ = 0: 112 vs 108)
and S408 (tied at both endpoints): the pair is endpoint-degenerate but interior-resolved —
the deformation sees structure the endpoint spectra cannot. The AP's maximizer path:
(1,14) → (3,14) → (23,107) → (14,65) → (11,51) → (23,78) → (31,112) → **(21,55) →
(13,34)** — it exits the denominator-14 Farey cell between θ = 1/8 and 1/4 and enters
Fibonacci-land at θ = 7/8 (21/55 and 13/34 are consecutive-Fibonacci fractions). Two
named phase marks: the **Farey exit ≈ 1/4** and the **golden entry ≈ 7/8**; the middle
band is neither world — and it is exactly where the twins separate. Follow-ups: exact
transition thetas (binary-search the jumps); whether interior-θ separation of endpoint-
tied pairs is GENERAL (a resolution tool for tight-locus questions: perturb θ, separate,
classify); the Markov-import question (S408) now has a concrete staging ground.

## Cross-links

S405-F2 (the shared mod-6 gate — mechanism now found) · S408 W1/W2 (deformation;
box-move hole) · THM-1291/1292 + ghost loop (the deleted penultimate convergent) ·
death-star gate tower (the trichotomy is the GW-shape row of the gate table, now with
mechanism) · mac-mini THM-1065 (the 6|(n−1) criterion, now explained) · script + out.
