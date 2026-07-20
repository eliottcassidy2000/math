# The regularity bridge, the ghost loop, and the sheet tower

**Instance:** opus-2026-07-19-S407 (owner: see more connections between previous threads;
long session, own hypotheses, investigated concurrently). **HYP-8005.** Script + frozen
out: `lrc_connections_tournament_slack2_ghost_opus_S407.py`.

## A. The observer-tournament regularity bridge — CONFIRMED (r = −0.926)

Hypothesis (mine, from THM-640's Paley frame + S399's "the tournament lens returns at the
wall"): at a family's binding time t*, canon's half-turn comparator (THM-373/381) on the
14 positions (observer at 0) yields a tournament whose 3-cycle count c3 anti-correlates
with M — tightness should force near-regularity (max c3 = 112 at n = 14).

Data (13 families, exact): every near-floor family (M < 0.08) sits at c3 ∈ {108, 112};
**AP13, the deep well, ladder3, ladder5, {1..12,14}, {1..12,15} all attain exactly
112/112**; loose families sit at 42–92 (clusters 42–52, primes13 42, randoms 68–92).
Pearson corr(M, c3) = **−0.926**. Two structure notes: (i) the AP's 112 is *forced* —
its positions k/14 are perfectly uniform, so the comparator gives the rotational
(circulant) tournament, which is regular (one-line proof); (ii) **the invariant separates
the two tight families**: AP = 112, GW = 108 — the first observable I know that
distinguishes the extremal pair at their common value 1/14. Honest scope: a 13-family
correlation and one forced instance, not a theorem; the conjecture-form worth proving is
"M ≤ f(112 − c3(t*))" or the qualitative "tight ⟹ c3 ≥ 108". Triage-compliant: c3(t*)
is translation-SENSITIVE (positions at t* move under V ↦ V + c), so axis 1 does not bar
it — this is the tournament-side floor detector the S403 Singer probe failed to find.

## B. Slack-2 single-far: reconfirmed empty at N ≤ 14; the real open row is multi-far

Exhaustive N = 6..14, x ≤ 300 (gridmax prune + exact scan): **zero first-gap-window hits
at all** — independently reconfirming death-star's census (members only at N = 19, 25,
31, all slack-1) with a different instrument. Consequence for the gate table: the
"slack-2 gate" question is VACUOUS in single-far (settled by census + this); the live
form of my hypothesis moves to **multi-far families: does any family with ≥ 2 far
elements realize a slack-2 rung D/((N+1)D − 2) inside the first gap?** — unexplored,
now the named open row.

## C. The ghost loop closes: the excluded element is the penultimate convergent

Verified 3/3: δ(ghost) = 1 at the maximizer — DW: δ(13) = 1 at 14/183; ladder3:
δ(12) = 1 at 17/41; F₄(31): δ(30) = 1 at 55/127. This is F1's identity re-read:
(N−1)·a ≡ −1 (mod Q) ⟺ δ(N−1) = 1 ⟺ N−1 is the **penultimate convergent denominator**
of the maximizer (δ = 1 characterizes it). One sentence now holds the whole structure:

> **The extremal family EXCLUDES the penultimate convergent of its own maximizer (the
> −1-ray generator), pays the K-duty for that exclusion (THM-1292's ghost gap M/K), and
> its determinant D is the CF remainder of the convergent it DOES include (THM-1291).**

Extremality = forbidding your own best rational approximation. The June forced-cover
obstruction, the July ghost channel, the CF active-leg law, and the ghost-packing
ceiling are four limbs of this sentence.

## D. The sheet tower (F5's rhyme → proved note)

Two-line identities: for the shift τ ↦ τ + 1/2^m: ‖v(τ + 1/2^m)‖ = ‖vτ + v/2^m‖, which
EQUALS ‖vτ‖ whenever 2^m | v, and shifts the phase by v/2^m otherwise. So the dyadic
shifts fix the 2^m-divisible sub-family's distances and move everyone else — the
THM-760/761 sheet dodge at c = 2^m is the m-th layer; mac-mini-S123's τ ↔ 1−τ palindrome
is the m = 1 involution face (‖v(1−τ)‖ = ‖vτ‖ always — reflection — with the central
column at τ = 1/2 fixed and even speeds' local structure preserved); the
witness-component parity fact and the t = 1/2^m witness criterion (⟺ no speed ≡ 0 mod
2^m, m ≤ 3 at threshold 1/14) are its fixed-point faces. F5 is no longer a rhyme: the
"2-adic mirror stack" IS the dyadic sheet tower, stated. What remains genuinely open is
whether Wall B's palindrome layer admits a depth-2 refinement (τ ↦ τ + 1/4 acting on
codex's blocker words) — one named session.

## Cross-links

THM-640 + THM-373/381 (the comparator; A) · S403 Singer negative (what A replaces) ·
death-star census THM-1255/1256 (B) · HYP-7985-F1, THM-1291, THM-1292 (C — the loop) ·
THM-760/761, mac-mini-S123, kps-S259 (D) · script + frozen out (S407).
