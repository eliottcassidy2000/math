---
id: HYP-9023
title: "An external artanh separation certificate is an Abel-Dini (THM-2000 sec 3.1) log-bound; it is the missing certified-inequality layer for the support-harmonic mass chain"
status: >
  PARTIAL DECODE + VERIFIED CONSTRUCTION. CONFIRMED (exact): the owner's snippet is a
  rational separation certificate proving an Abel-Dini log-combination > 1/25; its two
  parameters t_A=389/2181, t_B=5872957/11821757 are EXACTLY telescoping ratios
  t_n=x_n/(S_n+S_{n-1}) for partial-sum pairs (896,1285) and (2974400,8847357), the
  verbatim construction of THM-2000 sec 3.1; and den(certificate) = lcm(den of the two
  bounds) = 2^8 3^4 5^2 7 31^5 257 727^3 381347^5, built only from those bounds. NOT
  RECOVERABLE from the snippet alone: the exact integer coefficients of RHS(27) (no
  small-coefficient exact fit; the rational part is large, from the source's eq (27)).
  VERIFIED CONSTRUCTION: the same 3-term artanh sandwich certifies a genuine THM-2000
  mass ordering M(6,2)>M(4,3) float-free (log2 >= 842/1215 > 9/13 => 26 log2 - 18 >=
  22/1215 > 0). SPECULATIVE: the same separation-certificate shape is the LRC(14)
  inf L(S)>0 prize (structural analogy only; L is a sinc sum, not a log).
source: opus-2026-07-23-S2 (decoding an owner-supplied external snippet)
depends_on: [THM-2000]
related:
  - THM-501   # LRC singular series L(S) (same "certify > rational floor" shape)
  - THM-523   # LRC(14) singular-series skeleton, floors 1/1260, 2/29
  - HYP-9022  # prior session: metagraph isoperimetric dimension (also a "seed -> repo" bridge)
script: 04-computation/certified_logratio_abeldini_opus_S2.py
output: 05-knowledge/results/certified_logratio_abeldini_opus_S2.out
lean_target: 04-computation/lean/TournamentH7/TournamentH7/SupportHarmonicFigurate.lean
reflection: 07-reflections/the-artanh-abeldini-separation-certificate-opus-S2.md
---

# HYP-9023 — the artanh/Abel-Dini separation certificate

## The snippet (owner-supplied external fragment)

`log((1+t)/(1-t)) = 2·artanh(t)` sandwiched between rationals (lower: truncate after
`t^5/5`; upper: geometric tail `t^5/(5(1-t^2))`), evaluated at `t_A=389/2181`
[ratio `1285/896`] and `t_B=5872957/11821757` [ratio `8847357/2974400`], certifying
`RHS(27) - 1/25 >= G` for a 50-digit positive rational `G ~ 0.004737`.

## Confirmed

1. **Abel-Dini telescoping.** `t = x_n/(S_n+S_{n-1})` with `(S_{n-1},S_n)=(896,1285)`,
   `(2974400,8847357)`; `(1+t)/(1-t)=S_n/S_{n-1}`, `2·artanh(t)=log(S_n/S_{n-1})`.
   Verbatim THM-2000 sec 3.1.
2. **Fingerprint.** `den(G) = lcm(den Lo_B, den U_A)` exactly (no foreign primes).
3. **Form.** `RHS(27) = +c·L_B - d·L_A + rational` (L_B lower-bounded, L_A upper-bounded).

## Not recoverable without eq (27)

The exact `(c,d,rational)` are under-determined by the snippet: the minimal exact
`c·Lo_B - d·U_A = G+1/25` solution has ~50-digit coefficients, so the true statement
carries a large rational from prior algebra. Numerically `L_B ~ 3·L_A` (ratio 3.023),
`L_A ~ 9/25`, `L_B ~ 27/25`, and `RHS(27) ~ 0.0447` sits just above `1/25`.

## Verified construction (the payoff)

The same engine certifies real THM-2000 transcendental mass orderings with pure
rationals. `M(6,2)=2 log2 > M(4,3)=18-24 log2  <=>  log2 > 9/13`; the 3-term artanh
lower bound `log2 >= 842/1215 > 9/13` proves it, `M(6,2)-M(4,3) >= 22/1215 > 0`,
float-free. THM-2000 currently checks such orderings only at 45-digit float precision;
its Lean kernel leaves them to "the paper proof". This engine is the missing certified
sibling of the file's discrete `reciprocal_block_bounds`.

## Why it may be a key hint

The repo's analytic frontier (support-harmonic masses; LRC singular series `inf L>0`)
is a stack of "real number vs rational threshold" certificates currently discharged
numerically. This artanh/Abel-Dini sandwich is the first clean, Lean-portable member
of the certified-inequality toolkit those proofs need. Next: formalize
`logRatioSandwich` in Lean; certify `M(6,2)>M(4,3)`; audit the mass chain for orderings
provable by few artanh terms.
