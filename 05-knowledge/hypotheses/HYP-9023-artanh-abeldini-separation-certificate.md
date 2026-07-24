---
id: HYP-9023
title: "An external artanh separation certificate is an Abel-Dini (THM-2000 sec 3.1) log-bound; it is the missing certified-inequality layer for the support-harmonic mass chain"
status: >
  PARTIAL DECODE + VERIFIED CONSTRUCTION. CONFIRMED (exact): the owner's snippet is a
  rational separation certificate proving an Abel-Dini log-combination > 1/25; its two
  parameters t_A=389/2181, t_B=5872957/11821757 are EXACTLY telescoping ratios
  t_n=x_n/(S_n+S_{n-1}) for partial-sum pairs (896,1285) and (2974400,8847357), the
  verbatim construction of THM-2000 sec 3.1; and den(certificate) = lcm(den of the two
  bounds) = 2^8 3^4 5^2 7 31^5 257 727^3 381347^5, built only from those bounds. COEFFICIENT
  RESOLVED (opus-S4, was "not recoverable"): a height cliff (only d=+1 gives c-height 10^3.8
  vs ~10^52 for every other integer d) UNIQUELY pins (c,d,r)=(2457/6592, 1, 0), so
  RHS(27) = (2457/6592) log(8847357/2974400) - log(1285/896), no rational part; the earlier
  "31^5 381347^5 in den(r)" was an artifact of the t^5/5 truncation order, not eq (27).
  VERIFIED CONSTRUCTION: the same 3-term artanh sandwich certifies a genuine THM-2000
  mass ordering M(6,2)>M(4,3) float-free (log2 >= 842/1215 > 9/13 => 26 log2 - 18 >=
  22/1215 > 0). HOME (opus-S4/klein-S404/mac-mini-S169): the log-energy-beats-floor family,
  LRC(14) wider-gap n=13 leading (2457=3*S2(AP{1..13}), 6592=2^6*S1(GW{1..11,13,24}),
  1/25~1/(2n-1)); RULED OUT this session -- figurate-mass ordering (no M(s,d) difference
  matches; both logs below the mass band), BOTH candidate LRC papers (Bedert 2511.16636 and
  2604.23906, by direct fetch), and standard binomial Pade (integers are non-hypergeometric).
  Irrationality-MEASURE not fully killed. See 07-reflections/artanh-two-log-form-pinned-and-homed-opus-S4.md.
source: opus-2026-07-23-S2..S4 (decoding an owner-supplied external snippet)
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

## Partially recovered by p-adic pinning (opus-S3, `snippet_padic_decode_opus_S3.py`)

The denominators have **isolated primes**: `{7,257,727}` live only on the `U_A` (L_A)
side, `{31,381347}` only on the `Lo_B` (L_B) side. Clearing to `L=lcm` and reducing
`NX = c·M_B − d·M_A + rL` modulo each isolated prime pins the coefficients IF the
rational part `r` is prime to that prime:

- **`d = 1` EXACTLY and robustly** (coefficient of `L_A = log(1285/896)`). Verified:
  `den(r)` contains no `7,257,727` for any tested `c`, so the `U_A` side contributes
  with coefficient exactly 1 and no residue.
- **`c` (coefficient of `L_B`) is NOT pinnable**: `den(r)` always contains
  `31⁵·381347⁵ = (S_B+S_{B−1})⁵`. So **the rational part `r` inherits the Abel–Dini
  denominator `11821757 = S_B + S_{B−1}`** — the signature of an Abel–Dini *gap* term
  (`x_n/S_n` or `x_n/(S_n+S_{n−1})`). This corrupts the p-adic residue for `c`.

So the snippet fixes `d=1` and the Abel–Dini structure but leaves `(c, r)`
under-determined (one equation, two unknowns). `RHS(27) = c·L_B − L_A + r`, with `r`
an Abel–Dini rational carrying `11821757` in its denominator. Numerically `L_B ~ 3·L_A`
(ratio 3.023), `L_A ~ 9/25`, `L_B ~ 27/25`, `RHS(27) ~ 0.0447` sits just above `1/25`.
Finishing the decode needs eq (27)'s surrounding text (the source is external/lost).

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
