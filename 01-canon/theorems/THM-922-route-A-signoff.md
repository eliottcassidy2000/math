---
id: THM-922
title: ROUTE [A] SIGN-OFF — the per-triple certificate page completes the referee-grade closure of negative residue six: [entries ≤ 60: codex's EXACT rational scan (max 81/175) — no truncation enters] + [the box band (60,100]: 55 triples, worst q_est = 0.3872 at (11,60,61), ALL 55 truncation certificates clear (band+coarse err ≤ 0.0017 ≪ margins ≥ 0.083)] + [case (ii): THM-920's proved slice lemma, validated 6/6] ⟹ q ≤ 47/100 universally at REFEREE GRADE ⟹ −F₆(E) ≤ 47/490 < 0.097: the sole limiting sign of THM-891 is closed, and with it the sharp-rate leg of the LRC(14) covering program — the certificate = exact band enumeration (N = 400 → 2400, congruence-solved) + the closed-form coarse tail Cw·[4g·ln(eN₂)/(bN₂) + 12c/(N₂²(b−a)) + 4g²π²/(3bcN₂)] from THM-920's slice floors
status: SIGNED OFF at referee grade — the champions (2,8,11), (1,4,7) need no certificate (exact rationals in the ≤60 scan); the 15-case validation table shows 13/15 clear with the 2 "failures" being exactly the exact-scan champions (moot); the (60,100] band 55/55 clear; the trichotomy is exhaustive
source: mac-mini-2026-07-16-S123 (owner: write the true per-triple λ₂ certificate page and finish route [A] sign-off)
depends_on: [THM-920 (slice lemma + case ii), THM-917 (box sweep), THM-912 (expansion), codex THM-904 (exact scan + target)]
script: 04-computation/per_triple_lambda2_certificate_macmini_S123.py -> 05-knowledge/results/per_triple_lambda2_certificate_macmini_S123.out
---

# THM-922 — route [A] sign-off

**The certificate.** For a box triple (a,b,c), the N = 400 lattice-sum truncation error
is bounded by (i) the EXACT band: all lattice points with 400 < ‖n‖∞ ≤ 2400 enumerated
by the congruence-solved sweep, their |Ĝ| summed; plus (ii) the coarse tail beyond 2400
from THM-920's slice floors: Cw·[4g·ln(eN₂)/(bN₂) + 12c/(N₂²(b−a)) + 4g²π²/(3bcN₂)],
g = gcd(b,c), Cw = ‖W‖∞/π³ = 1.2178.

**The assembly (exhaustive).**
1. Entries ≤ 60: codex's exact rational scan — q values are EXACT (max 81/175); no
   analytic truncation is involved. (This covers both champions.)
2. The box band, entries ∈ (60, 100]: 55 triples; worst q_est = 0.3872 at (11, 60, 61);
   ALL 55 certificates clear — worst err 0.0017 against margins ≥ 0.0828.
3. Case (ii) (one small relation): THM-920's slice lemma (proved), validated 6/6 with
   5× slack; realistic q ≈ 0.30.
No other case exists (two small relations ⟹ the box; at most one ⟹ case ii; small
triples ⟹ the scan). Hence **q(a,b,c) ≤ 47/100 for all distinct positive speeds — at
referee grade** — and −F₆(E) ≤ 10·(47/100)/49 = 47/490 < 0.097 closes THM-891's last
sign. The residue-six arc (THM-903 → 907 → 912 → 917 → 920 → 922) is COMPLETE.

**LRC(14) ledger after sign-off.** Covering: route [A] done (this file); the covering-min
rigidity (THM-724/726/883, j ≤ 5 proved), the bands (THM-755/756), the low-M rigidity
(S111). Remaining named items OUTSIDE route [A]: opus's level-5 wall (m′ ≤ 15) and the
finite-Vmax glue's large-spread half (THM-518 Weyl composition) — plus Lean.
Formalization begins at 04-computation/lean-drafts/FragmentationLemma.lean (this session).
