# Analogy and reduction atlas: typed analogies, the Legendre tournament, and problem classification by logical form and machine size

2026-10-25, klein (worktree `codex/session-tournament-twist-unify`, third directive). Exploration
packet; statuses per note. **EXACT / PROVED** statements are elementary and proved inline;
**FINITE-EXACT** receipts agree normally and under `-O`; **CITED** external data are dated and
sourced; **OPEN** items are HYP-9375 and the ledger's open rows. No audit was requested or performed.

| Note | Content | Status |
|---|---|---|
| [ATLAS.md](ATLAS.md) | Cauchy–Schwarz ⟺ triangle inequality exactly when the norm is an inner-product norm; central limit theorem's compactness step is Tychonoff for `[0,1]^N`, the Gaussian limit requires characteristic-function identification; quadratic reciprocity makes the Legendre relation a tournament on primes `≡ 3 (mod 4)` (the numeral four is a coincidence); Möbius band = associated real sign bundle of square-root monodromy, spiral = log monodromy (argument principle), and the hostile `C_p` | EXACT / ANALOGY / COINCIDENCE, typed |
| [QR_TOURNAMENT.md](QR_TOURNAMENT.md) | Census of the Legendre tournament: all four 4-vertex classes, frequencies within `0.002` of random at 2,000 and 20,000 primes, not regular | FINITE-EXACT; HYP-9375 OPEN |
| [PROBLEM_CLASSIFICATION.md](PROBLEM_CLASSIFICATION.md) | The NP-completeness analogy repaired: machine encodings give level (`Pi^0_1`, `Pi^0_2`, decidable) and size (BB states), not a completeness class; external ledger (Goldbach BB(25)/BB(21), RH BB(744)…BB(95), Sloane BB(6,4), Erdős ternary BB(15), weak Collatz BB(14), odd multiply perfect BB(23,3), Antihydra BB(6)); repository ledger (LRC(14) decidable in principle; JC(2), `H >= disc`, black-core `Pi^0_1`; Collatz, Mahler, Rule 30 `Pi^0_2`); where tournaments actually enter | CITED + elementary |
| [LEDGER_15_25_MULTIPERFECT.md](LEDGER_15_25_MULTIPERFECT.md) | The fifteens (supersingular primes, Sloane's exponent, Erdős's 15 states), the twenty-fives (`II_{25,1}`, Goldbach's 25 states, the two regular groups of order 25, …) and multiply perfect numbers, each typed numerology-versus-structure | CITED / FINITE-EXACT |

Receipts: `qr-normal.jsonl`, `multiperfect-normal.json` and their `-optimized` twins; `SHA256SUMS`.

## Inheritance and board

Anchor: the classification ledger. Niche: the Legendre tournament. Wildcard: the typed analogies.
Closest proved mechanisms: THM-640 (Paley bridge, prime scope), THM-307 (QR decomposition of the
count), THM-2228/THM-4754 (Collatz-type shifts and tournaments), MISTAKE-214/232 (node repetition,
composite Paley). Canonical hostiles: `l^1` for Cauchy–Schwarz, stable laws for the central limit
theorem, `C_p` for the winding proof, the residue table versus orbit count for "four". Least-used
sidecars: the parallelogram law, the fixed-point equation of renormalization, the monodromy group of
a cover, the arithmetical level of a conjecture.

## Reproduction

```
python3 qr_tournament.py 120 2000      # exact census at 120 primes, sampled at 2,000 (stdlib, ~2 min)
python3 multiperfect.py                 # k-perfect numbers to 2,000,000 (stdlib)
```

Correction and extension: [orientation, monodromy and four-state transport](../orientation_monodromy_rule30_20261026/README.md). The limit-identification and cover/bundle types were repaired on 2026-10-26.
