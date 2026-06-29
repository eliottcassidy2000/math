---
id: HYP-3547
title: 28 = T(7) = C(8,2) = dim so(8) = the 2nd even perfect number is the arc-count of LRC(14)'s Forcade order 8; apex 7 carries OCTONION structure (QR{1,2,4} = Fano line = octonion rule = B_0(Paley T_7) = the 7 sectors); the 7-14-21-28 tower = Im(O)/G2/so(7)/so(8) with order 14 = dim G2; and apex 7 is one of only TWO primes (3,7) simultaneously Mersenne + Heegner + 3-mod-4 -- so its three arithmetic properties ARE the three proof pillars (2-adic descent, Fejer-Bochner SOS, Borsuk-Ulam), singling out LRC(6),(14)
status: (1) THEOREM (perfect = T(Mersenne prime), elementary); (2) SOLID structural (QR{1,2,4}=Fano=octonion=B_0(T_7), verified); (4) VERIFIED number-theoretic characterization (Mersenne ∩ Heegner ∩ 3mod4 = {3,7}); (3) dimension tower + cut⊕cycle (real) + G2/so resonance (numerology flagged). The three-pillars MAPPING is an organizing insight, not a proof.
source: mac-mini-2026-06-29-S14
related:
  - HYP-3546   # perfect numbers = Forcade arc-hypercube dims (28 = C(8,2))
  - THM-580    # 2-adic parity descent (the MERSENNE pillar: 14=2*7 peels to apex-7 face)
  - THM-581    # Borsuk-Ulam witness / saddle index (the 3-mod-4 pillar: (p-1)/2=3 odd)
  - THM-586    # Paley T_7 arithmetic (the apex DRT carrying QR/Fano/octonion)
  - THM-448    # Mersenne doubling B_0(T_7)=T_3
  - HYP-3535   # dihedral synthesis (Q(sqrt(-7)) Heegner, S75e Fejer-Bochner SOS = the HEEGNER pillar)
external: Euclid-Euler; Heegner numbers; octonions/Fano plane/G2; Forcade 1973
results:
  - 04-computation/twentyeight_octonion_heegner_macmini_20260629.py
  - 05-knowledge/results/twentyeight_octonion_heegner_macmini_20260629.out
---

# HYP-3547 -- 28, the octonion apex, and why LRC(14) is uniquely poised

## (1) THEOREM: even perfect numbers are the triangular numbers of Mersenne primes
`T(2^p-1) = (2^p-1)2^p/2 = 2^{p-1}(2^p-1)` = the Euclid-Euler even perfect number when `2^p-1` is
prime. So `6 = T(3)`, `28 = T(7)`, `496 = T(31)`, `8128 = T(127)`. Since `T(M) = C(M+1,2)`, each
perfect number is the **arc-count of the tournament on `M+1 = 2^p` vertices** (the Forcade order,
HYP-3546). For the LRC apex prime `7`: `28 = T(7) = C(8,2)` -- the 2nd even perfect number is the
arc-hypercube dimension of LRC(14)'s Forcade order 8. The Mersenne-prime apices `3,7,31,127` give LRC
orders `6,14,62,254`.

## (2) The apex 7 carries OCTONION structure (the solid anchor)
The quadratic residues mod 7 are `QR = {1,2,4}` -- simultaneously: the Paley `T_7` arc rule; the
out-neighborhood `B_0(T_7)` (the Mersenne doubling, THM-448, verified); a **Fano-plane line** (the 7
lines are the translates of `{1,2,4}` mod 7, = PG(2,2)); and the **octonion multiplication rule** (the
7 associative triples of `Im(O) = <e_1..e_7>` are exactly the Fano lines). So **LRC(14)'s 7 inner
sectors are the 7 Fano points = the 7 imaginary octonion units**, and the apex tournament's arc rule
IS the octonion product. This is structural, not numerical: the difference set `{1,2,4}` genuinely
equals the Fano/octonion incidence.

## (3) The 7-14-21-28 dimension tower (part real, part resonance)
`7, 14, 21, 28 = 7*{1,2,3,4} = dim Im(O), dim G2, dim so(7), dim so(8)`. Two of the steps are a REAL
tournament decomposition:
- `28 = 21 + 7`: the arc space of an 8-tournament `= cycle space (21 tiles = C(7,2)) (+) cut space
  (7 scores = base-path arcs)` -- the GF(2) cut⊕cycle split -- mirrors `so(8) = so(7) (+) R^7`.
- `21 = 14 + 7`: `so(7) = G2 (+) (7-coset)`.
The order `14 = dim G2 = Aut(O)` and `28 = dim so(8)` (triality) are resonances (flagged: not shown to
be group actions, only dimensions). But the cut⊕cycle = `so(7)⊕R^7` correspondence is genuine and ties
the half-tiling model to the reductive split.

## (4) Apex 7 = Mersenne AND Heegner AND 3-mod-4: the three pillars
The Heegner primes (`Q(sqrt(-p))` class number 1, p odd) are `{3,7,11,19,43,67,163}` (all `= 3 mod 4`);
the Mersenne primes are `{3,7,31,127,...}` (all `= 3 mod 4`). Their intersection is **exactly `{3,7}`**,
giving LRC orders `{6,14}`. For apex `p in {3,7}` the three arithmetic properties are precisely the
project's three proof pillars:

| property of apex `p` | what it provides | proof pillar |
|---|---|---|
| `p = 3 mod 4` (Paley exists) | saddle index `(p-1)/2` is ODD (p=7 -> 3) -> free `Z_2`, Borsuk-Ulam | **THM-581** (witness, odd degree) |
| `p` Heegner (`Q(sqrt(-p))` h=1) | the GENTLEST cyclotomic; totally-real `Q(cos 2pi/p)` Fejer-Bochner SOS minorant | **S75e SOS** (HYP-3535, the floor) |
| `p` Mersenne (`2p = 2^?...`, perfect arcs) | `2p = 2*p` peels one free `Z_2` to the all-odd apex-`p` face; perfect arc-count | **THM-580** (2-adic parity descent) |

So **LRC(14) is, after the PROVED LRC(6), the unique small case whose apex prime supplies all three
proof tools at once** -- and the tools are bespoke to `14` precisely because `7` is the only prime
besides `3` that is Mersenne, Heegner, and `3 mod 4` together. This reframes "why 14 is special" from
numerology to mechanism: the three pillars exist BECAUSE of the apex's three arithmetic properties, and
`28 = T(7)` is their shared signature (perfect = Mersenne descent; `Q(sqrt(-7))` = Heegner SOS; `(7-1)/2=3`
odd = Borsuk-Ulam). The next Mersenne apices `31,127` keep two pillars (Mersenne + 3-mod-4) but LOSE
Heegner -- so the bespoke `14`-strategy does not transfer, consistent with LRC's known hardness above 7.

## Extra 14/28 facts (cataloged)
`phi(14) = 6` (= the 6 inner sectors = the 6 witnesses = units mod 14 `{1,3,5,9,11,13}`, sum `42 = 3*14`).
`28 = 2+3+5+7+11` (sum of the first 5 primes). `28` proper divisors `{1,2,4,7,14}` sum to `28` and are
exactly the apex order-tower (powers of 2 and their 7-multiples). `14 = 2*7`, `28 = 4*7 = 2*14`.

## What it buys
Not a proof step, but an organizing principle: the floor/cap closure SHOULD use all three pillars (they
are simultaneously available only at `p=7`), and the `Q(sqrt(-7))`-Heegner SOS (S75e) is the linchpin the
Mersenne/Borsuk-Ulam structure alone cannot replace. Engineering-side: the apex-7 octonion/Fano structure
gives a canonical 7-symbol code (the 7 sectors = Fano points) for the LRC danger-comb alphabet.
