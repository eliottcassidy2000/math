# The polysemous constants: which numbers wear many hats, and which conflations break proofs

*klein-2026-06-29-S7. Asked, after the "two 7's" (apex prime vs b_1^-(5)=7, HYP-3563), to catalog the OTHER numbers this project sees in multiple guises that proofs must keep apart. The general lesson, then the table.*

## The lesson from 7, generalized

A number can occur in an **arithmetic** role (a prime, totient, modulus, factor) and in a
**dimensional/combinatorial** role (a count, a Betti number, a binomial, a dimension). When the same
number appears in both, there are exactly three cases, and telling them apart is the whole game:

- **BRIDGE** — the two roles are the *same* number for a *provable* reason, valid at all scales. SAFE and
  often powerful: use it.
- **TRAP** — the two roles coincide at one value only, by accident. Conflating them invents false
  structure. The diagnostic is **persistence**: a bridge survives `n -> n+1`; a trap holds at one `n` and
  dies. (`b_1^-(5)=7` died at `n=7`: `1772 = 2^2·443`. Trap.)
- **HOMONYM** — one symbol, several genuinely different objects (an overloaded *index*, or several
  order-2 maps). Not a coincidence — a notation hazard. Fix: track all of them explicitly (MISTAKE-033).

Run the persistence test before believing any numeric coincidence is structure.

## The table

| number | the roles it plays | verdict |
|---|---|---|
| **2** | (a) complement / reversal involution `R = σ` (`t↦1−t`); (b) the 2-adic **doubling** descent `δ` (`t↦2t`, THM-580); (c) the field GF(2) / cut⊕cycle; (d) the base of `2^{n-1}=E[H]`; (e) **three** distinct involutions — arc reversal `T^op`, path reversal, tile-bit complement (MISTAKE-033) | **HOMONYM (most dangerous).** `σ` and `δ` are the "two order-2 structures of the circle" and are *not* interchangeable (reflection wants an invariant; doubling computes). MISTAKE-033 was conflating the three involutions. Always name which 2. |
| **n** (the index) | (a) tournament vertex count (THM-584/587/588, every metagraph statement); (b) the LRC modulus `n=14`; (c) a Paley prime `p`; (d) a Forcade/Mersenne order `2^p`; (e) the apex prime of LRC(`2p`) | **HOMONYM.** "at n=5" (my b_1^- work, a tournament size) is a different `n` from "LRC(14)" (a modulus) and from apex `7`. The owner's own "H_1^- at n=5 / apex 7" mixed two indices. Always say *which* n. |
| **7** | (a) the **apex prime** of LRC(14): `Z_7`, Paley `T_7`, QR `{1,2,4}` = Fano line = octonion triple (HYP-3547/THM-586) — arithmetic, deep; (b) `b_1^-(5)=7`, `V(E_5)=7` — dimensional, coincidental | **BRIDGE (a) + TRAP (b).** The octonions live in `Z_7` (a), NOT in the n=5 metagraph homology (b). HYP-3563 killed the (b) route. Don't route octonion/descent work through `b_1^-`. |
| **14** | (a) the LRC modulus / runner count `=2·7`; (b) the **14-sheet count** `N_R` (THM-579, the floor 2nd moment); (c) `dim G_2` (HYP-3547) | **(a)+(b) BRIDGE** (the 14 sheets *are* `Z/14`, the modulus — same 14). **(c) NUMEROLOGY** (flagged in HYP-3547): `dim G_2 = 14` is suggestive but not load-bearing; do not build a proof on it. |
| **28** | (a) `2·14`; (b) the 2nd **perfect** number `= T(7) = C(8,2)` = arc-count of the Forcade order-8 tournament (HYP-3546, **proven** Euclid–Euler-through-tournaments); (c) `dim so(8)` (HYP-3547) | **(b) BRIDGE (a proven theorem!)** — perfect = `T(Mersenne)` = arc-hypercube dimension; safe and powerful. **(c) NUMEROLOGY** (flagged). The bridge is the arithmetic↔dimensional identity *with a proof*; the `so(8)` reading is decoration. |
| **6** | (a) `φ(14)=6` — the **inner sectors / units mod 14 / witnesses** (HYP-3538, the cap dimension; arithmetic); (b) `C(4,2)=6` — the n=4 **arc count** (dimensional); (c) `T(3)=6` — the **first perfect** number `=` LRC(6) | **Three unrelated 6's.** The cap's 6 (a) is a totient; the n=4 arc-cube's 6 (b) is a binomial; LRC(6)'s 6 (c) is `2·3`. No bridge among them — a TRAP if cross-identified (e.g., do not read the 6 cap sectors as n=4 arcs). |
| **3** | (a) the prime `3` — LRC(6) apex, `3≡3 mod 4` (Paley/Borsuk–Ulam, THM-581 saddle index `(p−1)/2`); (b) **3-cycles** = cyclicity = the metagraph Fiedler mode (THM-588); (c) the Fano **line size** 3 = octonion-triple length; (d) the "3 proof pillars" (HYP-3547) | **Mostly HOMONYM.** The prime 3 (a) and the 3-cycle (b) are different 3's (one arithmetic, one a cycle length); the Fano-line 3 (c) is the octonion-triple size, tied to apex 7 not to (a) or (b). |
| **red-herring class** (e.g. **89 = F_11**, the primes **17, 443**) | numbers that *appear* in computations (binding denominators `89`; factors of `b_1^-`: `119=7·17`, `1772=2^2·443`) and tempt a structural story (Fibonacci/badly-approximable) | **TRAP by default.** HYP-3551 explicitly flagged `89=F_11` as a red herring (tightness = core density, not bad approximability). Factorization primes like `17, 443` are *outputs*, not structure — never reverse-engineer meaning from a single factorization. |

## How to use this (proof hygiene)

1. **Name the role.** Write "the apex prime 7" or "the totient `φ(14)=6`" or "the 14 sheets", never the bare
   number, when two readings are in play.
2. **Persistence test before promotion.** Before treating a numeric coincidence as structure, check it at
   the next `n`. Bridges persist (28 = `T(Mp)` for every Mersenne `Mp`); traps die (`b_1^-=7` only at 5).
3. **Numerology stays flagged.** `dim G_2 = 14`, `dim so(8) = 28` are *mnemonics*, not load-bearing; a
   proof may not depend on them (HYP-3547 already flags this).
4. **Track every order-2 map.** Reflection `σ`, doubling `δ`, arc-reversal, path-reversal, tile-complement
   are distinct (MISTAKE-033). The safe bridges among them (`σ` = antipodal of the arc-cube, THM-584) are
   proven; the rest must be composed explicitly.

The deepest safe identities in the project are exactly the **proven** arithmetic↔dimensional bridges —
`28 = T(7) = C(8,2)` (HYP-3546), `P_n(-1) = SC(n)` (THM-587, the antipodal Euler number), `χ(T_p)=p`
(THM-129). Everything that merely *looks* like one of those should be run through the persistence test
before it touches a proof.

See MISTAKE-033, [[two-order-two-structures-parity-and-descent]] (the two 2's), HYP-3563 (the 7 trap),
HYP-3547 (28/14 bridges-and-numerology), HYP-3551 (the 89 red herring), THM-587 (the `P_n(-1)=SC` bridge).
