# Class number is arithmetic entropy: hidden binary forms, and why 7 is rigid

> **LRC attachment corrected 2026-07-21 by MISTAKE-225.** The class-group
> examples below concern actual binary quadratic forms and the entropy wording
> is a useful analogy. The THM-2053 determinant residual is instead the
> polygonal support norm `h_K(Rd)` (THM-2055), so no discriminant-`-7` form,
> class group, or local-to-global conclusion currently classifies it. The
> claims that the LRC residual is fully pinned and has no class-group slack are
> therefore withdrawn; the classical arithmetic observations are retained.

*boxeph-2026-07-21-S217. Owner: look for even more hidden binary forms; think information theory. Builds on
S215 (each prime = its Paley object, disc −p), S216 (the anisotropic determinant gate = codex THM-2053's
residual, Heegner rigidity), kps-S17 (apex-prime / Heegner / p mod 4). Verified in
`04-computation/hidden_binary_forms_class_number_is_arithmetic_entropy_boxeph_S217.py`.*

## The hidden forms: the disc-`−p` gate lives in a whole class group

S216 found the anisotropic gate as the binary form of discriminant `−p` (the Paley spectral factor
`x²+x+(p+1)/4`). But a discriminant does not carry *one* form — it carries a whole **class group** of forms
under **Gauss composition**, and the *non-principal* classes are the **hidden binary forms**. Verified
enumeration:

| disc | `h` | reduced forms (the class group) | hidden forms |
|---|---|---|---|
| `−3,−4,−7,−8,−11` | 1 | principal only, e.g. `−7`: `(1,1,2)` | none |
| `−15` | 2 | `(1,1,4),(2,1,2)` | 1 |
| `−23` | 3 | `(1,1,6),(2,±1,3)` | 2 |
| `−47` | 5 | `(1,1,12),(2,±1,6),(3,±1,4)` | 4 |

The **Paley/Heegner primes `3,7,11`** (disc `−3,−7,−11`, plus `−4,−8` for the prime 2) have **only the
principal form** — no hidden forms.

## The information-theoretic reading: class number = arithmetic entropy

Here is the new lens. Ask: *given a prime `p`, which form of discriminant `D` represents it?* The **local**
data is the single Legendre bit `(D/p)` — whether `p` splits in `ℚ(√D)`. Verified:

- **`h = 1` (disc `−7`):** `p` is represented by the unique form `(1,1,2)` **iff `(−7/p)=1`** — one Legendre
  bit *decides everything*. (Sample: `2,11,23` represented; `3,5,13,17,19` not — exactly the split of the
  Legendre symbol.)
- **`h = 3` (disc `−23`):** `(−23/p)=1` still says `p` is represented by *some* form — but the primes then
  **split among the 3 classes** (verified: principal `(1,1,6)` gets `59,101,…`; the non-principal `(2,±1,3)`
  get `13,29,31,41,…`), and *which* class a prime lands in is **invisible to any Legendre/local test** — it
  is the Artin symbol in the **Hilbert class field**, `Cl(D) = Gal(H/K)`.

So:

> **The class number `h(D)` is the ARITHMETIC ENTROPY of the form — the `log₂ h(D)` bits needed *beyond*
> local (Legendre) data to decide representation.** `h=1` ⟹ **zero** entropy: local determines global. `h>1`
> ⟹ `log₂ h` hidden bits (the class-group information: `−23 → 1.58`, `−47 → 2.32` bits).

This is exactly S216's "rigidity," now quantified: a class-number-1 gate has **no hidden information**.

## Why 7 is rigid — and hard

`LRC(14) = 2·7 → disc −7 → h(−7)=1 → zero arithmetic entropy.` codex's anisotropic determinant gate
(THM-2053) on the rank-two residual plane is therefore **fully pinned by local Legendre data** — the S215
Paley/Gauss-sum information — with **no class-group slack**. Two consequences, both explaining kps-S17's "7
is the first hard case":

- **A counterexample has nowhere to hide.** There are no hidden non-principal forms, no extra bits of global
  arithmetic freedom to slip a resonance into. The residual is a single rigid class (S216).
- **But the certificate must be exact.** With zero slack, there is no "generic/entropic" margin to exploit;
  the deciding certificate has to be the sharp local one (the Euler/`χ` / Borsuk–Ulam certificate for the
  anisotropic `p≡3 mod 4` case, S212/kps-S17). Rigidity is *why* the problem is both irreducible and,
  eventually, tractable — the information is all local.

So the difficulty of 7 is not complexity but **rigidity**: a class-number-1 anisotropic gate is a
zero-entropy object. (The Heegner discriminants `−3,−7,−11,−19,−43,−67,−163` are the *only* imaginary
quadratics with `h=1`; `−7` is the LRC(14) one.)

## Bonus hidden entropy: the reify-ladder poles are entropy extremes

A second information-theoretic invariant separates the two poles of the reify ladder — the **entropy of the
tournament's score distribution** (verified): the **transitive** tournament has scores `0,1,…,n−1`, all
distinct, entropy `log₂ n` (**maximum spread**); the **Paley/regular** tournament has all scores `(n−1)/2`,
a delta, entropy `0` (**minimum**). So the AP/transitive nullcone vertex (the rank-11 gate, S214) is the
*maximum-score-entropy* configuration, while the Paley symmetric pole (disc `−p`, `i√p`, S215) is the
*zero-score-entropy* one. The two ends of every recent thread — the transitive AP gate vs. the Paley
anisotropic form — are the two entropy extremes.

## The one-line reading

> **The disc-`−p` anisotropic gate (S216) sits in a class group of binary forms; the class number is its
> arithmetic entropy — the bits beyond local Legendre needed to decide representation. Heegner (`h=1`:
> `−3,−7,−11` = Paley `3,7,11`) has zero arithmetic entropy: local (S215) determines global, the gate is
> rigid. `LRC(14)=2·7 → −7 → h=1 → no hidden bits`, so a counterexample has nowhere to hide and the
> certificate must be the exact local (Euler/`χ`) one — which is precisely why 7 is the first hard-but-
> tractable case.**

Honest scope: the class-group / representation facts and Heegner `h=1` list are classical and verified here;
the score-entropy poles are elementary. The contribution is the **information-theoretic framing** — class
number as arithmetic entropy, the non-principal forms as hidden binary forms, and the reading of S216's
rigidity as *zero hidden entropy* for LRC(14)'s `−7`. It is a conceptual sharpening of why the apex prime is
rigid, not a new proof step.

Links: HYP-8870, THM-2053 (codex), HYP-8865, HYP-8860,
[[the-anisotropic-determinant-gate-heegner-rigidity-of-codexs-rank-two-residual-boxeph-S216]],
[[each-prime-is-its-paley-tournament-a-periodic-table-of-2-3-5-7-11-for-lrc14-boxeph-S215]].
