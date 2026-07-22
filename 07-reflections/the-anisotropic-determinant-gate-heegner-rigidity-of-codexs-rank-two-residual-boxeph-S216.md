# The anisotropic determinant gate: Heegner rigidity of codex's rank-two residual

*boxeph-2026-07-21-S216. Owner: keep pushing for LRC leverage; think "anisotropic determinant gate."
Builds on incoming codex THM-2053 (the anisotropic rank-two terminal, PROVED) + THM-2052 (rank-11 descent);
kind-pasteur S17 (apex-prime / Heegner / p mod 4); my S212 (χ Euler branch / HYP-8845), S214 (rank-11
AP-core), S215 (each prime = its Paley object, disc −p). Verified in
`04-computation/anisotropic_determinant_gate_heegner_residual_boxeph_S216.py`.*

## The gate already exists — it is codex's anisotropic terminal

codex THM-2053 (PROVED) says every rank-two relation plane `U` (a `Z`-basis `u,z`, columns
`c_i=(u_i,z_i)`) has the **anisotropic terminal**

> `max_i |a z_i − b u_i| ≤ (a²+b²)/91  ⟹  the direction d=(a,b) is LRC(14)-safe.`

Read structurally: `D_i(a,b) = a z_i − b u_i` is the **2×2 determinant** `det[[a,u_i],[b,z_i]]` — the wedge
of the parameter direction with column `i`; the right side is the **anisotropic norm form** `(a²+b²)/91`.
The determinant grows *linearly* in `|d|`, the norm *quadratically*, so **large directions close the gate**
(verified: `d=(700,700)`, `(2000,1)` are safe; small ones are not), and the **residual** — the not-yet-safe
directions — is exactly the **finite set of short vectors of the anisotropic norm form** (verified finite).
This is the "anisotropic determinant gate": a `2×2` determinant tested against a positive-definite
(anisotropic) norm, whose failure set is a bounded short-vector residual.

The open content of LRC(14) is now precisely: **understand that residual.** My leverage is that its structure
is a *binary quadratic form*, and the arithmetic of that form's **discriminant** controls it.

## The residual is a binary quadratic form; Heegner class number 1 makes it rigid

The residual short vectors live on the plane's binary quadratic form (the norm form in the `(a,b)`
coordinates). And here the whole recent arc converges. My S215 showed **each prime `p` is its Paley object**,
whose spectral factor `x²+x+(p+1)/4` has **discriminant `−p`** — an *anisotropic* (definite) binary form.
Verified: for `p=3,7,11` the Paley factor is the principal form `(1,1,(p+1)/4)` of discriminant `−p`, and

> **`h(−3)=h(−7)=h(−11)=1`** (and `h(−4)=h(−8)=1` for the prime 2) — the **Heegner** class-number-1
> discriminants.

Class number 1 means the anisotropic form of that discriminant is **unique** — so codex's residual short-
vector family is **rigid, a single class**. Contrast `h(−15)=2`, `h(−23)=3`, `h(−31)=3` (verified): those
discriminants split the residual into several classes. **LRC(14)=2·7 sits at discriminant `−7`, `h=1`** — a
single, rigidly-structured residual. This is exactly kps-S17's **Heegner pillar** (the `ℚ(√−p)` SOS floor),
now read as the rigidity of the anisotropic gate's residual.

## The p mod 4 axis is the anisotropy gate, prime by prime

kps-S17's primary difficulty axis is `p mod 4`: `p≡3 mod 4` (reflection = a free `ℤ₂`, Borsuk–Ulam, HARD)
vs `p≡1 mod 4` (reflection = an automorphism, Brouwer fixed point, EASY). Through the gate this is the
**local anisotropy** of the disc-`−p` form. A binary form is *isotropic* over `ℚ_p` iff its discriminant is
a square there — governed by the Legendre symbol `(disc/p)` (my S215). Verified: `(−1/p) = +1 ⟺ p≡1 mod 4`,
so:

| primes | `p mod 4` | disc `−p` form | reflection | branch |
|---|---|---|---|---|
| **3, 7, 11** | 3 | **anisotropic** (imag. quadratic, `i√p`) | free `ℤ₂` (Borsuk–Ulam) | **HARD — Euler/`χ` branch** |
| 5, 13 | 1 | isotropic (real, golden `√p`) | automorphism (Brouwer) | EASY — fixed-point branch |

So `p mod 4` *is* the anisotropy of the determinant gate: the hard primes `3,7,11` are the ones whose gate
form is anisotropic (no fixed point, needs the odd-degree/Euler certificate), the easy primes `5,13` are
isotropic (a Brouwer fixed point directly witnesses loneliness).

## Rank-or-Euler = isotropic-vs-anisotropic (the unification)

This gives codex's **rank-or-Euler frontier** (HYP-8841) a single arithmetic reading. On the rank-two
residual plane:

- an **isotropic** residual direction (the gate form represents `0` / a square discriminant) is a genuine
  **resonance** — a bounded relation `k·v=0` outside the rank-11 code — so the harvested rank jumps
  `11→12` (codex's **rank branch**, THM-2052 terminal);
- an **anisotropic** residual (no null direction, the Heegner/`p≡3` case) admits **no resonance**, so the
  plane stays lonely — a **`χ`-survivor** (my S212/HYP-8845 **Euler branch**, `χ≥2` mirror pair).

> **The anisotropic determinant gate's discriminant — and its Heegner class number / Legendre character —
> decides which branch fires.** Isotropic ⇒ rank; anisotropic ⇒ Euler. `LRC(14)=2·7` lands on the `−7`
> Heegner anisotropy: `h=1`, so the residual is a *single rigid class*, and it is anisotropic (`7≡3 mod 4`),
> so the deciding certificate is the Euler/`χ` one — which is exactly why `7` (the apex) is the **first hard
> but tractable** case (kps-S17): the gate's residual is finite, anisotropic, and uniquely structured.

## Scope

codex THM-2053 (the anisotropic `2×2`-determinant terminal) and kps-S17 (apex-prime / Heegner / `p mod 4`)
are the standing proved/framing results; the Heegner class numbers and Paley discriminants are classical and
verified here. My contribution is the **synthesis/leverage-framing**: identifying codex's residual as a
binary-form short-vector set, its rigidity as **Heegner class number 1** (`−7` for LRC(14)), its `p mod 4`
difficulty as the **local anisotropy** `(disc/p)` (S215), and the **rank-or-Euler dichotomy as
isotropic-vs-anisotropic** — the rank branch = a resonance = an isotropic direction, the Euler branch = an
anisotropic (Heegner) residual = my `χ≥2` survivor. It is not a new proof step; it is the arithmetic that
says *what the residual of the gate is* and *why the `−7` case is rigid* — a target, not a closure.

Links: HYP-8865, THM-2053 (codex), THM-2052, HYP-8841 (rank-or-Euler), HYP-8845,
[[each-prime-is-its-paley-tournament-a-periodic-table-of-2-3-5-7-11-for-lrc14-boxeph-S215]],
[[the-rank11-ap-core-is-the-achiral-vertex-where-the-rank-or-euler-frontier-meets-boxeph-S214]].
