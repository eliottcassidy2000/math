---
source: opus-2026-06-02-S567 (remote-control)
status: SYNTHESIS — tournaments as simplex(mesh) AND polygon(dihedral); the regular/dihedral face is odd-m only; LRC worry-set parity dichotomy
tags: [LRC, simplex, polygon, dihedral, parity, regular-tournament, round, n+2-recursion, ModeB, n14, trinity, S565, S566]
---

# Tournaments are simplex and polygon: the dihedral face is "every other n", and the LRC worry-set inherits it

**Prompt (user):** tournaments are both simplices (dimensions) and regular polygons
(dihedral groups) — sometimes acting like just their outside, sometimes like the
whole mesh. Dihedral groups occur every other n; this is related to n+2 recursion.

This is right, and each clause is checkable against the project's structures.

## 1. The two faces (the trinity, grounded)

A tournament on `n` vertices IS both (the-polygon-simplex-staircase-trinity):

- **SIMPLEX / "the whole mesh":** an orientation of `K_n` = the 1-skeleton of the
  `(n-1)`-simplex; all `C(n,2)` arcs, all `C(n,3)` potential 3-cycles, …. The
  general tournament needs the whole mesh — counted by **A000568**.
- **POLYGON / "just the outside":** the permutohedron/zonotope point; the **round**
  tournaments are exactly those reconstructible from their *outside* — a cyclic
  **gap-necklace** on the circle (each out-set a clockwise arc), counted by
  **A000016** (S565 orbit side). A round tournament *acts like its boundary*; a
  general one *acts like the full mesh*.

So "outside vs whole mesh" = **round (A000016) vs general (A000568)** — the
dual-Burnside orbit/fix split (S565).

## 2. The dihedral face is ODD-m only ("every other n") — verified

A **regular** tournament (every out-degree `(m-1)/2`, the maximally symmetric
rotational one with `Aut = C_m` and, with reversal `i↦-i`, the **dihedral `D_m`**)
exists **iff `m` is odd**:

| m | 2 | 3 | 4 | 5 | 6 | 7 | … | 13 |
|---|---|---|---|---|---|---|---|---|
| regular/dihedral exists? | – | ✓ D₃ | – | ✓ D₅ | – | ✓ D₇ | … | ✓ D₁₃ |

`(m-1)/2 ∈ ℤ ⟺ m odd`. So the regular-polygon / dihedral face genuinely appears at
**every other vertex count** — exactly your claim. Even `m` has *no* regular
tournament: it is forced into the irregular "simplex" mode (split out-degrees
`{d, d+1}`).

## 3. The LRC worry-set inherits the parity (the payoff)

The loneliest configuration is the regular `n`-gon (roots of unity, geometric `D_n`)
for **all** `n`. But as a **tournament** the runner block lives on `m = n-1`
vertices, and the dihedral/regular face is available **iff `m` odd ⟺ `n` even**:

> **EVEN `n` (the doubled primes 10, 14, 22, …):** `m = n-1` odd ⟹ the runner block
> sits on the **polygon/dihedral face** — the near-regular rotational encirclement
> (S566), with the apex `(q,q)` reflection symmetry (S547). The worry-set is
> **polygon-shaped**.
>
> **ODD `n`:** `m` even ⟹ **no** regular runner tournament; tightness is forced
> into the **simplex/irregular** mode.

Honest precision (out-degrees of the AP block at `t=1/n`): it is the **`n`-gon
minus the observer vertex** (n-1 of n equally-spaced points), so it is *near*-
regular (`{5:6, 6:7}` at n=14, `{2:3,3:3}` at n=7) with antipodal ties, and its
exact symmetry is the **apex reflection** (S547), not full `D_{n-1}`. The
"dihedral" lives on the full `n`-gon (observer + runners, `D_n`), marked at the
observer down to the apex axis. So the right statement is: **even-`n` LRC tightness
is the observer-marked regular `n`-gon (`D_n` broken to the apex reflection);
odd-`n` tightness has no such polygon.** This is *why* the project's hard cases are
the even doubled primes — they are where the polygon/dihedral face exists.

## 4. The n+2 recursion = the parity carrier

`n ↦ n+2` is the parity-preserving step (the project's **Mode B**, `n→n-2`,
Cayley–Dickson descent / meta-graph recursion). The dihedral face, being an
odd-`m`/even-`n` phenomenon, **recurs along `n+2`**: the regular-polygon tournaments
form the sub-sequence `m = 3,5,7,…` (`D_3,D_5,D_7,…`), and the LRC polygon-tight
worry instances form the even-`n` sub-sequence `…,10,12,14,…` within which the
doubled primes `2q` sit. So "dihedral every other n" and "n+2 recursion" are the
same parity skeleton: the simplex face is generic and grows every step; the
polygon/dihedral face switches on and off with parity and is carried by `n+2`.

## 5. What it buys for LRC

- **A parity split of the conjecture:** even-`n` (polygon/dihedral, doubled-prime)
  tightness vs odd-`n` (simplex/irregular) tightness are *structurally different*
  worry-sets. The repo's even-fold, `(q,q)`, mod-2×mod-7 gears, and dual-Burnside
  fix-side all live on the **even-`n` polygon face**.
- **A clean statement of the extremiser:** the even-`n` tight config is the
  observer-marked regular `n`-gon (`D_n`/apex-reflection); the conjecture says this
  most-symmetric polygon still leaves the observer a `≥2/n` gap. The simplex face
  (generic tournaments) is automatically lonely (S564 IGNORE / S565 orbit side).
- **Recursion target:** climb the `n+2` ladder of even-`n` polygon-tight configs
  (Mode B) — the doubled primes are the prime-channel-clean rungs (S546).

## 6. Honest status

Grounded: the simplex/polygon (round=A000016 vs general=A000568) duality (trinity,
S565), and **regular/dihedral ⟺ odd m** (verified). New synthesis: the **LRC
parity dichotomy** — even-`n` tightness is the observer-marked regular `n`-gon
(polygon/dihedral face), odd-`n` tightness is simplex/irregular — and its
identification with the `n+2`/Mode-B parity recursion. The precise symmetry is the
apex reflection (S547), not full `D_{n-1}` (the observer-marking breaks it). Not a
proof; a structural organising principle tying parity, the dihedral face, and the
doubled-prime worry-set together.

**Artifacts:** verification inline (`regular ⟺ odd m`; AP-block out-degrees).
Builds on the trinity, S565 (dual Burnside), S566 (strong/regular encirclement),
S547 ((q,q)/apex), Mode B. New: **HYP-2090**.
