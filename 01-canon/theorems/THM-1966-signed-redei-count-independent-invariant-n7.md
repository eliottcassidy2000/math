# THM-1966 — The signed Rédei count |R| is spectral (n≤5), (spectrum,H)-measurable (n=6), and a GENUINELY INDEPENDENT invariant from n=7

**Status:** PROVED for `n ≤ 6` by exhaustive labeled census; the `n=7` independence is PROVED by
an explicit verified witness. Decoupling *counts* at `n=7` are sampled (lower bounds).
**Author:** mac-mini-2026-07-21-S160 (cont.).
**Places** the signed Rédei count `R` (THM-1936) on the moment/co-spectral ladder of **THM-1780**
(where `H` leaves the spectral floor at `n=6`).

## Object

`|R(T)| = |Σ_{Ham paths π} sgn(π)|`, the signed Rédei count (THM-1936, join-multiplicative;
a tournament invariant since relabeling only flips the sign of `R`). `H(T)` = # Hamiltonian
paths (Rédei). Group labeled tournaments by their **exact characteristic polynomial**
(= trace-moment vector, by Newton) — the *co-spectral classes* of THM-1780.

## Statement — the ladder placement of |R|

1. **`n ≤ 5`: `|R|` is SPECTRAL.** It is constant on every co-spectral class — a function of the
   trace moments, exactly like `H` (THM-1780). [exhaustive: `n=4` 3 classes, `n=5` 9 classes,
   0 `|R|`-splits.]

2. **`n = 6`: `|R|` leaves the ladder *with* `H`, coupled, but adds nothing beyond `(spectrum,H)`.**
   `|R|` splits exactly the **same 3 co-spectral classes** that `H` splits — *perfect coupling*:
   no co-spectral class splits in one but not the other (0 `H`-only, 0 `R`-only splits). Yet
   `|R|` is **`(spectrum,H)`-measurable**: `(spectrum,H,|R|)` and `(spectrum,H)` distinguish the
   *same* 32 of the 56 iso classes, so `|R|` carries **no information beyond `spectrum+H`** at
   `n=6`. On the THM-1780 witness `(0,0,12,12,10,48)`: `(H,|R|) = (13,11)` and `(17,9)` — `|R|`
   separates it too, with reversed ordering.

3. **`n = 7`: `|R|` DECOUPLES from `(spectrum,H)` — a genuinely independent invariant.** There
   exist **co-spectral tournaments with the same `H` but different `|R|`**. Explicit verified
   witness — two non-isomorphic `7`-tournaments with
   ```
     char poly  x⁷ − 9x⁴ − 12x³ − 16x² − 8x − 1     (moments (0,0,27,48,80,291,763)),
     both  H = 81,   but  |R| = 1  vs  |R| = 17.
   ```
   (Exact integer char polys equal; `H` equal; `|R|` differ; canonical forms differ.) Hence
   `(spectrum, H, |R|)` **strictly refines** `(spectrum, H)`: `|R|` is a genuinely new tournament
   invariant from `n = 7`, independent of both the spectrum and the Rédei count. The `n=6`
   coupling was a small-`n` coincidence; at `n=7` `H` and `|R|` split co-spectral classes
   *independently* (≥10 classes same-`H`/different-`|R|`, ≥9 same-`|R|`/different-`H`, sampled).

## Witness adjacency matrices (rows, `1`=arc out)

```
T0 (H=81, |R|=1):  0001111 1001110 1100011 0010010 0011011 0000001 0101000
T1 (H=81, |R|=17): 0010001 1011101 0000111 1010011 1001010 1100001 0000100
```

## Consequences

- **`(H,|R|)` is a strictly finer tournament code than either `H` or the spectrum.** At `n=6`
  it distinguishes **31** of 56 iso classes, vs `H` alone **19** and the full spectrum **28** —
  two combinatorial Hamiltonian-path counts beat the linear-algebra spectrum.
- **`|R|` is the SIGNED partner of `H`.** Both are Hamiltonian-path statistics, both #P-flavored,
  both leave the spectral floor at `n=6`; coupled there, independent from `n=7`. Equivalently the
  2-D invariant `(#even-sign, #odd-sign Ham paths) = ((H+|R|)/2, (H−|R|)/2)` carries strictly
  more than `H`. This *extends THM-1780's ladder to the signed count* and answers the S160
  highest-leverage handoff: **`|R|` is genuinely new information beyond `H`.**
- **max|R| maximizer (handoff 2, partial):** at `n=7` the maximum `|R| = 147` is achieved by the
  **QR(7) Paley (regular) tournament** (`H = 189`). The `max|R|` sequence `3,3,15,15,147` is
  **not** the double factorial (`7!! = 105 ≠ 147`).
- **Strong-atom |R|-spectrum (handoff 1):** strong tournaments realize `|R|` in
  `n=3:{3}, n=4:{1}, n=5:{3,5,7,11,15}, n=6:{1,3,5,7,9,11,13}` — note strong-`6` caps at `13`
  while the decomposable `5▷1` reaches `15` (THM-1936); characterizing the strong-atom spectrum
  is the residual open question.

→ THM-1780 (H leaves the ladder at n=6); THM-1936 (R join-multiplicative); THM-1870 (cycle
counts spectral); THM-1862 (order-join, boxeph); Rédei's theorem. Fits the tournament-invariant
zoo/ladder (klein-S399, kps THM-1780/1870).
