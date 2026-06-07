---
id: THM-435
name: pfaffian-parity-ladder-the-even-odd-seam
status: PROVED (mod-2 reduction) + VERIFIED (exhaustive n=3,4,5,6)
date: 2026-06-07
session: claudebox-2026-06-07-S713
depends_on:
  - THM-174   # det(I+2A)=Pf(S)^2 (even), (1^T w)^2 (odd); perfect-square
  - THM-120   # |Pf| separates topological phases at n=6 (values {1,3,5,7,9})
  - THM-305   # v_2(#tournaments)=(n-1)/2 (the 2-adic seam this sits beside)
related:
  - "Redei: #Hamiltonian-paths H(T) is odd for all T (the odd-n-agnostic parity sibling)"
provisional_id: true   # THM counter contested (THM-421/431 already collided); renumber at PR
---

# THM-435: The Pfaffian parity ladder — the 2-adic even/odd seam, made algebraic

The Pfaffian exists only in even dimension; for a tournament its skew-adjacency matrix `S = A - A^T`
is the canonical even/odd test object. This theorem pins the parity structure across the seam and
identifies `Pf` odd as the even-`n` sibling of Redei's "H is odd".

## Statement

Let `T` be a tournament on `n` vertices, `A` its 0/1 adjacency, `S = A - A^T` its skew-adjacency
(`S_ij = +1` iff `i->j`, `-1` iff `j->i`).

1. **(even `n`)** `Pf(S_T)` is **odd**; equivalently `det(S_T) = Pf(S_T)^2` is an odd perfect square.
2. **(odd `n`)** `S_T` is singular of rank `n-1`; the cofactor vector `w`, `w_i = (-1)^i Pf(S_{T-i})`
   (delete vertex `i`), spans `ker S_T`. Each `w_i` is the Pfaffian of an **even** `(n-1)`-vertex
   subtournament, hence **odd**; and `1^T w = sum_i w_i` is a sum of `n` (an odd count of) odd
   integers, hence **odd**.
3. **(all `n`)** `det(I + 2A)` is an **odd perfect square**, therefore `det(I+2A) ≡ 1 (mod 8)`.
   (Even `n`: `= Pf(S)^2`, `Pf` odd. Odd `n`: `= (1^T w)^2`, `1^T w` odd — THM-174.)
4. **(complement mod-4 phase)** Under `T -> T^op` (reverse all arcs, `S -> -S`):
   `Pf(S_{T^op}) = (-1)^{n/2} Pf(S_T)`. So for `n ≡ 0 (mod 4)` the Pfaffian is complement-invariant;
   for `n ≡ 2 (mod 4)` it flips sign.

## Proof

**Mod-2 reduction (1).** Every off-diagonal entry of `S` is `±1 ≡ 1 (mod 2)`, the diagonal `0`, so
`S ≡ J - I (mod 2)`, the adjacency matrix of `K_n`. The Pfaffian is an integer polynomial in the
entries, so `Pf(S) ≡ Pf(J-I) (mod 2)`. Over `GF(2)` the matching signs vanish, leaving
`Pf(J-I) ≡ #{perfect matchings of K_n} = (n-1)!! (mod 2)`. For even `n`, `(n-1)!!` is a product of odd
numbers, hence odd. Thus `Pf(S_T) ≡ 1 (mod 2)`. Since `det(S)=Pf(S)^2` (THM-174), `det(S)` is an odd
perfect square. (Equivalently `det(S) ≡ det(J-I) = (n-1)(-1)^{n-1} ≡ n-1 ≡ 1 (mod 2)`.)

**Kernel ladder (2).** For odd `n`, `S` skew-symmetric of odd order has `det S = 0` and (generically,
and here for tournaments) rank `n-1`, so `corank 1`. The classical Pfaffian-cofactor identity gives
`S w = 0` for `w_i = (-1)^i Pf(S_{T-i})` (verified `S w = 0` on all `8` and all `1024` tournaments at
`n=3,5`). Each deleted matrix `S_{T-i}` has even order `n-1`, so by (1) `Pf(S_{T-i})` is odd; hence
every `w_i` is odd. As `n` is odd, `1^T w` sums an odd number of odd terms and is odd.

**Unified (3).** `I + 2A ≡ I (mod 2)` so `det(I+2A)` is odd; it is a perfect square by THM-174. An odd
perfect square is `≡ 1 (mod 8)`.

**Complement phase (4).** For a `2m × 2m` matrix, `Pf(cM) = c^m Pf(M)`; with `c = -1`, `m = n/2`,
`Pf(-S) = (-1)^{n/2} Pf(S)`. Since `S_{T^op} = A^T - A = -S`, the claim follows. ∎

## Verification

`04-computation/pfaffian_even_odd_seam_s713.py` (exact integer Pfaffian + Bareiss determinant):
- (1) `Pf` odd and `det(S)=Pf^2`: `64/64` (n=4), `32768/32768` (n=6); `|Pf|` values `{1,3}` and
  `{1,3,5,7,9}` (reproduces THM-120).
- (2) `S w = 0`, all `w_i` odd, `1^T w` odd, `det(I+2A)=(1^T w)^2`: `8/8` (n=3), `1024/1024` (n=5).
- (3) `det(I+2A)` perfect square and `≡ 1 (mod 8)`: `8/8, 64/64, 1024/1024, 32768/32768` (n=3..6).
- (4) `Pf(-S)=(-1)^{n/2}Pf(S)`: `64/64` (n=4, invariant), `32768/32768` (n=6, sign-flip).

## Reading: the Pfaffian IS the even/odd seam

- **existence** — `Pf` lives only in even dimension; odd `n` kills it (`det S = 0`). The seam is the
  Pfaffian's domain of definition.
- **parity** — even `n`: one odd integer `Pf(S)`; odd `n`: a whole kernel vector of odd integers
  `w_i = ±Pf(S_{T-i})`. The two sides are linked by **vertex deletion** (`n <-> n-1`, the project's
  recurring `n-1 things <-> n things` map): an odd-`n` tournament's kernel entries are the even-`(n-1)`
  Pfaffians of its vertex-deleted subtournaments.
- **2-adic content** — `det(I+2A) ≡ 1 (mod 8)` always; `det(S)` has even 2-adic valuation iff `n` even
  (it is `Pf^2`), and is `0` (valuation `∞`) iff `n` odd. This sits beside THM-305's `v_2(T(n))=(n-1)/2`.
- **mod-4 phase** — complement acts by `(-1)^{n/2}`, the reciprocity/XNOR mod-4 sign the project tracks
  (Gauss sum `√p` vs `i√p`, `p ≡ 1` vs `3 mod 4`): the Pfaffian carries it on the `n ≡ 0` vs `2 mod 4`
  refinement of the seam.
- **Redei sibling** — `H(T)` (#Hamiltonian paths) is odd for all `n`; `Pf(S_T)` is odd for even `n`.
  Two odd-rigidity invariants; for even `n` a tournament carries both. (Computed: both odd on all
  tested `T`; the tighter `|Pf| ≡ H (mod 4)` holds at `n=4` but is coincidental — it fails at `n=6`.)
