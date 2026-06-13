---
source: opus-2026-06-03-S586 (remote-control)
status: SYNTHESIS (inspirational) — triangular numbers as the bridge of addition/multiplication & odd/even; the +2 ladder = odd gnomons (→squares), the ×2 doubling = geometric; the worry-set modulus 2n-1 = √(8·C(n,2)+1) ties the triangular pair-count to the mod-(2n-1) worry structure
tags: [LRC, triangular, addition, multiplication, odd-even, plus2, times2, staircase, simplex, 2n-1, worry-set, n14]
---

# The triangular bridge: how +2 and ×2, odd and even, add and multiply

**Prompt (user):** understand odd/even, multiplication/addition, the +2 and ×2
recursions; see how they create the triangular numbers; poke for inspiration.

The triangular numbers are the **fixed point where addition and multiplication meet**,
and the two LRC recursions (`+2` ladder, `×2` doubling) are their two faces. Poking
around turned up a clean identity that ties the triangular *pair-count* directly to the
worry-set's modulus `2n-1`.

## 1. Triangular = the add/mult bridge (with odd/even inside)

```
T_k  =  1 + 2 + ··· + k          (ADDITION — the staircase rows / the antiderivative of ℤ)
     =  k(k+1)/2                  (MULTIPLICATION — consecutive product; one factor EVEN, /2)
```
The `/2` is exactly where **odd/even** lives: of the two consecutive `k, k+1`, one is
even and supplies the `2`. So a triangular number is *(odd)·(even)/2* — the parity
interplay made a number. Two more bridges:
```
T_k + T_{k-1} = k²              (two consecutive triangles = a square)
8·T_k + 1     = (2k+1)²         (octupled triangle + 1 = an ODD square)   ← the key
```
`8 = 2³` (three doublings); `2k+1` odd. So **triangular ↔ odd-square** via three `×2`'s.

## 2. The +2 ladder = ODD gnomons (→ squares)

A tournament on `n` vertices has `C(n,2) = T_{n-1}` arcs (the simplex/staircase/
permutohedron trinity). The LRC ladder `n→n+2` adds
```
C(n+1,2) − C(n-1,2) = 2n+1          (an ODD gnomon)
```
arcs per step — `7,11,15,19,23,27,…` So the **additive `+2` recursion is the world of
odd gnomons**, whose running sums are squares (`T_k+T_{k-1}=k²`). The tournament arc
count *is* the triangular staircase, built by odd increments.

## 3. The ×2 doubling = the geometric/2-adic face

The doubling `n→2n` acts geometrically: `T_{2k} = 4T_k − k = 3T_k + T_{k-1}`. This is the
multiplicative/2-adic face (S579/S580): the prime-`q` triangle lifts to `2q` (`T_6=21 →`
the n=14 structure), and the dynamical-rigidity fragmentation (S585) is *this* face
failing while the `+2`/odd-gnomon face keeps working.

## 4. The inspiration: `2n-1 = √(8·C(n,2)+1)` — the worry-set modulus is the odd-square face

`8·C(n,2) + 1 = 8·T_{n-1} + 1 = 4n(n-1)+1 = (2n-1)²`. So
```
2n - 1  =  √( 8 · (pair count) + 1 )      — the ODD-SQUARE ROOT of the octupled triangular pairs.
```
For `n=14`: `8·91 + 1 = 729 = 27²`, and **`2n-1 = 27` is exactly the modulus of my 64
self-converse worry-classes** (the antipodal transversals mod `27`, S570). So the
worry-set's modular home `ℤ/(2n-1)` is *not arbitrary* — it is the odd-square root of the
triangular pair-count, with the `8 = 2³` supplying the three doublings between the
*additive* pair-staircase and the *multiplicative* modular shell. The triangular bridge
(add↔mult via `8T+1=odd²`) **is** the bridge between the pair/pinch structure (additive,
`C(n,2)`) and the worry-set's modular structure (multiplicative, `mod 2n-1`).

## 5. The n=14 = 2·7 triangular tower

`7 · T_k = 7, 21, 42, 70, 105, 147` — the prime `7` times the triangulars:
- `21 = T_6 = C(7,2)` — arcs on the prime-`7` (the 2-adic kernel, S579);
- `42 = 2·21 = C_5` — triangulations of the **heptagon**;
- `7 − 21 + 42 = 28 = T_7` (Euler χ); `7 + 21 + 42 = 70 = 7·T_4`.

So `n=14`'s triangular signature is the **doubling of the prime-7 triangle** `21→42`,
exactly the `×2` face acting on the additive pair-count of the kernel.

## 6. The 2×2 (the unifying picture)

```
                 ODD                         EVEN
ADD (+2 ladder)  odd gnomons 2n+1   →  squares   |   even arithmetic progressions
MUL (×2 double)  odd part (the kernel q)          |  powers of 2 (the 2-adic tower)
                              \                   /
                               TRIANGULAR  T_k = k(k+1)/2 = (odd²−1)/8
                               = C(n,2) pairs ;  2n-1 = √(8T_{n-1}+1)
```
- The **+2/additive/odd** face = the staircase, the pinch-pair count, the gnomons.
- The **×2/multiplicative/2-adic** face = the doubling tower, the dynamical rigidity.
- They **meet at the triangular pair-count `C(n,2)`**, and the bridge `8T+1=(2n-1)²`
  carries the additive pairs to the multiplicative worry-modulus.

## 7. Creative hypothesis

> **H (triangular worry-modulus).** The worry-set's modular shell `ℤ/(2n-1)` is the
> odd-square face of the triangular pair-count: `2n-1 = √(8·C(n,2)+1)`. The two LRC
> recursions are the additive (`+2`, odd gnomons, the staircase pairs) and multiplicative
> (`×2`, the 2-adic doubling) faces of this triangular bridge. The `n=2q` obstruction is
> where the multiplicative face fragments (S585) while the additive face (the `2n-1`
> modular transversal structure) stays intact — so the worry-set's *static* modular
> rigidity (the 64 classes mod 27) survives even as the *dynamical* doubling fragments.
> **Test:** check whether the worry-class count `2^((n-2)/2)` and the divisor structure of
> `2n-1` (composite at `n=14`, `27=3³`) carry the `n=14` residual — i.e. the composite
> `2n-1` is the additive-face shadow of the same `2·7` seam.

## 8. Honest status

- **Exact/verified:** `C(n,2)=T_{n-1}`; `+2` adds the odd gnomon `2n+1`; `T_k+T_{k-1}=k²`;
  **`8·C(n,2)+1=(2n-1)²`** (so `2n-1=√(8T_{n-1}+1)`, `n=14: 27²=729`); the `7·T_k` tower.
- **Inspiration (rigorous link):** the worry-set modulus `2n-1` *is* the odd-square root
  of the octupled triangular pair-count — the additive pair structure and the
  multiplicative worry-modulus are two faces of one triangular bridge.
- **Hypothesis H:** designed, partially-grounded (the identity is exact; the
  "additive-face survives / multiplicative fragments at `2q`" reading is structural,
  untested).

**Artifacts:** `04-computation/lrc_triangular_recursion_s586.py`→`05-knowledge/results/
lrc_triangular_recursion_s586.out`. Builds on the staircase/simplex trinity, `7-21-42`,
S570 (mod 2n-1 worry-classes), S579/S580 (doubling tower), S585 (dynamical rigidity).
New: **HYP-2128**.
