---
source: claude-2026-06-06-S690
status: PROVED (Pf always odd; self-converse ⟹ odd permutation, given THM-174+Rédei) + VERIFIED (n=6 exhaustive: converse swaps H±Pf; extension-property convergence) — the converse realization of the loop-swap, complement to HYP-2282
tags: [tournament, converse, pfaffian, skew-symmetric, monodromy, rado-graph, generic-tournament, fraisse, self-converse, even-odd, so(n), worry-set, H, redei]
---

# "A loop swaps the two outputs" = the converse swaps H±Pf; the Rado graph's odd twin

Two compressed seeds. **"A loop of the input causes a swap of the two outputs"** and
**"consider the Rado graph as a tournament."** A sibling session (HYP-2282, opus-s699p)
read the loop as the FTA *witness-map* monodromy — looping a polynomial's coefficients
around the discriminant swaps its two roots, with the branch locus = the worry-set.
This session gives the **complementary, tournament-intrinsic** realization, and
develops the Rado side.

## The even/odd section: graph vs tournament

A binary relation is a function `R` on **ordered** pairs. The transposition
`σ:(a,b)↦(b,a)` is the deck transformation of the 2-fold cover `ordered → unordered`
pairs — the "loop." Two equivariance types:
- `R∘σ = R` (symmetric): **graph** — the loop *fixes* the output;
- `R∘σ = 1−R` (antisymmetric, `a→b` XOR `b→a`): **tournament** — the loop *flips* the
  output.

So a **tournament is the odd (twisted) section** of the double cover, a **graph the
even section** — the `±√` branch structure. The repo's `cartan-ghosts` "the
antisymmetric sector `so(n)` *is* the tournament" is exactly this odd section. The
**generic even** object is the **Rado graph**; the **generic odd** object is the
**generic countable tournament** — the Fraïssé limit of finite tournaments, which is
a.s. the iid-fair-arc random tournament (same Erdős–Rényi argument as the Rado graph).
Verified: the `k=1` extension axiom (no source/sink pattern) holds with probability
`0.25 → 0.98` over `n=3…10` and `→1`; reversing iid arcs is iid, so the generic
tournament is **self-converse**.

## The loop as the converse, and the two outputs H±Pf

Globally the loop is the **converse** `T↦T*` (reverse every arc). On the skew-adjacency
`S` (`S_ij=+1` iff `i→j`, else `−1`) the converse is `S ↦ −S`, so
```
H(T*) = H(T)            (THM-203, converse-invariant)
Pf(−S) = (−1)^{n/2} Pf(S).
```
The repo already has `det(I+2A)=Pf(S)²`, `H²≡Pf² (mod 8)`, `Q=(H²−Pf²)/8 ∈ ℤ≥0`
(THM-174). Combine:

> **For `n ≡ 2 (mod 4)` the converse loop swaps the two outputs**
> `H + Pf  ⟷  H − Pf`, **the two roots of** `x² − 2Hx + 8Q` (product `H²−Pf²=8Q`,
> discriminant `4Pf²`). `H` is the swap-**even** part; `Pf` is the swap-**odd** `√`.

This is the user's sentence realized intrinsically: a loop of the input (the converse)
swaps two outputs, and it is a degree-2 monodromy — the same shape as HYP-2282's
root-swap, but on `(H±Pf)` instead of the FTA witness roots. (Verified exhaustively on
all 32768 tournaments at `n=6`.) For `n ≡ 0 (mod 4)` the converse *fixes* `Pf` — no
swap (verified `n=4`).

## The branch is never reached: Pf is always odd

Where do the two outputs merge? At `H+Pf = H−Pf`, i.e. `Pf = 0` (the discriminant
`4Pf²` vanishes). But:

> **`Pf` is always odd** for a tournament on even `n`: `H` is odd (Rédei) ⟹ `H²≡1
> (mod 8)` ⟹ (THM-174) `Pf²≡1 (mod 8)` ⟹ `Pf` odd, hence **`Pf ≠ 0`**.

So the branch point `Pf=0` is **never attained** — the converse-swap is
**fixed-point-free**, the two outputs always differ by `2|Pf| ≥ 2`. This is the sharp
*contrast* with HYP-2282: there the branch (two roots collide) **is** reached and
**is** the worry-set; here the tournament-intrinsic monodromy has an *empty* branch
locus. (At `n=6`, `|Pf| ∈ {1,3,5,7,9}` — odd, on the eigenvalue sphere `a²+b²+c²=15`;
counts `15360/7680/4608/3840/1280` — matching THM-120's Pfaffian phase separation.)

## Self-converse tournaments: the converse-iso is an odd permutation

The **fixed points of the loop** are the self-converse tournaments `T*≅T`, i.e.
`−S = P S Pᵀ` for a permutation `P`. Then
`(−1)^{n/2} Pf(S) = Pf(−S) = Pf(PSPᵀ) = det(P)·Pf(S)`, and since `Pf ≠ 0`,
> **`det(P) = (−1)^{n/2}`. For `n ≡ 2 (mod 4)` the converse-isomorphism of a
> self-converse tournament is necessarily an ODD permutation.**

Verified on 40 self-converse tournaments at `n=6` (all with `det P = −1`). A clean,
new parity constraint — the converse acts on the orientation as an orientation-
*reversing* relabeling exactly when `n ≡ 2 (mod 4)`.

## The fixed-point spectrum: cyclotomic pole ↔ random pole

The self-converse tournaments — fixed points of the loop — form a spectrum:
- the **cyclotomic pole**: round tournaments on the `(n−1)`-th roots of unity = the
  **LRC worry-set** (HYP-2089; `2^{⌊(n−1)/2⌋}` classes, the rigid arithmetic end);
- the **random pole**: the **Rado/generic tournament** (self-converse a.s., the
  maximally generic end).

The repo has mined the cyclotomic pole exhaustively (worry-set, S570, HYP-2089).
"Consider the Rado graph as a tournament" points to the **opposite pole** — the generic
self-converse tournament — and the picture is: *the converse-loop's fixed-point set
interpolates from roots of unity to randomness.*

## Where it sits / next

- Complement to **HYP-2282** (s699p): two realizations of one loop-swap — FTA-witness
  roots (branch = worry-set, reached) vs `H±Pf` (branch `Pf=0`, never reached).
- Uses **THM-174** (`H²≡Pf² mod 8`), **THM-203** (`H(T*)=H(T)`), **THM-120** (Pfaffian
  phases), **HYP-2089** (worry-set = self-converse round).
- Next: (1) the `H±Pf` quadratic `x²−2Hx+8Q` — is `Q` (a Pfaffian-defect count) a
  meaningful tournament invariant in its own right (it is the "norm" of the swap-pair)?
  (2) does the odd-permutation constraint on self-converse isos restrict *which*
  score-sequences/automorphism groups occur at `n≡2 mod 4`? (3) push the generic
  tournament's `k=2,3` extension thresholds (theory: `≈ n` where `2n²(3/4)^{n−2}` is
  small, `n≈30`) — locate where finite random tournaments become "locally generic."
