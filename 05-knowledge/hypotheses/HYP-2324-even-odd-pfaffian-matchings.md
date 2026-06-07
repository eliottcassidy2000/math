# HYP-2324 — Even/odd and the Pfaffian: det = Pf², the Pfaffian is a signed matching-sum, and it is odd because (2n−1)‼ is odd

**Session:** S646
**Status:** CONFIRMED (det=Pf² and oddness formalized at the available scale; structure verified)
**Provenance forward:** math-lean `Math/Tournaments/PfaffianParity.lean` (sorry-free)
**Prompt:** understand even/odd and the Pfaffian. Develops S645 (tournament discriminant = `det(skew)`).

---

## 0. The Pfaffian, and why it is the even/odd object

The **Pfaffian** of a `2n×2n` skew-symmetric matrix is a signed sum over the **perfect matchings** of
`{1,…,2n}`:
```
  Pf(M) = Σ_{matchings π} sgn(π) · ∏_{(i,j)∈π} Mᵢⱼ ,        det M = Pf(M)².
```
It exists **only for even order** — a perfect matching of `2n+1` points does not exist — so the Pfaffian
is the *parity object*: even `n` ⟹ Pf defined, `det = Pf²` a square; odd `n` ⟹ no matching, `Pf = 0`,
`det = 0`. This is the S645 master switch, now with its mechanism.

---

## 1. Why the tournament Pfaffian is ODD (the S645 finding, explained)

The number of perfect matchings of `2n` points is the **double factorial**
```
  (2n−1)‼ = (2n−1)(2n−3)⋯3·1   = a product of ODD numbers ⟹ ODD.
```
So for a `±1` (tournament) skew matrix, `Pf(M)` is a signed sum of **`(2n−1)‼` many `±1` terms** — an
**odd** number of **odd** terms — hence **odd**. Therefore `det M = Pf(M)²` is an **odd square**. This is
*exactly* the S645 observation ("the tournament Pfaffian is always odd"), now derived from first
principles: the oddness is the oddness of the matching count.

Verified (`pfaffian_even_odd_matchings_s646.py`): `(2n−1)‼ = 1,3,15,105,945,10395` for `2n=2..12`, all
odd; every `±1` Pfaffian odd; `det = Pf²` exact on all tournaments through `n = 6` (`Pf ∈ {1,3,5,7,9}`).

**Formalized (math-lean, sorry-free): `Math/Tournaments/PfaffianParity.lean`**
- `skew_two_det_eq_pf_sq : (!![0,a;-a,0]).det = a²` — `det = Pf²` at `n = 2` (`Pf = a`).
- `odd_doubleFactorial : Odd ((2n+1)‼)` — the matching count is odd (the engine of Pfaffian oddness).
- (`n = 4`: `det = (af−be+cd)²` over the `3 = 3‼` matchings — verified computationally; Mathlib lacks
  `det_fin_four`/a Pfaffian, flagged.)

---

## 2. The matching count is the rising factorial of 1/2 (ties S644)

```
  (2n−1)‼ = 2ⁿ · (1/2)^(n)   where (1/2)^(n) = (1/2)(3/2)⋯((2n−1)/2) is the RISING factorial of 1/2.
```
Verified exactly (`= 1,3,15,105,…`). So the number of terms in the Pfaffian — the size of the even/odd
combinatorics — is the **rising factorial of `1/2`** flagged at the end of S644: the `1/2` is the **apex
/ `σ`-fixed point `n/2`**, and `(1/2)^(n) ~ (2n)!/(4ⁿ n!) ~ √(πn)`-scaled is the **`√π` / Wallis** constant
of the tournament tiling (CLAUDE.md fiber fraction `f(n) = (1/2)_{n−2}/(n−2)!`). So the Pfaffian's term
count, the tiling density, and the apex half-integer are one rising-`1/2`-factorial.

---

## 3. The even/odd master switch (the full picture)

| `n` | perfect matching? | `Pf` | `det` | `rank M` |
|---|---|---|---|---|
| **even** | yes (`(2n−1)‼` of them) | odd `≠ 0` | odd square | full `= n` |
| **odd** | **no** | `0` | `0` | `n − 1` (deficient) |

And the structural fact behind the rank column: **a skew-symmetric matrix has even rank, always**
(verified: `rank = n` (even `n`) or `n−1` (odd `n`) — both even). So "the discriminant vanishes," "the
rank is deficient by 1," "no perfect matching exists," and "`n` is odd" are **four names for one event**.
The Pfaffian lives only on the even side; the `n → n−2` Pfaffian recursion (Mode B, S645) stays inside
one parity class, while `n → n−1` (Mode A) flips it.

---

## 4. Synthesis: even/odd is the sign, the Pfaffian is the √
- `det = Pf²` factors the determinant — a sum over **all of `Sₙ`** (S644 rising side) — as the **square**
  of the Pfaffian — a sum over **matchings** weighted by `sgn` (the even/odd of `Aₙ`, S643/S644). The
  Pfaffian is the genuine square root that "remembers the sign"; its square forgets it.
- Everything-odd (`H` Rédei, `|Aut|` S643, `Pf` here) is the tournament's signature of the **sign-kernel**
  world; the even/odd of `n` is the master clock; the cube-root `3` sets the discriminant ceiling
  `3^{n−2}` (S645). Even/odd, the Pfaffian, the sign, and the rising-`1/2`-factorial are one structure.

## 5. New threads / handoffs
- **Build a Pfaffian in Mathlib** (genuine gap): then formalize `det = Pf²` (even `n`) and `Pf(±1) odd`.
- **Skew matrices have even rank** — formalize (Mathlib may have a partial result via `Matrix.rank`).
- **`(2n−1)‼ = 2ⁿ(1/2)^(n)`** — formalize the double-factorial/rising-`1/2` identity (ties S644).
- Is the maximal Pfaffian `3^{(n−2)/2}` (S645) the **Hadamard/Pfaffian bound** attained by the
  doubly-regular (Paley-type) tournament? (cube-root ceiling, S637/8/40/42/43.)
