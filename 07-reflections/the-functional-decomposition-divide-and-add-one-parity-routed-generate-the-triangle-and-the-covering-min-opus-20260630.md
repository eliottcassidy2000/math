# The project's functional skeleton: with a(x)=x/2 (DIVISION = descent) and b(x)=x+1 (ADDITION = observer), parity-routed f (the odd one of {x,x+1}) and g (half the even one = ⌈x/2⌉) give f·g = T_x EXACTLY — "everything is the triangle" = the parity-routed product; recursing g IS the 2-adic descent (g halves to 1, f is the odd core at each level); the dyadic group ⟨a,b⟩ is the Stern-Brocot/Farey tree of the observer's escape [0;n−1,n]; and the covering-min denominator Φ₆(n)=n(n−1)+1 = 2·T_{n−1}+1 = b applied to the doubled triangle (the lcm-killer pronic) — so divide builds the descent, add-one builds the observer, parity routes, and the product is the triangle

*opus-2026-06-30. Owner: a(x)=x/2, b(x)=x+1; f=b on evens/id on odds, g=a on evens/a∘b on odds; note f·g =
the triangular numbers; think everything recursively in this functional decomposition — even/odd,
positive/negative, addition/division. It is the skeleton: two operations (÷2, +1), routed by parity.*

## f·g = T_x, the parity-routed factorization (verified)
- `a(x)=x/2` (DIVISION), `b(x)=x+1` (ADDITION).
- `f(x) = b(x) if x even else x` `=` **the ODD one of `{x, x+1}`**.
- `g(x) = a(x) if x even else a(b(x))` `= ⌈x/2⌉ =` **HALF the EVEN one of `{x, x+1}`**.
- **`f(x)·g(x) = x(x+1)/2 = T_x`** for all x (verified 1..40). The `2` in `T_x` is absorbed by whichever of
  `x, x+1` is even — **parity routes the division**. The triangle is `(odd part)·(halved even part)`.

## Recursing g IS the 2-adic descent (the renormalization)
`g = ⌈x/2⌉` applied repeatedly halves to `1`; `f` reads the odd core at each level:
> `x=14: 14→7→4→2→1`, odd cores `f = 15,7,5,3,1`. Depth `~log₂ x`.
This is the **LRC descent** (`S → E/2`, peel the odd core, halve the even part) and the **binary expansion**
at once — `a` is the halving, `f` the odd core, the chain the descent levels. The cusp→doublet
renormalization, the staircase's 2-adic structure, and the binary tree are one recursion: **apply `g`.**

## The dyadic group ⟨a,b⟩ = the Stern-Brocot/Farey tree = the observer's escape
`a` (÷2) and `b` (+1) generate the **dyadic affine maps** — the Stern-Brocot / binary tree of the rationals.
This is the SAME tree as the observer's covering-min escape `[0; n−1, n]` (the continued fraction whose
partial quotients are the two tightest blockings). So **the functional decomposition's tree IS the
observer's Farey escape** — divide-and-add-one navigates exactly the Stern-Brocot path the covering-min sits
on.

## The covering-min IS the doubled triangle plus the observer's +1
The cleanest fusion: `Φ₆(n) = n²−n+1 = n(n−1)+1 = **2·T_{n−1} + 1**`.
> The **pronic `n(n−1)`** (the construction's lcm-killer outlier `(n−1)n`) is `2·T_{n−1}` — **the doubled
> staircase**. The covering-min denominator is `b` applied to it: `Φ₆ = b(2·T_{n−1})` — the **observer's `+1`
> on the doubled triangle**. So `M = n/Φ₆ = n/(2·T_{n−1}+1)`: the LRC covering-min is the runner count over
> twice-the-staircase-plus-the-observer. The triangle (tournaments) and the covering-min (LRC) are joined by
> exactly `×2` (the antipodal doubling) `+1` (the observer) — `a⁻¹` then `b`.

## The three dualities (the operation group)
| duality | functional content | in the project |
|---|---|---|
| **even / odd** | parity = which of `{x,x+1}` carries the factor `2` (the `f`/`g` split) | the descent level; the odd core; `H=I(Ω,2)` parity |
| **addition / division** | `b (+1)` vs `a (÷2)` | the **observer baseline** vs the **2-adic descent** |
| **positive / negative** | `b (+1)` vs `b⁻¹ (−1)` | the **antipode**: killer `≡ −1 mod Φ₆`, complement = reversal, blocking pair `{±1}` |
> Two operations (`÷2`, `+1`) and their inverses (`×2`, `−1`), routed by parity, generate the whole
> structure. `a` builds the **descent** (divide down to the doublet); `b` builds the **observer** (the
> irreducible `+1` — Rédei's `inshat=1+2…`, the LRC escape's Farey hair, the antipodal killer); parity
> routes; and the **product is the triangle** (`f·g=T_x`, the staircase, "everything is the triangle").

## The project, recursively decomposed
> **DIVIDE (`a=÷2`) builds the DESCENT; ADD-ONE (`b=+1`) builds the OBSERVER; PARITY routes; the PRODUCT is
> the TRIANGLE.** Concretely: the staircase `δ` = `f·g = T_x`; the renormalization = recursive `g` (`a` on
> the even part); the observer's baseline = `b` (Rédei, LRC escape, Farey hair); the covering-min denominator
> `Φ₆ = 2T_{n−1}+1 = b(pronic)`; the Stern-Brocot tree `⟨a,b⟩` = the escape `[0;n−1,n]`. Every object is a
> word in `a, b` read with parity.

## Status
- **Verified (opus):** `f·g = T_x` (parity-routed, 1..40); `g`-recursion = the 2-adic descent (odd cores
  `f`, depth `~log₂`); `Φ₆(n) = 2·T_{n−1}+1 = b(pronic)` (the doubled triangle + observer `+1`); `⟨a,b⟩` =
  the Stern-Brocot tree of `[0;n−1,n]`.
- **The skeleton:** two operations — `a=÷2` (descent) and `b=+1` (observer) — routed by parity, generate the
  triangle (`f·g`), the descent (recursive `g`), the observer (`b`), and the covering-min (`Φ₆=2T_{n−1}+1`).
  The three dualities (even/odd, ±, ÷/+) are the operation group.
- **Open (unchanged):** the realizability node; but it now sits literally on the `⟨a,b⟩` Stern-Brocot tree
  (the escape `[0;n−1,n]`), and the covering-min target is `n/(2T_{n−1}+1)` — runner-count over doubled-
  staircase-plus-observer.

Related: the-observer-abstraction + the-observer-on-the-tournament-side (the `b=+1` observer), the-covering-
min-as-a-function-of-n (Φ₆), everything-is-the-triangle (the staircase = `f·g`), the-descent-product-is-
renormalization (recursive `g`); OPEN-Q-108/039.
