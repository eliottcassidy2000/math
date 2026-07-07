# First-gap emptiness is non-monotonic in N — it is arithmetic, not window-width

**opus-2026-07-06-S118.** Asked to understand the LRC for other numbers of runners and leverage
it for proof progress, I mapped the first-gap emptiness across the number of speeds `N` and found
a fact that **corrects the fleet's "window narrows ⟹ empty" narrative**: emptiness is *not
monotone* in `N`, so it cannot be driven by the window width. It is arithmetic.

Throughout, `N` = number of nonzero speeds; the LRC bound is `1/(N+1)`; the first gap is the
window `(1/(N+1), 2/(2N+1))`; our target LRC(14) is `N = 12` (AP `{1,…,12}`, window `(1/13, 2/25)`).

## The data (M computed exactly, `b ≤ 2·max`, so the maximizer is captured — S109 lever)

| `N` speeds | window | first gap nonempty? | witness |
|---|---|---|---|
| 6 | (1/7, 2/13) | **yes** | `{1,5,6,11,16,17}` → `5/33` |
| 7 | (1/8, 2/15) | **yes** | `{1,2,3,4,5,7,18}` → `3/23` |
| 12 | (1/13, 2/25) | **no** (fleet-verified, ~9k families) | — |
| 13 | (1/14, 2/27) | **yes** (verified this session) | `{1,…,11,13,36}` → `3/41` |

The `N = 13` member is real: `M({1,…,11,13,36}) = 3/41 ∈ (1/14, 2/27)`, and `3/41` is exactly the
`N = 13` **mediant** `3/(3N+2)`.

## The consequence: emptiness is not metric

The window width is `1/((N+1)(2N+1)) ~ 1/(2N²)` — **strictly decreasing** in `N`. If emptiness
were a consequence of the window narrowing (the Selberg-width / ladder-density picture), then a
narrower window would be *more* likely empty, and emptiness would be **monotone**: once empty at
`N = 12`, empty for all larger `N`. But `N = 13` has a *narrower* window and is **nonempty**.

So **the window width is necessary-but-not-sufficient** (the same shape as every structural lens
before it): narrowing bounds the *complexity/depth* of achievable in-window values (fewer of
them), but it does **not** decide whether even the shallowest — the mediant — survives. That is
decided by the **arithmetic of `q = 3N+2`** (the mediant denominator).

## What the arithmetic looks like

The mediant far element is *universal*: in both achievers, the binding pair is `{5, 3(N-1)}`,
summing to `3N+2 = q` (since `q - 5 = 3(N-1)`). What varies with `N` is whether a near-tight base
completes this to `M = 3/q` exactly:

- **`N = 7` (`q = 23`) and `N = 13` (`q = 41`): mediant achievable.**
- **`N = 5` (`q = 17`, near-exhaustive), `N = 9` (`q = 29`), `N = 12` (`q = 38 = 2·19`): not.**

Primality is *not* the criterion (`N = 5`, `q = 17` prime, is unachievable). The achievable cases
so far are `N ≡ 1 (mod 6)` (`7, 13`), i.e. `6 | (N-1)`, so `18 | 3(N-1) =` the far element — but
two data points is a lead, not a law. The point that *is* solid: whether the mediant is
achievable is an **N-specific arithmetic property of `q = 3N+2`**, decoupled from the window width.

## Proof implications (the leverage)

1. **Guard-rail (like S111/S114).** Any argument that closes `N = 12` *via the window width or a
   uniform Selberg-width bound alone* is incomplete: it would falsely predict `N = 13` empty. The
   `N = 13` member is the standing counterexample. The width is real but not decisive.
2. **The right route is arithmetic — Fan–Sun's gcd template.** `N = 12` emptiness must come from
   the *factorization* `q = 38 = 2·19`: the mediant residue-covering system at `q = 38` is
   obstructed where the prime-`q`-ish cases `23, 41` are feasible. This is exactly the flavor of
   Fan–Sun's `n = 4` gap-emptiness proof (gcd/divisibility case analysis, S116). The obstruction
   is likely a CRT/covering-system infeasibility mod the factors `2` and `19`.
3. **New obligation (O-arith).** Characterize the arithmetic condition on `q = 3N+2` under which
   the mediant is achievable, and prove `q = 38` fails it. Combined with the complexity bound
   (S117: the mediant is the shallowest, deeper values need higher complexity), this closes the
   bounded case at `N = 12` by an *arithmetic* obstruction, not a metric one.

## Caveat recorded

Constructing gap members is genuinely hard, and weak searches false-report "empty": my first two
base-generators reported `N = 6` empty (it is not — `5/33` exists). Only a diff-`d` AP base with
borders and the correct `3(N-1)` far element recovered the known members and found `N = 13`. Any
cross-`N` emptiness claim must use these shapes or it is a false negative — the fleet's stronger
sweeps remain authoritative for the `N = 12` emptiness itself. What this session establishes with
confidence is the *nonemptiness at `N = 13`* and therefore the *non-monotonicity*, which is enough
to relocate the `N = 12` obstruction from metric to arithmetic.
