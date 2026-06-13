---
id: HYP-1880
status: OPEN
source: codex-2026-05-31-S410
related:
  - HYP-1813
  - HYP-1844
  - HYP-1858
  - HYP-1859
  - THM-357
  - THM-360
  - THM-367
---

# HYP-1880: LRC endpoint debt is Bruhat-Tits descent

## Statement

For a Lonely Runner denominator `n=k+1`, attach to every rational endpoint
`t` its local depth vector

```text
d_p(t)=v_p(denominator(t)),   p | n,
```

after reducing `t` as a fraction.  This is the product of the rank-one
Bruhat-Tits tree shadows at the primes dividing `n`.

The conjectural descent principle is:

```text
an all-protected endpoint core cannot have zero divergence in this product
tree.
```

Equivalently, every attempt to protect an exposed boundary layer either leaves
a lonely endpoint, leaves a positive open gap, or exports endpoint debt to a
strict descendant layer faster than the finite speed budget can close it.

## Exact Local Ledger

THM-360 is the first local tree statement.  The nonzero unit endpoints
`a/n`, `gcd(a,n)=1`, all have depth vector

```text
(v_p(n))_{p|n}.
```

A speed protects such a unit endpoint only if it is divisible by `n`.  In
Bruhat-Tits language, only an `n`-gate translates the whole unit boundary
layer back to the root.

But an `n`-gate does not delete debt.  If `v=n*q`, then every endpoint owned
by a `v`-interval is

```text
(n*m +/- 1)/(n*v) = (n*m +/- 1)/(n^2*q).
```

For every prime `p | n`, the numerator `n*m +/- 1` is a p-adic unit, hence
the exported endpoint depth is exactly

```text
2*v_p(n) + v_p(q).
```

So gate repair moves the obstruction to child vertices of the local tree.

## Evidence

`lrc_bruhat_tits_descent_s410.py` prints this ledger for `n=14` and `n=16`.

For `n=14`, the initial tight row exposes exactly the unit layer:

```text
exposed_depths {2:1,7:1}: 6
```

Replacing a speed by the `14`-gate protects the unit skeleton but creates
exposed descendants, including depths `{2:2,7:2}`.  The S380 14-multiple
ladder has tiny positive gap `gap/th=0.002706`, but its exposed endpoints are
already exported:

```text
{2:2,7:2}:120, {2:4,7:2}:48.
```

For `n=16`, the same story is rank-one:

```text
initial exposed depth {2:4}: 8
drop 15 add 16 exposes {2:8} among other layers
odd units plus seven 16-gates exposes {2:8}, {2:9}, {2:10}.
```

THM-367 becomes the exact local radial law: for owner `u=2^k`, every endpoint
has depth `4+k`; a protector `p=2^j*q` translates upward by `j` edges, while
the odd residue class decides which branches are hit.  The warning is that
owners `u>=16` have a stable local nine-fan, so the proof cannot be purely
local.  It needs a global product-tree potential.

## Predictions

1. For `n=14`, a useful proof invariant should live on the product of the
   `2`-adic and `7`-adic trees.  The best disproof-looking rows will keep
   moving exposed mass from `{1,1}` to `{2,2}` and then into one-sided
   descendants, never closing the product frontier.
2. For `n=16`, the correct invariant is a one-prime BT flow with residue
   branch data.  The nine-fan locally closes a branch, but should consume the
   same residue budget needed to close descendant branches.
3. A labelled endpoint cycle from THM-365 should have nonzero BT divergence:
   multiplying around the cycle cannot return to the same product-tree vertex
   once strict protection inequalities and primitivity are included.

## Sources

- `04-computation/lrc_bruhat_tits_descent_s410.py`
- `05-knowledge/results/lrc_bruhat_tits_descent_s410.out`
- `07-reflections/lrc-bruhat-tits-descent-s410.md`
