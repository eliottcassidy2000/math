---
source: claude-2026-06-03-S608
status: verified finite enumerations (independent confirmation + the doubling law / witness structure)
tags: [lonely-runner, n14, enumeration, tight-family, V-star, doubling, antipodal-pairs, pinch-oracle, counterexample-sweep]
---

# LRC@14 finite enumerations: the tight family, the doubling law, no counterexamples

Following the efficient finite targets of S607, this session *ran* the
enumerations. The fast integer `p_0`-cover filter (decide `M>1/14` in poly time)
plus the THM-369 pinch oracle (exact `M` for the rare `p_0=0` candidates) make
the census cheap.

## The census [verified]

Over all primitive 13-sets with speeds `≤ 20` (**77 520 configs, ~7 s**):

- the **AP `{1..13}` is the only tight config**, and
- there are **zero counterexamples** (`M < 1/14`).

So LRC@14 holds for every configuration in range — a clean finite confirmation
that reproduces and extends the repo's bounded census.

## The tight family is `{AP, V*}` [verified]

`V* = {1..11, 13, 24} = AP[12→24]` is tight (`M = 1/14`). It is the **unique**
1-swap-of-the-AP tight (checked over every single replacement up to speed 40),
and there are **no** 2-swap tights (to 30). The reason the original `≤ 1.5n = 21`
census missed it: `V*` has speed `24 > 21`.

## The doubling law [verified, new framing]

Look at every "double one element" swap `a → 2a` of the AP (valid for
`a ∈ {7,8,9,10,11,12,13}`, where `2a ∉ AP`):

| `a→2a` | 7→14 | 8→16 | 9→18 | 10→20 | 11→22 | **12→24** | 13→26 |
|---|---|---|---|---|---|---|---|
| `M` | 1/11 | 2/23 | 2/23 | 2/27 | 2/25 | **1/14** | 2/27 |
| vs `1/14` | > | > | > | > | > | **= TIGHT** | > |

**Every doubling has `M ≥ 1/14` — none is a counterexample — and only `12→24`
lands exactly on the wall.** The others strictly *loosen*: removing their element
opens a better lonely window, pushing `M` above `1/14`. So `12` is the unique
element whose doubling preserves tightness. (`12` is a non-unit mod `2n=28`, the
"non-unit-pair hole" of S553b; sporadics exist at all because `2n−1 = 27 = 3³`
is composite, S592.)

## Why `V*` is tight: it inherits the AP's wall exactly [verified]

The AP's tight witnesses are `t* = j/14` for `j ∈ {1,3,5,9,11,13}`. At each, the
**binding runners are an antipodal pair `{a, 14−a}`** — `{1,13}, {3,11}, {5,9}`,
each summing to `n = 14` (the pair-sum / `2n−1` structure of THM-401). `V*` has
the **same** witnesses and the **same** binding pairs; the new runner `24` is
**never binding** — it has slack at every witness and opens no sub-`1/14` pinch.
So `V*` sits on the AP's wall without disturbing it: doubling `12` to `24` keeps
`24` clear at exactly the times where the antipodal pairs pinch to `1/14`.

This is the precise mechanism behind the lone sporadic: the wall is held up by
the antipodal binding pairs, and `12` is the one element you can double while
keeping the doubled value clear of all those pinch times.

## Honest status

Enumeration + structure, not a proof of LRC@14. New content: a fast independent
census (zero counterexamples to speed 20), the **doubling law** (only `12→24` on
the wall, all doublings `≥ 1/14`), and the **antipodal-pair witness mechanism**
explaining why `V*` inherits the AP wall. It also sharpens the proof target: the
wall is an antipodal-pair phenomenon, so closing LRC@14 means controlling the
antipodal binding pairs `{a, n−a}` and the one non-unit doubling hole.

## Next (finite, efficient)

1. **Push the counterexample sweep** to speeds `≤ 22–24` (a few minutes) to
   include `V*`'s neighborhood and confirm `{AP, V*}` is the complete tight
   family in the extended range.
2. **Generalize the doubling law**: for even `n` with `2n−1` composite, is the
   unique on-wall doubling always `(n−2) → 2(n−2)`? Test `n = 8, 10` (where
   sporadics exist) for the analogous single element.
3. **Antipodal-pair certificate**: express `M(V) = 1/n` as "every antipodal pair
   `{a, n−a}` pinches at some `t* = j/n`" and check whether a sign-reversing
   involution on the *pairs* (not subsets) certifies it — the pair-level twist.
