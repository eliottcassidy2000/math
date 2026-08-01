---
source: kind-pasteur-2026-07-24-S159 (Opus 4.8)
status: RESULT (upgrades kps-S158). The FC(3) cyclic-weight leak map is defined over Z (integer moment table
  T(A,P,Q)=A!P!Q!*g, free parameters, NO omega), so its Jacobian rank can be computed EXACTLY modulo a prime --
  fast, no precision loss. This carries the transversality/capacity scan to full J>=P for D=4 and D=5: rank=P at
  D=2,3,4,5 (all). Capacity=P CONFIRMED through D=5; the cyclic-weight route to an FC(3) counterexample is closed
  through D=5. Strong evidence FC(3) holds in the cyclic-weight family (contra "may be false").
tags: [factorial-conjecture, S3-symmetry, invariant-theory, transversality, modular-arithmetic, capacity]
related: [kps-S156, kps-S157, kps-S158]
supersedes: [kps-S158 partial D=4,5]
---

# FC(3): modular rank confirms transversality through D=5

## 1. The unlock: the leak map is over Z, so rank mod p is exact
kps-S158 built the invariant reduction `L(s^A L1^P L1bar^Q)=A!P!Q!*g(A,P,Q)` (integer `g`-recurrence) but was
blocked at D>=4 by the factorial-scale cancellation (`~degree!`) which forces huge mpmath precision. Key
realization: in `(A,P,Q)` (DFT) coordinates the entire leak map is **defined over `Z`** -- the moment table `T`
is an integer, the parameters `c_i` are free, and **no `omega` appears** (the `omega`-weight is just the residue
rule `2P+Q mod 3`). So the generic leak-Jacobian rank equals its rank **modulo a prime** (for good `p`), which is
computed with fast exact modular arithmetic -- no precision, no giant integers. (Used `p=2^31-1`.)

## 2. Transversality scan -- COMPLETE for D<=5
Leak-Jacobian `J_{ji}=[c_i] d L(f^{3j})/(3j)! = L(f^{3j-1}M_i)`, random `c in F_p`, rank by Gaussian elimination:
| `D` | params `P` | leaks `J` | rank mod `p` | verdict |
|---|---|---|---|---|
| 2 | 2 | 3 | **2 = P** | transversal, capacity=P |
| 3 | 5 | 6 | **5 = P** | transversal, capacity=P |
| 4 | 10 | 12 | **10 = P** | **transversal, capacity=P CONFIRMED (J>=P)** |
| 5 | 17 | 19 | **17 = P** | **transversal, capacity=P CONFIRMED (J>=P)** |

> **Capacity = `P` at every degree `D<=5`.** `J>=P` for all four, so this is the full test (not the partial
> J<P of kps-S158). The `P` moment-leaks are independent; the `(P+1)`-th is generically forced nonzero. **The
> cyclic-weight tower terminates through D=5 -- no finite cyclic-weight `f` closes it.**

## 3. Reading (honest, unchanged shape)
- **The natural roots-of-unity / cyclic-weight route to an FC(3) counterexample is closed through D=5.** This is
  evidence FC(3) is **true** in that family, cutting against "FC(3) may be false" -- the obvious construction
  fails, transversally, at five consecutive degrees.
- **Caveats (still binding):** (i) *cyclic-weight ansatz* -- a counterexample could live outside it (but the
  roots-of-unity structure makes it the natural home; a non-cyclic-weight `f` would have to kill the `3 !| k`
  moments too, strictly harder). (ii) Generic full rank rules out a *family* and pins capacity, but not an
  *isolated special-point* counterexample; such a point needs a **non-character (Kontsevich-Zagier
  period-rigidity) coincidence** beyond the Conway-Jones vanishings -- the irreducible hard core (kps-S154). (iii)
  Confirmed to D=5; the method (modular, over `Z`) scales further with more time (`f^{3j}` dict still the cost).

## 4. What this pins down
The FC(3) question, restricted to the natural family, is now sharp: **it is exactly whether the moment-leak map
is EVERYWHERE transversal or has a special non-transversal point.** Five degrees of everywhere-full-rank is real
evidence for the former (FC(3) true in-family). A genuine counterexample must be an isolated KZ-rigidity failure
-- the same period-theoretic object the whole factorial/Jacobian bridge (kps-S154) points at. Numerics can raise
the degree; only period theory can settle the isolated point.

Files: `/tmp/{invred,modrank}.py`. Upgrades kps-S158 (D=4,5 now full, via modular rank over Z).
