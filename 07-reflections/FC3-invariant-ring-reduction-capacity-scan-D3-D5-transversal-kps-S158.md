---
source: kind-pasteur-2026-07-24-S158 (Opus 4.8)
status: RESULT (invariant-ring reduction + capacity scan). Implements the S3-invariant reduction (a single
  integer recurrence for the moment table L(s^A L1^P L1bar^Q)=A!P!Q!*g(A,P,Q)), collapsing the degree-9j
  three-variable expansion of L(f^{3j}) to dict-powering in (A,P,Q). Pushes the FC(3) cyclic-weight capacity scan
  to D=3 (FULLY confirmed transversal, capacity=P=5) and D=4,5 (partial: first 6 leaks independent). NO
  rank-deficiency / overshoot at any degree -> the transversality conjecture (kps-S157) holds where testable, so
  the cyclic-weight route to an FC(3) counterexample is closed; consistent with FC(3) true in this family.
tags: [factorial-conjecture, S3-symmetry, invariant-theory, transversality, capacity, exact-arithmetic]
related: [kps-S154, kps-S156, kps-S157]
---

# FC(3): the invariant-ring reduction, and a capacity scan to D=3-5

## 1. The reduction (the computational unlock)
`L(f^{3j})` was blocked by expanding `f^{3j}` in `X1,X2,X3` (degree `9j`). Work instead in `(s,L1,L1bar)` DFT
coordinates (`s=sum X_i`, `L1=sum omega^i X_i`), where the factorial functional has an integer moment table:
> `L(s^A L1^P L1bar^Q) = A! P! Q! * g(A,P,Q)`, `g(a,b,c)=[u^a v^b w^c] prod_j 1/(1-l_j)`,
> and from `(1-l1)(1-l2)(1-l3) = 1-3u+3(u^2-vw)-(u^3+v^3+w^3-3uvw)` (using `prod(u+omega^{j}v+omega^{-j}w)
> = u^3+v^3+w^3-3uvw`), `g` obeys the **all-integer recurrence**
> `g(a,b,c)=[a=b=c=0]+3g(a-1,b,c)-3g(a-2,b,c)+3g(a,b-1,c-1)+g(a-3,b,c)+g(a,b-3,c)+g(a,b,c-3)-3g(a-1,b-1,c-1)`.
`g` is a real integer, nonzero only on `2b+c≡0 (mod 3)` (the C3 selection rule). Validated: `T(0,3j,0)=(3j)!`
(the seed theorem, kps-S156) for `j=1,2,3`. Now `L(f^{3j})` = dict-power `f^{3j}` in `(A,P,Q)` + table lookup --
no `9j`-degree 3-variable object.

## 2. The diagnostic: leak-Jacobian rank = transversality
kps-S157 recast FC(3)-in-the-cyclic-weight-family as: are the moment-leaks `f -> (L(f^3),L(f^6),...)`
**transversal** (capacity = #params `P`)? Equivalent local test: the leak Jacobian
`J_{ji} = d/dc_i [L(f^{3j})/(3j)!] = L(f^{3j-1} M_i)/(3j-1)!` at a random nonzero point should have **full rank
`P`**. Full rank `P` => the leak map is an immersion => solutions of `j=1..P` are isolated and the next leak is
generically nonzero (tower terminates). Rank **deficiency** `< min(J,P)` would be the overshoot / FC(3)-false
signal (a positive-dimensional family of near-counterexamples). Computed at high precision (`mpmath dps=90`, to
resolve the `~degree!` factorial cancellation).

## 3. Results
| `D` | params `P` | leaks tested `J` | Jacobian rank | verdict |
|---|---|---|---|---|
| 2 | 2 | 3 | **2 = P** | transversal, capacity 2 CONFIRMED (matches exact kps-S157) |
| 3 | 5 | 6 | **5 = P** | **transversal, capacity 5 CONFIRMED, tower TERMINATES** |
| 4 | 10 | 6 | 6 = min(J,P) | first 6 leaks independent (partial; `J<P`, `f^{29}` infeasible) |
| 5 | 17 | 6 | 6 = min(J,P) | first 6 leaks independent (partial) |

> **No rank-deficiency at any degree.** D=2 and D=3 fully confirm transversality (capacity = `P`); D=4,5 show the
> first 6 leaks independent (no overshoot in reach). The cyclic-weight family's leaks behave transversally
> everywhere tested.

## 4. Reading (honest)
- **Within the cyclic-weight (C3-eigenvector) ansatz, the tower terminates**: `P` parameters kill `P` leaks and the
  next is generically forced nonzero. So **the natural roots-of-unity route to an FC(3) counterexample is closed**
  through D=3 (fully) and D=4,5 (first 6 leaks). This is evidence FC(3) is **true** in this family, contra the
  "FC(3) may be false" guess -- at least, the obvious construction does not work.
- **Caveats (do not overclaim):** (i) this is the *cyclic-weight* ansatz -- a counterexample could live outside
  it (though the roots-of-unity structure makes it the natural home). (ii) Full rank `P` rules out a *family* of
  counterexamples and pins capacity generically, but does **not** exclude an *isolated* special-point
  counterexample (rank is a generic/first-order invariant). Such a point would need a **non-character (KZ
  period-rigidity) coincidence** beyond the Conway-Jones vanishings of kps-S157 -- exactly the hard core.
  (iii) D=4,5 are only partial (`J<P`).

## 5. Status / next
Delivered: the S3-invariant reduction (integer `g`-recurrence) and a transversality scan confirming capacity=`P`
for D=2,3 and no overshoot for D=4,5. The reduction makes higher `j`/`D` reachable but the `f^{3j}` dict still
grows; the clean next step is to push D=4,5 to `J>=P` by computing the leak map **in the two reflection
invariants `theta_2,theta_3` directly** (a genuinely 2-variable object) rather than `(A,P,Q)`. If transversality
persists, it becomes a real conjecture that FC(3) holds on the cyclic-weight family; a rank drop would be a
concrete counterexample lead.

Files: `/tmp/{invred,rank}.py`. Builds on kps-S156/S157.
