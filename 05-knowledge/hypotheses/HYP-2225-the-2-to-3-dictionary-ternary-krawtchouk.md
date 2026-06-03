# HYP-2225 — The 2→3 dictionary: ternary Krawtchouk, signed depth, and the 3-structure of n=14

**Session:** claudebox-2026-06-03-S622. **Frame:** the Delsarte-LP consolidation (HYP-2215); the user: "look for
more concepts like the diagonal Krawtchouk-positive certificate but instead of being 2, they're 3."
**Threads:** HYP-2210 (binary Krawtchouk), HYP-2205 (flow-shell crosses), HYP-2200 (Helly-3), HYP-2175 (Collatz ×3).

## The whole toolkit lives at q=2 — here is its q=3 twin
Our depth weight enumerator is the BINARY (q=2, forbidden/safe) Hamming scheme, and the diagonal certificate uses
the binary Krawtchouk. The "instead of 2, they're 3" request = the ternary (q=3) parallel toolkit.

### The q-ary Krawtchouk (verified + formalized)
`Kq(q,k,n,x) = ∑_j (−1)^j (q−1)^{k−j} C(x,j) C(n−x,k−j)`. `q=2` recovers binary `K`; `q=3` is the **ternary**
Krawtchouk (factor `2^{k−j}`). Generating function `∑_k Kq z^k = (1−z)^x (1+(q−1)z)^{n−x}` (q=3: `(1−z)^x(1+2z)^{n−x}`,
verified). Baseline `Kq(q,k,n,0) = (q−1)^k C(n,k)` (q=3: `2^k C(n,k)`).
**Formalized (Krawtchouk/Basic.lean):** `Kq`, `Kq_two_eq_K`, `Kq_zero_index`, `Kq_at_zero`.

### The signed (ternary) depth refinement
Each runner near the origin is in one of THREE states, not two: `+` (in `(0,δ)`), `−` (in `(−δ,0)`), or `0` (safe).
The signed depth `(n₊, n₋, n₀)` refines the binary depth `n₊+n₋`; lonely ⟺ `(n₊,n₋)=(0,0)`. Since the `+` and `−`
arcs are disjoint, `(1−x⁺)(1−x⁻)=1−x⁺−x⁻=1−x`, so `p₀` is unchanged — but the ternary weight enumerator carries
MORE moments (n₊, n₋ separately) and the σ-symmetry `n₊ ↔ n₋` (HYP-2205), so the **ternary Delsarte LP has strictly
more constraints than the binary one** — a route to a feasible dual where the binary diagonal was vacuous.

## The 2→3 dictionary (more "3" concepts)
| q=2 | q=3 |
|---|---|
| binary Krawtchouk `K_k` | ternary Krawtchouk `Kq(3,·)` (factor `2^{k−j}`) |
| depth: forbidden/safe (2 states) | signed depth: `+ / 0 / −` (3 states) |
| Helly number 2 (intervals on a LINE) | Helly number 3 (arcs on a CIRCLE) — HYP-2200 |
| pairwise overlap `S₂` | triple overlap `S₃` = 3-term coincidence `a+b=c` — HYP-2195 |
| 2-adic doubling `÷2` | 3-adic / Collatz `×3` — HYP-2175 |
| antipodal `⟨−1⟩`, order 2 (the σ-pair / cross) | doubling `⟨2⟩` mod 7, order 3 (the 3-cycle `{v,2v,4v}`) |
| apex / 2-block at `n=2q` | 3-block at `n=3q` (ternary apex) |

## The 3-structure of n=14 (verified)
`n=14 = 2·7`, and **2 has multiplicative order 3 mod 7**, so the doubling map `⟨2⟩` has orbits the two 3-cycles
`{1,2,4}` and `{3,6,5}` on `(ℤ/7)*`. The antipodal σ-pairs `{1,6},{2,5},{3,4}` (order 2) sit *inside* these 3-cycles.
So n=14 carries BOTH a 2-structure (σ, the cross arms) and a 3-structure (`⟨2⟩`, the doubling triples) — the
2-adic/3-adic seam in one config. The ternary refinement is adapted to the 3-cycles; the σ-cross to the pairs.

## Why this could help LRC(14)
The binary Delsarte LP's diagonal duals are vacuous at the gap (HYP-2215). The ternary refinement keeps `n₊, n₋`
separate (the binary depth conflates them), adding the σ-symmetry and the per-arm constraints — a strictly larger
constraint set. The off-diagonal certificate (HYP-2215's target) may be a TERNARY Krawtchouk-positive dual,
adapted to the `⟨2⟩` 3-cycles, that the binary scheme cannot express.

## Formalized (math-lean, sorry-free) — `Math/Krawtchouk/Basic.lean`
`Kq` (q-ary Krawtchouk), `Kq_two_eq_K` (q=2 → binary), `Kq_zero_index`, `Kq_at_zero` (the q-ary baseline).

## Open
- The ternary signed-depth weight enumerator and its `Kq(3,·)` transform; the ternary Delsarte LP and a feasible
  off-diagonal dual for LRC(14).
- The q-ary genfun `∑_k Kq z^k = (1−z)^x(1+(q−1)z)^{n−x}` (formalize).
- The 3-block at `n=3q` (ternary apex) — does the apex-sheaf story have a q=3 analogue?
