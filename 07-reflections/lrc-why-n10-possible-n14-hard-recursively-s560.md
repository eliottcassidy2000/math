---
source: opus-2026-06-02-S560 (remote-control)
status: SYNTHESIS — recursive account of why n=10 is proven but n=14 is the wall; refines S559 (apex is shared, NOT the distinguisher)
tags: [LRC, n10, n14, recursion, finite-checking, scale, apex, doubled-prime, S558, S559]
---

# Why n=10 is possible but n=14 is hard — recursively

**Prompt (user):** understand recursively what makes n=10 possible but n=14 hard.

**The trap this dissolves.** n=10 = 2·5 and n=14 = 2·7 are *both* doubled primes,
so both have the apex zero-divisor `q` and both fail the polynomial-method shortcut
(S559). So the apex / `2q` structure — the thing S559 spotlighted — **cannot** be
what separates them. The separation is elsewhere, and it is almost entirely
**scale**.

## The recursion of the proof, layer by layer

LRC(n), `n=k+1`, is attacked by the finite-checking pipeline. Its difficulty
decomposes recursively; each layer either bottoms out or passes the cost down:

```
LRC(k)
 └─ (Tao) reduce to integer speeds; (Malikiosis–Santos–Schymura) bound the
    speed product:  ∏ uᵢ < B_k = (C(k+1,2)^{k-1}/k)^k          ← the "haystack"
     └─ rule out the haystack via the DIVISIBILITY SIEVE: for primes p, show
        every residue tuple mod p is eventually proper
         └─ per prime, LIFT by multipliers c = 2,3,…,k+1
             └─ the c=(k+1) lift is the deepest/most expensive layer; it is
                where the TIGHT tuple (1,…,k) and its relatives live
                 └─ BASE CASE for this layer:
                     • if k+1 is an ODD PRIME → polynomial method settles the
                       tight tuple with NO computation (Prop 4.1; works because
                       ℤ_{k+1} is a field — no speed ≡0, so no apex obstruction)
                     • if k+1 = 2q → the apex (speed q) is a zero-divisor; the
                       polynomial method FAILS (S559); the c=(k+1) lift must be
                       computed in full
```

So the recursion bottoms out at the **c=(k+1) lift of the tight tuple**, and there
are exactly **two ways to survive it**:

- **(A) Compute it directly** — feasible only while `k` is small.
- **(B) Skip it via the polynomial method** — available only when `k+1` is an odd
  prime.

## The haystack explodes (the only thing that actually changes)

`log₁₀ B_k`, and its per-step growth, for k=7..13:

| k | n | k+1 | log₁₀ B_k | Δ (per step) | route used |
|---|---|---|---|---|---|
| 7 | 8 | 2³ | 55 | — | (A) small |
| 8 | 9 | 3² | 80 | +25 | (A) |
| 9 | **10** | **2·5** | **110** | +30 | **(A) small — DONE** |
| 10 | 11 | 11 prime | 147 | +36 | (B) prime |
| 11 | 12 | 2²·3 | 189 | +42 | (A) still feasible |
| 12 | 13 | 13 prime | 237 | +48 | (B) prime rescues large lift |
| 13 | **14** | **2·7** | **291** | +54 | **NEITHER → WALL** |

The haystack grows from `10¹¹⁰` (n=10) to `10²⁹¹` (n=14), and the *increment*
itself grows (+30 → +54 per step). Route (A) silently ran out somewhere around
k=12–13.

## The crisp answer

- **n=10 (k=9, 2·5):** route (B) is *blocked* (apex zero-divisor 5, exactly as at
  n=14) — but `k=9` is small, so route (A) closes it. **Proven.**
- **n=11, n=13 (k+1 = 11, 13 prime):** route (B) is *open*, so even though they are
  larger than n=10 they are settled by the polynomial shortcut. The primes bought
  room and **masked the fact that route (A) was approaching its ceiling.**
- **n=12 (k=11, 2²·3):** route (B) blocked, but `k=11` is *still just small enough*
  for route (A).
- **n=14 (k=13, 2·7):** the **first** doubled-prime case where route (A) has become
  infeasible (`B₁₃ ≈ 10²⁹¹`, the c=14 lift ≫ the ~40-CPU-days that k=12 cost) **and**
  route (B) is unavailable (14 not an odd prime → apex zero-divisor 7 blocks
  Prop 4.1). It falls in the gap left when the prime sequence (11, 13) steps over to
  the next composite (14).

> **n=14 is not hard because of any new structure. It is the first
> `2q` (route-B-blocked) case large enough that the direct finite check
> (route A) finally fails — with no neighbouring prime to rescue it.**

## The counterintuitive refinement (n=14 is structurally *easier*)

Under the S559 reduction, the per-tuple residual rate (apex excluded, parity-
matched — the hardest class) is

```
q = 3   5   7   11    (n = 6  10  14  22)
res= .62 .41 .28 .14   ← DECREASING
```

So a *single* hard tuple at n=14 (q=7) is **easier** to correct than at n=10
(q=5): more units in ℤ_q^× means more freedom to dodge the forbidden residues,
and the apex repair works for every q. **The structural difficulty per tuple goes
down as q grows.** Everything that makes n=14 hard is the *number* of tuples
(the `B_k` haystack), not their individual hardness. n=10 and n=14 are
structurally isomorphic (apex, mod-2-forced + mod-q-free reduction, even-fold to a
proven base); only the scale differs.

## Recursive corollary for the repo's program

This relocates the target precisely. Since the per-tuple problem *eases* with q
while the count explodes, the productive attacks on n=14 are the two that fight
scale, not structure:

1. **Restore route (B) for `k+1 = 2q`** — the S559 program: an algebraic
   (polynomial-method) substitute that handles the tight tuple and clears the
   ratio-spread residual *without* the c=14 enumeration. The apex is already
   handled; the residual is a single `ℤ_q^×` ratio-cover condition that *shrinks*
   with q — encouraging.
2. **A structural (non-computational) proof** via the pinch/shield route
   (THM-396, HYP-2061) or even-fold, which is `k`-uniform and so indifferent to
   the haystack growth.

Both are scale-beating; the brute c=14 lift is not. The even-fold also bottoms out
recursively the same way for both n (even half of n=14 → proven LRC(13); of n=10 →
proven LRC(9)), confirming the only real difference is the size of the residual
check.

**Artifacts:** `04-computation/lrc_n10_vs_n14_cost_s560` (inline) +
`05-knowledge/results/lrc_n10_vs_n14_cost_s560.out`. Builds on S558 (method map),
S559 (apex/2q reduction), arXiv:2604.23906. No new HYP (synthesis; refines
HYP-2063).
