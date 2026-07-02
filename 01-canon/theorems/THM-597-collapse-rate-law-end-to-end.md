# THM-597: The collapse-rate law, end to end (general argmax formula + machine-checkable certificates)

**Status:** PROVED (elementary; exhaustively verified on the known tight locus q = 5, 6, 8, 14)
**Author:** mac-mini-2026-07-01-S97 (HYP-3853)
**Verification:** `04-computation/lrc_collapse_law_endtoend_macmini_S97.py` (+ `.out`): formula = direct last-cell slope on all 9 known tight sets including both beaters; all witness denominators divisible by q; all sets clean.
**Role:** the rigorous chain the owner asked for — "unit-residue rigidity + gap-sum formula ⟹ the collapse-rate law, machine-checkable end to end." Combines THM-593(A) (Lean-verified) with THM-592's co-area frame and closes klein's HYP-3835 formula as the clean case of a general unconditional law.

## Setting

`S ⊂ ℤ_{>0}` finite, primitive, **q-tight**: `M(S) = max_t min_{v∈S} ||vt|| = 1/q`, with no element `≡ 0 (mod q)`. Write `f_S(t) = min_v ||vt||` and `Λ_S(r) = |{t : f_S(t) ≥ r}|`.

## Statement

**(a) Argmax structure.** The argmax `A = {t : f_S(t) = 1/q}` is finite, and every `t* ∈ A` is a rational whose reduced denominator is a multiple of `q`. The `φ(q)` unit fractions `{a/q : gcd(a,q)=1}` all lie in `A` (THM-593 Corollary A2). Elements of `A` with denominator `qs`, `s ≥ 2`, are **deep witnesses** (none observed on the known tight locus; their absence per set is a finite check).

*Proof.* If `f_S(t*) = 1/q`, a binding `v` has `|v t* − c| = 1/q` for some `c ∈ ℤ`, so with `t* = p/q'` reduced, `|vp − cq'| · q = q'` in ℤ, giving `q | q'`. Finiteness: `f_S` is piecewise linear with finitely many pieces; its level set at the maximum is a finite set of points (no flat pieces at the max: a flat maximal piece would need two binding constraints with slope 0, impossible since slopes are `±v ≠ 0`; a maximal open interval of constancy of `f_S` at value `1/q` cannot occur because on any interval free of breakpoints `f_S` is affine with slope `±v`). ∎

**(b) Local width.** For `r < 1/q` above the last breakpoint `ρ_last` of `Λ_S`, the lonely set is a disjoint union of one interval per `t* ∈ A`, of width
```
(1/q − r) · ( 1/v₊(t*) + 1/v₋(t*) ),
```
where `v₊(t*)` (resp. `v₋`) is the **fastest** runner whose binding constraint decreases to the right (resp. left) of `t*` — i.e. the fastest `v` with `v t* ≡ −1/q` (resp. `+1/q`) `(mod 1)`.

*Proof.* Near `t*`, `f_S = min` over the binding constraints (tents `|v(t−t*) ∓ 1/q|` with slopes `±v`) and the non-binding runners, which stay `> r` in a neighborhood for `r` close to `1/q` (finitely many, continuity). On the right of `t*` the constraint envelope decreases at the rate of the fastest right-decreasing binding runner, so `{f ≥ r}` extends exactly `(1/q − r)/v₊`; similarly left. Components at distinct argmax points are disjoint once `r` exceeds every merge radius — i.e. beyond `ρ_last` (THM-592(i): finitely many breakpoints). Every surviving component shrinks (nested, THM-592) to a point of `A`, and every point of `A` seeds one: the correspondence is exact. ∎

**(c) The law.** On `(ρ_last, 1/q)`:
```
Λ_S(r) = c(S) · (1 − q r),   c(S) = (1/q) · Σ_{t* ∈ A} ( 1/v₊(t*) + 1/v₋(t*) ).
```

**(d) Unit-fraction floor and the clean case.** At each unit fraction `a/q`, the binding sides are the residue classes `±a^{−1} (mod q)` and the binding runners are the fastest members `v_max(±a^{−1})` (present by THM-593(A)). Hence
```
c(S) ≥ (2/q) · Σ_{u ∈ (ℤ/q)ˣ} 1/v_max(u),
```
with equality iff `A` = the unit fractions (**clean**) — recovering THM-593(B) and klein's HYP-3835 two-term form `c(S) = (1/q)Σ_k [1/max(class(+k^{−1})) + 1/max(class(−k^{−1}))]` exactly on clean sets.

**(e) Machine-checkable certificates.** For a given `S`: `A` is exactly the set of binding-grid rationals (denominators `v+w`, `w−v`, `2v`) at which the exact minimum equals `1/q` — a finite rational computation; `v±` are direct residue reads; hence `c(S)` is an exact per-set rational certificate, and the collapse law over the classified tight locus (HYP-3750 + THM-593(A3): drops are non-units) is a finite exact verification. Together with the Lean-verified THM-593(A) (`TournamentH7.LRCUnitResidue`), the chain [tightness ⟹ unit witnesses] → [argmax certificate] → [width formula] → [c(S)] has no analytic gaps: the only non-elementary input is the piecewise-linearity of `Λ_S` (THM-592(i), elementary real analysis).

## Verified instances (all exact)

| q | S | \|A\| | clean | c(S) |
|---|---|-----|-------|------|
| 5 | {1,2,3,4} | 4 | ✓ | 5/6 |
| 5 | {1,3,4,7} | 4 | ✓ | 29/42 (beater) |
| 6 | {1,2,3,4,5} / {1,3,4,5,9} | 2 | ✓ | 2/5 (both) |
| 8 | {1..7} / GW {1,2,3,4,5,7,12} | 4 | ✓ | 44/105 (both) |
| 8 | {1,4,5,6,7,11,13} | 4 | ✓ | 328/1001 (beater) |
| 14 | AP {1..13} / GW {1..11,13,24} | 6 | ✓ | 1666/6435 (both) |

## Remarks

- The deep-witness branch (`s ≥ 2`) is empty on the known locus; THM-593's addendum conjecture ("tight ⟹ no multiple of q") and the "clean" question are *separate* finiteness questions, but neither is needed: the law (c) is unconditional, and (e) certifies each set.
- Via the S96 MT-slice identity, `c(AP_q) = 2q·Σ_{p+q'=q} 1/(pq'(p+q'))`: the certificates of (e) are level-`q` Mordell–Tornheim data.

-> THM-592, THM-593 (+ Lean), THM-594, HYP-3750, HYP-3834/3835 (klein), HYP-3852, HYP-3853, OPEN-Q-108.
