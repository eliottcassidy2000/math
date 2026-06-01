# The Lonely-Runner gap is a balanced apex-pair problem; the tight locus is pairs summing to a multiple of n (S521)

*claudebox-2026-06-01-S521, long creative session. A new structural theorem for the
Lonely-Runner gap function, with proof, plus its corollaries: the tight locus, the
pairwise-sum moduli, and the apex/fiber-only connection. Candidate canon theorems
(assign THM ids at merge). Builds on THM-384/387, the apex reflection, and the
marked-tournament reduction.*

## The gap function

For distinct positive integers `v_1,...,v_m` (n = m+1 runners, threshold 1/n),
`M(v) = sup_t min_i ||v_i t||`.  LRC(n) is the statement `M(v) >= 1/n`.

## Theorem A (balanced apex-pair formula) — PROVED

> Unless all `v_i` are odd (then `M = 1/2`, attained at `t = 1/2`, every runner
> antipodal to the observer), `M(v)` is attained at a time `t* = k/(v_i+v_j)` for
> some pair `i,j` and integer `k`, at which `v_i, v_j` lie at equal distance
> `M` from the observer **on opposite sides**.  Equivalently
> ```
> M(v) = max over pairs (i<j), k in [1, v_i+v_j):
>          { d = ||k v_i/(v_i+v_j)|| = ||k v_j/(v_i+v_j)|| :
>            ||k w/(v_i+v_j)|| >= d for all other speeds w }.
> ```

**Proof.**  `D(t) = min_i ||v_i t||` is continuous, 1-periodic, piecewise linear
(a lower envelope of tents `||v_i t||` of slope `±v_i`), so its maximum is
attained.  At a maximiser `t*` let `A` be the active set (runners achieving the
min).  If `|A| = 1`, say runner `i`, then near `t*` `D = ||v_i t||`, whose only
interior maximum is a tent peak `||v_i t*|| = 1/2`; but then `D(t*) = 1/2` forces
every other runner `>= 1/2`, i.e. all active — so either `|A| >= 2` or every
runner is at distance `1/2` (all `v_l t*` half-integers), the all-odd case
`t* = 1/2`.  Otherwise take two active runners `i,j`.  For `t*` to be a maximum
the two tents must cross with **opposite slopes** (else `D` increases on one
side).  Write the increasing tent as `v_i t* = d in (0,1/2)` (`||·|| = d`,
slope `+`) and the decreasing tent as `v_j t* = 1-d in (1/2,1)` (`||·|| = 1 -
v_j t* = d`, slope `-`).  Adding: `(v_i+v_j) t* = 1`, so `t* = k/(v_i+v_j)` (the
`k`-th such time), and `d = ||k v_i/(v_i+v_j)||`.  The two runners sit at `d` and
`1-d`: equal distance, opposite sides.  All other runners are `>= d` by
definition of the min.  ∎

The binding pair `(i,j)` are the observer's two **circular neighbours** at `t*` —
the source/sink of the seam, i.e. the **apex tile** of the staircase.  So all of
`M` lives at the apex (fiber-only LRC), confirming the structural role of the
single source–sink arc.

Verified: `M_balanced = M_true` with 0 mismatches over 1326 random sets;
`t*` denominator divides `v_i+v_j` of a binding pair, 0 violations over 2651 sets.

## Theorem B (tight locus) — PROVED

> If `M(v) = 1/n` then some pair of speeds sums to a **multiple of n**:
> `n | (v_i + v_j)` for the binding pair.  Equivalently, every exactly-tight
> (extremal) speed set has a pair of speeds summing to `0 (mod n)`.

**Proof.**  By Theorem A, `M = ||k v_i/q||` with `q = v_i+v_j`.  Writing
`r = k v_i mod q`, `M = min(r, q-r)/q`.  `M = 1/n` gives `min(r,q-r) = q/n`, an
integer, so `n | q = v_i + v_j`.  ∎

Verified (n = 5,6,7): every tight set has such a pair; sets with **no** pair
`≡ 0 (mod n)` have `M` strictly above `1/n` (min `M = 3/13, 2/11, 2/11`).

## Corollary / reduction (heuristic, partially proved)

The extremal locus of LRC(n) — the only place `M` can equal `1/n` — is the family
of speed sets having a pair summing to a multiple of `n`.  For `n = 14` these are
the sets with two speeds summing to `14, 28, 42, ...`.  The arithmetic-progression
extremiser `{1,...,n-1}` has many such pairs (`j + (n-j) = n`).

**Honest caveat.**  Theorem B characterises *exact* tightness `M = 1/n`.  It does
**not** by itself rule out a sub-tight counterexample `M < 1/n` for a no-pair set
(the M-formula does not force `n | q` when `M < 1/n`).  Empirically no-pair sets
have `M > 1/n` with room, so the *conjecture* "no pair `≡ 0 (mod n)` ⇒ `M > 1/n`"
is strongly supported and would reduce LRC(n) to the pair-summing class — but a
quantitative lower bound on `M` for no-pair sets is still needed to prove it.

## Methodology refinement: the canonical moduli are pairwise sums

Theorem A says the **deepest** lonely time is `t = k/(v_i+v_j)`.  So the canonical
moduli for the residue/multiplicative-walk attack are the **pairwise sums**
`q = v_i + v_j` (not arbitrary `q`).  At such `q`, winning means a multiplier `k`
with `k·{speeds} mod q` avoiding the forbidden band `F_q`, with the binding pair
balanced (`k v_i ≡ -k v_j`).  This sharpens the search from "all denominators" to
the `<= C(m,2)` pairwise sums, and ties the multiplicative walk to the additive
structure of the speed set (the sums), which is exactly the addition/multiplication
meeting point the program keeps circling.

## Connections found this session

- **Three-distance theorem** (agent computation): for AP speeds the gap multiset
  of `{0, t, 2t, ..., (n-1)t}` takes at most 3 values (Steinhaus); for general
  speeds up to `n` distinct gaps.  The AP tightness is the 3-distance regular case
  with binding pair `(1, n-1)` summing to `n`.
- **Apex / fiber-only LRC**: Theorem A localises the whole gap function to the
  observer's two neighbours = the apex seam, confirming the marked-tournament
  reduction (runner sub-tournament observer-blind; LRC content at the apex).

## Seeds

a. **Prove the no-pair reduction.**  Show `M(v) > 1/n` whenever no pair sums to
   `0 (mod n)`, by lower-bounding the balanced distance at a well-chosen pair-sum
   modulus.  This would reduce LRC(n) to the pair-summing class — a genuine
   shrinkage.
b. **Tight-locus structure.**  Among pair-summing sets, characterise which are
   actually tight (`M = 1/n`) vs. which still have slack — Theorem B gives the
   necessary condition `n | v_i+v_j`; find the sufficient mod-n pattern (the
   regular-polygon residue structure) and count the tight sets.
c. **n=14 focus.**  Restrict the whole machine to sets with a pair summing to 14
   (the tight locus); the binding pair sits at `t = k/14` (the regular polygon),
   and the residual is the mod-14 balance of the other 11 speeds — exactly the
   fully-covered / apex residual from earlier S521 reflections, now reached from
   the gap-function side.
