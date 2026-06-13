# THM-407 — Twisted-involution (scaling × time-reversal) reduction of the LRC additive residual; n=14 collapses 13 shells → 3

**Status:** PROVED (Lemma A elementary; the reduction Corollary follows; the n=14 stratification
verified exactly). **Source:** opus-2026-06-03-S599g (remote-control). Realises Transfer T2/T8
(HYP-2160) as a concrete n=14 improvement, via the laminar-flow / criss-cross picture.
**Setup:** speeds `S={v_1,…,v_{n−1}}⊂ℤ_{>0}`, gap `1/n`, `M(S)=max_t min_i ‖v_i t‖`; LRC(n) ⟺
`M(S)≥1/n` for all `S`. Pair-sum modulus `2n−1` (THM-401); **shells** = the pairs `{a,−a}` in
`ℤ/(2n−1)`, `a=1,…,n−1`.

---

## The flow picture (where this lives)

Draw the trajectories `y=v_i t (mod 1)`: a **criss-cross lattice** of lines of slope `v_i`, with
crossings at `t=m/(v_i±v_j)` (the arrangement vertices = the candidate lonely times, S599d). The
danger bands `‖v_i t‖<1/n` are horizontal strips; a **lonely time** is a vertical line threading
the **laminar channels** between all strips. Two maps preserve this entire picture's tightness:

## Lemma A (scaling & time-reversal invariance of `M`) — PROVED

> For every integer `c≥1` and every sign, `M(cS) = M(S)` and `M(−S) = M(S)`.

**Proof.** `M(−S)=M(S)` since `‖(−v)t‖=‖vt‖`. For scaling: as `t` runs over `[0,1)`, `u:=ct mod 1`
runs over all of `[0,1)` (`c` times around), and `min_i‖(cv_i)t‖ = min_i‖v_i u‖`; taking the sup,
`M(cS)=sup_u min_i‖v_i u‖=M(S)`. ∎ **Verified** (`lrc_twisted_involution_flow_shells_s599g.py`):
`M=M(2S)=M(3S)=1/6` for the AP and the sporadic `(1,3,4,5,9)`; `2/5` for a loose config.

Geometrically: **time-reversal `t↦−t`** is an honest involution of the laminar flow (fixing
`t=0,½`); **doubling `S↦2S`** is the *twist* (the renormalisation that rescales the criss-cross
lattice by 2). Together they generate the tightness-preserving group.

## Corollary (twisted-involution shell reduction) — PROVED

> On the cover-relevant data — the residues mod `2n−1` (THM-401/S574) — scaling-by-2 and
> negation act as the multiplicative group `G = ⟨2,−1⟩ ≤ (ℤ/(2n−1))^×`. By Lemma A every
> `G`-orbit of shells is a single tightness class, so the worry-set's **additive residual is
> determined by one shell per `G`-orbit.** The number of cases is `#{G\text{-orbits on } ℤ/(2n−1)\setminus\{0\}}`,
> not `n−1`.

The reflections `(−1)·2^k ∈ G` are the **twisted involutions on the flow shells** (glide-reflections
of the criss-cross lattice); the orbits are their necklaces.

## The n=14 collapse (verified)

`2n−1 = 27 = 3³`, and `2` is a primitive root mod `27` (`ord_{27}(2)=18=φ(27)`). The `G`-orbits:

```
 orbit (size 18):  units {1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26}  → 9 shells   gcd=1
 orbit (size  6):  {3,6,12,15,21,24}                                        → 3 shells   gcd=3
 orbit (size  2):  {9,18}                                                   → 1 shell    gcd=9
```

> **The 13 shells of the n=14 additive residual collapse to 3 `G`-orbit representatives,
> stratified by `gcd(a,27)∈{1,3,9}` = the prime-3 tower `3⁰,3¹,3²`.** Verifying the worry-set on
> the additive face needs **3 representative shells**, one per `3`-adic stratum — precisely the
> sporadic prime-3 decomposition of S592, now derived as the orbit space of the twisted-involution
> group.

## Why n=14 is the frontier (structural reading)

Across even `n` (verified `n=6..20`):

```
 n        :  6   8   10  12  14  16  18  20
 2n−1     : 11  15   19  23  27  31  35  39
 #orbits  :  1   3    1   1   3   3   3   3
```

When `2n−1` is **prime with 2 primitive** (`n=6,10,12`: `11,19,23`), `G` is transitive on units and
the additive face is a **single** case — essentially closed. **`n=14` is the first even `n` whose
`2n−1` is an odd prime *power* (`3³`)**, so doubling can no longer be shell-transitive: the `3`-adic
strata `{3⁰,3¹,3²}` are `G`-invariant and irreducible. The arithmetic obstruction that makes
`n=14` the open frontier is exactly this failure of the twisted-involution group to act
transitively — the prime-3 fragmentation (companion to the prime-2 fragmentation of THM-404).

## Consequence / improvement

- The LRC(14) additive (`2n−1`) residual is reduced to **3 shell representatives** (gcd `1,3,9`).
  Combined with THM-404 (prime-2) and the solved prime-7 (`ℚ(ζ_14)=ℚ(ζ_7)`), the n=14 worry-set
  check is now: prime-2 doubling face × **3** prime-3 shell-strata.
- The `gcd=9` stratum is a *single* shell `{9,18}` (orbit size 2) — the smallest, most rigid core;
  the `gcd=3` stratum is 3 shells; the unit stratum is 9 shells but one orbit. The hard core is the
  `gcd∈{3,9}` (non-unit) strata — the multiples of 3, i.e. the `3`-divisible resonances.

## Honest status

- **PROVED:** Lemma A (`M(cS)=M(−S)=M(S)`); the `G=⟨2,−1⟩` shell reduction; the n=14 orbit
  structure `13→3` stratified by `gcd∈{1,3,9}`; the even-`n` orbit table.
- **Verified:** `M`-invariance on AP/sporadic/loose; the orbit computation `n=6..20`.
- **Open (the residual after reduction):** prove the 3 representative shells each admit the clock
  witness (loneliness) — i.e. close the prime-3 strata; and combine with the prime-2 face for a
  full LRC(14). The reduction makes this a **3-case** check on the additive face, not 13.

**Artifacts:** `04-computation/lrc_twisted_involution_flow_shells_s599g.py` (+`.out`). Builds on
THM-401 (`2n−1` modulus), THM-404 (prime-2 doubling), S592 (prime-3 sporadics), S599d (criss-cross
arrangement vertices), HYP-2160 (Transfer T2 involution). Companion reflection
`lrc-twisted-involutions-laminar-flow-shells-s599.md`. New: **HYP-2162**.
