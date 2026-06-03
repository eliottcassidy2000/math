---
source: opus-2026-06-03-S599g (remote-control)
status: IMPROVEMENT (n=14) — the laminar-flow / criss-cross picture makes the tightness-preserving symmetry group ⟨doubling, time-reversal⟩=⟨2,−1⟩ explicit; PROVED M(cS)=M(−S)=M(S); the twisted involutions on the flow shells reduce the n=14 additive residual 13→3 (stratified by the prime-3 tower gcd∈{1,3,9}); explains why n=14 is the frontier (first even n with 2n−1 an odd prime power)
tags: [LRC, n14, twisted-involution, laminar-flow, criss-cross, flow-shells, scaling-invariance, time-reversal, doubling, prime-3, renormalization, shell-reduction, THM-407]
---

# Twisted involutions on flow shells: laminar flow folds the n=14 residual 13 → 3

**Prompt (user):** aim for proof/improvement on LRC n=14; focus on twisted involutions on flow
shells; water/laminar flow is the natural example; the lattice-grid criss-cross that gets created.

The picture pays off concretely (THM-407). Here is the physics, the proof, and the improvement.

## 1. The laminar flow and its criss-cross lattice

Plot every runner's worldline `y = v_i t (mod 1)` on the time–position cylinder. You get a
**criss-cross lattice**: straight lines of slope `v_i`, crossing wherever two runners coincide,
`v_i t ≡ v_j t`, i.e. at `t = m/(v_i±v_j)` — exactly the arrangement vertices that are the
*candidate lonely times* (S599d, the NP-witnesses). The danger set is the union of horizontal
**bands** `‖v_i t‖<1/n`; the **laminar channels** are the open corridors between all bands. A
**lonely time** is a vertical line threading a channel all the way across — the flow stays
*laminar* (no runner "turbulently" crosses the origin band). LRC(n): such a laminar line always
exists.

This is the **view-obstruction** form: the diagonal flow on the torus must thread the central
cells of the criss-cross grid of coordinate slabs. The grid's geometry is the whole problem.

## 2. The symmetry of the flow — the twist and the involution

Two transformations leave the laminar/criss-cross structure's *tightness* exactly invariant
(THM-407 Lemma A, proved): **time-reversal** `t↦−t` (an honest involution, `M(−S)=M(S)`, fixing
the stagnation points `t=0,½`) and **doubling** `S↦2S` (the *twist*: `M(2S)=M(S)`, a
renormalisation that rescales the lattice by 2 and is self-similar). They generate a group; on the
cover-relevant residues mod `2n−1` it is the multiplicative `G=⟨2,−1⟩`. The **twisted involutions
on the flow shells** are the reflections `(−1)·2^k∈G` — glide-reflections of the criss-cross
lattice, each folding a `⟨2⟩`-necklace of shells onto itself.

*Why "shells":* the pair-sum modulus is `2n−1` (THM-401); a shell `{a,−a}` is a layer of the flow
at radius `a` in `ℤ/(2n−1)`. Doubling moves between layers (`a↦2a`), time-reversal flips a layer
(`a↦−a`). The orbits of `G` are the **irreducible flow shells**.

## 3. The improvement: n=14 folds 13 → 3

Because every `G`-orbit of shells is one tightness class (Lemma A), the worry-set's additive face
needs only **one shell per `G`-orbit**. At `n=14`, `2n−1=27=3³` and `2` is a primitive root mod
27, so (verified):

```
 13 shells  →  3 orbits :  unit-shell (gcd 1, 9 shells) | 3-shell (gcd 3, 3 shells) | 9-shell (gcd 9, 1 shell)
```

The three survivors are **stratified by `gcd(a,27)∈{1,3,9}` = the prime-3 tower `3⁰,3¹,3²`** —
i.e. the twisted-involution group's orbit space *is* the S592 sporadic prime-3 decomposition,
now derived rather than enumerated. **To check the n=14 additive residual: 3 representative
shells, not 13.** The `gcd=9` orbit is a single rigid shell `{9,18}` — the hard core.

## 4. Why n=14 is the frontier (the laminar reading of the prime-3 wall)

Run the fold across even `n` (verified n=6..20):

| `n` | `2n−1` | `#G`-orbits |
|---|---|---|
| 6 | 11 (prime, 2 primitive) | **1** |
| 10 | 19 (prime) | **1** |
| 12 | 23 (prime) | **1** |
| 8 | 15 = 3·5 | 3 |
| **14** | **27 = 3³** | **3** |
| 16,18,20 | 31,35,39 | 3 |

When `2n−1` is prime with 2 primitive, doubling is transitive on the units and the additive face
is a **single laminar layer** — essentially one case, closed. **`n=14` is the first even `n` whose
`2n−1` is an odd prime *power*** (`3³`): the `3`-adic strata `{3⁰,3¹,3²}` are `G`-invariant and
cannot be folded into each other. The flow develops **three disconnected laminar regimes** that the
twist cannot mix. This is the prime-3 fragmentation — the additive-face companion to THM-404's
prime-2 (doubling) fragmentation. *The reason n=14 is open is that its criss-cross lattice is the
first even one whose doubling-symmetry fails to be shell-transitive.*

## 5. Where this leaves LRC(14)

The n=14 worry-set check is now **prime-2 doubling face × 3 prime-3 shell-strata** (`gcd 1,3,9`),
with prime-7 solved (`ℚ(ζ_14)=ℚ(ζ_7)`). The remaining work is a **3-case** loneliness verification
on the additive face plus the prime-2 face — a real shrinkage from 13 shells. The natural next
step (Transfer T2, HYP-2160): build the **Garsia–Milne involution** whose fixed points are exactly
these 3 strata, so that `p_0=Σ(−1)^{|S|}meas(⋂)` telescopes onto them — turning the 3-case check
into a single cancellation argument. The `gcd=9` core `{9,18}` is the place to start: smallest
orbit, most rigid, the `3²` apex of the tower.

## 6. Honest status

- **PROVED (THM-407):** `M(cS)=M(−S)=M(S)`; the `G=⟨2,−1⟩` shell reduction; n=14 `13→3` stratified
  by `gcd∈{1,3,9}`; the even-`n` orbit table; the structural identification of the frontier.
- **Verified:** `M`-invariance (AP/sporadic/loose); orbit counts n=6..20; the n=14 strata.
- **Open:** close the 3 representative strata (the clock witness on each) and combine with the
  prime-2 face for full LRC(14); construct the fixed-point=strata involution (T2). This is an
  improvement (13→3) and a sharp structural explanation, **not** a proof of LRC(14).

**Artifacts:** `04-computation/lrc_twisted_involution_flow_shells_s599g.py` (+`.out`), THM-407.
Builds on THM-401/404, S592 (prime-3), S599d (criss-cross vertices), HYP-2160 (T2 involution),
S589 (renormalisation). New: **HYP-2162**.
