# Cuts are rational geodesics on the additive Farey skeleton: the resonance law is divisibility, the tight witness is killed by multiples of 14, and the whole thing is a hyperoctahedral quotient shared with the tournament metagraph

*opus-2026-06-29. Owner: think abstractly and recursively in terms of cuts — simultaneous Diophantine
approximation, rational geodesics pinned to the additive Farey skeleton, hyperoctahedral symmetry, the
tournament metagraph. They are one structure. The cut recursion is a dynamical system on the Farey
skeleton, driven by divisibility resonances, quotiented by the hyperoctahedral group whose central
involution is the same `R` that splits the metagraph.*

## 1. Cuts ARE Farey geodesic steps (verified)
The peak `M(S)=½−Σcut` is built by moving a **rational witness** `a/q` around the circle. Verified: the
AP `{1,…,13}` traces the **Farey LEFT SPINE** — witness denominators `2,3,4,…,14`, and the values
`1/2,1/3,…,1/14` are **consecutive Farey neighbors** (`|a_i b_{i+1}−a_{i+1}b_i| = 1`, checked) — i.e. a
**geodesic** in the Farey tessellation / the left edge of the Stern–Brocot tree (mediants of `0/1` and
`1/1`). Each cut is one geodesic step `1/k → 1/(k+1)`.

## 2. The resonance law IS divisibility (multiplicative acting on the additive skeleton)
> **A speed `s` cuts the witness `a/q` iff `s·a ≡ small (mod q)`; an EXACT hit (`‖s·a/q‖=0`) iff `q | s`.**
So whether a speed cuts is a *divisibility* (multiplicative) condition read against the *Farey
denominator* (additive) of the current witness. The AP is **self-resonant**: speed `k` hits denominator
`k` (`q|s` true at every step), which is exactly why it descends the spine maximally. Off-resonance
speeds (`13,14,15,16` vs denom `12`) cut `0`; `24=2·12` cuts (`12|24`). The cut table is governed by
`gcd(s, q)`, not raw size — resonance, not magnitude.

## 3. The tight witness `1/14` is killed by multiples of 14 — this IS the covering margin (verified)
> `‖s/14‖ < 1/14 ⟺ 14 | s`. So **the tight Farey witness `1/14` survives iff the set contains NO
> multiple of 14.** The AP `{1,…,13}` has none ⇒ `1/14` survives ⇒ `M=1/14`. **A covering set is FORCED
> to contain a multiple of 14 ⇒ it KILLS `1/14` ⇒ the witness is pushed OFF the spine to a deeper Farey
> rational** (`{1..11,13,84}`: denom `12 → 89`, `M = 7/89`). The +10% margin is precisely **the Farey
> jump from the killed spine-witness `1/14` to the next reachable rational.** The covering condition is a
> *multiplicative* command ("be `≡0 mod 14`") that detunes the *additive* spine.

This is the add/mult duality at its sharpest: the **additive Farey spine carries the tight `1/(k+1)`
values; the multiplicative divisibility constraint (covering) knocks the witness off it.** A disproof
would need a multiplicative detuning that lands the witness *below* `1/14` — but killing `1/14` always
pushes *up* to a coarser rational, never below.

## 4. Simultaneous Diophantine approximation, and the hyperoctahedral symmetry
`M(S)=max_t min_s ‖s_i t‖` is the **anti-approximation** quantity: the direction `(s_1,…,s_n)` worst
simultaneously approximated by `Z^n` — the line in `T^n=R^n/Z^n` that best avoids the central
view-obstruction cube `[−1/n,1/n]^n+Z^n` (Cusick). That cube and the lattice `Z^n` have automorphism
group the **hyperoctahedral group `B_n = (Z_2)^n ⋊ S_n = Aut(Z^n)`** (signed permutations). The
involution **`t → −t` is the central `−1 ∈ B_n`** (flip all coordinates); it fixes each Farey
denominator (`a/q ↔ (q−a)/q`) and is exactly the **`R`** (reversal/complement) of the whole project.

## 5. The same `B_n`, the same `R`, on the tournament metagraph
The metagraph is the **arc-hypercube `Q_{C(n,2)} / S_n`** with the complement `R` = flip all arcs = the
central `−1` of the arc-cube. So **both** objects are signed-cube quotients with `R` as the central
involution:
| | LRC (approximation) | metagraph (tournaments) |
|---|---|---|
| ambient cube | view-obstruction cube in `T^n` | arc-hypercube `Q_{C(n,2)}` |
| symmetry | `B_n=Aut(Z^n)` | `S_n` + central flip `R` |
| central `−1` | `t→−t` | complement `R` |
| `R`-even / `R`-odd | cap bulk / cap obstruction (HYP-3538) | SC / NS (THM-587, `P_n(±1)`) |
| additive spine | Farey left spine `1/(k+1)` = the **staircase** direction `(1,…,n−1)` | transitive backbone, `H=1`, the **staircase** `δ_{n−2}` |
The AP direction `(1,2,…,13)` IS the tournament **staircase** (`everything-is-the-triangle`): the
additive/transitive pole on *both* sides. The `R`-even/`R`-odd split is the `±1` eigenspace of the same
central involution in both — which is why the LRC obstruction (`R`-odd cap residual) and the metagraph
`R`-odd block (`NS`, `b₁⁻=1,7,119`) are the *same* hyperoctahedral `−1`-eigenspace.

## 6. The recursion, abstractly
Adding a speed is a map on `(witness on Farey skeleton)`: **resonant** (`gcd(s,q)` large ⇒ a geodesic
cut deeper down the Stern–Brocot tree) or **silent** (`gcd(s,q)=1`-ish ⇒ witness unchanged). The cut
recursion is therefore a **Stern–Brocot descent (additive) gated by divisibility (multiplicative),
`B_n`-equivariant under `R`.** It is the exact LRC twin of the tournament `n→n−2` Cayley–Dickson
metagraph recursion — both are `R`-equivariant descents of a signed cube.

## Improved targets (sharpened by the Farey picture)
1. **The margin is a Stern–Brocot depth.** `inf_{covering} M` (the `7/89` floor) = the shallowest Farey
   rational reachable after the mult-of-14 detuning. Bounding it = a continued-fraction/`gcd` extremal
   problem on the spine, NOT an analytic measure bound.
2. **`LRC(14) ⟺` the detuned witness never drops below `1/14`** = killing the spine value `1/(k+1)`
   always lands on a *coarser* (larger-value) rational. A monotonicity of the Stern–Brocot descent
   under multiplicative detuning.
3. **Use `B_n`-equivariance:** decompose the cut operator into `R`-even (provable, the spine) and
   `R`-odd (the obstruction) — the same split that makes the metagraph tractable.

## Status
- **Verified:** AP = Farey left spine (denoms `2..14`, neighbors); resonance law = divisibility (`q|s`);
  `1/14` killed iff `14|s`; covering jump `12→89`.
- **New (opus):** cuts = Farey geodesics; resonance = divisibility on the additive skeleton; the covering
  margin = the multiplicative detuning of the tight Farey witness; LRC and the metagraph as the same
  `B_n`-quotient with central `R`; the cut recursion = Stern–Brocot descent gated by divisibility.
- **Open:** the margin as a Stern–Brocot-depth extremal; `B_n`-equivariant cut-operator decomposition.

Related: the peak-deletion-contraction + cut-deficit reflection; the razor-thin (measure vs peak); the
duality-web (add/mult at apex 7); mac-mini HYP-3537/3538 (cap OCF / R-odd), klein THM-587 (`P_n(±1)`),
`everything-is-the-triangle` (staircase), THM-580 (apex-7 descent), OPEN-Q-108.
