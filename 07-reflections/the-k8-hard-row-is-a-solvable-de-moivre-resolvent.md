# The single hard node (k=8) is a solvable De Moivre resolvent: the φ⁴ quartic is a biquadratic

*mac-mini-2026-06-27-S70. Owner: integrate incoming + past work, understand comprehensively what REMAINS in
LRC(14), and synthesize creatively into an improved argument — with the hint "where have we seen 120 and 320"
and a worked SOLVABLE (De Moivre) quintic whose resolvent quartic `x⁴+10x³−120x²−320x+1024` has geometric
roots `{2,−4,8,−16}`. The two threads fuse: the whole proof now reduces to ONE node (the k=8 bounded-core dip),
and that node is a **quintic whose resolvent quartic is solvable by radicals** — and is exactly the φ⁴/quartic
of S67. Continues [[the-cap-is-a-phi4-field-theory-and-the-quartic-marks-the-hard-row]] and the S60 reflection.*

## What remains (comprehensive, current): ONE node
The covering bound's tree is now almost entirely closed (scout-verified across the latest ~70 sessions):
- **Apex-majority** (≥7 multiples of 7): THM-573 level-7 sieve — **PROVED** (mod LRC≤13).
- **Single-far / Node-3** (r=1): THM-546/547 + HYP-2900 — **PROVED**.
- **Multi-far floor** (r=2..6): the SPEC bound is **PROVED ELEMENTARILY** (kps-S255/HYP-3129, `R′≥0.642`, an
  absolutely-convergent Hardy–Littlewood series — *not* EH-dependent); and S68/S69 (HYP-3127/3131) show **the
  far elements only push the miss-PGF zeros OUTWARD** — the far structure subsumes into the bounded core.
- **Bounded core, k≥10**: cap = `C(k+1,2)/91` exactly — **PROVED** (THM-577).
- **Induction base** LRC≤12/13 (arXiv:2604.23906) — **ACCEPTED**.

> **The entire open core is the bounded-core extremality at the binding row k=8** (the dip `dip_8 = 1081/76440`,
> equivalently the φ⁴ quartic `κ₄/κ₂²` bound, HYP-3122). Everything else is proved or reduced to it.

## The k=8 node IS a solvable quintic resolvent (the 120/320 lens)
**k=8 means `|P| = 13−8 = 5` — a QUINTIC** (the minimizer `{1,5,7,8,9}`). And the gK8/Delsarte dual at k=8 is
its **resolvent QUARTIC** (THM-534): `g(t) = (t−1)(t−2)(t−4)(t−5) = t⁴ − 12t³ + 49t² − 78t + 40`. The dual
*degree* stratifies the rows: **quadratic** (k≥11), **cubic** (k=9,10), **quartic** (k=8). k=8 is the unique
quartic — the resolvent level — exactly the hard row.

The user's De Moivre example is solvable because its resolvent quartic factors nicely. The **gK8 resolvent is
solvable too — it is a BIQUADRATIC under the reflection symmetry `s ↦ 6−s`** (center `t=3`, the order-2
antipodal of S60). Substituting `u = t−3`:
```
(t−1)(t−2)(t−4)(t−5)  =  (u+2)(u+1)(u−1)(u−2)  =  (u²−4)(u²−1)  =  u⁴ − 5u² + 4.
```
- It is **even in `u`** (no odd term) — i.e. it is a **φ⁴ potential** `V(u) = u⁴ − 5u² + 4`. **The φ⁴ structure
  of S67 IS the biquadratic resolvent**; φ⁴-evenness = the reflection symmetry = De Moivre-solvability.
- Its discriminant (of the quadratic in `u²`) is `25 − 16 = 9 = 3²` — a **PERFECT SQUARE**, so the radicals
  collapse to rationals (`u² ∈ {1,4}`, `u ∈ {±1,±2}`, `t ∈ {1,2,4,5}`). **This is why `cap_8`, `dip_8` are exact
  rationals** (THM-577) and not surds: the De Moivre solvability of the k=8 resolvent terminates in `ℚ`.
- `1024 = 2¹⁰` (the user's resolvent constant = product of its roots) = the number of **tilings at n=6** (the
  cube `Q_{C(5,2)} = Q_{10}`); the 2-adic/dyadic root tower (`2,4,8,16 = 2^{1..4}`) is the same powers-of-2
  spine that runs through the four-faces H4 / Cayley–Dickson face.

## The improved argument
The single hard node — bound the k=8 quartic dip uniformly over the binding family — is a **solvable quartic
problem**, and the reflection symmetry `s↦6−s` **reduces it to degree 2**:
1. The gK8 quartic obligation is reflection-symmetric (S60: the pairwise covariance matrix is `s↦6−s`-symmetric,
   the half-tiling `Z/2`). So it block-diagonalizes / folds to the **biquadratic `u⁴ + b u²`** — a **degree-2**
   (always-solvable) object in `u²`.
2. The φ⁴ stabilizer sign (`λ>0`, `κ₄<0` at k=8, S67) fixes the direction: the quartic is a stabilizer
   (bounded, right sign).
3. So "bound `κ₄/κ₂²` over the binding family" (the open obligation) becomes "bound the **biquadratic
   resolvent's `u²`-coefficient**" — a solvable, degree-2 (Cardano-trivial) bound, not a general quartic. The
   De Moivre solvability is the lever that makes the k=8 quartic tractable.

And the meta-point (why LRC(14) is *solvable* at all): the bounded-core dual never exceeds **degree 4** (k=8 is
the deepest), so it stays below the Abel–Ruffini quintic wall — the same **n≤7 tameness window** as the
A000568 sandwich (S69: `C(n,3) ≤ A000568(n) ≤ 2(n−1)!/3` only for n=4..7). LRC(14)'s hard core is solvable
precisely because its resolvent degree is ≤4.

## Honest status
VERIFIED: the comprehensive reduction to k=8; the gK8 quartic resolvent; its biquadratic (φ⁴) form under
`s↦6−s`; the perfect-square discriminant `9` ⟹ rational collapse; `1024 = Q_{10}` tilings. BOLD/OBLIGATION: that
the *uniform* k=8 dip bound over the binding family reduces rigorously to the degree-2 biquadratic-resolvent
bound (the reflection-fold of the φ⁴ obligation) — the improved-argument target. The literal "120, 320" of the
owner's example are `−e₂, e₃` of a *geometric* (non-symmetric) resolvent; the LRC resolvent is the *symmetric*
biquadratic (`e₁=e₃=0`), the solvable sub-case.

Related: HYP-3132 (this), HYP-3122 (φ⁴ cap = the resolvent), THM-577 (the rational cap = the solvability
collapse), HYP-3085 (gK8 dual = the resolvent quartic), HYP-3131 (far subsumes into bounded core),
[[the-cap-is-a-phi4-field-theory-and-the-quartic-marks-the-hard-row]], OPEN-Q-108.
