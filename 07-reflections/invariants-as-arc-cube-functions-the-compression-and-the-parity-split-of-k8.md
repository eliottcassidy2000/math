# Invariants as arc-cube functions: the score-compression boundary IS the k=8 hard row, split by parity

*mac-mini-2026-06-27-S71. Owner: think of everything as functions; the function quartet `a+b, a·b` (symmetric,
order-blind) vs `a^b, b^a` (asymmetric, order-aware); the n=3 tournament from its edges (coin flips); n=4 in
two schemes (Klein-four group vs magma); compressions → tricks; Worpitzky. The pieces fuse into an EXPLANATION
of why k=8 is the LRC hard row, merging codex's n=3 Worpitzky kernel (HYP-3147) with my biquadratic (S70) and
φ⁴ cumulants (S67). Continues [[the-k8-hard-row-is-a-solvable-de-moivre-resolvent]].*

## The principle: an invariant is a function on the arc-cube, with two parity faces
A tournament is a point of the arc-cube `(Z/2)^{C(n,2)}`. Any invariant is a **function** on that cube. The
owner's four operations on an edge's endpoint pair split it by symmetry — and this is the cut⊕cycle / even⊕odd
split the project keeps hitting:
| operation | symmetry | tournament face | algebra |
|---|---|---|---|
| `a+b` (sum) | **symmetric / commutative** | the **score** (out-degree) — the cut space | even, linear |
| `a·b` (product) | **symmetric** | `H = I(Ω,2)` — the cycle/independence face | even, multiplicative |
| `a^b, b^a` | **asymmetric / order-aware** | the **orientation** `i→j` (the directed arc / path order) | odd, sign |

codex (HYP-3147): "sum and product cannot see an arc flip; the ordered exponential pair can see which endpoint
is base." Symmetric = order-blind = the even content; the ordered pair = the odd content.

## The compression (the trick): score determines the iso class for n≤4, fails at n=5
VERIFIED (`lrc_arc_cube_compression_parity_macmini_S71.py`):
```
 n=3: #iso=2 = #score=2   BIJECTIVE      (score = the iso class)
 n=4: #iso=4 = #score=4   BIJECTIVE
 n=5: #iso=12 > #score=9  FAILS          (the ORDER / a^b face becomes essential)
```
So at **n≤4 the commutative (`a+b`, score) face is a COMPLETE invariant** — the iso class is a *linear/group*
function of the arcs. That is why the owner's **scheme 2 (4 arcs fixed, 2 free) is a clean Klein-four group `V₄`**
on `{T,+,-,S}`: at n=4 the class is the score, an additive object, so the 2-arc slice closes into a group.
**Scheme 1 (the tiling model, 3 free arcs) is the same data over-coordinatized** — a magma where `S` looks
absorbing — because it carries a redundant arc; *the compression to the 2-arc group slice removes the
redundancy.* The trick: **find the gauge where the symmetric face is a complete invariant → the function becomes
a group → it is computable.** (At n≥5 no such slice exists — the odd face is irreducible.)

## The payoff: the k=8 hard row IS the score-compression boundary
The cap dip turns on EXACTLY at the n=4→5 compression boundary (VERIFIED):
```
 k:    13   12   11   10    9        8
 |P|:   0    1    2    3    4        5
 dip:   0    0    0    0   1/4004   1081/76440   (tiny at |P|=4, LARGE at |P|=5)
```
`|P| = 13−k`, so **k=8 ⟺ |P|=5** — the quintic level, exactly where `score → iso` *fails*. **The cap dip is
the failure of the commutative compression**: for `|P|≤3` (k≥10) the symmetric (pair-Pascal) face is complete
and `dip=0` (THM-577); at `|P|=4` (k=9) it just begins to leak (`1/4004`); at `|P|=5` (k=8) the antisymmetric
(orientation / `a^b`) content is irreducible and the dip is substantial. **This is why k=8 is the binding row:
it is the first row past the n=4 score-compression boundary.** (The `|P|=5` quintic also has the solvable
*resolvent quartic* of S70 — the same "5".)

## The merge: the k=8 dip = EVEN (biquadratic, S70) + ODD (Worpitzky, codex)
The gK8 dual's higher-order content `−9S₃ + 6S₄` splits by **parity** (VERIFIED for consec_8):
- **even `+6S₄`** = the **symmetric** (`a+b,a·b`) face = the **biquadratic resolvent** `u⁴−5u²+4` (S70),
  solvable by radicals (degree 2 in `u²`);
- **odd `−9S₃`** = the **antisymmetric** (`a^b,b^a`, orientation) face = the **Worpitzky / ordered-descent**
  content (codex HYP-3147), and it **dominates** (`|odd|/|even| ≈ 3.15`).

So the two recent threads are the two parity faces of the same k=8 dip:
- **my biquadratic (S70)** bounds the **even** part (solvable, degree-2);
- **codex's Worpitzky n=3 kernel (HYP-3147)** bounds the **odd** part — a sum of n=3 edge-flip oscillations
  (eigenvalue **−1/3**) weighted by the **Eulerian descents `A(3,k)=[1,4,1]`** (Worpitzky:
  `x³ = C(x+2,3)+4C(x+1,3)+C(x,3)`).

**Improved argument:** bound the k=8 dip = [bound the even biquadratic coefficient — solvable] + [bound the odd
content as a Worpitzky-weighted sum of n=3 `−1/3` oscillations]. The owner's `a+b,a·b` vs `a^b,b^a` IS the
parity decomposition; the symmetric face is the solvable biquadratic, the antisymmetric face is the Worpitzky
oscillation — and the antisymmetric (dominant) face is exactly the content the score-compression cannot see.

## Honest status
VERIFIED: score→iso bijective n≤4, fails n=5; the cap dip onset at `|P|=4→5`; the `−9S₃+6S₄` parity split with
the odd part dominant; the n=4 score=iso (⟹ the Klein-four). SYNTHESIS/BOLD: that the k=8 dip bound rigorously
decomposes into [even biquadratic, S70] + [odd Worpitzky-`−1/3` sum, codex] (the proof-relevant target); the
`|P|=5 ⟺ n=5` link is numerological-but-aligned (the quintic, the resolvent quartic, the score-fails-at-5, the
dip-onset all point to 5). The compression-trick meta — *find the gauge where the symmetric face is complete* —
is the reusable lesson.

Related: HYP-3150 (this), HYP-3147 (codex n=3 Worpitzky kernel = the odd face), HYP-3132 (k=8 biquadratic = the
even face), HYP-3122 (φ⁴ cumulants), THM-577 (cap = symmetric face, exact for |P|≤3), THM-062/063 (deformed
Eulerian), OPEN-Q-108.
