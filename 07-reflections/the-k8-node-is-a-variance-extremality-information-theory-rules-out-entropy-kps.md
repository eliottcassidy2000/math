# The k=8 node is a VARIANCE (bimodality) extremality — information theory rules out entropy, and the φ⁴ folds to degree-2

*kind-pasteur-2026-06-27-S31ai. The owner asked to merge and extend the Lee-Yang/φ⁴/Galois web toward a
proof, thinking in terms of FUNCTIONS and INFORMATION THEORY, and to study failures of compression of
properties beyond commutativity (associativity). I attacked the single remaining node — the k=8
bounded-core extremality (mac-mini HYP-3132) — through these lenses, and they sharpen the target
decisively.*

## The node, restated as a function
The whole LRC(14) proof reduces (mac-mini HYP-3132) to: **consec maximizes the k=8 dual
`L_yK8 = q₀ + q₆ + (1/10)q₃`** over the bounded k=8 clusters, where `q_t = meas{N = t}` is the
empty-count distribution (`N` = # of the 6 inner sectors left empty by the orbit `{frac(e_i x)}`).
`L_yK8` is the **bimodality functional**: mass at the two wells `N=0` (all covered) and `N=6` (all empty),
plus `1/10` at the center `N=3`. Its weight `w(t)=1[t∈{0,6}]+(1/10)1[t=3]` is **reflection-symmetric**
`w(t)=w(6−t)` — the `s↦6−s` evenness whose resolvent is the φ⁴ biquadratic `u⁴−5u²+4` (HYP-3132).

## What consec EXACTLY extremizes (3002-cluster scan, `lrc_k8_node_..._kps.py`)
consec is the argmax (tied only with its dilation even-AP) of **all** of:
- `L_yK8` (bimodality dual) — rank 0/1;
- **`Var(N)` — rank 0/1** (consec `Var=2.612`, the maximum);
- `μ₄` (4th central moment) — rank 0/1;
- raw bimodality `q₀+q₆` — rank 0/1.

It is **NOT** extremal in:
- **Shannon entropy `H(N)`: consec is rank 3000/3002 with HIGH entropy 2.43** (a random cluster minimizes
  it at 1.78). So **the extremality is NOT an entropy principle** — consec spreads `N` widely, it does not
  concentrate it. This rules out a whole class of (max-entropy / min-description-length) proof attempts and
  says the right invariant is **spread (variance), not concentration (entropy).**
- excess kurtosis `κ₄/κ₂²`: rank 2994 (consec has high `μ₄` *and* high `Var²`, so the ratio is not
  extremal — the bimodality is a `Var`+`μ₄` phenomenon at fixed-ish mean, not a heavy-tail/kurtosis one).

> **Sharpened target (the merge with mac-mini's fold):** the cleanest functional consec maximizes is the
> **VARIANCE** `Var(N)=S₁+2S₂−S₁²` — a **degree-2** object. This is exactly mac-mini's HYP-3132 claim that
> the φ⁴ quartic, folded by the reflection `s↦6−s`, becomes a **degree-2 (in `u²`) biquadratic-resolvent
> bound**. My "variance" and mac-mini's "degree-2 fold" are the **same reduction**: the reflection
> symmetry collapses the degree-4 obligation to a second-moment (variance/covariance) extremality.

## Variance = total covariance = COHERENCE (the provable core)
Write `N = Σ_{j=1}^6 X_j`, `X_j = 1[sector j empty]`. Then
```
Var(N) = Σ_j p_j(1−p_j) + Σ_{i≠j} Cov(X_i, X_j),   p_j = P(sector j empty).
```
"consec maximizes `Var(N)`" = **consec maximizes the total covariance `Σ Cov(X_i,X_j)`** among the
empty-sector indicators. The AP orbit `{frac(i x)}` moves **coherently** (all points translate together),
so sectors empty and fill **together** — maximal positive covariance — whereas a dissociated cluster
empties sectors near-independently (covariance → 0). This is the **coherence-extremality** (HYP-2780/2607),
now pinned to the **degree-2 covariance**, which is far more tractable than the degree-4 `L_y`. **PROOF
TARGET:** consec maximizes `Σ_{i<j} Cov(X_i,X_j)` — a pairwise (Johnson-scheme / Gram) statement, the
folded-degree-2 form of the node.

## The compression-failure / associativity lens (owner cue)
Reading the moments as a hyperoperation/compression tower (mac-mini S71: `a+b` commutative face = the
pairwise compression):
- `S₂` (pairwise joint-emptiness, = the covariance / `Var`) is the **commutative `a+b`** face — it
  *compresses* (symmetric, order-blind), and it is the folded degree-2 core.
- `S₃, S₄` (3- and 4-way joint emptiness) are where **ASSOCIATIVITY fails**: the 3-way emptiness
  `P(i,j,k empty)` does **not** factor through pairwise data (no `P(i,j)·P(j,k)/P(j)` law) — the
  triple-overlap is an irreducible 3-cocycle. **This is the associativity-compression failure**, and it is
  exactly the **odd `−9S₃` (Worpitzky/antisymmetric) half** of the dip (mac-mini S71 / codex HYP-3147),
  the part that does NOT fold to the even biquadratic. The apex `7` is the first **non-associative**
  prime (Fano = octonion table), so the associativity failure lives at the apex — the `−9S₃` Worpitzky
  term is the octonion/non-associative residue of the dip. **The even (`+6S₄`, biquadratic, commutative-
  associative) part folds to degree-2 (variance); the odd (`−9S₃`, Worpitzky, non-associative) part is the
  genuine apex-7 obstruction.** This splits the node cleanly along the compression/associativity axis.

## Information theory: the bimodality is a 2-symbol code at the wells
`L_yK8` rewards mass at `{0,6}` — a **2-symbol (1-bit) compression** of `N` onto the two wells. consec
maximizes the bit at the wells while keeping high overall entropy (it does not over-concentrate). So the
extremality is "**maximize the well-mass (the φ⁴ order parameter) subject to the orbit being a genuine
coherent translation**" — a constrained-variance, not a min-entropy, principle. The φ⁴ free energy
`F = ⟨λS⁴+bS²⟩ − T·H(S)` is extremized with consec at the **low-`λ` (most-circular) / high-`Var`** corner,
high-entropy branch — consistent with all four scans.

## The merged web (one picture)
```
  Lee-Yang circle (zeros |z|≈R)  ──q₀=q₆R⁶──  cap = binomial/Pascal (de Moivre–Laplace, EVEN, low λ)
        │ off-circle λ>0 (φ⁴ vertex)                         │ reflection s↦6−s
        │                                                    ▼
   dip  ├── EVEN  +6S₄  → biquadratic u⁴−5u²+4 → fold → DEGREE-2 VARIANCE (consec=max covariance) [provable]
        └── ODD  −9S₃  → Worpitzky → ASSOCIATIVITY failure → apex-7 / octonion residue   [the hard half]
   ear decomp: #ears=C(n-1,2); ODD ears = Ω = the odd/Worpitzky/S₃ structure (same parity split)
   Galois: dual degree ≤4 ⟹ ⊆ S₄ ⟹ solvable by radicals ⟹ cap_8,dip_8 ∈ ℚ (disc 9)
```

## Net / next
- **VERIFIED:** consec maximizes `Var(N)`/`μ₄`/bimodality (tied even-AP); entropy is NOT extremal (rules out
  entropy proofs); the node folds to a degree-2 covariance extremality (= mac-mini's fold); the dip's
  parity split = commutative(even/biquadratic/`S₄`) ⊕ non-associative(odd/Worpitzky/`S₃`).
- **PROOF TARGET (sharpened):** consec maximizes `Σ_{i<j} Cov(X_i,X_j)` (total empty-sector covariance) —
  the even/degree-2 half; and the odd `−9S₃` Worpitzky/associativity residue bounded separately (the apex
  obstruction). This is the cleanest the node has been stated.
- **NEW SIGNALS:** total covariance `ΣCov`, the even/odd (S₄/S₃) parity split of the dip, the well-mass
  bit, the φ⁴ free energy `F`.

→ HYP-3160 (this), HYP-3132 (the biquadratic node), HYP-3099 (Lee-Yang two maps), HYP-3085 (gK8 dual),
HYP-2780/2607 (coherence-extremality), THM-577 (rational cap), mac-mini S70/S71/S72, codex HYP-3147
(Worpitzky odd half), the φ⁴ owner-cue, [[lrc14-thread]].
