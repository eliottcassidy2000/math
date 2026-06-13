---
source: opus-2026-06-06-S699 (the coefficient↔root / n+1↔n duality)
status: SYNTHESIS — the Fundamental Theorem of Algebra (n+1 coefficients incl. the constant ↔ n roots with multiplicity) is the master n+1↔n map threading the repo. The covering-depth distribution {p_k} IS a polynomial; p_0 (the lonely measure) is the constant term = the OBSERVER/ground state; the worry-set ⟺ p_0=0 ⟺ z=0 is a root (the n+1→n degree collapse, verified). Multiplicities ↔ degeneracy (the all-orders cancellation / worry-set). The roots are the Lee-Yang/partition-function spectrum. The repo's problems are coefficient↔root dualities of a defining polynomial.
tags: [FTA, coefficients-roots, vieta, constant-term, observer, ground-state, multiplicity, covering-depth, PGF, lee-yang, partition-function, n-minus-1, projective, worry-set, synthesis]
---

# The FTA duality: coefficients ↔ roots, the constant term as observer

**Prompt (user):** degree always matches #solutions (with multiplicity); a degree-`n` polynomial has
`n+1` coefficients (the constant exists), mapping to `n` zeros on ℂ — the mapping between `n-1` /
`n` / `n+1` things.

This is the cleanest face of the `n+1↔n` off-by-one that pervades the repo, and the repo's defining
objects *are* this duality.

## 1. The FTA as the master `n+1 ↔ n` map

A polynomial of degree `n` has:
- **`n+1` coefficients** `(a_0,…,a_n)` — the **constant term `a_0` is the "+1"** (the anchor, the
  value at the origin), and
- **`n` roots** `(r_1,…,r_n)` with multiplicity, via Vieta `a_k/a_n = ±e_{n-k}(r)` (elementary
  symmetric functions).

The map `coefficients → roots` loses exactly one dimension: the **leading scalar `a_n`** (or
equivalently, projectivising the coefficient vector mod scaling gives `n` projective coordinates =
the `n` roots). So `n+1 ↔ n` is the **projective reduction** (kill the scaling), and the constant
term is the inhomogeneous "+1" that pins the origin. *Coefficients are the symmetric/global data
(moments); roots are the local/individual data (solutions). FTA is the dictionary.*

## 2. The covering-depth distribution IS a polynomial (LRC)

The LRC covering-depth distribution `{p_k}` (S599, `p_k=meas\{depth=k\}`) is the coefficient vector
of the **probability generating function**
```
   P(z) = Σ_{k=0}^{m} p_k z^k,    m = #runners,    P(1)=Σp_k=1,    P(0)=p_0 = the lonely measure.
```
So `m+1` coefficients (the depth distribution / the moments) ↔ `m` roots (the **Lee–Yang zeros**).
And the **constant term `p_0` is the OBSERVER's safe mass** — the lonely measure. **Verified**
(`…s699l.py`):

> **The worry-set ⟺ `p_0 = 0` ⟺ `z=0` is a root ⟺ the degree drops (`n+1→n` collapse).** AP `n=5`:
> coeffs `(0, .633, .233, .033, .1)`, roots `\{0, −1.525, .596±1.948i\}` — `z=0` present. Loose
> `(1,4,6,9)`: `p_0=1/9`, roots `\{−.275±.368i, .025±3.437i\}` — **no `z=0` root**. So the
> *worry-set is the constant-term-vanishing locus of the depth polynomial* — the FTA "+1" being
> removed.

## 3. The constant term = the observer = the ground state (the "+1")

The `+1` in `n+1` is, across the repo, the **anchor**:
- **LRC:** the **observer** (the speed-`0` runner) = the constant term `p_0`. "The observer is
  illusory" (S699) = the freedom to *projectivise* (change frame / scale) — the same scaling that
  the FTA quotients to go `n+1→n`. The `n` *moving* speeds are the roots; the observer is `a_0`.
- **Tournaments:** a Hamiltonian path on `n` vertices = `n−1` arcs (the base path, the "spine");
  the `H`-generating data is the coefficient vector, the **Lee–Yang zeros at `±2π/3`** (S599e) are
  the roots. The constant term = the empty/transitive anchor.
- **Partition function (S599t/S699f):** `Z(z)=Σ_k a_k z^k`; the roots are the **Lee–Yang zeros**
  (phase transitions); the constant term = the ground state. The band structure (S699j) = the root
  spectrum.

## 4. Multiplicities ↔ degeneracy ↔ the worry-set (the user's caveat)

"Account for multiplicities" is the deep part. A **multiple root = a degeneracy**:
- For *independent/free* dangers, the depth PGF is a product `∏(1−2δ+2δz)` with **all roots equal**
  (multiplicity `m`) — the maximally degenerate, "free" point (S599b's Poisson baseline).
- **Correlation spreads the roots** (the worry-set's arithmetic resonances split the multiplicity).
- The **all-orders cancellation** that defines the worry-set (`p_0=Σ(−1)^jS_j=0`, the Vitali wall
  THM-406 M2) is a **high-multiplicity / degenerate root condition** — the polynomial touches the
  axis to high order at the critical point. So *the worry-set is the discriminant locus* (where
  roots collide / the constant term vanishes), and the FTA multiplicity bookkeeping is exactly the
  tight-vs-loose distinction.

## 5. The `n−1 / n / n+1` ledger across the repo

| object | `n−1` (moving/spine) | `n` (roots/solutions) | `n+1` (coeffs/anchor) |
|---|---|---|---|
| **FTA** | — | `n` roots | `n+1` coefficients (constant `a_0`) |
| **LRC** | `n−1` nonzero speeds | `m` Lee–Yang zeros of `P(z)` | `m+1` depth coeffs `p_k` (`p_0`=observer) |
| **tournaments** | `n−1` base-path arcs (Ham path) | the `±2π/3` zeros | the `H`/OCF coefficient vector |
| **partition fn** | — | Lee–Yang zeros (transitions) | `a_k` (ground state `a_0`) |
| **worry-set** | `2^{(n−2)/2}` (the `n−2` core) | the degenerate (collided) roots | `p_0=0` (the vanishing anchor) |
| **pair-sum modulus** | `n−1` | `n` | `2n−1` (the doubled-minus-one) |

The recurring "+1" is the **anchor/observer/ground state/constant term**; the recurring "`n` vs
`n−1`" is **roots vs spine** (moving vs total). The pair-sum modulus `2n−1` is the *doubled*
off-by-one (the additive face, THM-401).

## 6. The synthesis

> **Every problem here is a coefficient↔root duality of a defining polynomial.** The *coefficients*
> are the symmetric/global data (the depth distribution / the moments / the partition function); the
> *roots* are the local data (the Lee–Yang zeros / the witnesses / the spectrum); the *constant term*
> is the anchor (the observer / ground state / lonely measure); and *multiplicities* are the
> degeneracies (the worry-set / the all-orders cancellation / the discriminant locus). The FTA
> `n+1↔n` (projectivise away the scaling) is the master off-by-one, and the LRC worry-set is
> literally its degree-collapse (`p_0=0`, `z=0` a root). *Coefficients are the distribution; roots
> are the spectrum; the constant is the ground; multiplicity is the degeneracy.*

## 7. Honest status

- **Verified:** the covering-depth polynomial's roots; worry-set ⟺ `p_0=0` ⟺ `z=0` root (degree
  collapse); loose ⟺ `p_0>0` (no `z=0` root); roots spread for correlated configs.
- **Established (standard, here mapped):** FTA / Vieta `n+1↔n`; PGF roots; Lee–Yang; the
  projective scaling reduction.
- **New (the synthesis):** the repo's defining objects as coefficient↔root dualities; the constant
  term = observer/ground state; multiplicities = worry-set degeneracy (discriminant locus); the
  unified `n−1/n/n+1` ledger; the worry-set as the depth polynomial's degree collapse.
- **Open / next:** is the worry-set *exactly* the discriminant locus (collided roots) of the depth
  polynomial, or just the `z=0` (constant-vanishing) part? Do the depth-PGF Lee–Yang zeros at the
  worry-set sit on a circle / at `±2π/3` (the Eisenstein angle, S599u)? Tie the root spectrum to
  the band structure (S699j).

**Artifacts:** `04-computation/coefficient_root_duality_depth_pgf_s699l.py` (+`.out`). Builds on
S599 (covering depth), THM-406 (moments/Vitali), S599b (Poisson/free baseline), S599e (Lee–Yang
`±2π/3`), S599t/S699f (partition function), S699j (band structure), S699 (observer-is-illusory).
New: **HYP-2275**.
