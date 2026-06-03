---
source: opus-2026-06-03-S599f (remote-control)
status: CREATIVE MENU + one VERIFIED anchor — how the SOLVED face (Rédei/GF(2) involution on the permanent-shaped Ham-path count) transfers truths & speedups to LRC/Collatz/Riemann. VERIFIED: worry-set = self-converse round-tournament classes = 2^⌊(n-2)/2⌋ (64 at n=14) ⟹ an enumeration SPEEDUP (count tournaments, not configs) and a Rédei forced-parity certificate. Plus 8 further transfer mechanisms, each labelled by rigor.
tags: [redei, OCR, OCF, round-tournaments, worry-set, involution-principle, garsia-milne, GF2, determinant-permanent, circulant-spectrum, collatz, transfer, speedup, LRC, n14]
---

# Mining the solved face: how Rédei transfers to the unsolved siblings

**Prompt (user):** see how the solved face can let us access truths or computational speedups
in the others; spend a long time being creative.

The solved face is **Rédei**: the Hamiltonian-path count of a tournament is a *permanent-shaped
sum* whose all-orders cancellation is **resolved** — its parity is forced to `1` by a
sign-reversing involution (= a determinant over `GF(2)`, `⊕P` collapse). The solved face is a
**toolkit**, not a curiosity, and several of its tools cross over. I give one *verified* anchor
and then a long menu of transfer mechanisms, each labelled **[verified] / [rigorous tool] /
[concrete conjecture] / [speculative]**.

---

## ANCHOR — the worry-set IS a tournament-counting problem  **[verified]**

THM-402 mapped the LRC-tight set to **round (locally transitive) tournaments**. This session
made the count exact (`redei_transfer_round_tournaments_s599f.py`):

```
 m (runners) :   3    4    5    6
 round iso-classes:   2    3    6   11
 SELF-CONVERSE classes:  2    2    4    4     =  2^⌊(m−1)/2⌋
 round Ham-counts: {1,3}, {1,3,5}, {1,3,5,9,11,15}, {1,3,5,9,11,15,17,23,41,45}  — ALL ODD
```

`2^⌊(m−1)/2⌋` with `m=n−1` gives `2^⌊(n−2)/2⌋` = **64 at n=14** — exactly the repo's worry-set
count. So:

> **The LRC worry-set = the self-converse round-tournament iso-classes.** Finding the worry-set
> is therefore a *tournament-counting* problem (a small, structured iso-class enumeration the
> repo's transfer-matrix / canonical-form machinery does fast), **not** an exponential search
> over integer speed sets. **Speedup:** enumerate `2^⌊(n−2)/2⌋` tournament classes, not `~N^{n}`
> configs. **Truth:** every one has an *odd* Hamiltonian-path count (Rédei) — a forced-parity
> invariant attached to each worry-set element.

This anchor is what makes the rest credible: the solved face's objects *are literally* the hard
residual's objects.

---

## TRANSFER 1 — the forced-parity certificate (Rédei ⟹ an LRC obstruction)  **[concrete conjecture]**

Rédei forces `#HamPaths(T) ≡ 1 (mod 2)` for *every* tournament. A **counterexample** to LRC at
`n` would be a config whose round tournament `T(S)` fails to admit the clock witness. Program:
find a *loneliness functional* `L(T)` (built from the odd-cycle / Ham-path data of `T(S)`) such
that `M(S) < 1/n ⟹ L(T(S))` is **even**. Then Rédei's forced oddness *refutes the
counterexample*. The repo's **OCR / Odd-Cycle Collection Formula** is the natural `L`: it
already expresses tournament invariants as odd-cycle collections with a forced parity. **This is
using the solved cancellation (forced odd) to kill the unsolved one (no counterexample)** — the
single most valuable possible transfer. Status: the link `M(S)<1/n ⟹ L even` is unproven; the
*vehicle* (OCR + Rédei) is rigorous and in-repo.

## TRANSFER 2 — the involution principle for (★)  **[rigorous tool]**

Rédei's proof is a **sign-reversing involution** on path-pairs. Garsia–Milne's **involution
principle** is the general machine. Import it to `p_0 = Σ_S(−1)^{|S|}meas(⋂_S D_i)` `(★)`: seek
an involution `ι` on the subset lattice with `meas(⋂_S)=meas(⋂_{ι(S)})` and `|ι(S)|=|S|±1`, so
the alternating terms cancel in pairs, leaving `p_0 = Σ_{fixed points}`. The **fixed-point set
is the worry-set** (the un-cancellable core). This converts "evaluate/bound the all-orders sum"
into "exhibit the involution and its fixed points" — exactly how Rédei sidesteps the permanent.
The involution candidate: pair a danger-overlap `⋂_S` with its image under the doubling map
`t↦2t` (THM-404) or the antipode `t↦t+½` — both are measure-preserving on arcs. Status: the
tool is rigorous; constructing the *specific* measure-preserving fixed-point-=worry-set
involution is the concrete task.

## TRANSFER 3 — the mod-2 determinant speedup  **[rigorous tool / computational]**

`#HamPaths mod 2` is a **`GF(2)` determinant** (poly-time), even though the count is `#P`. The
general principle (`⊕P` via `det` over `GF(2)`): the *parity* of a permanent-shaped sum is often
a determinant. Transfers:
- **LRC:** the parity of the worry-set's lattice-point counts (the resonance sums `S_j mod 2`)
  may be a `GF(2)` determinant of the owner/slack two-block (S595) — a poly-time *parity
  certificate* for emptiness, complementing the bounded-CRT automaton.
- **Collatz:** the *number of cycles of bounded length mod 2* via a `GF(2)` linear system on the
  `2^E−3^k` residues — "no cycle" as a `GF(2)` kernel triviality (à la Rédei's forced parity).

## TRANSFER 4 — the circulant-spectrum speedup (→ THM-406)  **[rigorous tool]**

Round tournaments are **circulant** (rotational), so their adjacency/skew matrices have
**closed-form eigenvalues** (sums of roots of unity). Via THM-406 (`{p_k}` = spectral measure of
the depth operator), the worry-set's depth distribution *near the floor* is governed by the
circulant spectrum — computable by **FFT** in `O(n log n)` instead of the arrangement scan. The
solved face's algebra (circulant diagonalization, the same that makes Rédei's count tractable)
**speeds up the master-object computation** at the critical configs. Concretely: the `p_0`-vanishing
rate (the singleton exponent α=1, S599) is read off the smallest circulant eigenvalue's behaviour.

## TRANSFER 5 — the OCF as the resonance accountant  **[speculative]**

The repo's **OCF** organizes a tournament's parity invariant as a signed sum over odd-cycle
collections, with cancellations that *telescope* (the determinant-ization). The LRC resonance
overlaps `meas(⋂_S D_i)` are themselves an odd/even-cycle-weighted sum on the runner graph
(a 3-term-relation = a 3-cycle, S577). **Proposal:** run the OCF's cancellation bookkeeping on
`(★)` — the odd-cycle collections of the resonance graph should be the surviving (worry-set)
terms, the even ones cancelling. This would *derive* the all-orders cancellation structure from
the OCF rather than discovering it ad hoc.

## TRANSFER 6 — Collatz via the GF(2) cycle space  **[speculative]**

The repo's **tournament ↔ even-graph ↔ cycle-space (`GF(2)` homology)** ladder (CLAUDE.md) is a
solved transfer machine. The **Collatz reverse tree** has a cycle space; a nontrivial Collatz
cycle = a nonzero class in a `GF(2)` homology built from the `×3+1 / ÷2` parities. The solved
face says when such a homology *vanishes* (the even-graph is acyclic over `GF(2)`). **Proposal:**
encode the cycle equation `a₁(2^E−3^k)=S` as a `GF(2)` boundary condition and show the only
0-cycle is trivial — a homological restatement of "no nontrivial Collatz cycle," importing the
repo's even-graph vanishing criteria.

## TRANSFER 7 — the AP ↔ rotational-tournament extremality (a loneliness monotone)  **[concrete conjecture]**

The AP (the extremal worry-set element) maps to the **rotational tournament `R_m`** (positions
`{1,…,m}` evenly spread at `t*=1/n`); the data shows round Ham-counts range up to `45` at `m=6`,
the rotational/regular one extremal. **Conjecture:** `#HamPaths(T(S))` is a **monotone proxy for
loneliness** — more Ham-paths = more "spread/regular" = larger `M(S)`. If so, a *lower bound on
the Ham-path count* (which Rédei + round structure supply) gives a **lower bound on `M(S)`** —
i.e. a loneliness guarantee from a tournament count. Testable: regress `M(S)` against
`#HamPaths(T(S))` over configs. This would turn a *counting* theorem into a *metric* (geometric)
bound — the deepest kind of transfer.

## TRANSFER 8 — the non-natural-certificate meta-guidance (dodge the barrier)  **[rigorous, meta]**

THM-406 M2: no finite-moment / measure / "natural" property certifies the sign of `p_0` (the
Vitali wall = the natural-proofs barrier). The solved face shows the *form* of a certificate that
**evades** the barrier: an **algebraic involution / parity invariant**, not a measure bound.
So the meta-lesson, now rigorous, is **prescriptive**: the LRC and Collatz closures must be of
*Rédei type* (an exact cancellation certificate over the right ring), and any measure/energy
attempt is provably blind. This *rules out* a class of approaches and *names* the target.

## TRANSFER 9 — Riemann's explicit formula as a permanent→determinant collapse  **[speculative]**

RH's `ψ(x)−x = −Σ_ρ x^ρ/ρ` is an all-orders oscillatory sum whose sign control is the
zero-location. The solved-face template ("find the structure turning the permanent-sum into a
determinant") suggests reading the **Weil explicit formula as a trace/determinant** (Connes'
spectral interpretation / the would-be Hilbert–Pólya operator). The repo contribution is modest
but pointed: the *same involutive/parity bookkeeping* that collapses Rédei is the combinatorial
shadow of "the zeros are eigenvalues" — the resonance terms pair up except on a self-adjoint
(real-line) core. Speculative, but it places RH in the same template with a concrete instruction.

---

## The unifying instruction (what the solved face teaches)

Every entry is one sentence: **find the structure that turns the permanent-shaped all-orders sum
into a determinant-shaped (certifiable) one.** Rédei's answer is *an involution over `GF(2)`*;
the transfers ask whether LRC/Collatz/Riemann admit the analogous structure (an involution, a
`GF(2)` determinant, a circulant diagonalization, an OCF telescoping, a spectral trace). The
**verified anchor** (worry-set = self-converse round classes = `2^⌊(n−2)/2⌋`, all Rédei-odd) is
proof that the solved face's *objects* already coincide with the residual's, so these are not
analogies of convenience — they are the same combinatorial species, one member of which is
solved.

## Honest status

- **[verified]:** worry-set = self-converse round-tournament classes, count `2^⌊(m−1)/2⌋`
  (`2,2,4,4` for m=3..6 ⟹ 64 at n=14); all round tournaments have odd Ham-counts (Rédei) ⟹ the
  enumeration speedup and the per-element parity certificate are real.
- **[rigorous tool] (cross-over technique, application open):** involution principle (T2),
  `GF(2)`-determinant parity (T3), circulant-spectrum/FFT (T4), barrier-dodging form (T8).
- **[concrete conjecture] (testable next):** OCR forced-parity obstruction (T1); Ham-count as a
  loneliness monotone (T7).
- **[speculative] (directional):** OCF resonance accountant (T5), Collatz `GF(2)` homology (T6),
  RH explicit-formula trace (T9).
- **Not claimed:** any resolution. The deliverable is a *prioritised transfer menu* with a
  verified anchor and a clear next test (T1/T7).

**Artifacts:** `04-computation/redei_transfer_round_tournaments_s599f.py` (+`.out`), THM-402,
THM-406, companion `lrc-collatz-…-cancellation-family-s599.md`. Builds on THM-402 (round↔tight),
THM-404 (doubling involution), S577 (3-cycles), S595 (two-block), the repo Rédei/OCR/OCF/
even-graph machinery. New: **HYP-2160**.
