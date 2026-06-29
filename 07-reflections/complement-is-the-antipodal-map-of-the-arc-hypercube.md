# Complement is the antipodal map of the arc-hypercube — so the R-eigenspace split, the merged metagraph, and Borsuk–Ulam are one object

*klein-2026-06-29-S1. Synthesizing a Klein-four-group observation about the n=4 classes with the project's live organizing principle (HYP-3538: reversal R has ±1 eigenspaces; the −1 part is the signed obstruction). Companion to THM-584.*

## The one-line synthesis

A labeled tournament is a vertex of the hypercube `Q_d`, `d = C(n,2)`, one bit per arc. Three things
that the project has been treating as separate are then **literally the same map**:

> **complement** `T -> T^op` = **reverse every arc** = **flip every bit** = the **antipodal map**
> `x -> x ⊕ 1` of `Q_d` = the map **Borsuk–Ulam is the theorem about**.

Everything else falls out of this identification, because the antipodal map of `Q_d` acts as
`(-1)^k` on hypercube level `k`, and `S_n` (vertex relabeling) permutes arcs **within** a level. So
the `S_n`-quotient — the iso-class metagraph — inherits a clean grading, and the reversal `R` is just
**level parity**.

## Why this is the right frame for HYP-3538

The owner's principle says every two-index phenomenon is the `±1` eigenspace split of `R`, with the
`eps=+1` (R-even) part the SOS/Brouwer-provable bulk and the `eps=-1` (R-odd) part the signed
obstruction. THM-584 makes that split *explicit and computable on the metagraph*:

- **R-even = even hypercube levels** `k = 0, 2, 4, …`. Contains the Perron mode (level 0, the bulk).
- **R-odd = odd levels** `k = 1, 3, 5, …`. The signed residual.
- A hard fingerprint: R-even eigenvalues `≡ d (mod 4)`, R-odd `≡ d−2 (mod 4)` (verified n=3..6).

This is the metagraph-side twin of THM-583's witness-side result. There, the odd index `f`
(palindromic Hamiltonian-path count) is a count on the half-system, **lossless iff you keep `φ`** —
because `φ` *is* the `eps=-1` coordinate. Here, the **merged metagraph `G_n/Z_2` is the R-even
projection** (its size `V_merged = (A000568+SC)/2` is exactly `dim V_+`), so folding by complement
**discards the odd-level coordinate**. HYP-2685 warned that the half/complement quotient "destroys the
discarded side's independent bit data, so it is not a full encoding." THM-584 names that destroyed
bit: it is the **odd-level (R-odd) eigenspace coordinate**, supported precisely on the
non-self-complementary classes.

So the witness side (`f`, half-system, keep `φ`) and the cap side (`M = M_even ⊕ M_odd`) and the
metagraph side (`A = A_even ⊕ A_odd`, even/odd levels) are three faces of **one antipodal involution
on three different `S_n`-quotients of the arc-cube**.

## The Klein four-group is the seed

At n=4 the whole structure is visible on a 2×2 square (`Q_6 / S_4` = 4 classes = `V_4 = (Z_2)^2`):

- two bits with meaning — `x` = "source destroyed", `y` = "sink destroyed";
- `R` = complement = the **coordinate swap** `x ↔ y` (reversing arcs swaps sources and sinks);
- the **R-even coordinate** is `u = x+y = #boundary defects ∈ {0,1,2}` — the three merged classes
  `{T, [±], S}`, i.e. `G_4/Z_2`, i.e. the even-level grading;
- the **R-odd coordinate** is `w = x−y = source−sink ∈ {−1,0,+1}` — zero on the SC diagonal `{T,S}`,
  `±1` on the swapped NS pair `{+,−}`, and `w -> −w` under `R`.

`V_4 ⋊ ⟨swap⟩ = D_4` (dihedral of order 8): source-kill and sink-kill generate `V_4`; complement is
the reflection that swaps them. The project's recurring **dihedral** motif
([[the-dihedral-recursion-existence-is-even-witness-is-odd]]) and its **two order-2 structures**
([[two-order-two-structures-parity-and-descent]], reflection vs doubling) both sit here: reflection =
the `R`-swap; the two generators `x, y` are the two boundary involutions; and the n→n−2 doubling is
the descent between levels.

## A clean prediction that FAILED — and why the failure is informative

It was tempting to think the `eps=-1` coordinate is also what makes *compression* hard. Concretely:
axis-aligned face-compression of the tiling cube is information-tight at n=4,5 (you can drop down to
`⌈log_2 A000568⌉` free tiles) but acquires a **+1 overhead at n=6** (needs 7 free tiles where
information allows 6). Hypothesis (HYP-3539): the overhead is the cost of *resolving the R-odd
coordinate* — telling the two members of an NS pair apart.

**Refuted, exactly.** Covering only the complement-**merged** classes (R-even, where you never have to
distinguish `+` from `−`) *also* needs a 7-face at n=6. So the overhead lives in the **R-even bulk
itself** — it is a covering-geometry deficiency of the merged metagraph, not a parity obstruction.
This is the *consistent-with-HYP-3538* outcome, read the other way: the bulk is the R-even side, and
the bulk is exactly where the hardness turned out to be — the obstruction here is "Brouwer-flavored"
(packing/covering), not "Borsuk–Ulam-flavored" (signed). The `eps=-1` coordinate is small (dim 22 of
56 at n=6) and cheap to carry; the expensive thing is packing the `eps=+1` bulk.

## Where this points

1. **The metagraph spectrum is a level-multiplicity sequence.** `A`'s eigenvalue `d−2k` has
   multiplicity = `dim(S_n-invariants at level k of Q_d)`. n=6 R-even shows `−1` with multiplicity 13,
   `3` with 10, … These multiplicities are a combinatorial sequence worth identifying (likely an
   orbit-counting / Burnside expression per level). Open as HYP-3540.
2. **Borsuk–Ulam on `Q_d`.** Since complement = antipodal, any `S_n`-equivariant odd map / signed
   count on tournaments is a Borsuk–Ulam statement on `Q_d` modulo `S_n`. Rédei's "`H(T)` is odd" and
   the project's hunt for an *odd index* (the missing topological key in
   [[two-order-two-structures-parity-and-descent]]) should be looked for as a **nonzero degree of the
   antipodal action on a level-parity-graded complex** — the R-odd block is where it must live.
3. **The mod-4 law is a free audit.** Every metagraph eigenvalue computation across the project can be
   checked against `≡ d` (R-even) / `≡ d−2` (R-odd) mod 4.

The mathematics kept pointing past the Klein group: a 2×2 multiplication table turned out to be the
shadow of the antipodal map on a 15-dimensional cube, and the project's whole even/odd, Brouwer/
Borsuk–Ulam, SC/NS, merged/unmerged vocabulary turned out to be **the level parity of that one map.**

See [[the-pm-one-eigenspace-of-reversal-is-the-whole-split]] (HYP-3538), THM-584, THM-583,
[[merged-metagraph-invariants]], [[everything-is-the-triangle]]. Hypotheses: HYP-3539 (refuted),
HYP-3540 (open).
