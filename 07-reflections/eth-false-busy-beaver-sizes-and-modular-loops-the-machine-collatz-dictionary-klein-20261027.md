# "ETH is false", busy-beaver sizes, and modular loops: the machine–Collatz dictionary, what transfers, and what does not

*klein-2026-10-27 (second checkpoint). A REFLECTION, written to think rather than to compute, on the owner's request:
"how to leverage the fact that the exponential time hypothesis is false and how that relates with our work on Turing
machines and classifying famous open problems in number theory as halting questions; how a Turing machine of a given
size is isomorphic to the way the different integer series loop through modular arithmetic". The refutation of ETH is
the EXTERNAL PREPRINT CLAIM of [the external note](../05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md)
section 3 (unverified; no public confirmation found). Everything below is typed; nothing changes a status. The
classification it refers to is [PROBLEM_CLASSIFICATION.md](../05-knowledge/results/analogy_reduction_atlas_20261025/PROBLEM_CLASSIFICATION.md)
(level + machine size, not a completeness class) with its inherited sources (procgen atlas, THM-4622, camion–busch,
reflection S599).*

## 1. Two different `n`

The preprint's content is: counting the satisfying assignments of a bounded-width CNF on `n` Boolean variables costs
`2^{o(n)}` instead of `2^{n}`. Its `n` is the length of a binary certificate, and the verifier is a *local, sparse*
constraint system of size `O(n)`. The cryptids of the atlas (Goldbach as BB(25), Erdős ternary as BB(15), weak
Collatz as BB(14), Antihydra as BB(6)) are `Π⁰₁` sentences whose counterexample is an integer `N` found by
searching the integers themselves: the search space is unary in the bit-length `log N`, and the verifier of one
candidate (primality, digit patterns, an orbit) is a circuit of size polynomial in `log N` that is not local in
the bits. Three consequences, each a NO TRANSFER:

- Bounded search "is there a counterexample below `N`" is, in the bit-length, a `Σ₂`-type question (`∃ N' ≤ N ∀ p`),
  not an NP question; Toda's theorem does let counting reach the polynomial hierarchy, but its reductions multiply the
  number of variables polynomially, and `2^{o(polylog N)}` is larger than the brute-force `N`. Nothing becomes
  `N^{o(1)}`.
- The Cook–Levin tableau of a verifier running in time `t` has `Θ(t)` variables, so a subexponential #SAT counter
  helps only verifiers of *linear* size in the certificate; number-theoretic verifiers are polynomial in the
  bit-length and the certificate is the number.
- The levels do not move: a `Π⁰₁` sentence stays `Π⁰₁`, its busy-beaver size is a description-length coordinate,
  and no entry of the atlas cites ETH or SETH (grep: zero). ETH itself quantifies over all algorithms (in its
  standard form it is `Σ⁰₄`: `∃δ ∀A [A total and correct ⟹ ∃x slow]`), far above the `Π⁰₁`/`Π⁰₂` rows, and its
  refutation is the only resolution in this neighbourhood that comes as an **explicit witness with a proof** (an
  algorithm plus its analysis) rather than as a universal statement. That is a remark for the atlas, not a lever.

## 2. The dictionary: a machine of size `k` is a modular loop with `O(k)` branches

The owner's intuition is a theorem family. A generalized Collatz map is `g(x) = (a_i x + b_i)/d_i` for
`x ≡ i (mod m)`; Conway (1972) showed such maps simulate every computation through FRACTRAN (counters as prime
exponents, one fraction per instruction), and Kurtz–Simon (2007) that the generalized Collatz problem is
`Π⁰₂`-complete (CITED in the camion–busch note and the procgen atlas). Michel's busy-beaver work shows the converse
direction in the small: the champions and the cryptids *are* Collatz-like maps (Antihydra: `a ↦ ⌊3a/2⌋` split by
parity, the repository's own reading `⌊3a/2⌋ = −T(−a)` with `T` the Mahler carry map). The correspondences, each
exact up to constants:

| Turing machine | modular loop |
|---|---|
| a state (with head direction) | a residue class `i (mod m)` with its affine rule |
| a sweep across the tape (a macro-step) | one application of `x ↦ (a_i x + b_i)/d_i` |
| unary blocks / binary counters on the tape | the integer `x` (or a counter tuple packed as `∏ p_j^{c_j}`) |
| halting | reaching a residue class without a rule (the stop class) |
| cycler / translated cycler (regular deciders) | an eventually periodic residue sequence: a cylinder set `mod m^k` |
| bouncer / halting segment / FAR (counting deciders) | a branch word with a fixed Parikh pattern: a semilinear invariant |
| cryptid | a map whose stop question is a `Π⁰₁` number-theoretic sentence (Antihydra: `Σ(2v_i − 3) = −6`) |
| `k` states | `O(k)` branches (FRACTRAN ↔ Minsky ↔ Turing, linear overhead each way) |

So the atlas's size column has a twin: **the number of affine branches of the modular loop that encodes the
problem**, and Collatz itself is the two-branch member of a `Π⁰₂`-complete family. The sentence "a machine of size
`k` is isomorphic to how integer series loop through modular arithmetic" is exactly this table; what it is *not* is
an isomorphism of *particular* machines to *particular* sequences: a given `k`-state machine corresponds to one
map with one seed, and the "series" is its orbit.

## 3. Three levels of invariant, and where the ETH algorithm lives

Deciders for small machines, and arguments about modular loops, come in three strengths.

1. **Regular (automatic).** Closed tape languages for machines; for a map on integers, conditions on `x mod m^k`,
   i.e. the cylinder sets. Terras's theorem is the statement that for Collatz the `k`-step branch word is a
   bijective coding of `Z/2^k`: every regular condition on the seed is exactly a condition on a prefix of the orbit,
   and nothing more. This is *why* no CTL-type argument can settle Collatz: the regular invariants are the parity
   vectors, and they are a full shift.
2. **Semilinear (Parikh).** Linear conditions on the letter counts of the branch word: the Steiner–Eliahou cycle
   identity `T log 2 − d log 3`, Antihydra's halting criterion `Σ(2v_i − 3) = −6`, the repository's valuation-sum
   "prices" and "debts". The bouncer-type deciders live here. The interaction of a semilinear condition with the
   regular structure is the whole difficulty: which cylinders contain seeds whose orbit's Parikh vector hits the
   line. That is a *discrepancy* question (how the Archimedean size and the 2-adic residues correlate), the same
   shape as the LRC `disc_v` term and the "cancellation family" of reflection S599.
3. **Beyond.** The actual cryptids.

The ETH algorithm is a level-2 device over `F₂`: its engine is the *shared linear summary* (groups of variables
communicate through integer measurements `y = C x` of bounded size, enumerated and multiplied) and an exact count
recovered modulo `2^{ℓ}`. It is entirely on the 2-adic side. The number-theoretic halting questions sit at the
2-adic/Archimedean interface, which no `mod 2^{ℓ}` computation sees; this is the same boundary the repository met
as "the Fourier gap is transforms' blindness to the atom" and "conjugation moves content into the boundary term".

## 4. Two tool identities worth keeping (IDENTITY of tools, not of statements)

- **The replacement interpolant is a Bonferroni truncation.** The preprint's `f_I(z) = (−1)^{A_0} Σ_{j<R} (−2)^j
  C(z − A_0, j)` is the truncation at `R` terms of `(1 − 2)^{z−A_0} = (−1)^{z−A_0}`; LRC(14)'s depth-`D`
  Bonferroni certificate `B_D` is the truncation of `(1 − 1)^z = [z = 0]` at `D` terms (LEM-020: exact iff
  `maxmult ≤ D`). Same family (`t = −2` versus `t = −1`); the preprint pays for the window `R` by a binomial
  concentration bound and measurement rows where LRC pays by a multiplicity cap.
- **Modulus amplification is the Hensel idempotent.** Beigel–Tarui's `p_ℓ(Y) = (1−Y)^ℓ Σ_{j<ℓ} C(ℓ+j−1, j) Y^j ≡
  1_{2|y} (mod 2^ℓ)` holds because `Y^ℓ` vanishes `mod 2^ℓ` on even `y` and `(1−Y)^ℓ` on odd `y`: the two
  2-adic neighbourhoods, the same trick as the repository's lifting-the-exponent clocks (`2^{3^n} ≡ −1 mod 3^{n+1}`).
  Reading the exact count as a residue `mod 2^{n+1}` is an Ax–Katz-flavoured use of 2-adic valuations of point
  counts.

## 5. What would transfer if the preprint stands

Counting the solutions of sparse polynomial systems over `F₂` with `O(n)` variables and *high treewidth* in
`2^{o(n)}`: tournaments with forbidden local substructures counted in `2^{o(n²)}` over the `C(n,2)` arc bits (a
subexponential census of tournament Ramsey-type classes, practically useless at the orders where `gentourng` works);
zero counts of sparse Boolean functions; fixed points of high-treewidth sparse maps. The modular loops of number
theory are the opposite case: the carry chain of `3x+1` in binary has bounded treewidth, so cylinder counts are
already transfer-matrix cheap (Terras), and lattice dynamical systems (one-dimensional automata such as Rule 30 on a
ring) count their periodic points by transfer matrices in polynomial time; the preprint's bound is worse there.
**The gain is exactly on the expander-like side of the 2-adic world, which is the side the repository's problems do
not occupy.**

## 6. What to do next (no computation yet)

1. Add the branch-count coordinate to the classification ledger beside the busy-beaver size, with the FRACTRAN
   constants, and the remark that ETH is a `Σ⁰₄`-type sentence refuted by an explicit witness.
2. State the three invariant levels of section 3 as a definition and locate every atlas row by the strongest level
   at which its known partial results live (regular: Terras-type; semilinear: Steiner–Eliahou/Antihydra; beyond).
3. Cross-reference the two tool identities of section 4 from LEM-020 and from the 2-adic clock notes.
4. Keep the verdict: the refutation of ETH, true or not, is not a lever on any `Π⁰₁` or `Π⁰₂` entry of the atlas.
