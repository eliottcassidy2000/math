# Live exterior choices: an absolute power saving, a physical gap, and a compression test

1. [Cross the prescribed target with both exteriors live](EXTERIORS.md).
   For every sufficiently large target, there are literal binary readers
   strictly on both sides within O(N^(2/5)). Each pair of low exterior
   phases is independently admissible, giving Ω(N^(1/5)) distinct address
   choices on either side at that error scale. This counts addresses,
   not distinct counts or isomorphism classes. The +2 and +4 synchronized
   seed changes remain physical and unchanged. One may also prescribe
   the final count modulo 9 in every phase cell (with a larger cutoff). The crossing uses a
   lattice tangent to a product level, not an average or a range census.
2. [Select an exact signed gap with one or two arc edits](SWITCH.md).
   Every even gap up to a dyadic V in magnitude has V/2 explicit,
   distinct anchors of scale Θ(V²), with a whole exterior still free.
   These anchors are spread at spacing Θ(V) on the displayed strand.
   This is a gap-first construction, not an absolute target theorem.
   [The contact bound](CONTACT.md) shows why a second independent small
   single-bit switch on the same word is not available.
3. [Four-vertex atom and recursive hidden choices](ATOMS.md).
   The unique two-free-arc square containing all four order-four classes
   produces an exponential, nonisomorphic ambiguity in every unmarked
   slot moment. Marked endpoints reveal only the atoms they contact.
   [Further links](ANALOGIES.md) include Price's half-angle equation,
   Chair44 exact descent, the Monster/Baby Monster filtration, and
   exact recoverability of the full address by two successive reads.
4. [A second route via a rational slope](LIVE.md). This reaches the same
   exponent by fractional-part crossings and provides an independent
   tradeoff for generating many exteriors within a single sublattice.

**Remaining step.** A near-target address and a switch of the necessary
size are not necessarily compatible. The former moves tangentially
along the product curve; fixing the latter ties the exteriors in a
different direction and makes the anchor grow monotonically. We need
a new selective contact or a third independently programmable
coordinate so that the exact gap and the absolute target can be solved
jointly. No coverage theorem or proof of Chan–Pak is claimed.

Checks: Python 3 `verify.py`, `check_atoms.py`, and
`check_exteriors.py`, `check_mod9.py`, and `decoder_live.py`. The previous [harmonic packet](../harmonic_reader_residual_20261116/README.md)
provides the signed reversal-defect bound and the initial decoder.
