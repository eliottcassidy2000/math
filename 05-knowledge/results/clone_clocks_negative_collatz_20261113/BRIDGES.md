# Typed bridges to older repository work

klein, 2026-11-13. Connections are typed as AGENTS.md asks: source → target, map, preserved predicate, lost
data, sidecar, and the cheapest decisive test. Only the first bridge carries a mathematical transfer. The
others are recorded analogies with a test that was actually run. **No Jacobian, Hadwiger–Nelson, Hodge or
Collatz statement is claimed.**

## 1. The 2–3 axis (HYP-9339) — substantive

* **Source → target:** physical clone moves on Π (CLOCKS §3–4) → the multiplicative semigroup `⟨2,3⟩` acting on
  `H−1=60·2^k·3^j`.
* **Map:** dyadic growth `×2`, ternary growth `×3`, and the one-arc reversal `×3/2`. That last ratio is the
  Pythagorean fifth, so in log scale the reversal is the lattice step `(+1,−1)`.
* **Preserved:** `H−1` exactly, for every module pair. **Lost:** the module interiors, which are recovered by the
  P-form.
* **Sidecar:** HYP-9339's anchors and Størmer pairs. The shadowing exchange rates of COLLATZ §2 are exactly the
  balanced cycles `{−1}, {1,2}, {−5,−7,−10}` at `2→3, 4→3, 8→9`. The 462 ray sits on the Catalan pair `8|9`.
* **Test (run):** `verify.py` §1 and §3.

## 2. Planar Jacobian response fibres (THM-2230, inherited via the H17 history lane) — analogy with a run test

* **Source → target:** a two-module skeleton → its response tensor `c_{mm'}` in
  `H=Σc_{mm'}P_M(m)P_N(m')`.
* **Analogy:** THM-2230 computes a whole response fibre by target shears around a unit Jacobian. Here inserting
  the middle vertex `x0` acts on the tensor as the shear `(c_{12},c_{11},c_{32}) : (−1,1,30) → (0,1,60)`. It moves
  the dyadic defect and fixes the unit coefficient `c_{11}=+1`.
* **Preserved:** `c_{11}` under the `x0` shear. **Lost:** any Keller condition; none is transported.
* **Test (run):** the shear fixes `c_{11}=+1` inside the block family only. Across all two-clock lines `c_{11}`
  takes many odd values (Ψ has `−3`), and `−1` occurs in products and one genuine family. So the "unit
  coefficient" is not an invariant of two-clocks. The analogy's prediction failed twice, and it is kept only as
  the record of a failed transfer.

## 3. Hodge-type bidegrees — vocabulary only

The P-form is a finite Künneth-like decomposition into bidegrees `(m,m')`. Two-clocks are the tensors supported
in `{(1,1),(3,1),(3,2)}` (plus the `σ_2` term `(1,2)`). The C2-line exclusion concerns the support,
`(1,2)` absent, jointly with the value `−1` of the `(1,1)` component. This is a description of the observed data, not a transfer. The repository's earlier Hodge lane
(connecting homomorphisms, H17 history) is about lifting quotient loops, which these exact identities do not need.

## 4. Hadwiger–Nelson — shape coincidence, flagged as such

At the ternary root `u`, the two back-bundles of size two make two pairs of cyclic triangles. Each pair shares
the arc `a→u` (respectively `c→u`), and `p=q=2` is necessary for a two-clock in the block family. Two rhombi of
two triangles meeting at a vertex is the combinatorial shape of the Moser spindle. Here, though, `u` lies on the
shared edge rather than at a tip, and no map from colourings to path words is known. **Test (run):** bundle sizes
`(1,2), (2,1), (3,·)` destroy ternary shape (`scripts/block_family.py`). Nothing is inferred about the plane's
chromatic number.

## 5. Same-night H17 and multicore sessions

The [H17 exact-repair packet](../h17_exact_repairs_20261112/README.md) already generalized the H17 preparation
into reusable exact-target repairs and a closed recursion `C(q,s)` with dyadic factor `2^{k−1}`. Its target image
is sparse, of order `O(√X)`. The Π family is a second closed recursion, simpler and with a two-dimensional grid of
exact returns, but on the +1 anchor. Both packets land on the same boundary: exact reuse is achievable, while
target supply at the C2 anchor remains open.

## 6. The old dualities (rational/irrational, odd/even, positive/negative, addition/multiplication)

[The rational/irrational reflection](../../../07-reflections/the-rational-irrational-duality.md) (corrected
2026-11-13) frames proofs as a bounded census plus a sweep that certifies the census is non-vacuous. Each
duality meets this packet in a checkable place.

* **Finite certificate ⇒ all sizes (the census/sweep hinge).** A finite visit table, computed once with
  rigorous caps, certifies a two-clock for every `(j,k)` and every pair of modules (CLOCKS §2, §6). The finite
  object is what is verified; the exponential-polynomial identity is the sweep. The refuted constant-term
  hypothesis (CLOCKS §5) is the failure mode the reflection warns about: a census boundary was read as a law.
* **Odd/even.** Rédei oddness forces every clone clock to have even amplitude and an odd digit (Theorem 2 of
  CLOCKS). The even C2 digits `0, 2` (`3H`, `3H+4`) are exactly the ones clone clocks cannot produce.
* **Positive/negative.** The +1 and −1 anchors are exchanged by `x=−H`. The certified family Π runs the
  negative-integer Collatz moves physically, with one arc for the odd step (CLOCKS §4). The C2 line is the
  positive-side mirror, and it is the one still missing.
* **Addition/multiplication.** Ordinal sums multiply counts, so constant terms multiply too. This is what
  refuted the additive "constant never −1" reading: `(60·3^j−1)(2^k+1)` has constant −1 for a purely
  multiplicative reason. The surviving C2-line exclusion is stable under products, because odd digits make
  every product carry a nonzero pure-`2^k` term.
* **Rational/irrational.** Each clock family converges p-adically to a rational point, its anchor: dyadic rays
  2-adically to their offset, ternary clocks 3-adically to `−σ`. The irrational slope `log2/log3` governs only
  how the `(j,k)` lattice of values spreads on the log scale. Such families always have density zero, which is
  why they seed lines (CLOCKS §6) and do not cover the spectrum.

## 7. The attached briefs of the continuation (typed; owner-supplied, read as inspiration)

* **Wrap and split words** (Session23 report). `Wrap(U)=C3[K1,U,K1]` has `H=Π_2(U)+H(U)`, the canonical unit
  dyadic clock (`1+2^k`). The wrap law of CLOCKS §8.3 says that a dyadic module placed on a digit-1 clone clock
  always ends up in this *outer* position, `H−1=2^k·(α3^j−1)`. The split-word cores `3^t(1+2^{a+2t})` multiply
  by 3 through triangle letters. That is a different "three" from three-bin clone growth, matching the
  four-way distinction of the concurrent ternary-core grammar packet.
* **Controlled repair** (2026-10-19 research map: "preserve exactly the state needed by the next
  argument"). For clone growth that state is the complete visit table. The map's matching certificate
  `ν≥2b−E` and Session 197's matching number `ν_u` play the same role, a matching size bounding exchange
  capacity. Here `ν_u=2` caps the ternary support, and in every census digit-1 clock the crossing matching is
  realized by two cherries.
* **Metagraph briefs** (pure-blue atoms `K1, C3, C5`; `a_j=a_{j−1}+a_{j−3}`). No transfer: clock skeletons are
  near-transitive strong cores, and their counts are not the merged-node masses those briefs classify.
* **Quantum brief** (native moves, Paley frames). The one-arc Collatz move of CLOCKS §4 is a single-arc
  operation with an exact count law, `H↦(3H−1)/2`. That is the kind of named-move control the brief asks for, but
  no frame or measurement statement is made.
* **Square-sum panels** (prime-indexed graphs, "67 missing"). 67 is also a persistent exception of the
  near-join and C2-sweep spectra (THM-4865 (e)). The repository has already typed the 25-prime grouping as
  numerology; it is not redone here.

Two checks prompted by this list are recorded for completeness. `Π−v` is itself a digit-17 clock
(`H=60·3^m−17`, the most frequent digit), so the one-arc Collatz run of Π stops when the `v`-block has one vertex.
And THM-5133's hyperbolic kernel cannot be composed with a clock (CLOCKS §8.5).
