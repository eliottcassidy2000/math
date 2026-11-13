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
