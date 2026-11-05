# The charge-minus-two interface misses thirteen at width three

Session 197. This is an all-order bounded-output certificate for a particular rooted interface, not a global Hamiltonian-path exclusion.

## Statement and scope

Let a nonempty tournament have an ordering of backward span at most three, and let its designated root be the **last vertex of that ordering**. Write `h=H(T)`, `e=e_root(T)`, and `R=R_root(T)` for the one-path and ordered two-path-cover responses defined in `../ribbon-terminals/PROOF.md`.

**Certified theorem.** Among all such ordered rooted tournaments, with no bound on their order, the counts `h<=13` compatible with

`-2h+e+R=-2`

are exactly `3,5,9,11`. In particular count 13 cannot occur at this charge.

The statement includes arbitrary and strong tournaments. It allows every span-three ordering with its own last root; it does not assert an exclusion for a specified root that is not last in a span-three ordering. The literal concatenation compiler in `../ribbon-primitive/PROOF.md` always returns its root last, so the theorem applies to that compiler even if its three cells are replaced by arbitrary cells, or their choices depend on the full current state.

Count 13 itself is attainable at backward span two. For example the width-two word `00113` has count 13 and last-root profile `(13,2,18)`, whose charge is `-6`. Thus this is a failure of retaining the particular charge `-2`, not a scalar spectrum obstruction.

## 1. Exact frontier readout of the rooted responses

Use the independently reviewed path-forest transfer from `../ribbon/transfer.py`. Its states are unordered collections of directed path components represented by their ordered endpoint pairs. The retained labels are `0,1,2`, with label 2 the newest vertex; forgotten endpoints have label `-1`. A state may have at most one forgotten start. Its coordinate counts the corresponding spanning directed path forests.

For a nonempty word, the three responses are nonnegative sums of frontier coordinates:

- `H`: states having exactly one component;
- `e`: one-component states whose end is label 2;
- `R`: states having either one or two components and containing a component whose start is label 2.

For `R`, the component starting at the root is the nonempty second path `Q`. When there is one component, the first path `P` is empty; when there are two, the other component is `P`. This assignment is unique, so no extra ordering factor occurs. Such a cover has at most one forgotten start, namely the possible start of `P`; therefore the transfer discards none of the desired covers. Conversely every state counted by this readout gives exactly a cover of the prescribed type.

These formulas are checked against separate graph-level subset dynamic programming for every ordered width-three tournament through order six.

## 2. Why the signed charge is safe under saturation

Set `M=13` and `C=2M+1=27`. Clip every individual frontier coordinate at `C` after every transition. All transfer coefficients are nonnegative integers, so componentwise clipping commutes with every transition: for nonnegative integers `a_i,v_i`,

`min(C,sum a_i v_i)=min(C,sum a_i min(C,v_i))`.

Thus a clipped state is exactly the coordinatewise clipping of the true state for its word.

Discard a state if the sum of its clipped one-component coordinates exceeds `M`. This test is exact for `H<=M`: a clipped summand would itself be at least `C>M`, so every retained `H` summand is exact. Conversely a true `H>M` cannot be reduced to at most `M` by clipping. Since `e` is a sub-sum of the `H` coordinates, every retained `e` is also exact and satisfies the physical inequality `0<=e<=H<=M`.

**No uniform small bound on `R` is assumed.** The relevant conditional physical bound is instead obtained from the proposed charge itself:

`e+R=2H-2`, with `e>=0`, implies `R<=2H-2<=2M-2<C`.

For any nonnegative coordinate sum `R`, its clipped readout obeys

`R_clipped >= min(R,C)`.

Indeed if any summand is clipped, that one summand already contributes `C`; otherwise the entire sum is exact. Consequently a clipped readout can pass `e+R_clipped=2H-2` **only if every coordinate contributing to `R` is exact**. If clipping could change `R`, the clipped sum is at least `C`, larger than the entire right-hand side. Large positive `R` cannot cancel against a clipped negative contribution, because the only negative term, `-2H`, is exact and has magnitude at most `2M`.

The charge test therefore has neither false positives nor false negatives. We do not subtract two independently saturated unbounded quantities, and we do not claim that all `R` readouts in the closed set are exact.

## 3. Finite closure proves an all-order exclusion

There are 119 frontier coordinates at width three. Starting from the empty initial frontier, follow all eight mask letters with the above saturation, retaining precisely those states with `H<=13`. Breadth-first closure terminates at exactly **2,382** distinct vectors.

The certificate contains every vector in this set. The checker verifies that it contains the initial vector and that every retained transition from every listed vector returns to the listed set. It also independently regenerates the reachable closure and compares it with the certificate. Thus induction on word length includes every word whose count is at most 13. Pruning loses no later low-count word: Hamiltonian-path insertion injects paths of an induced tournament into its one-vertex extension, so `H` never decreases when a vertex is appended.

Scanning the exact-safe charge readout on the closed set gives only counts `3,5,9,11`. Literal witnesses are:

- width-three word `002`: `(h,e,R)=(3,1,3)`;
- `0003`: `(5,1,7)`;
- `04217`: `(9,1,15)`;
- `00036`: `(11,3,17)`.

All four profiles are also verified directly from their adjacency graphs. The finite object here is a closed state space, not an enumeration to a chosen graph order. The conclusion is therefore valid at every order.

## 4. Consequence for the primitive-state approach

The three-cell family in `../ribbon-primitive/PROOF.md` showed that primitive state and primitive observable continuation content can coexist with a sparse scalar image. The present theorem shows a stronger limitation of the particular preserved charge: **no** span-three last-root construction that always maintains `-2H+e+R=-2` can produce every permitted odd target. It misses the valid target 13 irrespective of gadget-size bounds or control rules.

This does not rule out changing charges, multiple interfaces, a different root-order convention, or a different construction. Those would require fresh physical transfer laws or separate analysis.

## Reproduction and dependencies

Run `python session197/ribbon-charge/verify.py`. The verifier is standard-library Python. It pins these reviewed or separately frozen dependencies:

- `../ribbon/transfer.py`, SHA-256 `2ef89174390c503494d4e88e81f1b7ee149f985b6746f2912430ae88e46a289d`;
- `../ribbon-terminals/verify.py`, SHA-256 `9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e`, used for independent graph profiles.

`closure.json` is the full vector certificate. `VALIDATION.json` is the deterministic checker receipt. The proof does not depend on the larger exploratory bound 101.
