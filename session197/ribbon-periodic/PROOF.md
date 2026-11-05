# Why fixed append blocks do not give binary induction at span 3

Session 197. This packet attacks an all-target construction strategy rather than enlarging a finite target census. It proves a periodic-growth dichotomy, rules out every fixed-block binary-affine append rule at backward span 3, and records two explicit counterexamples to tempting stronger claims. It does **not** settle the full span-3 spectrum.

The definition and mask convention are those in `../ribbon/PROOF.md`: all backward edges in the supplied vertex order have length at most `w`; a new mask specifies its outgoing edges to the previous `w` vertices. The transfer and monotonicity used below were proved in that packet. The path-merging inequality used below was proved in `../ribbon-compression/PROOF.md`.

## 1. A reset obstruction

Append `w` vertices that each lose to every older vertex, namely the word `0^w`. This preserves the Hamiltonian-path count `h`.

At the resulting frontier, all original vertices have been forgotten. No appended vertex points to any original vertex. Therefore any surviving path-cover state must have selected a single Hamiltonian path on the original vertices: two original components would retain two forgotten starts, and cannot be repaired. Once that old path is selected, the possible frontier choices on the appended sinks are independent of its internal identity. Consequently the resulting vector is

`h q_w`,

where `q_w` is one fixed nonnegative integer vector, independent of the original nonempty tournament.

Every later transfer is homogeneous with integer entries. Thus **every future count is divisible by h**. In particular a forward-sink reset followed only by appending vertices cannot give `2h-1` or `2h+1` when `h>1`.

The same argument applies to any proposed reset that maps all inputs of count `h` to `h` times a fixed integer response profile. A count-preserving normalization that erases every independent response also erases the possibility of unit offsets. This does not rule out resets retaining a genuinely independent response.

## 2. A false boundary-deletion proxy

A tempting surrogate for the full frontier is the collection of counts obtained by deleting subsets of the last three vertices. Their gcd is **not** a divisibility invariant for future extensions.

At width 3, the word `0001500` has `H=15`. Deleting the subsets of its last three vertices in binary-mask order gives

`15, 5, 15, 5, 15, 5, 15, 5`.

All are divisible by five. But appending mask `1` gives word `00015001` with `H=91`. Thus these eight deletion responses do not capture every future contribution. The counts 15 and 91 already show that the gcd of *all* continuation counts is one. No conclusion of eventual coverage follows from that gcd being one.

All nine displayed counts are checked by independent subset Hamiltonian-path dynamic programming, not just by the frontier transfer.

## 3. A general periodic-growth dichotomy

Let `P` be any fixed prefix and let `W` be a nonempty extension block of length `L`. Let `T_k` be the width-`w` tournament represented by `P W^k`. Assume **L>w**.

**Periodic-growth theorem.** Exactly one of the following alternatives applies:

1. The infinite periodic suffix has no directed triangle. Then `H(T_k)` is eventually constant.
2. The suffix contains a directed triangle. Then, for a constant `C` depending only on its location and on `W`,

   `H(T_k) >= 3^(k-C)` for all sufficiently large `k`.

These alternatives concern the actual graph, not the eigenvalues of an arbitrary matrix representation.

### Transitive suffix case

A triangle-free tournament is transitive. Consider the first `w` suffix vertices, a set `S`. Every suffix vertex more than `2w` positions after the prefix lies after all of `S` in the suffix's transitive ordering, because its edges from `S` are forced forward.

Let `U` be the initial segment of that transitive ordering ending at the last member of `S`. It contains all of `S`, is contained in the first `2w` suffix vertices, and is fixed once those vertices have arrived. Every remaining suffix vertex is after `U` in the transitive ordering. It also has natural index more than `w` beyond the old prefix, since `U` contains `S`, so every old-prefix vertex points to it.

Thus, after at least `2w` suffix vertices, the whole graph is an ordered join

`T[P union U] -> transitive remainder`.

Its Hamiltonian-path count is the fixed value `H(T[P union U])`. This proves eventual constancy without assuming the prefix is transitive or ignoring its hostile short backward edges.

### Triangle suffix case and the overlap issue

Fix a suffix triangle with labels `a<b<c`. Translating its vertices by `jL` preserves all its orientations. Any directed triangle in a span-`w` ordering has one of these forms:

- `a -> b -> c -> a`, with `c-a<=w`;
- `b -> a -> c -> b`, with `b-a<=w` and `c-b<=w`.

Hence its diameter is at most `2w`.

Since `L>w`, the first orientation has disjoint translates. In the second orientation, overlaps between translates can only occur at consecutive translations: `2L>2w`, so more distant translates cannot intersect. Neither short gap can equal `L`. Thus overlap is possible **only** when `c-a=L`, identifying the last vertex of one triangle with the first vertex of the next.

If no overlaps occur, `m` translates are disjoint directed triangles. The deterministic path-merging inequality gives at least `3^m` Hamiltonian paths on their union, and induced-subtournament monotonicity gives the same lower bound for the whole graph.

If overlaps occur, write `b=a+d`, with `d<=w` and `L-d<=w`. The union of `m` translates has the ordered vertices

`a, b, a+L, b+L, ..., a+(m-1)L, b+(m-1)L, a+mL`.

Every adjacent arc points backward, by the triangle orientation. Every **nonadjacent** pair has index distance at least `L>w`, so its orientation is forced forward. There are no unspecified cross-arcs here. The induced graph is exactly the width-1 chain `R_(2m+1)` from the companion packet.

Its counts obey `h_0=h_1=h_2=1` and `h_n=h_(n-1)+h_(n-2)+h_(n-3)`. Monotonicity of this sequence gives `h_(n+2)>=3h_n` for odd `n>=1`: check `n=1`, and for `n>=3` subtract the two recurrences to reduce the inequality to `h_(n-1)>=h_(n-3)`. Therefore `h_(2m+1)>=3^m`.

In either case, the first `k` blocks contain at least `k-C` translates for a fixed `C`; the initial triangle's offset only changes this constant. This proves the theorem.

The threshold matters. Repeating the single reversed-adjacency mask has the tribonacci growth of `R_n`, whose per-vertex growth is below two, so it would be false to assert a universal factor-three bound for all nonempty periodic blocks.

## 4. The short-block certificate at width 3

The general theorem handles every block of length at least four when `w=3`. To handle lengths one, two and three, `representation.json` supplies an exact integer output representation derived from the 119-state frontier:

- a `119 x 62` matrix `B`;
- eight `62 x 62` integer matrices `R_a`;
- a length-62 integer column `f`.

If `M_a` is the original frontier matrix and `F` its one-component acceptance column, the verifier checks every entry of

`M_a B = B R_a`, and `B f = F`.

Consequently, if the current frontier row is `v`, every extension word `u` has output

`v B R_u f`.

No minimality assertion or lattice-generation claim is needed for this argument. The representation can have signed coordinates; the graph counts remain those of the original nonnegative frontier.

The verifier checks `det(2I-R_W) != 0` for every word of lengths 1–3. It suffices to check cyclic rotation classes, since cyclically rotated products of square matrices have the same characteristic polynomial. There are exactly 8, 36 and 176 such classes, respectively: **220 classes covering all 584 nonempty words** of these lengths.

Every determinant is certified nonzero modulo at least one of the primes `3,5,7,11,13`. The number whose first successful modulus is each prime is:

- Length 1: `3,4,1,0,0`.
- Length 2: `31,5,0,0,0`.
- Length 3: `112,56,4,0,4`.

These are exact modular computations; no floating-point eigenvalue calculation is used. The verifier rebuilds the original frontier and verifies the intertwining identities before checking the determinants.

## 5. No repeatable fixed-block binary-affine rule

**Theorem.** Fix a span-3 prefix `P`, a nonempty append word `W`, and any constant `c`. The sequence

`H_k = H(P W^k)`

cannot be eventually nonconstant while satisfying

`H_(k+1) = 2 H_k + c`

for all sufficiently large `k`.

**Proof for length at least four.** Any nonconstant sequence with that recurrence has the form `A 2^k-c`, with `A>0` because counts are nonnegative. The periodic-growth theorem gives either eventual constancy or a lower bound `3^(k-C)`, contradicting that form.

**Proof for lengths one to three.** Let `A=R_W`. The exact certificate says two is not an eigenvalue of `A`. By Cayley–Hamilton, every output sequence `v B A^k f` satisfies the characteristic-polynomial recurrence. If a tail were `alpha 2^k-c` with `alpha != 0`, applying this recurrence would give

`alpha p(2) 2^k - c p(1) = 0`

for every sufficiently large `k`. Since `p(2) != 0`, this is impossible. Zero eigenvalues and finite transients do not affect this argument. QED.

Thus an indefinitely iterable binary construction cannot use the same fixed nonempty append word for its `2h+1` branch, or for a `2h-1` branch started at `h>1`, while keeping span three. This rules out a concrete all-target induction strategy.

**Scope.** It does not rule out words selected from the changing full response state, controllers with multiple states and longer cycles, modifications of old edges, a new vertex ordering after each step, or broader physical gluing constructions. It is not an obstruction to span-3 completeness itself.

## 6. A false stronger spectral conjecture, with a physical closed form

It would be tempting to strengthen the short-block certificate to “no span-3 block has eigenvalue two.” That statement is false.

Take `W=001476`, of length six. `eigenvalue2-counterexample.json` provides a nonzero integer vector `z` with

`R_W z = 2z`.

Moreover, the mode actually contributes to physical Hamiltonian-path counts. For every `k>=1`,

`H(W^k) = (7 * 75^k + 66 * 2^k) / 73`.

The initial values are `9, 543, 40461, 3034047, ...`. The verifier checks the initial two values and the exact row-vector identity that yields the recurrence

`H_(k+2) = 77 H_(k+1) - 150 H_k` for `k>=1`.

The displayed closed form follows. This is not merely an unobservable or unreachable matrix eigenvalue. The factor-75 contribution is exactly why the presence of a doubling mode does not provide a repeatable physical doubling construction. The periodic-growth theorem remains valid.

## 7. Reproduction and remaining question

Run from the repository root:

`python session197/ribbon-periodic/verify.py`

The verifier uses Python 3 only, with assertions enabled. It imports the audited frontier from `../ribbon/transfer.py`; no third-party package or symbolic algebra system is needed. It checks the representation, all 220 short-block classes, every nonempty frontier-basis reset row, the boundary-deletion counterexample, the explicit doubling eigenvector and physical mixed-exponential recurrence, all overlapping-triangle geometries tested in the recorded finite range, and hostile-prefix transitive-tail cases.

`VALIDATION.json` records the exact checks. The proofs, rather than the geometry samples, establish arbitrary-order and arbitrary-block validity.

The main open task remains a physically valid state-dependent realization principle that reaches all positive odd targets except the proposed exceptional values. Neither primitiveness of a small response list, a count-only reset, nor a fixed binary append block supplies that principle.
