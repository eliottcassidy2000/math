# Independent review: full binary-mirror bin operator

**PASS, no correction.** The all-order full-profile formula, its empty boundary, the positive delta/Delta form and the explicit ratio-collision example are valid. No total-nonnegativity or primitive-content assertion beyond the stated scope is inferred.

## Frozen binding

Reviewed source: `session197-orbit-wildcard/mirror-operator/`.

- `PROOF.md`: `6bb9cb49dc5ac2c3a3b69082b714f02d75c31d167c7bb174f1dcbc6d5a79fb9b`.
- `verify.cpp`: `aeb8ede9ddeb7ed403bc9ae45516b5110fd0366cd332056ab9f5acaa1c8cbb92`.
- `receipt.txt` and `receipt-debug.txt`: `0d6755fdff8e737a2cac778ae437704147cf4c793cb84208849bf8cf3f8f2bec`.

The source was copied unchanged to `source.cpp` and replayed in optimized and debug modes. Both receipts exactly reproduce the frozen output. The named universal context and endpoint reciprocity dependencies are unchanged; the reviewed binary-boundary proof supplies the inherited strongness, symmetry and Hamilton-cycle-count statements for the explicit family.

## All-order proof review

The only edge allowing a directed cycle to move backward among the ordered blocks I,T,O is the reversed edge O_a->I_a. Hence any cycle meeting several blocks uses that edge, and a vertex-disjoint collection contains at most one such cycle. A crossing cycle cannot have two T-runs because this would require a second use of that edge. Its decomposition into an I-path starting at I_a, an optional T-path, and an O-path ending at O_a is unique after anchoring at the reversed edge.

The parity filter has the correct sign. For side path total length ell and input-path size u, oddness means ell+u is odd. The factor 2x times its indicator is x[1-(-1)^(ell+u)]. Summing input subsets and their Hamilton paths against residual bin profiles gives f(x+1) and f(x-1), respectively. This works when the input is empty as well as nonempty. The unique invalid empty-input joining configuration has both side paths singleton and formal length2; its even parity removes it. All other empty-middle joins use an actual I->O edge. Thus the displayed converse construction is correctly restricted to odd total length.

For each side, residual parity gives sum_P (-1)^|P| B_(S-P)(x)=(-1)^|S| E_S(-x). The two sides square that sign, so the negative term is exactly -x E_S(-x)^2 f(x-1), with no extra parity sign. The no-crossing contribution B_S^2 f counts the three disjoint internal collections once. This proves the complete operator, not merely a positive-clone fit.

The scalar specialization retains -e^2 f(0), using the reviewed endpoint reciprocity sign. Thus f=1 gives h^2+q^2-e^2, whereas every nonempty input has f(0)=0. The zero boundary is neither omitted nor replaced by a positive-length extrapolation.

Splitting E^2 into its even and odd parts U,V yields xU delta+xV(Delta+2I). Therefore P=B_S^2+2xV, Q=xU and R=xV are exactly the stated nonnegative-coefficient polynomials. This is entrywise-positive operator data through the inherited calculus, not a statement that arbitrary resulting matrices are TN.

The explicit side S_k has exactly the stated backward-edge odd cycles and endpoint paths. Substitution and polynomial expansion give the Q_k formula and its B(2) expression. For Q_3 the profile is x^10+16x^8+64x^6+64x^4+16x^2, with evaluations161,10304,216009 at1,2,3. Its pair (H,B(2)) equals161 times the transitive-six pair, but its third evaluation is not proportional. Homogeneous linear pair-dependent coefficients therefore scale by161, preserving cancellation and changing integer content. The restriction to pair-dependent contexts is essential and correctly stated.

## Independent controls

`check.py` imports no source or research helper. Its actual bin profiles use directed odd cycles counted by least-vertex-anchored subset words, followed by a disjoint-cycle packing recursion. This differs from the source's Hamilton-path cover/falling-factorial computation. Endpoint paths are independently enumerated as literal directed words.

Normal and `-O` runs agree on:

- 65 rooted side contexts, including every labelled rooted side through order3 and deterministic sides of orders4 and5.
- 325 complete physical output polynomials, each with empty, singleton, transitive-two, transitive-three and cyclic-three inputs.
- 65 exact endpoint signs and the delta/Delta coefficient decomposition.
- 2,395 literal crossing-cycle words, with exact equality to the uniquely reconstructed side/middle-path triples.
- 145 excluded singleton/singleton/empty-middle length-two configurations.
- Four actual Q_k full profiles through order12, the closed B(2) expression, and the exact Q_3 ratio collision.

Independent checker SHA256: `d1fdf966fa8a13c04a413a9e133e39df5e8b3dbcf0755767a2be4078cbdcbb99`.
Independent normal/optimized receipt SHA256: `aad7ca42db47f20dfb366436f50dc5aa91b600ba09b04b2603a58993ffab906e`.
Source optimized/debug receipt SHA256: `0d6755fdff8e737a2cac778ae437704147cf4c793cb84208849bf8cf3f8f2bec`.

Run `python3 check.py` or `python3 -O check.py`. The source replay uses standard C++17 in `source.cpp`. These finite controls support the separately reviewed all-order proof; they are not used to extrapolate it.
