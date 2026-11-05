# Independent audit of the equal-response five-times window

Session 196, 2026-11-05. **PASS. No counterexample or unresolved proof gap found.** This audit supports the architecture-specific theorem, not unrestricted Hamiltonian-path spectrum coverage.

## Audited conclusion

From any physical rooted tournament with state `(h,d,h,d)`, allowing arbitrary finite rooted one-seam root-switch gadgets and interleaved stationary K,D moves, the exact scalar support through `5h` is

- `{h,3h,5h}` when `h=d`;
- `{h,2h-d,3h,4h-3d,5h}` when `h>d`.

In the second case physical structure forces `h>=3d`. Consequently, for `h>=3` and delta in `{2,4}`, a literal child `3h+delta` is reachable exactly when `h-3d=delta`; no fixed equal-response input reaches both children. These statements concern the specified marked input and operation architecture.

The final primary proof inspected was `../global-induction/PROOF.md`, SHA256 `ed89a4821d2e577cf3a92a673f4b61d160c5703e80ae582743f8e87bead9c5ac`. The final primary verifier inspected was `../global-induction/verify.py`, SHA256 `777d0599a81e886b93a819d10b25b90d0a3d5682e9d58093568b64f56d7986e4`.

## Analytical review

1. **Physical ratio.** Zero merge slack forces root outdegree at most one. The nontrivial case has exactly one out-neighbor a, and `e=d` implies no Hamiltonian path of the deleted graph ends at a. Thus a is in a nonterminal strong component. A rooted cover contains all earlier components in its first path and independently splits each later component between the two paths. The factorization in the proof is exact, and at least one later component supplies a factor of at least two. The singleton-birth formula then gives `h>=3d`. This handles weak as well as strong inputs. The endpoint fact used here can also be proved by taking a longest path ending at a prescribed vertex in a strong tournament: if it omits vertices, strongness provides an omitted vertex beating a path vertex, which can be inserted without changing the prescribed endpoint.
2. **Arbitrary gadgets.** Nontransitivity implies at least three Hamiltonian paths by insertion from a directed triangle. Before the first non-D operation, the gadget bound is at least `6h-3d>=5h`. Later, `R>=2h-d` gives at least `3R>=5h`. Hence nontransitive gadgets cannot contribute below the boundary, regardless of gadget order.
3. **Large transitive gadgets.** A source-marked transitive gadget of order at least two has the same scalar as A but has `e'=0` and `R'>=4h`. Any next non-D operation costs at least `H'+R'>=6h-d>5h`, including after intervening D moves. It is therefore terminal in the strict window and can be replaced by A for scalar reachability. Nonsource-marked transitive gadgets contribute only `3h` below the boundary, with no later non-D continuation remaining below it.
4. **Depth and boundary.** The four two-letter A/K states and all next-step scalar formulas check. Initial, internal and terminal D placements are covered. Reaching `5h` or more cannot later produce a smaller scalar. The boundary `5h` is supplied explicitly by DDA, so the proof does not need an unsupported claim of strict growth after hitting it.
5. **Degenerate cases.** The singleton signature is exactly A. The diagonal `h=d`, including `h=1`, is handled by integer-linear scaling of the physical singleton grammar and oddness of tournament Hamilton counts. Neutral K loops on diagonal states add no scalar. Terminal D moves preserve H; the closure algorithm records outputs with `H<=cap` even if `R>cap`, and only prunes their future exploration. Zero-endpoint states are handled without dividing by e or assuming it positive.

## Independent finite controls

The supplied primary verifier was rerun normally and under `python3 -O`. The two outputs agree byte-for-byte and equal its checked-in `receipt.json`.

The audit's separate C++ implementation uses permutations and legal cut positions, not the primary verifier's subset recurrence:

- Enumerated all **33,867 labelled tournaments through order six**, covering **202,013 roots**.
- Checked every physical inequality and `zero slack iff outdegree<=1` at every root. There are 38,615 zero-slack rooted inputs.
- Checked the ratio `h>=3d` for every nontrivial equal-response root.
- Found 12 distinct nontrivial `(h,d)` pairs: `(3,1),(5,1),(9,1),(9,3),(15,3),(15,5),(17,1),(17,5),(19,5),(23,5),(27,3),(37,5)`.
- Regenerated all **345 distinct eight-coordinate gadget signatures through order six**, checked the gadget inequalities and all transitive signature formulas, and verified that all 67 frozen signatures occur in this independent bank.

The separate Python saturation implementation adds every transitive rooted gadget through order twenty, making **534 distinct signatures**, and records terminal outputs before applying the response-coordinate exploration bound. It confirms the five-times support for:

- all 12 physical equal-response pairs found by the census, plus `(49,13)`;
- **5,151** integer parameter pairs with odd `d=1,3,...,101` and `h=3d+gap`, for even `gap=0,2,...,200`;
- nine diagonal parameters, including `h=1` and `h=10^50+1`.

The two independent closure runs, normal and optimized, have identical receipts. Parameter pairs outside the enumerated physical list are algebraic controls, not claims of physical realizability. These finite controls support the audit; the arbitrary-gadget and arbitrary-depth conclusion rests on the analytical proof.

## Reproduction

From this directory:

    g++ -O3 -std=c++17 census.cpp -o /tmp/s196_audit_census
    /tmp/s196_audit_census bank6.txt > census-receipt.txt
    python3 verify_closure.py > closure-receipt.json
    python3 -O verify_closure.py > closure-receipt.optimized.json
    cmp closure-receipt.json closure-receipt.optimized.json

All checks use explicit exceptions or nonzero return codes, not removable assertions. The only inherited files read independently were the named combined-growth and zero-slack proofs at commit `ee07a5aac9cdfd8d3834fa834df4fb357fcaf769`, through the GitHub connector. No dirty repository checkout was read or modified.
