# A balanced order14 atom

Status: VERIFIED explicit strong tournament and PROVED construction upper bound. No exhaustive order14 census, exact minimum order, or global optimality claim is made.

There exists a strong tournament T on14 vertices with

 H(T)=11390625=15^6.

An explicit witness is the adjacency code292187853813044375413825338. Bit k, starting at zero, belongs to the kth lexicographic pair i<j and equals one exactly when i->j. The full decimal integer is stored without floating-point conversion. `verify.py` reconstructs the orientation, checks strong connectivity, computes both start- and end-point subset tables, verifies their sums and checks the arc-incidence normalization. The independent suffix-recursion checker uses separately implemented adjacency and reachability routines.

Let m(h) denote the least order of any tournament with exactly h Hamiltonian vertex-word paths. Dominance joins multiply path counts and add orders. Therefore q copies of T yield count15^(6q) on14q vertices. For k=6q+s,0<=s<6, append s copies of an order-five H15 atom. This proves

 m(15^k)<=14q+5s,     lim_k m(15^k)/k<=7/3.

Existence of the limit follows from the inherited subadditivity argument. The new bound improves the prior13/5 upper bound, while leaving that earlier witness theorem valid. This does not prove m(15^6)=14. Larger atoms may improve the asymptotic rate again, and no positive global lower bound on the linear rate is claimed.

## Discovery and method boundary

The preceding balanced-atom-oracle packet proves the exact single-arc reversal formula H(T')−H(T)=D_vu−D_uv for u->v, where D_uv is the complementary-subset endpoint convolution. The same audited source was used here, with mt19937_64 seed20813, target11390625 and a cap of2000 restarts. It reached the exact target at zero-based restart1525, step42. An initial cap300 run had reached only a near miss; no exclusion claim followed from that finite failure.

Each selected arc difference was checked against a fresh DP recount. The discovery graph reverses single arcs and is not the repository's native path-complement metagraph. No native-metagraph degree or fibre consequence is claimed. The proof rests on the explicit independently verified witness, not a convergence guarantee for this search or a spanning-tree approximation. The attached preprint's main logarithmic theorem is not used.

## Reproduction

Run `python3 verify.py` and `python3 independent_check.py` with assertions enabled. Their outputs match the archived outputs; the main verifier deterministically regenerates RESULT.json.

Optional discovery replay:

    g++ -O3 -std=c++17 target_search.cpp -o /tmp/s208-target14
    /tmp/s208-target14 14 2000 11390625 > /tmp/target14.out 2> /tmp/target14.log

This replay is much slower than direct witness verification. The source is the unchanged audited oracle from commit386fa872bd9b0baab3c3c55d86c818d7acaf8342, `session208/balanced-atom-oracle/target_search.cpp`. All verification inputs are bundled; no network access is needed to verify the witness.
