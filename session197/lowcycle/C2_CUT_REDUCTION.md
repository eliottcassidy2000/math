# The remaining bridge-free C=2 problem is a bounded 2-strong problem

Session 197. **Partial all-order reduction, not a complete C=2 classification.** This uses the separating-vertex product proved in `C1_COMPLETE.md` and a published Hamiltonian-cycle lower bound.

## 1. A cut vertex costs at least four cycles when there is no strong bridge

Let T be strong and let r be a separating vertex. Write D_1,...,D_k for the ordered strong components of T-r, with k>=2. The exact product from `C1_COMPLETE.md` is

    C(T)=A * product_(i=2)^(k-1) H(D_i) * B,

    A=sum_(a in D_1, r->a) s_a(D_1),
    B=sum_(b in D_k, b->r) e_b(D_k).

Every summand and factor is a positive integer. If A=1, r has exactly one out-neighbor a in D_1, and r->a is the sole arc from the complement into D_1. It is a strong bridge. Dually, B=1 exposes a strong bridge b->r from D_k into its complement.

Therefore, if T has no strong bridge, both A and B are at least two, giving C(T)>=4. In particular:

**Proposition.** Every strong bridge-free tournament with C=2 or C=3 is 2-strong. Equivalently, every non-2-strong tournament with C=2 or C=3 has a strong bridge.

This is only a lower bound of four in the cut-vertex case; no sharpness claim is made. A bridge-free non-2-strong example with C=9 has rows

    [10,12,9,112,39,71,23].

The verifier checks all its arc deletions and its separating vertices directly.

## 2. A rigorous order bound for the residual atoms

The permitted primary source is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics 13 (2006), N3, [DOI](https://doi.org/10.37236/1141), [official paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf). Its **Corollary 2**, checked in the primary text, states that a 2-connected tournament of order n has at least

    beta^(n/32-1) Hamiltonian cycles,   beta=5^(1/3).

Here 2-connected is vertex-strong connectivity, the 2-strong condition used above. This is a published bound, not a new estimate. For C(T)=c it implies

    5^(n-32) <= c^96.

For n>=32, the exact integer comparisons

    5^41 <= 2^96 < 5^42,
    5^65 <= 3^96 < 5^66

give the following bounds. Orders below 32 are already within them.

- Every bridge-free C=2 tournament has order at most 73.
- Every bridge-free C=3 tournament has order at most 97.

The first result reduces the conjecture that the regular five-tournament is the sole bridge-free C=2 atom to a finite statement through order 73. This is not a feasible exhaustive-search claim, and the finite statement is **not proved here**.

## 3. Finite controls and what remains open

Complete labelled enumeration through order six finds exactly 24 2-strong five-tournaments, all isomorphic to the regular five-tournament and all with C=2. At order six the 2-strong minimum is four, so no new C=2 or C=3 atom occurs there. These small checks do not establish the claim through order 73.

The unique-cycle theorem in `C1_COMPLETE.md` is complete at all orders. The present C=2 reduction deliberately stops short of a parallel classification. In particular it does not assume that all arcs of a 2-strong tournament lie on Hamiltonian cycles: that assertion is false, as the six-vertex rows `[14,20,56,50,33,3]` show. They are 2-strong, have five Hamiltonian cycles, and their arc 2->5 belongs to none.

`verify_c2_cut_reduction.py` uses the previously frozen literal enumeration helpers from `verify_c1_complete.py`, checks every labelled tournament through order six, the two explicit scope-control graphs, and the exact integer bounds. Normal and optimized runs must agree. The product argument and published bound prove the all-order reductions; the enumeration is finite evidence only.
