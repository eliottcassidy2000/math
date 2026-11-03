# How far the shared-center fan extends

**Scope.** This theorem is for one explicit family and one prescribed order of edge flips. It does not exclude longer equal-jump chains from other sources or using different edges.

Fix n≥4. On the labelled transitive base 0<1<⋯<n−1, let Z_t reverse exactly the edges 0–i for 2≤i≤t+1, where 0≤t≤n−2. The step Z_t→Z_{t+1} reverses exactly one edge, 0–(t+2). Put G(X) for triangle substitution X→a→b→X. For 0≤t≤n−3, these give actual edge reversals between strong tournaments of order n+2.

Every directed odd cycle in Z_t consists of 0,1 and a nonempty odd-sized subset of {2,…,t+1}, in increasing order. Indeed all cycles use 0, its only possible first vertex on a cycle is 1, and only a reversed endpoint can return to 0. All cycles share 0, so at most one can occur in a cycle packing. The odd-cycle packing profile formula therefore gives

    F_Zt(r) = r^n + 2 ∑_{q odd} binom(t,q) r^{n−q−1}.

In particular the difference when the next edge is reversed is

    Δ_t(r)=2∑_{s even} binom(t,s) r^{n−s−2}.                  (1)

Since H(G(X))=F_X(1)+F_X(2), the jump J_t in the wrapper is

    J_0 = 2 + 2^(n−1),
    J_t = 2^t + 2^(n−2−t)(3^t+1),    1≤t≤n−3.              (2)

The two initial jumps are equal. Whenever t≥1 and both J_t and J_{t+1} are defined,

    J_{t+1}−J_t = 2^t + 2^(n−3−t)(3^t−1) > 0.            (3)

Thus in this specific ordered fan the equal initial pair never extends to a third equal jump. Examples: n=5 gives +18,+18,+24; n=8 gives +130,+130,+164. The latter formulas come from actual labelled graphs, not formal state vectors.

There is also an immediate global constraint on the *unwrapped equal-profile mechanism* when starting from any transitive tournament: three successive edge reversals each with F-difference 2r^d would produce path counts 1,3,5,7, contradicting the established exclusion of 7 in [S157](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session157_gaps_20261103/gaps/gaps-note.md). This does not prevent longer wrapped equal-jump chains, because equal wrapped jumps need not have the same full response or even the same unwrapped H difference.

## Finite labelled seven-vertex test

Start from the specific G(T₅), with original vertices 0,…,4 transitive and 0,…,4→5→6→0,…,4. Exhaust all 21 possible edge reversals at each state, keeping precisely those raising H by 18. Subset-and-endpoint dynamic programming computes H; identical labelled endpoint graphs are deduplicated. The numbers of distinct graphs at successive counts are

    H=33: 1;   H=51: 3;   H=69: 14;   H=87: 0.

No strongness test is imposed, so the no-third-step statement holds even if intermediate or final graphs are not strong. This is a reproducible finite fact for this exact starting graph and seven fixed vertices, not a claim about all starting graphs or other orders.

Run `python3 fan_check.py`, `python3 -O fan_check.py`, `python3 chain_bfs.py`, and `python3 -O chain_bfs.py`. The first checker compares slot-cover profiles via independent subset/path convolution and directly checks wrapped graph counts and jump formulas for small n; the second exhausts the entire specified labelled edge graph up to depth three. Exceptions remain active under optimized Python.
