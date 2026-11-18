# Exact mass transport on the actual merged metagraph

Status: all-order elementary consequences of inherited THM-796/THM-643; independent finite replay through order five. No new universal selector or complete adjacency recurrence is claimed.

## Carrier and inherited priority

Let X_n be tournaments on the fixed labels 0,...,n−1 containing the directed path 0→1→...→n−1. There are M_n=2^((n−1)(n−2)/2) such tilings. Let π_n map to tournament isomorphism classes modulo converse. The merged metagraph lines are pairs {t,κt}, where κ reverses every NONPATH arc. They are not vertex switching classes.

For a merged node u represented by T, put χ_u=1 if T is self-converse and 2 otherwise. The fibre weight is

    w_n(u)=χ_u H(T)/|Aut(T)|.

Automorphisms act freely on ordered Hamiltonian paths; relabel a path into the displayed order. This proves the ordinary H/Aut count and then the converse factor. These are inherited THM-643/796 laws. The three-sorted recursion and its failure to lump by merged nodes were already proved/tested in THM-796, especially `three_sorted_metagraph_recursion_codex_S9.out`. The present note exposes its probabilistic and endpoint coordinates rather than claiming another recursion discovery.

## Cross-order transport and exact endpoint reading

Append vertex n with the forced arc n−1→n. The remaining n−1 new arcs are free. Delete the last label to recover the unique parent. Define E_n(u,v) as the number of (n+1)-tilings whose parent is in u and whose merged node is v. Then

    Σ_v E_n(u,v)=2^(n−1)w_n(u),
    Σ_u E_n(u,v)=w_(n+1)(v).

These are integer equalities, not equality of unweighted numbers of classes.

There is an intrinsic endpoint formula. For a representative D of v, e_D(x) is its Hamiltonian-path count ending at x. Let R(v) contain D alone when self-converse, otherwise D and its converse. Then

    |Aut(D)| E_n(u,v)
       = Σ_(B in R(v)) Σ_(x : merged(B−x)=u) e_B(x).

Proof: classify Hamilton paths of each B by endpoint x; delete that endpoint. Divide by the free automorphism action on whole paths. No freeness on vertices is assumed. The total numerator is divisible by |Aut D| even if its individual endpoint terms are not.

Consequently the exact forward and backward kernels are

    K_n(u,v)=E_n(u,v)/(2^(n−1)w_n(u)),
    L_n(v,u)=E_n(u,v)/w_(n+1)(v).

Under μ_n(u)=w_n(u)/M_n, they obey Bayes transport

    μ_n(u)K_n(u,v)=μ_(n+1)(v)L_n(v,u).

This yields correct marginals and a Markov model obtained by refreshing uniformly within each fibre. It does NOT imply the literal projected growth process is Markov: past observations can select different presentations inside a fibre.

## A concrete continuation obstruction and what to retain

All of the following source tilings are the unique strong four-vertex class H=5. Encode forward arcs in order (01,02,03,12,13,23). Exhaustively append the eight possible new vertices with 3→4 forced:

* source 100101 gives child-count histogram {5:1,9:4,11:2,13:1};
* source 110101 gives {5:1,9:4,11:1,13:1,15:1}.

Thus even existence of a next H=15 differs between presentations of the SAME tournament class. A metagraph node, its H, its automorphism group, and its converse type do not tell whether this particular continuation admits 15. The retained last vertex matters. To exploit the merged transport for target selection, preserve an endpoint-marked class or finer continuation signature, rather than replacing it by its averaged row. Replay finds two nonlumpable parent merged fibres at the 4→5 transition.

## Information conservation and target mass

For a uniform tiling, Shannon entropy in natural logarithms decomposes exactly:

    log M_n = Entropy(μ_n)+E_(μ_n) log w_n.

Therefore

    Entropy(μ_(n+1))−Entropy(μ_n)
      =(n−1)log2 + E_(μ_n)log w_n − E_(μ_(n+1))log w_(n+1).

The new n−1 bits divide between identifying a merged node and residual Hamilton-path presentation freedom. This is elementary conditional entropy, not evidence of a particular modular form or physical energy law.

For a positive count N,

    μ_n{H=N}= N/M_n · Σ_(ordinary classes C:H(C)=N) 1/|Aut C|.

Positivity is equivalent to existence, but conserved total mass does not force this individual coefficient positive. This is the precise gap between an enumeration/distribution theorem and Chan–Pak. A large mass may concentrate at already attained counts.

## Why the n! switching identity is different

For a fixed labeled tournament T, vertex switching reverses all arcs across a cut. There are 2^(n−1) distinct switches. Every vertex ordering is Hamiltonian in exactly one switch: prescribe its first switch bit, and each consecutive arc uniquely prescribes the next bit. Therefore

    Σ_switches H(T^s)=n!.

This counts Hamiltonian orderings in ONE switching class. The tiling/metagraph identity instead sums H/Aut over ALL ordinary classes and equals M_n. Neither an orbit of cuts nor a path presentation can be substituted for the other. For example at order five the numbers are 120 (switching total), 16 (switches), 64 (tilings), 12 (ordinary classes), and 10 (merged nodes).

## 3,5,13,24: genuine connections and boundaries

The recovered metagraph source `metagraph_connection_corollaries_20261021/session24-connections/CONNECTIONS.md` gives a rigorous triangle-substitution injection of Hamilton-path orbits n→3n, preserving converse fixed/pair types. It also gives the decisive obstruction to transporting WIRING: κ(C3[C3]) has H=2983, whereas κ(C3)[C3] has H=27. Entry/exit contacts cease to be uniform; retain ports. The pure-blue strong atoms {1,C3,R5} belong to this exact carrier, not to the prime list13 or a dimension24 by numerology.

The independent repo Paley→Golay→Leech construction (`paley_golay_leech_moonshine_20261018.out/json`) yields vector families 97152,98304,1104 totaling196560, and the dimension identity196560+24+300=196884. This is an actual chain of different mathematical objects. It does not identify any24 blue self-lines,13-path source, or merged fibre with those vector families. Importing a code would require an explicit map preserving the marked extension signatures above. Counting coincidences alone cannot provide it.

## Reproduction

`python metagraph_mass_transport_replay.py`; also rerun under `python -O`. All assertions hold and outputs match. Exhaustive fixed-path tilings through order five, independent permutation canonicalization, subset-DP endpoint counts, ordinary/converse automorphism factors, all transport entries, and entropy identities. Results are deliberately small controls of the all-order proofs, not new coverage evidence.
