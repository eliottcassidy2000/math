# Globally least positive-unit ternary clone amplitude

## Theorem and scope

Among all rooted finite tournaments, the globally least A>0 for which transitive root cloning satisfies

    H_m=A*3^m+1

for every sufficiently large positive integer m is **72**. Every eventual identity of this form is exact for every positive m. No self-converse assumption is made.

An attaining strong source has rows

    [414,260,8,50,294,327,415,318,12], root4,

and complete positive root-run profile(217,216,72). It has counts217,649,1945,5833,17497 at multiplicities1 through5.

This strengthens the earlier marked self-converse minimum by removing the symmetry hypothesis. It does not classify every equality source, every affine graph operation, or the unrestricted Hamilton-path spectrum. The existing negative-unit minimum60 theorem remains unchanged.

## 1. Exact run reduction and integrality

Let a_j count words using every nonroot vertex once and j indistinguishable root symbols, with no adjacent root symbols and with every consecutive pair a tournament arc. Replacing the root by a transitive m-vertex module lifts a j-run word in j!S(m,j) ways: assign the m ordered clone vertices surjectively to the j ordered runs, preserving their increasing order within each run. Thus

    H_m=sum_j a_j j!S(m,j),  m>=1.

The largest nonzero index K supplies exponential base K with positive coefficient a_K. Hence an eventual positive-amplitude ternary identity forces K=3 and A=a_3 to be an integer. Independence of the exponential bases1,2,3 then gives

    a_1=3A+1,  a_2=3A,  a_3=A.

The identity is consequently exact for all positive m. Since H_1 is odd, A is even. A proposed A<72 is therefore at most70, and its source has H_1<=211.

The maximum run index is one more than the bipartite matching number from root out-neighbors to root in-neighbors. The direct proof is in the separately audited `../MINIMUM_TERNARY.md`: each internal fragment supplies a disjoint crossing edge, while a matching can be extended into the required fragment word by tournament path insertion and cutting. Thus only roots of matching rank two need be examined.

At rank two, the complete coefficient equations are equivalent to the two finite tests

    H_2=3H_1-2,  H_3=9H_1-8.

Indeed H_2=a_1+2a_2 and H_3=a_1+6a_2+6a_3. The tests determine a_2=H_1-1 and a_3=(H_1-1)/3, giving the stated identity. A larger matching rank cannot be hidden by cancellation of its strictly larger positive leading exponential.

## 2. All-order reduction to strong cores of order at most ten

Hamilton counts multiply over the ordered strong components. If the marked root is a singleton SCC, transitive cloning leaves the count constant, contradicting A>0. Otherwise its nontrivial SCC remains strong under cloning; all other components contribute a fixed odd product D.

Because the eventual identity is already exact at m=1,2, D divides both3A+1 and9A+1. Their gcd divides2, and both numbers are odd. Hence D=1. Every nontrivial strong tournament has at least three Hamilton paths, so all exterior SCCs are singletons. Removing them preserves the entire clone-count family.

The remaining strong root core has H<=211. Busch's published Theorem1 gives H^3>=5^(n-1) for nontrivial strong tournaments. The exact integer comparison is

    5^10=9,765,625 > 211^3=9,393,931.

Therefore the core has n<=10. This is an unrestricted all-order cutoff, not an inference from a finite search window. It applies equally to arbitrary singleton padding and to eventual identities.

The named source is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006)N3, Theorem1, DOI10.37236/1141, https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf. The same theorem was checked directly for the earlier negative-unit and symmetric packets. No assumption that strong order6 has minimum17 is used; its actual minimum is15.

## 3. Complete capped strong-core generation

Two elementary facts justify the finite bank.

**Induced-count monotonicity.** Given a Hamilton word of T-v, insert v at a deterministic valid position. If v beats the first vertex, prepend it; if the last vertex beats v, append it; otherwise a transition in their incidences supplies consecutive vertices x,y with x->v->y. Removing v recovers the original word, so this is an injection. Iteration gives H(S)<=H(T) for every nonempty induced subtournament S. Thus it is safe to discard an intermediate graph as soon as its count exceeds211.

**Strong deletion lemma.** Every strong tournament of order at least4 has a strong induced subtournament one vertex smaller. Choose a proper strong subset S of maximum cardinality; a cyclic triangle exists. If its complement has one vertex, there is nothing to prove. Otherwise no exterior vertex can be mixed on S, because adding such a vertex would give a larger proper strong subset. Partition the exterior into A with S->A and B with B->S. Strongness forces both to be nonempty and an edge a->b with a in A,b in B. The set S union{a,b} is strong. If there are at least three exterior vertices, it is a larger proper strong subset, a contradiction. If there are exactly two, delete any vertex of S instead. The remaining nonempty S-part, together with a,b, is strong through the cycle of blocks S-a-b-S, regardless of its internal connectivity. This gives the desired order n-1 strong graph.

Consequently every strong tournament can be built from C3 by repeatedly adding a vertex with at least one incoming and one outgoing edge. If its final H is at most211, every intermediate deletion graph also has H at most211. Enumerating every mixed orientation at each step, then applying an exact isomorphism key, is therefore complete under the cap.

The key first groups vertices by outdegree, then minimizes the relabelled adjacency tuple over every permutation within those groups. Since every isomorphism preserves outdegree, equality of these keys is equivalent to graph isomorphism. This is not a heuristic hash or a labelled-seed extension assumption.

The forward Hamilton DP uses exact nonnegative counts and rejects only when a subset count exceeds the cap. Such a subset is an induced tournament, so the preceding injection justifies the early rejection. Retained graphs are independently recounted without a cap in Python.

## 4. Complete finite certificate

The exact retained strong class counts at orders3 through10 are

    1,1,6,35,353,1074,644,174.

There are2,288 retained graphs. The extension census is:

- order4:6 candidates,6 accepted before isomorphism,1 class;
- order5:14 candidates,14 accepted,6 classes;
- order6:180 candidates,180 accepted,35 classes;
- order7:2170 candidates,2170 accepted,353 classes;
- order8:44478 candidates,6450 accepted,1074 classes;
- order9:272796 candidates,3612 accepted,644 classes;
- order10:328440 candidates,974 accepted,174 classes.

Among all retained graphs, exactly2,642 roots have H=1 modulo6 and crossing matching rank two. Every such root is tested at clone multiplicities2 and3. **None passes both positive-unit equations.** Therefore no A<72 is possible at any order.

The supplied verifier regenerates the bank from the triangle. Its Python stage independently checks every retained graph's exact canonical key, strongness and uncapped Hamilton count; every candidate root is independently matched by augmenting paths and checked using reverse-endpoint clone counts. It separately regenerates every labelled strong tournament through order6 and compares the resulting canonical bank. These are distinct from the C++ forward DP and subset-state matching used in generation.

## 5. Attainment and reproducibility

The displayed nine-vertex source is strong and has matching rank two. The verifier directly counts its actual clone graphs at m=1,...,5 and obtains

    217,649,1945,5833,17497.

The first three counts give(a_1,a_2,a_3)=(217,216,72), and the run-degree bound proves the identity72*3^m+1 for all positive m. This attains the global lower bound.

Run `python verify.py` and `python -O verify.py`; deterministic receipts must agree. The C++ source is `verify.cpp`. A separate unoptimized C++ replay with library checks produces the identical bank and search result. All files stay inside this session-specific scratch directory; no worktree or unrelated checkout is changed.

The source attainer was already independently reviewed in the peer's `session197-orbit-wildcard/positive-clock/PROOF.md`. The present new content is the unrestricted strong-core cutoff and complete global smaller-amplitude exclusion. It does not claim a new discovery of that source or a literature-priority result.
