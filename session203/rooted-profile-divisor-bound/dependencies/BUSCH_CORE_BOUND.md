# Unrestricted logarithmic cores and an exact numeric-target algorithm

Session197. **An elementary corollary of Busch's published strong-order bound, with exact finite controls. No novelty or literature-priority claim.** This is an unrestricted decision-and-witness algorithm, not a constructive proof that all odd targets other than7 and21 occur. Its running time is quasi-polynomial in the numerical target, not polynomial in its binary length.

## 1. Definitions and the published input

For a nonempty finite tournament T, H(T) counts directed vertex words using every vertex exactly once. A singleton has count1. The code counts an empty graph as0; an empty component product of1 below is only bookkeeping. Bit j of row i means i→j.

Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics13 (2006), N3, Theorem1 and Corollary1, prove the following sharp minimum f(n) over nontrivial strong tournaments:

- n≡0 modulo3: f(n)=3·5^((n−3)/3).
- n≡1 modulo3: f(n)=5^((n−1)/3).
- n≡2 modulo3: f(n)=9·5^((n−5)/3).

Here n≥3, since no tournament of order2 is strong. In particular H(T)^3≥5^(n−1).

Verified primary source: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf

DOI: https://doi.org/10.37236/1141

The paper also records Rédei's oddness of H and ordered strong-component factorization. We give the factorization argument below. The new statements here are consequences of that published input, not a new proof of Busch's theorem.

## 2. A uniform component inequality, including small orders

**Lemma.** Every nontrivial strong tournament S of order s satisfies

    3^s <= H(S)^3.                                      (1)

For s=3, S is a directed triangle and H=3, so equality holds. For s=4, a strong tournament contains a spanning directed cycle: choose a cyclic triangle, whose remaining vertex must be mixed on it, and insert that vertex at a suitable cycle edge. Cutting this4-cycle gives four distinct Hamilton paths. H is odd, so H≥5; hence H^3≥125>81=3^4.

For s≥5, Busch gives H^3≥5^(s−1). The inequality5^(s−1)≥3^s holds at s=5 because625≥243, and is preserved as s increases since the left side multiplies by5 and the right by3. This proves (1) at every permitted strong order, with no rounding assumption.

## 3. Removing irrelevant singleton components

The strong components of a tournament are totally ordered, with every edge from an earlier component to a later one. A directed path can never return to an earlier component. Every Hamilton word is therefore exactly the concatenation, in this component order, of one Hamilton word from each component. Consequently

    H(T) = product_i H(S_i).                            (2)

Singleton components contribute1. If H(T)=h>1, at least one component is nontrivial. Delete every singleton component and retain all vertices of the other components. The resulting nonempty induced tournament T0 has the same H=h. Write N=|T0|. Applying (1) to each retained component gives

    3^N <= h^3,
    N <= floor(3 log_3 h).                              (3)

The floor is outside the multiplication by3. Replacing this by3 floor(log_3 h) is false: a strong four-tournament has H=5, and no tournament of order at most3 has H=5.

The bound is computed exactly as

    B(h)=max {N>=0 : 3^N<=h^3}.

The implementation uses successive integer powers, never floating-point logarithms. The target h=1 is handled separately by a singleton witness. In that case B(1)=0, and deleting every singleton component would leave the empty graph, so no nonempty preservation statement is made for that reduction.

Ordered joins of k directed triangles have N=3k and h=3^k. Thus (3) is sharp as a bound on the size of a reduced tournament. This does NOT assert that those h require that many vertices in a smallest alternative witness.

## 4. Optional sharp minimum at each reduced order

Let g(N) be the minimum H among tournaments of order N having no singleton strong components. The permitted positive orders are N≥3. Then

    g(3q)=3^q,
    g(3q+1)=5·3^(q−1)       for q>=1,
    g(3q+2)=25·3^(q−2)      for q>=2,
    g(5)=9.                                             (4)

This gives a slightly stronger bound max{N>=3:g(N)<=h} for nonunit targets. Setting g(0)=1 is convenient for component-partition bookkeeping; no positive reduced tournament of order1 or2 exists.

**Proof.** First minimize the product of costs3,5,9 over partitions into sizes3,4,5. A size5 part together with a size3 part can be replaced by two size4 parts, lowering the cost27 to25. A size5 with a size4 can be replaced by three size3 parts, lowering45 to27. Two size5 parts can be replaced by sizes3,3,4, lowering81 to45. Therefore a size5 part occurs in an optimal partition only when the total size is5. Three size4 parts can be replaced by four size3 parts, lowering125 to81. For every other permitted size, the number of size4 parts is thus0,1,or2, uniquely fixed by the residue modulo3. This proves (4) for the partition problem.

Each published strong minimum f(s) is at least this partition value g(s). Equality holds at3,4,5. At6,7,8 the comparisons are15≥9,25≥15,45≥25. For s≥9, f(s)=5f(s−3), while g(s)=3g(s−3), so induction finishes the comparison. Replace each SCC of an arbitrary reduced tournament by a size3/4/5 partition of cost g(s)≤f(s)≤H(S). The union of those partitions has cost at least the minimum g(N), proving H(T)≥g(N).

Conversely, the partition is physically realized by ordered joins of directed triangles, strong four-tournaments with H=5, and the exceptional strong five-tournament with H=9. Explicit rows for these three seeds are

    [2,4,1], [6,12,8,1], [14,28,24,16,1].

The verifier counts them directly. This proves equality at every reduced order. It is an exact extremal corollary of the published strong minima, not a spectrum classification.

## 5. Direct unrestricted decision algorithm

For a positive integer target h:

1. If h=1, return the singleton. Even h are impossible by Rédei's theorem.
2. For each order1<=n<=B(h), enumerate all2^binom(n,2) labelled tournaments.
3. Count Hamilton paths by the subset endpoint recurrence. Return the first graph with H=h; if no graph matches, report that h is unattainable.

Completeness follows from the nonempty reduced witness in (3). Labelling causes duplicate work but cannot omit a witness. The endpoint recurrence is

    E({v},v)=1,
    E(A,v)=sum_(u in A−{v}, u→v) E(A−{v},u),
    H(T)=sum_v E(V(T),v).

All entries may be capped at h+1: addition is nonnegative, and min(h+1,x+y) is unchanged if each summand is first capped. The final result therefore equals H when H<=h and otherwise certifies overflow. This involves no subtraction or signed readout.

Write b=ceil(log_2(h+1)). Since B(h)=O(b), the graph enumeration has2^O(b^2) iterations. The DP requires O(n^2 2^n) additions of O(b)-bit capped integers per graph. The complete bit running time is

    2^O(b^2) = 2^O((log h)^2) = h^O(log h).              (5)

Processing one graph at a time uses2^O(b) bits of workspace, including the DP and graph. Integer powers and enumeration bookkeeping add only polynomial factors. Neither the logarithmic vertex count nor this exponential-in-b verification establishes unrestricted NP membership for binary targets.

## 6. More economical implementation: enumerate only strong factors

The supplied algorithm.py implements an equivalent method with a smaller enumeration order. Compute

    S(h)=max {n>=3 : f(n)<=h},

using exact integers, with S(h)=0 if the set is empty. Since f(n)^3≥5^(n−1), every such order satisfies n<=1+3 log_5 h.

Enumerate all labelled tournaments through S(h), retain the strong ones, and record one representative for every count c<=h. Test whether h is a product of recorded counts c>1. One can restrict every candidate factor and every intermediate product to divisors of h. A successful product gives an actual witness by ordered joining the stored component representatives; an unsuccessful finite search proves unattainability.

**Completeness.** Every nontrivial SCC of any h-witness has count c dividing h, hence c<=h and order at most S(h), so its representative is in the bank. Conversely any product of bank counts is realized by the indicated ordered join. Singleton components are irrelevant, except that h=1 uses a singleton witness.

The strong-factor algorithm still has bit time2^O(b^2) and space2^O(b). At most h distinct counts and at most h divisor states are retained, each of polynomial bit size; this is within the stated exponential space. It does not assume that any prime is realizable, that indecomposable counts are prime, or that a gcd-one interface gives coverage. Code returns either literal adjacency rows or None.

## 7. Prior sources and reproducible checks

The authorized S196 guide and its explicitly linked symmetry, spectral-moves and orbit-spectrum packets were inspected at commit cf3b44b9a442c3cff0f28cab924e73c4c00bb00c. Its linked S195 guide and logarithmic-growth/strong-structure packets were inspected at ee07a5aac9cdfd8d3834fa834df4fb357fcaf769. The logarithmic-growth result there concerns a specified rooted grammar; it is not the unrestricted bound above. The session197 ribbon-compression result gives a fixed-width logarithmic witness and a fixed-width NP conclusion by an additional constant-size transfer. This note does not export that NP conclusion to unbounded width. No claim that the present elementary corollary is absent from all other repository or literature sources is made.

Starting guide: https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/00-navigation/S196.md

Restricted prior bound: https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/logarithmic-growth/PROOF.md

Fixed-width contrast: https://github.com/eliottcassidy2000/math/blob/c5bf78bb711ae666ab398ca40a45cea3077d09a1/session197/ribbon-compression/PROOF.md

Run python verify.py and python -O verify.py; the receipts must agree byte for byte. The standard-library verifier uses explicit exception gates, not assert. It exhausts all labelled tournaments through order6, checks literal Hamilton words through order5, independently obtains SCCs by transitive closure, tests count-preserving singleton removal and all small minima, checks exact multiplicative partition minima through order200, and tests integer rounding for targets1..10000. A complete strong bank through6 decides all targets through23; the finite conclusion is justified at arbitrary graph order by f(7)=25. The known omissions7 and21 are consistency checks, not claimed as new omissions. Negative/zero targets, even targets, h=1 and the h=5 floor-placement counterexample are explicit controls.

Run python algorithm.py TARGET for the actual complete decision procedure. Its worst-case cost is large; no practical efficiency at large targets is promised. The full proposed {7,21} spectrum remains open here.
