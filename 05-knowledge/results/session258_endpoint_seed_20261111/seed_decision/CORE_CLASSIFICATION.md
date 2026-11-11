# Every positive unit-offset tail seed has exactly one marked strong core

**PROVED, with independent literal finite controls.** The structural proof does not use the strong-order floor. The finite classifier and logarithmic witness-order bound use the explicitly cited Busch floor already recorded in `PROOF.md`; that inherited floor is not re-audited here.

## Definitions and provenance

Let H(B) count directed Hamilton paths. Let K_ab(B) count ordered pairs (P,Q) of nonempty disjoint directed paths covering B with P ending at a and Q starting at b. For selected O⊆V(B), define

    W_B(O)=Σ_{a,b∈O}K_ab(B),
    U_B(O)=Σ_{b∈O}(s_b(B)+Σ_a K_ab(B)),

where s_b is the number of Hamilton paths starting at b. Thus U includes an empty first path. A **positive unit-offset tail seed** satisfies H(B)−W_B(O)=1 and U_B(O)>0. This is specifically the selected-set/transitive-tail construction of THM-5126 and HYP-9658, not an arbitrary doubling-vertex construction.

Inherited: order-join factorization of H (THM-1862-order-join-reduction-principle.md); the tail formula H=1+2^t U; source padding preserving H,W,U (`05-knowledge/results/metagraph_aspect_sequences_20261111/PROOFS.md`, Section 3); and HYP-9658's transitive-base restriction |O|≤1, slope 1. The present proof supplies the converse structural classification, the exact sink-padding law, the cross-component barrier, and the resulting finite strong-core search. A narrow check of these closest sources found no statement of this full converse; no exhaustive priority claim is made.

## 1. Cross-component selected pairs cost more than the whole count

Write the canonical SCC decomposition as an ordinal word

    B=S_1→S_2→...→S_m.

Take a∈A=S_i and b∈C=S_j with i<j. For a tournament D define Q_2(D) as its number of ordered covering pairs of directed paths allowing either path to be empty. Define

    U_a(A)=number of such covers of A with second path nonempty, starting at a;
    F_b(C)=number of such covers of C with first path nonempty, ending at b.

Let E be the product of H(S_l) over l<i or l>j. Exact path restrictions give

    K_ba(B)= E · U_a(A) · F_b(C) · ∏_{i<l<j} Q_2(S_l).       (1)

All vertices before A must lie in the first path, and all vertices after C must lie in the second. In A the second path starts at a; in C the first path ends at b. Every intermediate block may be split arbitrarily between the two paths, including an empty local restriction. Conversely all these choices concatenate uniquely because every inter-block arc points forward. This proves (1), including all empty-local-path cases.

For every Hamilton path of A, split immediately before a; this injects H(A) paths into U_a(A). Likewise split a Hamilton path of C immediately after b to obtain F_b(C)≥H(C). Finally Q_2(D)≥2H(D), from the two choices of placing a whole Hamilton path in one path and leaving the other empty.

If j>i+1, equation (1) therefore gives K_ba(B)≥2H(B). If j=i+1, it gives K_ba(B)≥H(B), while the opposite order has

    K_ab(B)=E · e_a(A) · s_b(C)>0.                         (2)

Endpoint positivity in a strong SCC follows by breaking a directed Hamilton cycle at the desired endpoint (Camion's standard theorem); singleton SCCs have endpoint count one. Hence in every case

    K_ab(B)+K_ba(B)>H(B).                                 (3)

If an intermediate SCC exists, K_ab=0, because no path can cover that intermediate component while ending in A or starting in C. This is consistent with the separate strict bound above.

## 2. Canonical strong-core form

**Theorem.** Every positive unit-offset tail seed is exactly of the form

    B=TT_a→S→TT_b,     a,b≥0,
    O⊆V(S),

where S is one SCC (possibly a singleton), all exterior vertices are unselected, and

    H(S)−W_S(O)=1,
    U_B(O)=2^b U_S(O)>0.                                  (4)

Conversely every such decorated core and every a,b≥0 give a positive unit-offset tail seed.

**Proof.** Since U>0, O is nonempty. If O meets two different SCCs, W contains the two nonnegative terms in (3), so W>H, contradicting W=H−1. Thus O lies in a single SCC S.

For every selected pair of endpoints inside S, earlier SCCs are forced wholly into the first path and later SCCs wholly into the second. Therefore

    W_B(O)=H(exterior) W_S(O),
    H(B)=H(exterior) H(S).

Their difference is one. Both factors in

    H(exterior) [H(S)−W_S(O)]=1

are integers, with H(exterior)>0, so H(exterior)=1 and H(S)−W_S(O)=1. Every exterior SCC has count one, hence is a singleton: a nontrivial strong tournament contains a triangle and has H≥3 by the odd-cycle formula.

For U, the second path starts in S, so all leading vertices are forced into the first path. Each of the b trailing transitive vertices may be assigned independently to either path, producing exactly 2^b assignments. Their induced order uniquely determines the local path, and concatenation preserves the second path's start. This proves U_B=2^bU_S. It also proves the converse and that leading padding is neutral. ∎

The singleton-core case is precisely a transitive B with one selected vertex. Its U is 2^b. It is not necessary for B itself to be strong; the theorem identifies the unique selected SCC and the exact harmless exterior padding.

## 3. Selected-start injection and the half-multiplier bound

**Lemma.** For every tournament B and every O,

    U_B(O)≥|O| H(B).                                      (5)

For each pair (Hamilton path P, selected vertex a), split P immediately before a. The resulting ordered cover has second start a. Concatenation recovers P, and its second start recovers a, so this is an injection. No symmetry or endpoint-positivity assumption is needed.

A nontrivial strong unit core has H(S)>1, and therefore |O|≥2, because W vanishes for |O|≤1. Consequently

    H(S)≤U_S(O)/2≤U_B(O)/2.                               (6)

More precisely H(S)≤U_S(O)/|O|, which can be used when enumerating a specified selected-set size.

## 4. Complete fixed-q classifier without automaton search

Let R(h) be the inherited upper bound on the order of a strong tournament with H≤h (the code's `bound(h)`), derived from the supplied Busch strong floor. For fixed q≥1, a tail seed with U=q exists iff either:

1. q is a power of two, supplied by a singleton core followed by v_2(q) unselected sink vertices; or
2. there is a marked nontrivial strong S with

       |S|≤R(floor(q/2)),
       H(S)−W_S(O)=1,
       q/U_S(O) is a nonnegative integral power of two.

The criterion is necessary by the structural theorem and (6), and sufficient by adding exactly b=log_2(q/U_S(O)) trailing singleton vertices. Leading padding can be omitted. This is an explicit finite classifier over bounded strong cores. The nine-state automaton remains correct but is no longer needed for decidability.

If a nontrivial fixed-q seed exists, there is one on at most

    R(floor(q/2)) + v_2(q)

base vertices. In the trivial case there are at most 1+v_2(q). Since the inherited floor grows exponentially with order, this is O(log q), substantially sharper than the earlier finite-state bound. It is a conditional-on-existence size bound, not an existence theorem for every q.

For odd q no trailing padding is possible: every exact-U=q positive tail seed with q>1 has a strong representative of order at most R(floor(q/2)). For a general odd **slope** c, U=c2^j allows unbounded j. Stripping exterior padding leaves the same odd part c, but does not bound the strong core order in terms of c. Therefore universal odd-slope seeds remain open.

## 5. Executable controls

`core_check.py` uses literal permutations and split locations, independently of the endpoint-DP formulas, to compute H, endpoint counts, and every ordered-cover K_ab. It checks all labelled tournaments of orders 1 through 5 and every selected subset:

- strict cross-component bound (3) for every pair in different SCCs;
- injection bound (5) for every selected subset;
- all positive unit seeds have exactly one selected SCC and singleton exterior;
- exact U=2^b U_core and unit offset of the extracted core;
- the strong-core classifier recovers the prior all-order q=1,...,13 decisions, and extends the complete fixed-q classification through 29.

For q≤29, every nontrivial permitted core has H≤14 and hence order at most five under the stated inherited floor (order-six strong floor 15). The exact positive q≤29 are {1,2,4,6,8,11,12,16,20,22,24,29}; every other q in this interval is excluded, with the same exact-q scope as `AUDIT.md`. Literal witnesses for c=5,U=20 and c=29,U=29 are included, with the first two tail totals independently counted as (41,81) and (59,117). These are witnesses and fixed-q minima where stated, not comparisons to the mixed-mechanism ray bank.

    python3 core_check.py > core_receipt.json
    python3 -O core_check.py > core_receipt_optimized.json
    cmp core_receipt.json core_receipt_optimized.json

Every check raises RuntimeError explicitly and remains active under optimization.

### Explicit checked core witnesses

Outgoing-bit rows and zero-based selected vertices:

- c=5: S=(20,1,2,7,14), O={1,2}. The exact data are H=9,W=8,U=20=5·2^2. Thus its tails give N=1+5·2^s for every s≥3; literal first totals are 41 and 81. The complete q≤29 exclusion of q=5 and q=10 proves that exponent j=2 is the earliest possible exact multiplier for this slope **within the selected-set tail family**. This says nothing about earlier individual ray values realized by other mechanisms.
- c=29: S=(24,1,3,6,14), O={1,2}. Here H=11,W=10,U=29. Its tails give N=1+29·2^s for all s≥1, starting at 59 and 117.

Both cores have five vertices. No core of order at most four has a multiplier that can be padded to 20 or 29, as the exhaustive profiles show; this is a minimum-core-order statement within this normal form, not a mixed-mechanism ray minimum.

Recorded replay: all 1,099 labelled bases, 33,866 selected subsets, 3,212 cross-component unordered pairs, and 2,885 positive unit seeds through order five passed. Normal and optimized JSON outputs are byte-identical.
