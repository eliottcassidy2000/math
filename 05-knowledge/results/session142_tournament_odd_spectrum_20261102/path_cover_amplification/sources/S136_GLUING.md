# Module-context tomography and finite spectral substitution algebras

**Status: PROVED (all-order statements below); FINITE-EXACT (explicit grammar closures).**
Universal realization of every positive odd Hamiltonian-path count except 7 and 21 remains OPEN.
No external priority claim is made.

## Inheritance and actual advance

- THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`, already proves the complete path-cover substitution kernel and its exact C3 diagonal. We specialize, not rediscover, that law.
- THM-3134, `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`, identifies the path-cover profile with an endpoint jet and gives the same-H masks 40/76 used below. That hostile is inherited.
- THM-3235, `01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md`, distinguishes a factor profile from the quotient as a substitution operator. We preserve that distinction.
- THM-4804, `01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md`, uses the no-immediate-re-entry lemma for modules. The two-exterior-vertices-per-return mechanism is inherited at lemma level; the exact context-equivalence filtration, decoder, sharp H floor, finite algebra, and all-depth grammar certificates are the present consequences.
- S133, `05-knowledge/results/session133_fractran_spectrum_refinement_20261101/core_reduction/README.md`, already records the transitive-module hostile H=2^n+1 and finite prescribed-spectrum decidability via a reduced-core bound. Neither the hostile family nor finite unrestricted decidability is claimed new here. The new finite algorithm respects a specified substitution grammar.
- S134, `05-knowledge/results/session134_constructive_spectrum_frontier_20261101/README.md`, provides marked fresh-block grammars and capped arithmetic. The current carrier is instead arbitrary homogeneous substitution over a finite quotient bank.

The new question is not whether a profile composes (known), but **which profile coordinates are necessary for which contexts, and when an infinite exact context space becomes finite at a fixed spectral cutoff**.

## 1. Objects and context coefficients

For a nonempty tournament T, p_T(c) counts covers of V(T) by c unordered nonempty vertex-disjoint directed paths; p_T(1)=H(T). Set p_T(c)=0 beyond its order.
A one-module context is a tournament Q with distinguished vertex x. Q[T] replaces x by T, every other Q vertex being a singleton; all exterior vertices see the entire block T homogeneously. Its exterior size is m=|Q|-1.

Let W_c(Q,x) count directed words in Q which use every vertex other than x exactly once and x exactly c times. Adjacent letters are distinct. Define a_c=c!W_c. Splitting an actual Hamiltonian path at maximal T runs gives the inherited specialization

    H(Q[T]) = sum_{c>=1} a_c(Q,x) p_T(c).                 (1)

The c actual path components are unordered, and c! bijects them to the ordered x occurrences. In particular a_1=H(Q). The coefficients are nonnegative integers and depend on Q,x, not T.

Between consecutive x occurrences there must be at least two exterior letters: a single letter y would require both x->y and y->x. These gaps are disjoint. Consequently

    a_c=0 if 2(c-1)>m.                                   (2)

This is a homogeneous-module statement, not a claim for arbitrary vertex sets or arbitrary apex attachments.

## 2. Exact bounded-context equivalence and minimal detectors

For m>=0, call T and U m-context-equivalent if H(Q[T])=H(Q[U]) for every one-module context with at most m exterior vertices. Then

    T ~_m U  iff  p_T(c)=p_U(c) for 1<=c<=floor(m/2)+1.  (3)

**Proof.** Sufficiency follows from (1)-(2). For necessity the zero-exterior context recovers p(1). For each 1<=s<=floor(m/2), use the actual context

    D_s[T] = C3[T, L_s, L_s],                            (4)

where L_s is transitive of order s. It has 2s exterior vertices. Its coefficient of p_T(s+1) is

    (s+1)! (s!)^2,                                      (5)

and all higher coefficients vanish. Indeed, in the inherited C3 kernel the only possibility with c_T=s+1 is c_2=c_3=s; its quotient-word multiplicity is one, and each L_s has exactly one all-singleton cover. Equation (5) is nonzero. Induction recovers p(2),...,p(floor(m/2)+1). QED.

Thus if the first difference between two profiles occurs at c, the **minimum exterior size of any detecting homogeneous context is exactly 2c-2**. For c=1 the identity context suffices. For c>1, (2) proves the lower bound and (4) with s=c-1 attains it. This conditional statement does not assert that actual tournament pairs first differ at every conceivable c.

Consequently complete profile equality is equivalent to exact H equality in all one-module contexts. This is a minimal sufficient invariant in the operational sense, not an isomorphism classifier and not a claim that profiles determine the tournament as an outer quotient.

The first decoder rows, in p(1),p(2),... coordinates, are

    D_1: (3,2)
    D_2: (11,34,24)
    D_3: (51,458,1296,864).

The recursive decoder divides the residual D_s response by (5); for actual profiles this division is exact. All entries are positive integers, but inverse decoding involves subtraction and is not a positivity/realizability theorem.

## 3. Interaction criterion and sharp module-size floor

Call x **active** if it lies in a nontrivial strong component of Q, equivalently on a directed cycle. Then

    a_2(Q,x)>0 iff x is active.                          (6)

If x is inactive all a_c for c>=2 vanish.

**Proof.** A repeated-x word contains a directed return, so makes x cyclic. Conversely a cyclic vertex in a tournament belongs to a directed triangle: choose a shortest directed cycle through x; any length >=4 can be shortened using the chord between x and the vertex two steps later, or the appropriate first transition from out-neighbors to in-neighbors along the cycle. Write the triangle x->a->b->x. Replace x by two twins x1->x2. The directed path (x1,a,b,x2) has separated twins. Insert every remaining quotient vertex, one at a time, using the elementary tournament path-insertion lemma. Insertion preserves the original subsequence, so the twins never become adjacent. Collapsing their names gives a word counted by W_2. Thus a_2=2W_2>=2. QED.

For n=|T|,

    p_T(2) = (1/2) sum_{empty != A proper subset V(T)}
                       H(T[A]) H(T[V(T)\A])
             >= 2^(n-1)-1.                              (7)

For n=1 the sum and p_T(2) are zero. Every induced nonempty tournament has at least one Hamiltonian path, giving the inequality. A transitive T attains it. This does not assert that equality in (7) alone forces transitivity (C3 is a hostile at n=3).

If x is active then H(Q)>=3: Q contains a directed triangle and insertion injects its three paths into those of Q. Equations (1), (6), and (7) yield

    H(Q[T]) >= 3H(T)+2p_T(2) >= 2^n+1.                  (8)

This is sharp for every n>=1:

    H(C3[L_n,1,1]) = 3+2(2^(n-1)-1)=2^n+1.              (9)

Hence every nonempty proper homogeneous module M in a strong tournament S satisfies

    |M| <= floor(log_2(H(S)-1)).                         (10)

Indeed collapse M; the quotient remains strong, so the collapsed vertex is active. This bound is module-specific, not a replacement for the inherited global minimum-H theorem for strong tournaments. It is sharp in (9).

## 4. Infinite exact states versus finite capped states

Every transitive L_n has H=1, yet (9) distinguishes all of them using just two exterior vertices. Therefore no finite exact-H context quotient exists, even restricted to the H=1 fibre. Source/sink padding is H-neutral for ordinal concatenation but is not context-neutral under future homogeneous substitutions. A reduced-core certificate is not an exact substitution state.

Fix an integer B>=1 and put cap(z)=min(z,B+1), beta=B+1. Let

    N_B=min{n>=1 : 2^n+1>B}.                             (11)

For n>=N_B every active one-module context saturates at beta by (8); every inactive context returns cap(H(Q)H(T)). Thus large modules need only capped H. For n<N_B retain the complete profile (or its entrywise cap) and the small-order tag. Define the finite signature

    sigma_B(T)= small(n, cap(p_T(1)),...,cap(p_T(n)))  if n<N_B,
                large(cap(H(T)))                    if n>=N_B.       (12)

Two identical signatures have identical capped H in every homogeneous context. For small inputs, use (1) and nonnegative integer capped arithmetic; for large inputs, use the active/inactive rule above. The tag is essential. The state is sufficient, not claimed to be the minimal contextual quotient. The threshold (11) is sharp uniformly because of (9). In the H=1 fibre the capped contextual quotient has **exactly N_B classes**: every H=1 tournament is transitive (a directed triangle would force H>=3), the orders1,...,N_B-1 have distinct unsaturated responses in (9), and all orders>=N_B share the large(H=1) state. This minimality statement is confined to the H=1 fibre.

There are at most

    sum_{n=1}^{N_B-1} (B+2)^n + (B+1)                    (13)

formal states; actual tournament signatures form a smaller subset. In particular the true capped-output contextual equivalence relation has finite index. The bound is deliberately coarse.

## 5. A finite algebra for every substitution grammar

The signature (12) composes under arbitrary tournament substitution, not merely under bounded-size contexts. Here is an explicit transition for any fixed quotient Q and nonempty child states.

1. If every child is small and their total order n is <N_B, compute the full output profile using the nonnegative integer THM-3121 kernel, with arithmetic capped at beta. Output small(n,profile).
2. Otherwise the output is large. If any large child occupies an active vertex of Q, return large(beta). To justify (8) when other children are non-singletons, first expand every other block: the chosen large block is still one homogeneous module at a cyclic quotient vertex.
3. If no large child is active, every such child is a singleton strong component of Q. Factor H across Q's ordered strong components. Each large singleton component contributes its capped H; a component whose children are all small is evaluated from their capped profiles by THM-3121. Multiply these contributions with capped arithmetic. Return large(the product).

These cases are exhaustive. They do not need the unknown exact order/profile of a large child. Capped addition and multiplication are exact for nonnegative integer expressions. Equivalently, replacement of one child at a time preserves all contextual capped H by Section 4; the explicit transition additionally computes the small output sidecar.

For any specified finite bank of quotient operations and finite seed bank, start from seed signatures and repeatedly apply every operation to all tuples of retained signatures. The finite closure terminates. By induction every retained state has an actual construction witness (soundness); induction on expression depth puts every constructible tournament in the closure (completeness). Thus the algorithm decides the grammar's **entire unbounded-depth spectrum up to B**, not only a depth-limited search. It says nothing about tournaments outside the specified grammar.

## 6. Exact all-depth grammar certificate at B=127

`module_contexts.py` implements the theorem for the seed K1 and operations J(A,B)=A->B and C(A,B,C)=C3[A,B,C]. Here N_B=7. For J, the complete profile transfer is also checked directly:

    p_{A->B}(d)=sum_{a,b,t: a+b-t=d}
       p_A(a)p_B(b) binom(a,t)binom(b,t)t! .              (14)

Choose t components on either side and match them; each matched pair concatenates. The other components remain separate.

The implementation retains exact small profiles, a finite refinement of (12). Its closure has 58 states; small profile counts at orders 1,...,6 are (1,1,2,3,6,13). Its exact attained values <=127 are

    1,3,5,9,11,15,17,23,25,27,33,37,45,47,51,55,
    65,69,75,81,85,91,93,95,99,107,111,115,121,123,125.

In particular 13 is impossible in this grammar at every order and depth. It is not a universal spectral hole: the five-vertex tournament with adjacency bitmasks

    (18,8,3,5,14)

has H=13, independently checked by subset DP. Its pair-bit mask is41 (pairs lexicographic, bit1 means low index -> high index). The grammar therefore cannot settle the full-spectrum conjecture, even though it has unbounded depth and block size.

The machine output also records a second closure whose seed bank contains every tournament through order five, with the same J/C operations. Its scope is still exactly that seed/operation grammar. This second closure has65 states, small counts(1,1,2,3,8,15), and adds exactly13,39,117 through127;19 is still absent at every depth. These are completed finite fixed-point computations, not predictions beyond127.

## 7. Reproduction and hostile boundaries

Run from any directory:

    python3 -B module_contexts.py
    python3 -B -O module_contexts.py

The primary audit checks all 1,099 labelled tournaments through order five by two profile engines: canonical-subset partitions using subset Hamilton DP, versus literal permutations and cuts. It checks all 1,099 labelled one-hole contexts through order five against seven distinct factor profiles through order four, giving 7,693 direct vertex-level substitution checks. The active criterion is checked independently using transitive closure; 23,079 threshold-state checks use six cutoffs. Five triangular decoder rows and eight padding hostiles are included. Normal and optimized output match byte-for-byte. Every recorded attainable value is additionally reconstructed from its expression witness and checked by vertex-level Hamilton DP:31 witnesses in the first grammar and34 in the second, maximum order12.

Inherited masks40/76 on five vertices both have H=15, but profiles (15,39,33,10,1) and (15,45,35,10,1), and triangle responses123/135. This is the first-derivative loss already proved in THM-3134, retained as a hostile control rather than claimed anew.

**Loss ledger.** Source: actual homogeneous factor tournament. Target: all capped or exact H responses under substitution. Map: maximal-run splitting. Kept: the required path-cover prefix, or the finite signature at fixed B. Lost: labels, isomorphism type, quotient-operator data, and nonhomogeneous attachment information. Boundary: increasing B requires a new cutoff/signature; finite grammar completeness does not imply unrestricted spectrum completeness. The pure H scalar is inadequate before capping, and arbitrary marks cannot be erased.

## 8. A spectral resource: minimum substitution width

Define w(H) as the least w for which some tournament with Hamiltonian count H can be built from K1 using arbitrary-depth substitution by quotient tournaments of order at most w. This is a scalar realization resource: it minimizes over tournaments with the same H, rather than assigning a width to one particular witness. It is finite for every attained H because that witness itself may be used as a quotient. No bounded-width realization theorem for the entire spectrum is asserted.

Every tournament on at most four vertices belongs to the J/C3 grammar. At order four the complete list is

    L4, C3->1, 1->C3, C3[L2,1,1].                       (15)

One can see this by choosing a directed triangle. The fourth vertex sees it uniformly, giving the second or third case, or nonuniformly, giving the last case; absence of a triangle gives L4. Substitution is associative, so allowing all quotients of order at most four is exactly the J/C3 grammar. The literal finite audit checks all64 labelled order-four tournaments against these four recipes, with multiplicities(24,8,8,24).

Since the completed J/C3 closure misses13 and19, but both have width-five witnesses,

    w(13)=w(19)=5.                                       (16)

For13 use the order-five mask41 above. For19, let Q be the order-five tournament of pair-bit mask10, and replace its vertex1 (zero-based labels) by L2. The resulting adjacency masks are

    (40,5,1,6,15,30),                                    (17)

and literal permutation counting gives H=19. The quotient itself has H(Q)=9.

This also gives a concrete **seed/operator separation**. The grammar with every tournament of order at most five supplied as a seed, but only J/C3 available as operations, still cannot attain19 at any depth. Allowing the one order-five mask10 as an outer quotient does attain19. A tournament supplied as a factor is not the same resource as the same tournament supplied as a substitution operator. This is the operational distinction already emphasized by THM-3235, here instantiated as an exact spectrum obstruction.

## 9. Complete width-five grammar at cutoff127; width(31)=6

**FINITE-EXACT, with inherited Moon–Busch pruning.** The supplementary `width5_probe.py` replaces a seed-only enlargement by all quotient operations of order at most five. It imports the standard-library companion `module_contexts.py`; it does not require repository computation code.

A proper nontrivial homogeneous module in an order-five tournament has size2,3,or4. Contracting that module leaves a quotient of order4,3,or2. The factor and quotient therefore both belong to (15). Thus every nonprime order-five tournament is already an operation expressible by J/C3. This direct argument avoids any need for an unproved all-order decomposition theorem.

The exhaustive1024 labelled order-five tournaments are tested for proper modules and grouped by literal permutation isomorphism. There are exactly264 prime labelled tournaments in three classes:

    pair-bit mask          10       41       76
    labelled multiplicity 120      120       24.         (18)

All three are strong. Hence J, C3, and these three order-five quotient operations generate exactly the width-at-most-five class; this is not merely a selected bank of convenient five-vertex quotients.

For evaluation at B127, retain the same finite signatures with N_B=7. The supplementary search prunes strong output total order at10 using the **inherited/cited** exact strong minimum from S133 `core_reduction/THEORY.md`:

    f(3k)=3*5^(k-1),  f(3k+1)=5^k,  f(3k+2)=9*5^(k-1).

In particular f(11)=225>127 and the displayed floor is nondecreasing. Substitution into a strong quotient is strong, so any C3 or prime-five output on at least11 vertices is necessarily saturated. Active large children are saturated separately by (8). Ordinal joins remain unpruned and use exact capped products. This is the sole external mathematical import needed for the supplementary size pruning; the earlier58/65-state closures do not use it.

For each prime quotient the script independently generates directed quotient-word multiplicities by last-letter dynamic programming, then contracts them against the factorial profiles using the inherited THM-3121 kernel. Child tuples are generated by total order<=10, rather than a blind Cartesian product. Each full small output (order<7) is also checked directly on its constructed graph.

The completed fixed point has86 states, with small profile counts

    (1,1,2,3,8,21) at orders1,...,6.

It attains49 odd values through127. Its complete missing odd list is

    7,21,31,35,43,53,67,71,77,83,89,101,103,113,127.     (19)

Every one of the49 positive values has an explicit adjacency witness in the output and a direct subset-DP Hamilton-count check; maximum witness order11. The graph reconstruction deletes only singleton SCCs (an inherited H-preserving operation), and re-pads to order7 when a large-state representative is required; it does not contract modules. These normalized witnesses remain width<=5: the class generated by all quotient tournaments of order<=5 is hereditary. Induct on an expression; deleting vertices deletes empty children, replaces the outer quotient by its induced subquotient (still of order<=5), and recursively takes induced subtournaments in surviving children. Adding the transitive prefix is again an allowed ordinal join. Thus graph normalization preserves both H and membership in the certified width class.

Finally the order-six tournament of pair-bit mask146 has adjacency masks

    (36,17,2,7,13,30),                                   (20)

and both literal permutation counting and subset DP give H=31. Equation (19) excludes width<=5, while (20), used as a quotient on six singleton children, gives width<=6. Therefore

    w(31)=6.                                            (21)

This is an exact scalar realization-resource threshold. Nothing here proves that w(H) is globally bounded or unbounded, that any fixed width has zero-density spectrum, or that all positive odd H except7/21 occur. Those questions remain OPEN. The finite state count grows with B; finite closure at one cutoff cannot be silently upgraded to a uniform-in-B assertion.
