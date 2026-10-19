# Tight suffixes give polynomial reconfiguration diameter for graceful upper-tail fibers

2026-10-19. **PROVED, using a cited strict-rank transversal theorem.** Independently [reviewed](review/INDEPENDENT_REVIEW.md). The result improves a sufficient reconfiguration diameter; it does not prove the Graceful Tree Conjecture or one-color connectivity.

## 1. The inherited mechanism and the missing decomposition

The repository's [earlier upper-tail proof](../graceful_upper_tail_two_color_connectivity_20261019/UPPER_TAIL_TWO_COLOR_CONNECTIVITY.md) proves that independent transversals of classes A_1,...,A_s satisfying

    rank(A_i) >= s-i+1

are connected by singles and closure-preserving double replacements in at most 2^s-1 moves. It then obtains graceful upper-tail width-two diameter at most 2^(r-1).

The subsequent [three-color release proof](../graceful_three_color_release_20261019/THREE_COLOR_RELEASE_DICHOTOMY.md) supplies the maximal tight-color-kernel reduction and imports Fon-Der-Flaass's strict-rank connectivity theorem. The missing combination is that the rank staircase forces **every tight color subset to be a suffix**. This turns the possibly complicated tight-kernel structure into a chain, so each successive block has a full-rank buffer and strict-rank remaining classes.

The closest proved mechanism is therefore tight-kernel contraction plus a connected buffer. The canonical hostile is the three-color cyclic matching obstruction, which has no full-rank color. The corrected near miss is requiring every initially correct target coordinate to stay fixed: the repository's five-state rank-staircase example disproves that. The present algorithm may disturb a correct buffer within its block; only the already completed inner blocks are retained.

Concept board: rank staircase; tight subsets; contraction by a fixed flat; strict-rank single-coordinate paths; buffer lifts; arithmetic residue chains. No tournament structure is introduced.

## 2. The cited strict-rank theorem

**Fon-Der-Flaass (1996), Theorem 1 and the distance bound following Theorem 2.** For n classes V_1,...,V_n in a matroid, if

    rank(union_{i in I} V_i) > |I|

for every nonempty I, their independent transversals are connected by single-coordinate replacements with distance at most

    g(n)=n²-n+1       (n>=1).

Set g(0)=0 for the empty system. The theorem allows non-simple matroids and distinguishable parallel elements. We use no stronger statement.

Primary source: D. Fon-Der-Flaass, *Local equivalence of transversals in matroids*, Electronic Journal of Combinatorics3(1),R24(1996), DOI https://doi.org/10.37236/1248; publication page https://www.combinatorics.org/ojs/index.php/eljc/article/view/v3i1r24 . The strict-rank statement is on PDF page2, and the bound on page3. On2026-10-19 the primary PDF's indexed text explicitly reproduced the bound `2+4+...+2(n-1)+1=n²-n+1`; a prior repository source audit records the precise theorem hypotheses. The [independent reviewer](review/INDEPENDENT_REVIEW.md) successfully retrieved and checked the [primary PDF indexed text](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v3i1r24/pdf/), including the exact hypotheses and bound. No saved PDF image is asserted.

The proof in the source constructs local replacements. The present note establishes a short-path theorem and a blockwise construction relative to that strict-rank path routine. Its finite checker may use BFS for that routine; it does **not** claim that a polynomial-time oracle implementation has been supplied or benchmarked.

## 3. Tight sets are suffixes

Let M be a finite matroid and let pairwise disjoint classes A_1,...,A_s satisfy the staircase. A state chooses one actual element per class, with the chosen set independent. Write A_I=union_{i in I} A_i. Call I tight if rank(A_I)=|I|.

The staircase guarantees a state exists: choose elements from A_s,A_(s-1),...,A_1 in that order. At stage i, the s-i previously chosen elements cannot span a class of rank at least s-i+1.

**Lemma1.** Every nonempty tight subset is a suffix {i,i+1,...,s}.

**Proof.** Put i=min I. Then

    rank(A_I) >= rank(A_i) >= s-i+1 >= |I|.

Equality rank(A_I)=|I| forces |I|=s-i+1; since I is contained in that suffix, it equals the suffix. QED.

Consequently all tight subsets are nested. Their ranks are determined by only s suffix-union rank tests. Include the empty suffix {s+1,...,s}. Let K be the largest nonempty tight suffix, if one exists; otherwise K is empty. Its flat H=cl(A_K) is spanned by the selected K coordinates in EVERY state.

The remaining classes outside K satisfy strict rank in M/H. Indeed

    rank_(M/H)(A_J)=rank_M(A_(K union J))-|K| >= |J|,

and equality would make a larger tight subset than K. The inequality follows because any full state supplies an independent selection on K union J.

## 4. Consecutive tight suffixes produce strict-rank blocks

Take consecutive tight suffixes

    K_inner={j,...,s} subset K_outer={i,...,s},

where no suffix strictly between them is tight; j=s+1 is allowed. Put B={i,...,j-1}, b=j-i, and contract H_inner=cl(A_K_inner). Restrict to A_B if desired.

**Lemma2.** In this quotient:

1. rank(A_B)=b;
2. rank(A_i)=b, so its first color is full rank in the block;
3. every nonempty J subset B\{i} has rank(A_J)>|J|.

**Proof.** The difference of the two tight ranks gives1. Since contraction lowers a set's rank by at most rank(H_inner)=s-j+1,

    rank_(M/H_inner)(A_i) >= (s-i+1)-(s-j+1)=b.

Together with1 this gives2. Independent transversals give the weak rank inequalities in3. Equality for J would make K_inner union J a tight set in M. By Lemma1 it would be a suffix strictly between the two displayed suffixes, impossible. QED.

Every block selection in a state is therefore a basis of the fixed block flat in the quotient. Changing it to another such basis preserves the flat spanned by all K_outer coordinates; any outer coordinates remain independent.

## 5. Quadratic reconfiguration theorem

Define

    c(1)=1;
    c(b)=(b-1)²-(b-1)+1=b²-3b+3   for b>=2.

**Theorem3.** Any two independent transversals of a rank-staircase system on s colors are connected using singles and closure-preserving changes of at most two coordinates. If its consecutive tight-suffix blocks have sizes b_1,...,b_t and its prefix outside the largest tight suffix has a colors, a path has length at most

    g(a)+sum_j c(b_j) <= g(s)=s²-s+1   (s>=1).

For s=0 the unique state needs no moves.

**Proof.** First process tight blocks from the innermost suffix outwards, keeping completed inner coordinates at their target values. A size-one block is changed directly to its target in at most one move; its closure is fixed.

For a block of size b>=2, delete its first-color coordinate temporarily. Lemma2 gives the strict-rank hypothesis for the remaining b-1 classes in M/H_inner. The cited theorem supplies a path of at most g(b-1) single-coordinate changes between the current and target partial selections.

Lift each partial step as follows. Its new partial selection has rank b-1. The first color has rank b, so choose an element of that color outside the new partial closure. Simultaneously change the affected partial coordinate and, if necessary, the first-color coordinate. The resulting block is a basis of the same fixed block flat. This is a legal move on at most two coordinates, and after adjoining the unchanged inner and outer coordinates its full selected-set closure is unchanged.

On the final partial step, choose the exact target first-color element, which is available because the new partial selection is then the target. Thus no extra final buffer step is needed when the partial path has positive length. If it has zero length, change only the buffer to target. The bound is max(1,g(b-1))=c(b).

After the tight kernel has been placed at target, the a outside classes satisfy strict rank in M/H. Their g(a)-step single-coordinate path lifts with the entire kernel fixed. These singles need not preserve the selected-set closure; the theorem only requires that double replacements do so.

For the final bound, c(b)<=g(b) for every b>=1. For positive x,y,

    g(x+y)-g(x)-g(y)=2xy-1>=1.

Since a+sum b_j=s, repeated merging bounds the sum by g(s). QED.

This proves a polynomial bound under precisely the earlier rank-staircase hypothesis. It does not require an added target-monotonicity hypothesis or the strict-rank condition on the original full system.

## 6. Every graceful upper-tail fiber has quadratic diameter

Fix labels0,...,m and any acyclic choice F of one edge of every difference color r+1,...,m. Keep F fixed. Contract it in the graphic matroid. The quotient has r+1 vertices and rank r, and the free class A_d consists of every actual color-d edge, d=1,...,r. Loops may be discarded; parallel realizations remain distinct.

The full distance-d graph on0,...,m has exactly d components, its residue chains modulo d. Adding F cannot increase that count. Therefore

    rank(A_d)>=r+1-d.

Color1 is connected, so the union of all free classes has rank r: the entire color set is tight. There is no outside prefix a in Theorem3.

**Corollary4.** Every such upper-tail completion fiber is connected under moves changing at most two free-color edges, with diameter at most

    sum_j c(b_j),

where b_j are its tight-suffix block sizes. In particular, the uniform bounds are

    0                         if r=0;
    1                         if r=1 or r=2;
    r²-3r+3                  if r>=3.

The block sum may be sharper than the uniform bound. The r<=2 bound follows directly by replacing all differing free edges at once. For r>=3, any partition of r has sum c(b)<=max(r,c(r))=c(r): merge all parts of size>=2, then merge any singleton parts into that part, using c(x+y)>=c(x)+c(y) except at x=y=1. If all parts are singletons the sum is r<=c(r).

Only the path-length bound is being improved: width-two connectivity was already proved. The arithmetic theorem still changes the abstract tree type and does not supply a graceful labeling of every abstract tree. The factorial-star family in the preceding packet has exact diameter r-1, so this uniform quadratic upper bound is not asserted optimal.

## 7. Decisive boundaries and algorithmic status

The connection is from arithmetic difference classes to a rank-staircase matroid system, then to its ordered chain of tight flats. It preserves independent transversals and exact allowed color coordinates. Contracting a tight flat forgets its internal basis, retained as an actual fixed-coordinate sidecar during the lifted moves.

The generic three-color cyclic assignment obstruction has class ranks2,2,2 in rank3, so it fails the first staircase requirement. It cannot be imported as a counterexample here. The existing target-preservation hostile survives: a block buffer may have to leave its already-correct value.

Computing the tight-suffix decomposition requires s rank tests once suffix unions are available. Buffer lifts require rank or independence tests on actual elements. The remaining strict-rank path routine is the cited constructive theorem. The finite verification included here uses complete state graphs to test that interface and the lifts on bounded systems. No claim is made that the provided BFS checker is a polynomial-time implementation for arbitrary inputs. Turning the source proof into a carefully tested polynomial-oracle implementation is a separate next step; the present universal polynomial diameter does not rest on that implementation.

## 8. Bounded exact controls

Run from this directory:

    python3 -B verify_polynomial_staircase.py
    python3 -B -O verify_polynomial_staircase.py

The two receipts agree. The checker performs29,727 explicit checks and validates2,158 full moves. It tests200 seeded binary linear-matroid systems, including loops and distinct parallel realizations, with1,000 endpoint routes. It independently enumerates all gracefully labeled interval trees through labels0,...,7, groups them by every fixed upper tail, and checks821 routes across585 nonempty fibers. For each system it directly tests every color subset to confirm the suffix classification, all strict-rank conditions of invoked blocks, exact target restoration, full independence, double-move closure preservation, and the stated bounds. It also retains the five-state target-preservation hostile and the cyclic three-color obstruction outside the hypotheses.

These finite controls test the decomposition and lift implementation; they do not establish unbounded quantifiers by extrapolation. The proof and the cited strict-rank theorem establish those quantifiers. Receipt: `controls.json`.
