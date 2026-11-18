# An invisible ordinal reorder exposed by selective recontact

Status: PROVED elementary composition identities and a complete ±2 classification for the stated operation. FINITE-EXACT tests below. Not a universal Chan–Pak selector. No priority claim; the identities recover familiar slot-cover mechanisms in a useful combined form.

## Inheritance

* `05-knowledge/results/arithmetic_seams_20260921_tournaments.md`, sections 4–6: fixed-path interval versus endpoint families, odd-cycle overlaps, and the warning that substitution need not multiply Hamilton counts.
* `01-canon/theorems/THM-160-overlap-weight-complement-formula.md`: complementary subtournaments are the sidecar lost by scalar cycle-overlap counts.
* `05-knowledge/results/procgen_smallgraph_20260926_small_graphs_arithmetic.md`: the small-graph arithmetic viewpoint and distinction between numerical balance and physical endpoint placement.

The niche combination here is not another seed census: retain the full slot polynomial under a neutral reorder, then expose exactly the information a selective contact reads. The local reference checkout is old and dirty; these are provenance, not claims of current novelty.

## Definitions and the neutral operation

Let H(T) count directed Hamilton paths, with H(empty)=1. Let P_q(T) count ordered q-tuples of directed paths, allowing empty paths, whose vertex sets partition V(T). In particular

P_2(T) = sum_{S subset V(T)} H(T[S]) H(T[V(T)\S]).

For the ordinal sum A→C, all cross-arcs point from A to C. Every path has at most one segment in each block. Cutting the q paths at that seam, or concatenating their corresponding segments, gives the bijection

P_q(A→C)=P_q(A)P_q(C).

Consequently arbitrary permutations of ordinal blocks preserve **every P_q**, not merely H and P_2. This does not preserve labeled endpoint data, nor does it assert isomorphism of the resulting tournaments. For example a singleton→directed triangle and triangle→singleton are distinct converse classes, with identical all-slot counts.

## The selective-contact compiler

Take nonempty tournaments A,C, a singleton b, and form A→b→C. Add x with x→A, b→x and x→C. Denote this genuine tournament F(A,C). Write a=H(A), c=H(C), u=P_2(A), v=P_2(C). Then

H(F(A,C))=c(u+a),
H(F(C,A))=a(v+c),
H(F(A,C))−H(F(C,A))=cu−av.                    (1)

Proof: C is a terminal ordinal block of F(A,C), so contributes c. In the remaining A+b+x tournament, paths either start at x (a choices), have all A before b→x (a choices), or pass b→x and then a nonempty suffix path in A. The third class is an ordered two-path cover of A whose second slot is nonempty; it has u−a choices. Thus the inner count is 2a+(u−a)=u+a. The reverse order proves the second identity.

The unrecontacted sources A→b→C and C→b→A have identical all-slot counts. Yet the new count sees the determinant of the two profiles (H,P_2). Therefore an all-slot scalar invariant does not suffice for arbitrary selective recontact. The missing coordinate is which physical block supplies a specified contact region.

## Complete classification of the ±2 case

For every tournament T with n>=2, P_2(T) is divisible by 4. Partition the subsets S into complementary pairs. Each product H(T[S])H(T[V\S]) is odd by tournament path parity (including empty sets). There are 2^(n−1) pairs, an even number, so their sum is even and P_2, twice that sum, is divisible by four.

If |A|,|C|>=2, both terms cu,av in (1) are divisible by four. Hence (1) cannot equal ±2.

If A is a singleton, then a=1,u=2, and the absolute determinant is

v−2c = sum_{empty≠S≠V(C)} H(C[S])H(C[V\S]) >= 2^|C|−2.

Every summand is positive. For |C|=1 this is zero. For |C|=2 it is exactly two. For |C|>=3 it is at least six. Thus |H(F(A,C))−H(F(C,A))|=2 **if and only if** the block sizes are {1,2}; every two-vertex tournament is transitive. The two resulting counts are precisely 3 and 5. Empty modules were excluded; allowing one gives an odd determinant and hence no additional ±2 case.

This is an all-order classification, not an inference from the finite scan.

## Explicit first pair and the constructive survivor

Use vertices A, then b, then C, then x. With A=[0] and C=[0,1] in outgoing-bitmask notation, F(A,C) has row masks

[14, 28, 0, 4, 13],

and F(C,A) has row masks

[12, 13, 24, 0, 11].

Their counts are 3 and 5. The replay prints these matrices directly, so the label convention is independently checkable.

For A=TT_r,C=TT_s, H=1 and P_2(TT_r)=2^r. Formula (1) gives

H(F(TT_r,TT_s))=2^r+1,
H(F(TT_s,TT_r))=2^s+1.

Thus the same fixed source size r+s+2 can recontact to two Fermat-block-type counts. The first +4 pair uses {r,s}={2,3}, giving 5 and 9 at order seven. It is a physical recontact, not a single-arc reversal. This recovers a useful mechanism, but not new density or prime coverage.

The exact prime obstruction is stronger than congruence: the output factors as c(u+a). Since u+a>=3, a prime output requires c=1, meaning C is transitive. Thus this particular design cannot use a nontrivial terminal Hamilton-count reservoir to place a prescribed prime. The next productive variation must destroy that terminal ordinal cut, for example by one selective return contact from C into A or into the new reader. Merely enlarging A and C leaves the factorization unchanged.

## A second independent constraint: source-order parity modulo three

Using the repository slot specialization P_q(T)=G_T(2q,q), where a spanning odd-cycle/singleton packing with k cycles and s singletons contributes x^k y^s, one gets

P_2(T) ≡ (−1)^|T| H(T) (mod 3).

Indeed k+s≡|T| mod2. Modulo three the P_2 monomial is (−1)^s and the H monomial is (−1)^k, differing by (−1)^|T|. Therefore an odd-order A forces P_2(A)+H(A)≡0 mod3. Any prime output other than 3 of F(A,C) requires C transitive **and A of even order**. The claimed slot specialization is an inherited input, not reproved here; the ±2 classification above is independent of it. This recovers the same parity/mod-three filter that earlier endpoint and floor families encounter, now as a direct filter on a physical source compiler.

## Falsifiable next bridge

Keep the neutral ordinal reorder but add a second contact which removes C as a terminal module. Search for pairs of resulting contact patterns whose difference is the determinant above **plus an independently controlled mixed-cover term of residue 2 modulo 4**. A working bridge must simultaneously (i) break the terminal factor c, (ii) retain the absolute center, and (iii) keep the contact-region labels needed for a repeat. This is a specified open construction target, not a claimed theorem.

## Verification

`python niche_bridge/replay.py` enumerates all labeled tournaments through order five (1099 total), checks P_2 divisibility and all profiles, tests literal recontacts of all profile representatives whose total core order is at most eight, and tests the TT_r/TT_s formula with r+s<=8. The separate exploratory `probe6.py` scanned all labeled tournaments through order six and found 39 distinct (H,P_2) profiles; the only determinant-two pair was (1,2),(1,4). This evidence is not needed in the proof.
