# A growing induced sevenfold lift and its inverse contact budget

Use outgoing masks, with bits numbered from zero. Begin with R₀=(38,20,32,7,77,26,47). Replace its vertex 6 by any nonempty tournament B, uniformly relative to the other six vertices, to form R. Add two vertices q,d: q beats precisely old vertices 0 and 1 and loses to all other old vertices; d beats only q. Then

    H(R)=34P₂(B)+5H(B),      H(R+q+d)=7(34P₂(B)+5H(B)).

Here P₂ is the number of ordered two-path covers, with empty slots allowed. The old R is literally induced, with no old reversals. Both skeletons are strong (checked directly on their finite quotients, which implies strength after any nonempty substitution). The entire R is not a module in the larger tournament; q cuts across it.

Vertex 6 has exactly one in-neighbor in each quotient, so any path intersects the B block in at most two runs: each entry except an initial run requires that sole predecessor. The one-run patterns multiply H(B), and the two-run patterns multiply P₂(B)−2H(B). With B replaced by a singleton the two counts are 73 and 511. With B replaced by a transitive pair they are 141 and 987. A transitive pair has H=1 and P₂−2H=2. Thus the pattern counts are (73,34) and (511,238), proving the identities. This is not merely a nonuniform extension of a fixed seed: both variable coordinates can change.

Examples of prime quotients include 73, 277, 1093 and 4357 (transitive B of size 1,3,5,7). A count-five strong block followed by five ordinal singletons gives (H,P₂)=(5,1024), hence the prime 34841 with sevenfold companion 243887. No claim of infinitely many prime quotients on any one ray.

## General architecture

Fix a deletion R, with J=H(R) and marked set Q of its vertices. Add a return block B with R→B→q. The new q beats precisely Q in R. For R define W as the number of ordered path covers (X,Z) of R with Z nonempty and starting in Q; X may be empty. Define U as the subcount with X nonempty and also ending in Q. Write h=H(B), P=P₂(B). Then the new count is

    (P−h)W + h(J+W−U).

Indeed if q ends the path the contribution is Jh. Otherwise q precedes the tagged start of Z. If it is reached from X or starts the path, the term is h(W−U) (the return block comes wholly last). If it is reached from B, the first B-run is nonempty and precedes q, with an optional B-run at the end, giving (P−h)W. For Q empty the new count is hJ and cannot be 7J. Henceforth Q is nonempty.

Every Hamilton path of R can be cut just before any chosen tag, giving W≥|Q|J. Cutting only at the first tag gives W−U≥J. Also P≥(|B|+1)h by splitting each full path at each cut. Thus H(new)≥[(P−h)|Q|+2h]J. If h≥3, then |B|≥3 and this exceeds 7J. Hence h=1. By the positive odd-cycle-pack formula (or a directed triangle and induced extension monotonicity), this is a transitive B, say of size k, with P=2^k. The inequality forces k≤2. If k=2 then |Q|=1 and U=0, so exact equality would require 4W=6J, impossible since J is odd. Thus any exact 7J lift here must add exactly q and a single vertex d. Its necessary and sufficient test is

    2W−U=6J.

The permitted number of tags is at most **four**. To see this, the bound above first gives |Q|≤5. If |Q|=5, equality forces W=5J and W−U=J. The first equality says that every W-cover concatenates: all J cuts at each tag already occur in W. The second says that the tags are contiguous in every full path: each tag-run start counts in W−U. This implies Q is homogeneous. Otherwise an outside vertex of mixed adjacency gives a directed tag→outside→tag path, extendible by single-vertex insertions preserving the relative order. To see extension, insert each missing vertex at the start if it beats the first, at the end if it loses to the last, or at a forward transition where it loses to one vertex and beats the next.

Fix the exterior prefix/suffix of any full path around its Q-run. By homogeneity, every ordered nonempty two-cover of Q can now be embedded as (prefix+X,Z+suffix) in W. The only possible failure of concatenation is between the internal X and Z; all W-covers concatenate. Therefore every nonempty two-cover of Q joins in its indicated order, giving P₂(Q)=(5+1)H(Q). But P₂ of every tournament on at least two vertices is divisible by four, as each unordered complementary subset pair contributes 2 times an odd number and there are an even number of such pairs. This contradicts 6H(Q)≡2 (mod 4).

Equivalently write A=W−|Q|J and D=(W−U)−J, two nonnegative integer excesses. The exact lift test is A+D=(5−|Q|)J. These are two coupled cover defects, not freely chosen digits.

## Prime inverse: the only two shapes

Let J=p be prime, p≥5. The only nontrivial SCC C of R has count p. Any prefix preceding the first tagged SCC is an ordinal prefix even in the extended tournament and can be trimmed if trivial. If the prefix contained C the remaining count would be 7, impossible. If a tagged singleton preceded C, its W contribution alone would be at least P₂(C): take the tag as the start of Z and distribute C over X and Z, placing earlier singleton vertices in X and assigning other later vertices to Z. More formally in ordinal SCC products a tag in component F contributes H(prefix)·W_F(tag)·P₂(suffix). Since |C|≥4, P₂(C)≥5p, with strict inequality by divisibility by four if |C|=4. This contradicts W≤5p. Thus we may delete an untagged singleton prefix and take C to be the first SCC and tagged.

An untagged trailing singleton suffix of length s multiplies W by 2^s, leaves U unchanged, and leaves J unchanged. (For U all suffix vertices must go to Z, since X ends at a tag.) At length s≥1 and |Q|≥2, this forces 2W−U≥(2^{s+1}−1)|Q|p+p>6p, using W₀−U₀≥p before that suffix. For one tag U=0 and 2^{s+1}W₀=6p is impossible. Hence the final SCC is tagged.

There cannot be three or more trailing singleton SCCs after C: a single tag in C would already contribute at least 8p. With exactly two singleton SCCs, the tag in C and the last singleton contribute at least 4p+p, exhausting the W≤5p budget. The middle singleton must be untagged, so all full paths have at least two tag runs, contradicting W−U=p. With exactly one singleton after C, it must be tagged. Three or more tags in C contribute at least 6p before the final tag; two tags in C force W≥5p and hence one tag run on every full path. That would make Q homogeneous, impossible for a set meeting only part of C and its following singleton: an omitted vertex of C sees the singleton outward, so must beat all included vertices in C, contradicting strength of C.

After stripping the irrelevant prefix there are therefore precisely these possible shapes:
1. R=C with 1–4 tagged vertices. The test is 2W−U=6p. For one tag this is the exact weighted-start condition W=3p.
2. R=C→{y} with exactly two tags x∈C and y. Set w=W_C({x}) and e equal to the number of full C paths ending at x. Then W=2w+p, U=w+e. To count U, either X ends at y and Z starts at x (w ways), or X ends at x and Z is the singleton y (e ways). The test is 3w−e=4p.

These are iff tests within the specified architecture, not universal exclusions on 7p. In particular shape 1 holds for the variable B family above whenever its quotient is prime. The challenging part is furnishing the profile at the prescribed p.

## A second rail sharing the variable block

Use base rows (116,57,34,181,36,0,62,119), substitute B at vertex 7, and mark 1 and 5 for the two-vertex extension. The old count is 22P₂(B)+5H(B) and the new count is seven times this. The old tournament has a terminal singleton (vertex 5) and an otherwise strong part; the new tournament is strong. The one-entry proof is the same, using the literal counts 49 and 93 for the old graph and 343 and 651 for its lift. The case B singleton gives a literal induced 49→343; B=TT₃ gives 181→1267 and, in the first scaffold, 277→1939.

These two independent paired lifts share the same B. A prescribed pair of ordinary counts (p₁,p₂) determines its required profile uniquely:

    P₂(B)=(p₁−p₂)/12,     H(B)=(17p₂−11p₁)/30.

These equations are necessary and sufficient exactly when their recovered pair is physical; in particular p₁−p₂ is divisible by 48 if |B|≥2. This is a decoder for the required profile, not a profile selector.

For a nontransitive example, the first scaffold with the eight-vertex H223 tile followed by three singleton SCCs uses (H,P₂)=(223,32768) and yields the prime 1115227 and its sevenfold companion 7806589. The H599 tile with the same padding gives prime 2231219 and 15618533. The accompanying script checks the tile profiles, arithmetic primality, literal counts on smaller substitutions, and the needed one-entry and strength conditions.

## Why ordinal growth cannot preserve the same exact contact lift

Write D=2W−U−6J. In each marked block also write E,S for counts of full paths ending/starting at a tag, and Z for the end-tag analogue of W (the tagged path is the first slot, and the other slot may be empty). Under R→A the exact concatenation formulas are

    W'=W_R P_A+J_R W_A,
    U'=J_A U_R+J_R U_A+E_R S_A+W_R Z_A.

For U, the two internal cases contribute when both marked endpoints are in the same block; if they are split across blocks, either each slot is wholly in one block or both slots cross the cut. This gives the last two terms, including cases where a slot's untagged segment is empty. Consequently

    D(R→A)−J_A D(R)
       = W_R[2(P_A−J_A)−Z_A]+J_R(2W_A−U_A)−E_R S_A.

For nonempty old tags and nonempty A the first term is strictly positive since Z_A≤P_A−J_A. The last two terms together are nonnegative: if A has any tags then 2W_A−U_A≥W_A≥J_A, while E_R≤J_R and S_A≤J_A; if A has none both vanish. Similarly prepending a tagged nonempty A to nonempty R gives correction W_A(2P_R−Z_R)−J_R U_A−E_A S_R≥W_A(P_R−J_R)>0 to J_A D(R). An untagged prefix scales D exactly. Thus an exact lift cannot be renewed by an ordinal suffix with unchanged old tags, even if the new suffix gets tags; the successful variable substitution is internal, not an ordinal decoration of the marked source.
