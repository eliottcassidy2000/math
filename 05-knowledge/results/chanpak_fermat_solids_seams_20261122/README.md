# Which rotations can survive at a prescribed perfect-number neighbor?

**PROVED.** For every exponent \(p\ge2\), the upper neighbor \(E_p+1\) of \(E_p=2^{p-1}(2^p-1)\) excludes every Fermat prime from its automorphism support. At \(p=61\) **every** tournament with this upper count has trivial automorphism group, regardless of its order. At \(p=31\) the group can only be trivial or cyclic of order 29, with a single 29-cycle in the latter case. The first composite Fermat number sharing a factor with any upper target is \(F_{18}\), detected by the Lucas–Lehmer doubling map started at \(-3/2\). A sharp \(2\sum q_i\)-vertex construction gives independent rotations of prescribed prime orders in modularly prime tournaments. The results select the necessary geometry of possible corrections. Absolute one-arc placement at exponent 61 is open here. [The actual exponent-31 edges](../chanpak_e31_physical_20261122/README.md)

## A guard on every prime order, independent of vertex count

Write \(H(T)\) for the number of directed Hamiltonian paths. Every subgroup of \(\operatorname{Aut}(T)\) acts freely on the ordered vertex lists of these paths, so its order divides \(H(T)\). Let \(q\) be a prime dividing the group order. It is odd, since an involution exchanges some two vertices and reverses their arc. By Cauchy's theorem an automorphism of order \(q\) exists. Every nontrivial orbit of this element induces a vertex-transitive regular tournament on \(q\) points. [LEM-003, free action on paths](../../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md)

Set \(C=q(q^2-1)/24\) and \(D=(q^2-1)/8\). Such an orbit contains \(C\) directed triangles, since each transitive triple has a unique source and each vertex beats \((q-1)/2\) others. Transitivity gives exactly \(D=3C/q\) directed triangles through each point. After choosing \(j\) disjoint triangles we exclude at most \(3jD\) possible next triangles. For \(1\le t\le\lceil q/9\rceil\) the product below therefore bounds the number of **ordered** disjoint \(t\)-tuples from below. The odd-cycle packing formula counts each unordered such collection with weight \(2^t\). Define the integer
\[
 L(q)=1+\sum_{t=1}^{\lceil q/9\rceil}
         \left\lceil \frac{2^t}{t!}\prod_{j=0}^{t-1}(C-3jD)\right\rceil .
\]
If the element has \(m\) nontrivial orbits, taking any collection independently in each shows
\[
                       H(T)\ge L(q)^m.
\]
The term 1 permits an empty collection. Every factor \(C-3jD=(q^2-1)(q-9j)/24\) increases throughout its positive range, so \(L\) increases across odd primes. This gives a finite target selector using only integer operations. It grows superexponentially: taking \(t=\lfloor q/18\rfloor\) gives a term at least \((C/t)^t\) because \(C-3jD\ge C/2\) there. [THM-002, odd-cycle packing formula](../../../01-canon/theorems/THM-002-ocf.md)

For an elementary coarser bound, a maximal packing of triangles in the orbit leaves a transitive subtournament. Its source beats the rest, so its size is at most \((q+1)/2\). At least \(\lceil(q-1)/6\rceil\) triangles were packed. Every subcollection contributes and their weights sum to at least \(3^{\lceil(q-1)/6\rceil}\) per orbit.

**Two exact applications.** At \(p=61\), \(L(89)>E_{61}+1\), and direct modular division by every prime through 83 yields no divisor. Thus the automorphism group must be trivial. At \(p=31\), \(L(53)>E_{31}+1\); the only dividing prime through 47 is 29, with valuation one. Also \(L(29)^2>E_{31}+1\), so a nontrivial 29-rotation acts on just 29 vertices. The remaining vertices are fixed; each sees the 29-cycle uniformly. Thus, if the tournament has more than 29 vertices, its unique nontrivial orbit is a proper module. The proof does not need the factorization of either large target.

## Constructible polygons and the two flanks

Let \(F_i=2^{2^i}+1\) for \(i\ge0\). The identity \(8(E_p+1)=(2^{p+1}-1)^2+7\) gives the inherited splitting restriction that every prime factor other than 7 of \(E_p+1\) belongs to \(1,2,4\bmod7\). Since \(2^i\equiv1,2\bmod3\) alternately, \(F_i\equiv3,5\bmod7\). These are nonsquares. Consequently no whole Fermat number divides the upper neighbor, even when the Fermat number is composite. In particular every odd straightedge-and-compass-constructible polygon order is coprime to the upper neighbor: the Gauss–Wantzel condition makes its prime factors Fermat primes. The stronger assertion of coprimality with all *composite* Fermat numbers does not follow from this argument. [Quadratic norm argument](../chanpak_perfect_wieferich_seams_20261122/PROOF.md) · [Wantzel's 1837 paper](https://www.numdam.org/item/JMPA_1837_1_2__366_0/) · [Earlier Fermat-five ladder](../fermat_five_pepin_collatz_20261025/FIVE_LADDER.md)

For odd \(p\ge3\), the lower factorization is
\[
                    E_p-1=(2^{p-1}-1)(2^p+1).
\]
For every prime \(q\mid F_i\) the order of 2 modulo \(q\) is exactly \(2^{i+1}\), since its \(2^i\)-th power is \(-1\). If \(i\ge1\), the second factor on the right cannot have such a divisor: it would force \(p\equiv2^i\pmod{2^{i+1}}\). The first factor contains the whole \(F_i\) precisely when \(2^{i+1}\mid p-1\); when this fails, it is coprime to \(F_i\). For \(i=0\), both factors contain 3. Hence among Fermat primes the possible lower-neighbor factors are exactly the initial indices \(i<v_2(p-1)\) for which \(F_i\) is prime. For \(i\ge1\) and \(d=2^{i+1}\mid p-1\), the exact multiplicity of each \(q\mid F_i\) is \(v_q(F_i)+v_q((p-1)/d)\) by lifting the exponent, as \(2^d-1=F_i(F_i-2)\). Thus exponent 31 admits just 3 on the lower side, and exponent 61 admits 3 and 5. This concerns possible group-order divisors, not sufficient realization conditions.

The rotational groups of the Platonic solids have orders 12, 24 or 60. Any nontrivial subgroup contains an element of order 2, 3 or 5. None can embed in the tournament automorphism group at \(E_p+1\). This excludes even the triangular and pentagonal cyclic rotations that a general tournament admits. In the six-label family \(4,2,5,3,7,13\), the sole possible upper prime is 7, and it occurs to the first power exactly when \(p\equiv2\bmod3\). An arbitrary subgroup supported on these labels is consequently trivial or cyclic of order 7. The norm interpretation and the Leech-lattice glue use different maps and carry different information. [Prior Platonic rotation classifications](../session190_constructive_links_20261104/groups/session190-groups/PROOF.md) · [Six Leech deep-hole types and glue](../chanpak_leech_perfect_phases_20261122/README.md)

## The first partial Fermat intersection and a Lucas–Lehmer orbit

**PROVED.** Put \(s_0=1\) and \(s_{i+1}=s_i^2-2^{2^i+1}\). For every prime power \(q^a\mid F_i\),
\[
   q^a\mid E_p+1\text{ for some }p\ge2
                   \quad\Longleftrightarrow\quad q^a\mid s_i.
\]
Thus \(\gcd(F_i,s_i)\) is the least common multiple of the individual overlaps \(\gcd(F_i,E_p+1)\) over all exponents. A single gcd replaces a scan of \(2^{i+1}\) phases. Normalizing \(b_i=s_i/2^{2^{i-1}}\) for \(i\ge1\) gives **exactly the Lucas–Lehmer doubling map** \(b_{i+1}=b_i^2-2\), here started from \(b_1=-3/2\). The target predicate is common divisibility by a Fermat factor. [Classical Mersenne Lucas–Lehmer iteration](https://www.mersenne.org/various/math.php)

For the proof set \(n=2^i\), \(d=2n\), and let \(\alpha,\beta\) be the roots of \(f(X)=X^2-X+2\). They satisfy \(\alpha\beta=2\) and \(s_i=\alpha^n+\beta^n\). Also \(f(2^p)=2(E_p+1)\), and 2 has exact order \(d\) modulo \(q\) as well as \(q^a\). If \(q^a\mid f(2^p)\), the roots mod \(q^a\) are \(2^p\) and \(2^{1-p}\). Their \(n\)-th powers are opposite signs since \(2^n=-1\) there. Thus \(q^a\mid s_i\).

Conversely, first work over the algebraic closure of \(\mathbb F_q\). Vanishing of \(s_i\) gives \(\alpha^n=-\beta^n\) and hence \(\alpha^{2n}=\beta^{2n}=-2^n=1\). The \(d\)-th roots of unity are exactly the powers of 2 in \(\mathbb F_q\). The roots of \(f\) are distinct, since its discriminant is \(-7\) and \(7\nmid F_i\). Both lift uniquely to roots modulo \(q^a\). Now \(q^a\mid s_i\) again gives \(\alpha^{2n}=1\) in the lifted ring. Every \(d\)-th root of unity mod \(q\) has a unique lift, since \(q\nmid d\). As 2 has order \(d\) modulo \(q^a\), these lifts are precisely its powers. Thus one exponent works and the two phase classes are \(p,1-p\pmod d\).

**FINITE-EXACT.** The first nontrivial overlap is at **\(F_{18}\)**. All gcds at indices \(0\) through \(17\) equal one, while
\[
                   \gcd(F_{18},s_{18})=13,631,489=q.
\]
This is prime. The two exponents modulo \(524288\) are \(96230\) and \(428059\). For example, the prime exponent \(9865243\) lies on the latter phase. This establishes divisibility for these exponent shapes without requiring the corresponding Mersenne number to be prime.

The phase continues to every power of \(q\). Indeed
\[
          (2^{524288}-1)/q\equiv725967\pmod q,
\]
so \(q\) is non-Wieferich in base two: its prime-to-\(q\) exponent quotient \((q-1)/524288\) preserves the valuation one. For a solution mod \(q^a\), adding \(k\cdot524288\cdot q^{a-1}\) to its exponent changes \(f(2^p)/q^a\) modulo \(q\) by \(k\cdot725967\cdot2^p(2^{p+1}-1)\). This coefficient is a unit, since \(f\) has simple roots. Exactly one lift succeeds above each root at every level. The two classes mod \(524288q\) at level \(q^2\) are \(3578626406374\) and \(3568199698459\), respectively. This gives a concrete distinction between an entire Fermat number, its prime factors, and their higher-power visibility at the upper neighbor.

## Edit distance from independent rotations

Suppose \(T\) admits automorphisms of prime orders with pairwise disjoint moved supports, with \(m_q\) chosen generators of order \(q\). For any tournament \(U\) on the same vertices of Hamilton count \(V\),
\[
 d_{\rm flips}(T,U)\ \ge\ \left\lceil\frac12\sum_q (m_q-v_q(V))_+\right\rceil .
\]
Indeed every rotation whose support avoids all changed endpoints still preserves \(U\); these rotations generate a direct product. At most \(v_q(V)\) chosen rotations of order \(q\) can survive untouched. Each edited arc touches at most two supports. This proof applies with any relabeling of either tournament. For a strong source with arbitrarily many protected rotations, substitute a directed triangle at each vertex of a strongly connected tournament on \(t\ge3\) vertices. The resulting strongly connected source requires \(\lceil t/2\rceil\) edits to reach *any* height coprime to 3.

**Sharp vertex budget for independent rotations.** Given any list of odd primes \(q_1,\ldots,q_t\) of length \(t\ge2\), the smallest order of a modularly prime tournament admitting order-\(q_i\) automorphisms with pairwise disjoint moved supports is exactly
\[
                              2\sum_{i=1}^t q_i.
\]
A module is a vertex subset seen uniformly by every outside point. A prime tournament has no module of size between two and its order minus one. Each support here is proper. A prime-order permutation with a single nontrivial orbit forces that orbit to be a module. Each generator therefore needs at least two full orbits, proving the lower bound.

For equality, form a block of size \(2q\) for each requested prime. On its two copies \(A,B\) of \(\mathbb Z/q\), orient \(A\) as the carousel with forward steps \(1,\ldots,(q-1)/2\), and \(B\) in the opposite direction. Across copies orient \(A_u\to B_v\) exactly when \(u=v\). Simultaneous translation of the labels preserves the tournament.

The block is prime. Two distinct vertices from one copy in a module force both matched vertices of the other copy to join it. If a module contains matched pairs on a nonempty proper label set \(I\), strong connectivity of the first carousel supplies \(k\notin I\) with \(A_k\to A_i\) for some \(i\in I\). This vertex loses to \(B_i\), so it too is forced, followed by its mate. Iterate to exhaust the block. For a mixed pair \(A_u,B_v\) at \(q\ge5\), there exists \(k\notin\{u,v\}\) beating \(u\) in the first carousel, since the indegree is at least two. Then \(A_k\) distinguishes the pair. At \(q=3\) the same proof applies unless \(v\) is the unique in-neighbor of \(u\). In that case the remaining third label \(k\) beats \(v\) in the first carousel, so \(B_k\) beats \(A_u\) and loses to \(B_v\). This deals with all mixed pairs.

Order the blocks. Between blocks \(i<j\) orient alike types (both \(A\) or both \(B\)) from \(i\) to \(j\), and unlike types in reverse. A module containing an interblock pair must contain a second vertex from the first member's copy: that copy sees the far member uniformly and has both arc directions at the near member. Its intersection with the first block is now a module of size at least two, so it contains the whole block. If instead two original points are in the same block, this same conclusion holds immediately. Every other point sees the two types of this block oppositely and must also join. The full tournament is prime and each independent block rotation survives. In particular for \(q_i=3\) these modularly prime sources require at least \(\lceil t/2\rceil\) edits to reach any count coprime to 3, on exactly \(6t\) points. This is unbounded even after restricting the lower source to have no proper module.

**FINITE-EXACT, a modularly prime strong source.** There is an 18-vertex example with three disjoint supports of order-3 rotations and no proper nontrivial module. Its vertices are six triples \(A_0,\ldots,A_5\), each cyclically ordered and internally directed as a triangle. Group the triples as \((A_0,A_1),(A_2,A_3),(A_4,A_5)\). In each paired block, a vertex of the first triple beats exactly its same-index mate in the second. Between triples belonging to different pairs, direct the entire bipartite set from higher to lower triple index, except for \(A_0\to A_2\) and \(A_0\to A_4\). Simultaneously rotating the two triples of any paired block gives a rotation supported on those six vertices. For each of the 153 vertex pairs, adjoining every outside vertex that distinguishes two members eventually recovers all 18 vertices. This proves primeness by the definition of a module. The literal path DP gives \(H=480496617\). Its upper scalar \(480496619\) lies in the certified spectrum and every single-arc edit misses it. [Corrected scalar prefix](../../../session203/coverage-recertification/README.md)

A protected *lower* host occurs at the actual exponent 61 and uniformly at \(p_s=60\cdot2^s+1\), \(s\ge0\), without postulating a Mersenne primality pattern. Factor
\[
 E_{p_s}-1=(2^{60}-1)\prod_{j=0}^{s-1}(2^{60\cdot2^j}+1)(2^{p_s}+1).
\]
In the first factor take two directed triangles and the seven-vertex carousel of height 175 (cyclic steps 1, 2, 3). The remaining factors of \((2^{60}-1)/(9\cdot175)\) are \(11,13,31,41,61,151,331,1321\), all individually attained. Each factor \(2^k+1\) is the count of the one-apex spike with just its top bit marked. Ordinally joining the witnesses realizes the lower target and supplies disjoint rotations of orders 3, 3, 7. The upper count excludes all three because \(p_s\equiv1\bmod3\). Thus this particular lower host, and any ordinal padding with singleton vertices, needs at least two flips before the upper height is possible. The modularly prime example above shows how the same obstruction operates without an ordinal decomposition. [Spike formula](../chanpak_perfect_wieferich_seams_20261122/PROOF.md)

Run `python3 -B 05-knowledge/results/chanpak_fermat_solids_seams_20261122/check.py` from the repository root for the integer thresholds, the tournament and its 153 module closures.

## Sparse asymmetric contacts at the absolute exponent-61 target

**FINITE-EXACT.** A tournament with transitive apex deletion and count \(E_{61}+1\) must have at least nine wins and nine losses incident to its apex. The bound allows arbitrary tournament order, adjacent contacts, and completely asymmetric mark positions. It supplements the earlier exhaustive radius test on centered symmetric masks. [Centered inverse](../chanpak_perfect_wieferich_seams_20261122/CENTERED-INVERSE.md)

Normalize the transitive spine word by stripping initial zeros and final ones, which contribute no inversions. With a nonempty nonconstant word, the largest marked exponent is \(x\) and position zero is unmarked. Even its inversion from \(x\) to zero contributes at least \(2^{x-1}\) when \(x\ge2\). Hence \(x\le121\) at the target. If \(x\le120\), at most eight marks have individual contributions of at most \(2^{120}+\cdots+2^{113}<2^{121}-2^{60}\), too small already before subtracting interactions. All remaining words to test have largest mark exactly 121.

After visiting some descending marks let \(A\) be the sum of their powers of two, \(h\) their current count, and \(v\) the smallest mark already present. Appending any smaller positive mark \(u\) changes the count by exactly
\[
                    d(u)=2^u-A/2^{u+1}-[u+1=v].
\]
The quotient is an integer; this is the weighted inversion formula including the distance-one correction. For \(1\le u<v\), \(d(u)\) is strictly increasing, even at the last step where the adjacency cost first appears. Hence the last unknown position is determined by binary search. When further positions below \(u\) remain, their total positive terms are at most \(2^u-2\). Their geometric negative terms have total magnitude at most \((A'+2^u)/2\), where \(A'=A+2^u\), since every partial address is below \(A'+2^u\) and \(\sum_{w\ge1}2^{-w-1}=1/2\). At most one additional adjacency correction occurs per remaining mark. These bounds justify both pruning steps in the checker. It visits the feasible prefix masks for each budget from one through eight, then inverts the final mark. All eight searches are empty. Reversing the spine and complementing all bits gives the converse tournament after relabeling, preserving the count and exchanging wins and losses. [Apex inversion formula](../chanpak_perfect_wieferich_seams_20261122/PROOF.md)

Reproduce with `c++ -O3 -std=c++17 sparse_upper.cpp -o /tmp/chanpak_sparse_upper && /tmp/chanpak_sparse_upper` from this directory. The positive control is the known three-mark word at exponent 8. The analogous absolute placement on general marked cores remains open.
