# Half-supported reroutings force an effective deficit-two cutoff

Status: new proof candidate; independent symbolic and finite verification received; final audit requested. Local six-vertex lemma is finite-exact and is supplied with an exhaustive integer verifier. No remote publication. 2026-10-20.

Let T be self-converse, A=Aut(T), and Delta=H(T)/|A| minus the number of blue path orbits. A directed Hamilton word P is blue exactly when its positional reversal r_P is an anti-automorphism. The induced converse/reversal action on path orbits fixes exactly the blue orbits, so Delta=2 means precisely one pair of nonfixed black orbits.

## 1. Half-supported even-order comparison

Take two Hamilton words P,Q. Write Q=qP, with q a permutation of vertices. Suppose the support S of q satisfies S intersect r_P(S)=empty. If q has even order, the following three events are all impossible:

1. P,Q are both blue. In that case r_Q=q r_P q^{-1} and r_P are anti-automorphisms, so r_Q r_P=q(r_P q^{-1}r_P) is an automorphism. Its two factors have disjoint support and the same order as q. Therefore their product has even order, impossible in a tournament automorphism group.
2. P,Q lie in the same path orbit. Then q is itself an automorphism, impossible because it has even order.
3. The path orbits of P,Q are interchanged by converse/reversal. Fix any anti-isomorphism sigma of T. Then Q=a sigma r_P P for some a in A, so q r_P=a sigma is an anti-automorphism. Its square q(r_P q r_P) is an automorphism with two disjoint factors of even order, again impossible.

The no-even-automorphism fact is elementary: a suitable power of an even-order automorphism is a nontrivial involution; each transposed pair has its arc reversed.

Consequently, for Delta=2, any graph of Hamilton words whose edges satisfy these half-supported even-order conditions is bipartite, colored blue versus black. An edge cannot have two blue endpoints by (1); its black endpoints would have to satisfy (2) or (3), since there are only the two converse-paired black orbits.

More generally, a clique of t such words gives Delta>=2(t-1): at most one word is blue, and the black words occupy different path orbits, none paired with each other.

## 2. The four-cycle prefix certificate

Suppose a four-vertex induced tournament S is strong, so it has a directed four-cycle C. Suppose a vertex v outside S is beaten by at least three vertices of S, and T-S has a directed Hamilton word beginning with v. If n>=8, then Delta(T)>=4.

Proof. Cut C immediately after any vertex of S that beats v, and append the common Hamilton word on T-S. There are at least three such words. The relative permutations between distinct four-cycle cuts are nontrivial powers of a four-cycle, hence have order two or four. Their supports lie within the first four positions, disjoint from their positional reversals when n>=8. Apply the clique bound with t=3. If all four vertices of S beat v, the stronger Delta>=6 follows. QED.

## 3. A finite local configuration lemma

FINITE-EXACT LEMMA: every strongly connected tournament on six vertices contains a strong induced four-vertex tournament S and a vertex v outside S that is beaten by at least three vertices of S.

The supplied verifier checks all 2^15 labeled tournaments on six vertices. It independently tests forward and backward reachability to define strongness. Exactly 22,320 inputs are strong; all have the desired witness. This is a complete finite proof by enumeration, not sampling. A separate 512-case verifier adjoins two vertices to the unique strong tournament on four vertices; 450 extensions are strong and all have an explicit witness. The two producer verifiers use different strongness tests (two-sided reachability versus all dominating cuts). An independent auditor has also reproduced both scans with different matrix/closure code.

By vertex-pancyclicity of strong tournaments, every strong tournament of order at least six has a strongly connected induced six-vertex subtournament: take any directed cycle of length six, whose induced tournament remains strong. Therefore the local lemma holds at every order at least six.

## 4. Delta=2 forces connectivity at most four

Let T be strong and self-converse with n>=8 and Delta=2. If its vertex connectivity were at least five, apply the local lemma to obtain S,v. Removing the four vertices of S leaves a strong tournament; Camion gives a Hamilton cycle there, whose cut beginning at v gives the required Hamilton word. Section 2 contradicts Delta=2. Therefore

    kappa(T)<=4.

Here kappa(T) is the minimum number of vertices whose deletion makes T nonstrong. A one-vertex tournament is considered strong. The n>=8 condition is needed for disjoint prefix/reflected-prefix supports; it is not silently applied at n=6 or7.

## 5. Rotation symmetry gives a connectivity lower bound

Let C be any Hamilton cycle in T with successor R. Let |A intersect <R>|=m, q=n/m, so that R^q is an automorphism of odd order m. Its q vertex orbits O_0,...,O_(q-1) have size m, and each induced tournament T[O_i] is regular. The cycle C contains a perfect directed matching from O_i to O_(i+1 mod q).

A regular tournament of odd order m has vertex connectivity at least ceil(m/3). Indeed, suppose deletion of s vertices leaves a nontrivial split X->Y. Put x=|X| and y=|Y|. The average outdegree of vertices in X is at least (x-1)/2+y, and that of vertices in Y is at most (y-1)/2+s. Regularity, with common degree (m-1)/2, implies y<=s and x<=s. Thus m=x+y+s<=3s.

Delete fewer than ceil(m/3) vertices from T. Each remaining orbit induces a strong tournament by the preceding argument. None is empty. Each perfect matching between consecutive orbits retains at least one edge because fewer than m vertices were deleted in total. Therefore the union remains strong, proving

    kappa(T)>=ceil(m/3).

This proof also handles q=1, where there is only the regular orbit. For m=1 the assertion is simply kappa(T)>=1, inherited from strongness.

### 5.1. Sharpening to ceil(m/2) when q>=2

In fact q>=2 gives

    kappa(T)>=ceil(m/2).

Suppose deletion of s<m/2 vertices leaves X->Y, with X,Y nonempty. No rotation orbit is empty because s<m. Call an orbit split if it meets both X and Y. For each split regular orbit O_i the preceding average-degree proof gives a_i=|X intersect O_i|<=s_i, b_i=|Y intersect O_i|<=s_i, hence s_i>=m/3. There can be at most one split orbit because two would require s>=2m/3>m/2.

If there is no split orbit, the cyclic list of orbit labels X,Y has a Y-to-X transition. Its perfect matching cannot have all m edges incident to the fewer than m deleted vertices, contradiction.

Suppose O_0 is the unique split orbit. The unsplit orbit labels O_1,...,O_(q-1) cannot have a Y-to-X transition for the same matching reason. Thus they consist of an X segment followed by a Y segment, with either segment possibly empty.

- If all these labels are X, the matching O_0->O_1 sends every Y vertex in O_0 into a deleted vertex of O_1. Thus b_0<=s_1. Also a_0<=s_0 by regularity. Hence m=a_0+b_0+s_0<=2s_0+s_1<=2s, contradiction.
- If all labels are Y, the matching O_(q-1)->O_0 gives a_0<=s_(q-1), and regularity gives b_0<=s_0. Again m<=2s.
- If both segments are nonempty, q>=3 and O_1 and O_(q-1) are distinct. The same matching constraints give b_0<=s_1 and a_0<=s_(q-1), so m=a_0+b_0+s_0<=s_0+s_1+s_(q-1)<=s, contradiction.

This proves the sharpening. Since m is odd, the condition s<ceil(m/2) is precisely s<m/2.

## 6. Effective cutoff and residual support

For Delta=2, the inherited cycle-cut quotient lemma gives q<=4. In addition m is odd and q divides n. For q>=2, sections 4 and5.1 give ceil(m/2)<=4, hence m<=8 and, because m is odd, m<=7. Therefore n=q m<=28.

For q=1, the weaker bound ceil(n/3)<=4 gives n<=11. The inherited odd-cycle congruence restriction excludes n=11, while the previously verified small orders exclude every such center through order nine. Thus q=1 cannot occur for any still-unknown center.

The previously verified small orders supply the unique S4 center and no other center through order nine. The inherited odd/even cycle-cut restrictions sharpen the possible new orders to:

- Odd: n in {15,21}, with q=3, m in {5,7}.
- n=2 mod4: n in {10,14}, with q=2, m in {5,7}.
- n=4 mod8: n in {12,20,28}, with q=4, m in {3,5,7}.

Thus D(z) has degree at most28 and support contained in {4,10,12,14,15,20,21,28}. This does NOT determine the seven unknown coefficients or prove D=z^4.

## 7. Boundaries and next obligations

S4 has Delta2 and connectivity1, so connectivity5 is genuinely a sufficient obstruction rather than an equivalence. The half-support condition cannot be omitted: regular C5 is pure-blue despite its Hamilton words having even-order relative permutations with overlapping reflected support. The complete order-six configuration scan is essential here; it is not replaced by an unsupported assertion about an arbitrary four-subset.

Next: independently verify the local finite lemma and symbolic group/degree proof; seek a small human case analysis and a strengthened κ threshold; exploit the finite q-cyclic residual presentations instead of enumerating all tournaments through28.

## 8. Primary source for cycle extension

J. W. Moon, "On subtournaments of a tournament", Canadian Mathematical Bulletin 9 (1966), Theorem 1, p.297, states that each vertex of a strong n-vertex tournament lies on a directed cycle of every length k=3,...,n. The k=6 case supplies the induced strong six-vertex subtournament; k=n supplies Hamiltonicity and arbitrary starting cuts. DOI https://doi.org/10.4153/CMB-1966-038-7 ; primary PDF https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf . The independent auditor verified the exact theorem statement and page.
