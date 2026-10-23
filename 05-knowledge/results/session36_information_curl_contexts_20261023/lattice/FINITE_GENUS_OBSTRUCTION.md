# No fixed root-shell theta depth determines Hamiltonian count

Session36. Status: PROVED elementary nonconstructive existence theorem. This is about bounded-genus coefficients restricted to the norm2 root shell of A_(n-1), not the full decorated lattice theta series. The initial proof allows arbitrary tournaments; §7 strengthens the result to strongly connected tournaments with identical labelled scores. No effective small witness or optimal depth bound is claimed.

## 1. The finite observation

For a tournament T on n vertices, let R_n={e_i-e_j:i!=j} be the roots of A_(n-1), with the usual Euclidean metric and alternating form omega_T(x,y)=x^T M_T y. For each j<=g, retain the full joint distribution of the symmetric Gram matrix and alternating pairing matrix over all ordered j-tuples of roots, allowing repetitions.

Denote the combined signature by F_g(T). This is exactly the norm2-root-shell portion of the joint theta coefficients at genera1,...,g from `METRIC_CURL_RECOVERY.md`. Labels of the original vertices are not retained; the whole distribution is invariant under relabelling.

THEOREM. For every fixed integer g>=1, there is an order n and two tournaments T,U of that same order such that

F_g(T)=F_g(U), but H(T)!=H(U).

In fact the construction below forces such collisions for every sufficiently large value of an integer B, and hence infinitely many n in an arithmetic progression. Thus no fixed root-shell observation depth can recover H at all orders. This complements the sufficient depth n-1 proved for each fixed n; it does not say that n-1 is optimal.

## 2. Polynomially many signatures at fixed depth

For roots x,y, the metric inner product belongs to {-2,-1,0,1,2}. Their alternating value is an integer between-4 and4, because x^T M y is a sum of four entries of M with signs. Therefore, for each j, there are at most

45^(j(j-1)/2)

possible joint bins. Diagonal metric entries are2 and alternating diagonal entries0; the off-diagonal entries above the diagonal determine the entire pair. Many of these45 choices cannot occur together, which only makes this upper bound looser.

Put D_g=sum_(j=1)^g 45^(j(j-1)/2). Every count is a nonnegative integer at most

|R_n|^j=(n(n-1))^j <= n^(2g).

Hence the number of possible complete signatures is at most

(n^(2g)+1)^D_g.

The crucial statement is polynomial growth in n for each fixed g, even though the exponent is enormous.

## 3. A supply of multiplicatively independent Hamiltonian counts

Let N_a be the tournament on a+2 vertices consisting of a transitive chain A={v_1,...,v_a} and two additional vertices b,c, with A->b->c->A. This is the nearly transitive tournament obtained from a transitive order by reversing its first-last arc, after relabelling. Its H formula is inherited in THM-1370, `THM-1370-h-spectrum-omits-7-21-all-n.md`, and THM-4744, `THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`; here is an independent direct proof of precisely the formula needed.

LEMMA. H(N_a)=2^a+1 for a>=1.

Proof. Vertex b has only c as an outneighbor, and c has only b as an inneighbor. There is one Hamiltonian path starting b (b,c, then all of A in its transitive order), and one starting c (c, all of A, b). Every other Hamiltonian path begins with a nonempty subset S of A in its unique increasing order, then b,c, then A\S in its increasing order. All nonempty S are allowed, including S=A. These2^a-1 paths plus the first two give2^a+1. QED.

For distinct odd primes p,q greater than3,

gcd(2^p+1,2^q+1)=3,
and v_3(2^p+1)=1.

Proof of the gcd: if d divides both odd numbers, then the multiplicative order of2 modulo d divides both2p and2q, hence divides2. Thus d divides2^2-1=3; both numbers are divisible by3. For the valuation, the binomial expansion of(-1+3)^p+1 is congruent to3p modulo9. Since p!=3, it is not divisible by9.

Consequently b_p=(2^p+1)/3 are integers greater than1, pairwise coprime, and coprime to3. Products of the numbers2^p+1 over a fixed collection of such primes, with arbitrary nonnegative multiplicities, determine those multiplicities uniquely: choose any prime divisor of b_p and read its exponent.

No prime-density estimate or primitive-divisor theorem is used; only infinitely many primes are needed.

## 4. Fixed-order family with too many H-values

Fix g and choose k=2gD_g+1 distinct primes p_1,...,p_k greater than3. Put

C=sum_(i=1)^k (p_i+2).

For each B>=1 and vector m=(m_1,...,m_k) in{0,...,B}^k, form a tournament T_m by taking an ordinal sum, in any fixed declared ordering, of:

- m_i copies of N_(p_i), and
- B-m_i transitive tournaments of order p_i+2,

for every i. All tournaments T_m have the same order n=BC. Hamiltonian paths multiply under ordinal sum: a directed spanning path must pass through the factors in their given order, and any collection of factor paths concatenates. Since a transitive factor has H=1,

H(T_m)=product_i (2^(p_i)+1)^(m_i).

By §3 these(B+1)^k values are pairwise distinct. Hence these are also(B+1)^k different tournament isomorphism classes.

On the other hand, there are at most((BC)^(2g)+1)^D_g possible signatures F_g. For B>=1,

((BC)^(2g)+1)^D_g <= 2^D_g C^(2gD_g) B^(2gD_g).

Since k=2gD_g+1, any integer B>2^D_g C^(2gD_g) makes B^k exceed this upper bound, and therefore also makes(B+1)^k exceed it. Pigeonhole supplies two different m,m' with the same F_g and different H. The displayed explicit threshold also proves infinitely many such orders n=BC. QED.

## 5. What this adds, and what it does not

- The genus3 Paley8 example shows an explicit short recovery step after universal genus2 blindness.
- This theorem prevents extrapolating that success to a universal fixed-genus root-shell formula for H.
- The full integral metric alternating pair still determines H, with genus n-1 root frames giving an upper bound on the needed depth at each n.
- Full genus2 or genus3 decorated lattice theta coefficients involve unbounded vector norms and are not covered by the finite-bin estimate. No conclusion is drawn about their completeness.
- The initial proof does not force strong witnesses, but §7 recovers a stronger fixed-score, strongly connected version using Moon interval chains.
- Counts of arbitrary bounded induced patterns have the same polynomial information budget; the precise root-shell theorem above avoids asserting a converse reconstruction of such decks.

## 6. Controls

`finite_genus_controls.py` independently counts Hamiltonian paths for N_a for1<=a<=12; verifies the gcd and3-valuation claims on explicitly listed primes; and checks injectivity of a finite multiplicity box. These are non-vacuity and arithmetic controls, not an enumeration of the enormous pigeonhole witnesses. The all-g theorem is the proof above.

## 7. Stronger version recovered from Moon's interval chains

THEOREM (strong, fixed-score horizon). The witnesses in §1 can be chosen strongly connected and with the same labelled score vector. Thus the information-depth obstruction is not an artifact of SCC multiplicativity, nor repaired by adding all vertex scores.

Recovered mechanism: THM-4744, `01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, §1 and proof1. Its product law is attributed there to Moon1972, not claimed as new. We use the general product law, not its separate finite inventory claims.

For a composition a=(a_1,...,a_s) with all a_i>=2 and sum N, the interval chain IC(a) has vertices0,...,N, every forward arc i->j (i<j) except that each interval from the preceding partial sum l to the next partial sum r has the endpoint arc reversed to r->l. It satisfies:

- IC(a) is strong, because each induced interval is a strong nearly transitive tournament and consecutive intervals share a vertex.
- H(IC(a))=product_i (2^(a_i-1)+1).
- Its labelled scores in the fixed spine order are s_0=N-1, s_t=N-t for1<=t<=N-1, and s_N=1, independent of the composition. Each shared block endpoint gains one outgoing arc from the preceding block and loses one from the following block.

For completeness, a proof of the nontrivial product identity is given in §8 below.

Choose the same k and primes as in §4. For each i provide B slots of total interval length p_i+1. In m_i of those slots put one part p_i+1, and in the remaining B-m_i slots put (p_i+1)/2 parts equal to2. The latter is integral because p_i is odd. Concatenate all these parts in a fixed order.

Every such chain has

n=1+B sum_i(p_i+1),

and the same labelled score vector. Its Hamiltonian count is

H=3^(sum_i[m_i+(B-m_i)(p_i+1)/2]) product_i b_(p_i)^(m_i).

The factors b_(p_i) are pairwise coprime and coprime to3, so again there are(B+1)^k different H-values. Put C'=1+sum_i(p_i+1); for B>=1 we have n<=BC'. The same signature bound and any B>2^D_g (C')^(2gD_g) force a repeated F_g with unequal H. Both witnesses are strong with identical labelled scores. QED.

## 8. Independent proof of the interval-chain path product

Induct on the number of parts. A one-part chain is N_(a_1-1), already counted in §3. For more parts, let A be the first interval, v its last vertex, and B the remaining chain including v as its first vertex. Set A'=A\{v}, B'=B\{v}; every cross-arc points A'->B'.

Every Hamiltonian path must therefore have the shape

X_A X_B v Y_A Y_B,

where X_A,Y_A partition A' and X_B,Y_B partition B', each run in its inherited order. This follows because a B'->A' transition is impossible on either side of v.

Such a path restricts to directed Hamiltonian paths X_A v Y_A in A and X_B v Y_B in B. Only two joins might fail:

1. If X_A is nonempty but its last vertex does not point to v, that last vertex must be the first spine vertex0, since v has only0 as an outneighbor within A. Vertex0 has only v as an inneighbor in the full chain, so it must be the first vertex of the entire path and X_A=(0). Then Y_A is nonempty (a_1>=2), but v cannot point to its first vertex because its only outneighbor in A' is the already-used0. Contradiction.

2. If Y_B is nonempty but v does not point to its first vertex, that first vertex must be w, the last vertex of the next interval, which is v's only inneighbor inside B. Therefore X_B is empty, since otherwise its last vertex would also have to be w. The interior of that next interval is nonempty. Consider its smallest-index interior vertex. It lies after w in Y_B and cannot have a predecessor: its possible inneighbors are v, vertices of A', and smaller-index interior vertices, all already before w or nonexistent. Contradiction.

Conversely, any pair of factor Hamiltonian paths P_1 v P_2 in A and Q_1 v Q_2 in B has the unique merge P_1 Q_1 v P_2 Q_2. Every cross join uses a forward A'->B' arc, and the two joins at v come from the factor paths; empty runs cause no difficulty. The decomposition and merge are inverse. Thus H(A union_v B)=H(A)H(B), and induction proves the formula.

The finite controls include all compositions with parts>=2 of totals2,...,12, checking strong connectivity, the labelled scores and product count with a separate subset Hamiltonian-path DP. They do not instantiate the huge horizon pigeonhole pair.

## 9. Eventual failure at every order, and polynomial observation budgets

COROLLARY. For each fixed g, every sufficiently large order n admits two strong tournaments with the same labelled scores and the same F_g but different H.

Proof. Keep k and the primes fixed, and put L=sum_i(p_i+1). For a target order n take

B=floor((n-1)/L)-1.

For sufficiently large n, B>=1. The residual r=n-1-BL lies between L and2L-1, so r>=2. Every integer r>=2 is a sum of2s and3s (use only2s when even, and one3 plus2s when odd). Append one fixed such residual composition to every chain in §7. It contributes a common positive factor to all H-values, so their(B+1)^k-way distinction remains. The score vector remains the universal interval-chain score vector at this new order n.

We have n<=(B+2)L<=3BL. Therefore the signature count is at most2^D_g(3L)^(2gD_g)B^(2gD_g). Any B>2^D_g(3L)^(2gD_g) gives the same strict pigeonhole inequality. B tends to infinity with n, proving the eventual every-order statement. QED.

More generally, the argument applies to any fixed finite collection of integer-valued observations whose possible joint signatures grow only polynomially with n. Choose k larger than that polynomial exponent and use the same interval-chain family. Thus finitely many bounded-degree pattern counts or bounded finite-shell joint coefficients cannot determine H uniformly. This last sentence concerns observations with a proved polynomial information budget; it does not include an unbounded theta series, all degrees simultaneously, arbitrary real encodings, or an observable that already encodes H.


## 10. A small explicit strong fixed-score hostile at genus2

FINITE-EXACT: IC(2,2,4) and IC(2,3,3) both have order9, the labelled score vector(7,7,6,5,4,3,2,1,1), c3=7 and c4=10, but H=81 and75 respectively. The full norm2 root-pair distributions agree, both by direct enumeration and by the general(c3,c4) characterization in `METRIC_CURL_RECOVERY.md`. Both tournaments are strong. `finite_genus_controls.json` gives their exact outneighbor masks and independently recomputed counts. This is not a minimality claim among all tournaments or all interval chains.
