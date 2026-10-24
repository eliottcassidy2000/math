# An all-prime-size family: exact ordinary, self-converse and merged counts

Status: PROVED for the matching-coupled cyclic-regular family defined here. The counts are not A000568 and do not enumerate all tournaments of order2p. Burnside's lemma and the involution generating function are classical; the point is the explicitly typed family and its reconstruction obstruction.

## 1. Intrinsic cells and exact group action

Let p be an odd prime, m=(p−1)/2, and let R_p be the cyclic regular tournament on F_p with steps1,...,m. Given a permutation q of F_p, form T_q from two copies A,B of R_p, with A_i→B_j exactly when j=q(i). Every other cross edge goes B→A.

A vertices have degree m+1; B vertices have degree m+p−1. These degrees differ, so an ordinary isomorphism preserves the two cells. Inside R_p, the outgoing arc i→i+d has exactly m−d common outneighbors for1≤d≤m. Thus the short arcs i→i+1 are intrinsic, and every automorphism is a rotation. Conversely every rotation is an automorphism. Therefore Aut(R_p)=C_p.

It follows that ordinary isomorphism is exactly the action

q ↦ t_b q t_(−a),     a,b in F_p,

where t_a(i)=i+a. The family is strongly connected: each internal cell is strongly connected and both cross directions occur.

## 2. Ordinary count

The identity pair fixes all p! permutations. A pair with exactly one nonzero translation fixes none. If a,b are both nonzero, the fixed equation q(i+a)=q(i)+b determines q from q(0), giving exactly p solutions. Burnside yields

N_p = (p! + p(p−1)^2)/p^2
    = ((p−1)! + (p−1)^2)/p.

The p−1 affine slopes each give an orbit of size p. Every nonaffine permutation has trivial rotation-pair stabilizer and orbit size p². In particular the affine carrier exhausts the case p=3, but not any p≥5. It is invalid to replace S_p by Aff(F_p) at higher cell sizes.

## 3. Converse, involutions, and merged count

Write r(i)=−i. Exchanging the two cells and reversing their internal orientations sends the matching to

q ↦ r q^(-1) r.

This involution normalizes the rotation-pair action and thus acts on its N_p orbits. Let F_p count its fixed ordinary orbits, exactly the self-converse types in this family, and let M_p count converse-merged types. Then M_p=(N_p+F_p)/2.

To compute F_p, apply Burnside to the enlarged group of order2p². A permutation q fixed by a converse-coset element satisfies

q=t_b r q^(-1) r t_(−a).

Set u=q t_a r. The equation is equivalent to u²=t_(b−a). If a=b, the number of solutions is I_p, the number of involutions in S_p. If a≠b, t_(b−a) is a p-cycle. Any square root commutes with its square, so belongs to the cyclic centralizer C_p; because p is odd, there is exactly one square root there. Summing over the converse coset gives p I_p+p(p−1). Hence

F_p = (I_p+p−1)/p,
M_p = (N_p+F_p)/2,
I_p = sum_(j=0..floor(p/2)) p!/(2^j j! (p−2j)!).

The formula for I_p chooses j disjoint transpositions and leaves the other symbols fixed. Its exponential generating function is

sum_(n≥0) I_n z^n/n! = exp(z+z²/2).

This is an exact S1/S2 decomposition of converse-fixed permutations, not a numerical analogy between representation theory and primes. Primality is used because every nonzero translation is one p-cycle; the displayed fixed-point counts need modification at composite cell sizes.

|p|ordinary N_p|self-converse F_p|merged M_p|
|---|---:|---:|---:|
|3|2|2|2|
|5|8|6|7|
|7|108|34|71|
|11|329900|3246|166573|

The cases3,5,7 are also exhausted directly. The11 row is a theorem-formula evaluation, not an exhaustive tournament census.

## 4. The scalar incidence repair necessarily stops

The rectangle contraction J(q) sums one sign over the p(p−1)/2 internal arcs of R_p. Thus J has at most p(p−1)/2+1 integer values. Any identification up to converse can reduce N_p types by at most a factor2.

For every odd prime p≥7,

N_p/2 > p(p−1)/2+1.

Indeed (p−1)!>p³ at p=7, and this inequality persists for every larger integer p by multiplying the factorial by p while the cubic ratio is smaller than p. Consequently N_p>p²>p(p−1)+2. Therefore even EXACT J cannot classify the merged types for any prime p≥7. The successful p=5 scalar repair is a genuine finite boundary, not the start of a universal one-statistic rule.

At p=7 the complete test gives108 ordinary and71 merged types, but only16 observed J values and46 pairs of exact (J,K) values, where K is the short-cycle contraction from the five-cell note. Thus even both exact contractions fail at the next prime. The restricted mod8/mod4 bank has only8 labels there. These are exact finite collision counts, not an assertion that every conceivable shared-incidence bank fails.

## 5. Asymptotic size and remaining obligation

Along odd primes, N_p ~ p!/p². Self-converse types form a vanishing fraction of this family. To see the latter without a delicate asymptotic theorem, positivity of the involution generating function gives, for any t>0,

I_p/p! <= exp(t+t²/2)/t^p.

Taking t=sqrt(p) shows p I_p/p! tends to zero. The explicit formulas then give F_p/N_p→0 and M_p~p!/(2p²).

This is a tractable all-size orbit enumeration and a lower bound on how much information a classifying signature must retain in this family. It is not a fast enumeration algorithm for all tournaments. The next reconstruction obligation is to retain enough of the actual permutation's rotation-double-orbit, rather than only its two aggregated orientation correlations.

## 6. A complete all-size shared-incidence carrier, rather than another scalar

The failure of J and K has a constructive repair within this matching-block family. Define the displacement word

d_i=q(i+1)−q(i) in F_p, for i modulo p.

Both internal short cycles are intrinsic. Thus each d_i records the actual cyclic displacement between two matching-owned images, not an independently matched histogram. A codomain translation leaves the word unchanged; a domain translation cyclically rotates it. Conversely, equal words up to cyclic rotation integrate to permutations differing only by those two translations. Therefore the displacement NECKLACE is a complete ordinary isomorphism signature for this family at every odd prime p.

Its realizability condition is exact: the total sum of the p displacements must be zero, and the p proper partial sums

0, d_0, d_0+d_1, ..., d_0+...+d_(p−2)

must be pairwise distinct in F_p. Necessity follows from q being a permutation. For sufficiency take those partial sums as q(0),...,q(p−1); distinctness covers F_p and the zero total closes the final step. No positivity, hidden matching existence, or phase consistency assumption remains in this restricted carrier.

For converse, first reconstruct q and form r q^(-1) r, then take its displacement necklace. One must not substitute a generic word reversal for this inverse-permutation operation. The smaller of these two necklaces is a complete merged signature. This reconstructs the matching-coupled tournament, not its Hamiltonian-path flag.

The carrier has p field entries instead of a p-by-p binary cross matrix. It is an exact compression of this matching family, not of arbitrary regular cross blocks. Its target is actual shared ownership; discarding the cyclic order returns to the scalar collision problem.

## 7. Automorphisms and metagraph fibre weights

The stabilizer in C_p×C_p is nontrivial exactly when q is affine. In that case it has order p; otherwise it is trivial. Thus precisely p−1 ordinary types have automorphism group C_p, and the remaining N_p−(p−1) types are asymmetric. Among the affine types, only slopes+1 and−1 are self-converse; all other slopes pair with their inverses.

Every tournament automorphism acts freely on its directed Hamiltonian paths: fixing one full ordered path fixes every vertex. Hence the number of fixed-path presentations in an ordinary class is H(T)/|Aut(T)|. This is the inherited metagraph fibre rule, now with its automorphism denominator classified for an all-prime-size family.

For the eight p=5 types in the preceding table, ordinary fibre sizes are1653,8245,8221,8201,8193,1635,1635,1649. The slope2/3 pair merges to one fibre of size3270; the other six fibres remain separate. These are weights of selected nodes in the full order10 metagraph, not a claim that these nodes form a closed wiring component or that their degrees determine their internal edge pairing.

## Reproduction

python3 check_prime_cell_growth.py
python3 -O check_prime_cell_growth.py

The checker exhausts all permutations for p=3,5,7; verifies ordinary and converse orbit formulas, shared-incidence collision counts, and evaluates the proved formulas at11,13,17. Separate independent audit enumerates the small quotients. All integer formulas are exact.
