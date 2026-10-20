# Second moment of the hollow-frame path count

Status: PROVED; independently reviewed. This is an add-only theoretical supplement, separate from the finite maximum census in `../HOLLOW_FRAME_BLUE_EXTREMA.md`.

This is a separate theoretical addendum to the exhaustive blue-maximum census. It does not claim a frame census beyond k=5. For k>=1, D is uniform on the 2^(k(k-1)) labelled external frames (A,b), and h(D) counts rho-symmetric directed Hamiltonian vertex words.

## Formula

Set R(z)=(1+z^2/4)/(1-z^2/4), interpreted as a formal series. Then

    E[h(D)^2]/k! = [z^k] R(z) sum_(m>=0) m! (z R(z))^m.

Equivalently,

    E[h(D)^2]
      = k! sum_(j=0,...,floor(k/2)) (k-2j)!/4^j
             * [t^j] ((1+t)/(1-t))^(k-2j+1).

The mean is E[h]=k!, so subtract (k!)^2 for the variance.

Exact variances for k=1,...,7 are

    0, 1, 6, 75, 1500, 89685/2, 1875825.

The values through k=5 agree with the complete external-frame census. The k=6,7 values are consequences of this formula, also checked through the much smaller permutation/sign-pair sum below; they are NOT frame-census results.

## Proof

A potential symmetric path is determined by a permutation p of the k pairs and a choice s_i in {0,1} for each pair, where s_i=0 means x_i and s_i=1 means y_i. Its final pair must choose x, so s_(p_k)=0. It requires k-1 distinct fair external bits. Write I_(p,s) for its indicator.

For a second candidate (q,t), let d_i=s_i XOR t_i. Constraints of the two paths can overlap only on unordered pair adjacencies shared by p and q. On a shared edge {i,j}:

- If d_i != d_j, the paths constrain different external bits and their probability ratio relative to independent constraints is 1.
- If d_i=d_j=0, that ratio is 2 when the paths traverse the pair adjacency in the same direction, and 0 otherwise.
- If d_i=d_j=1, it is 2 when they traverse it in opposite directions, and 0 otherwise.

This is because equal vertex choices preserve an arc constraint, complementary choices reverse it under rho, and changing exactly one endpoint switches between an A bit and a b bit.

Let epsilon_e=+1 for the same direction and -1 for opposite directions, and r_i=(-1)^d_i. The edge factor is

    w_e = 1 + epsilon_e (r_i+r_j)/2.

The joint probability is 2^(-2(k-1)) times the product of these factors.

Summing over all admissible s,t for fixed p,q makes d uniform on {0,1}^k when p_k != q_k. Indeed there are 2^(k-2) choices of s for each d. When p_k=q_k, d is uniform subject to d_(p_k)=0, since there are 2^(k-1) choices of s for each such d. Thus the signed-candidate contribution for p,q is the expectation of the product of w_e under independent fair r_i, except that r_(p_k)=1 when the final pair agrees.

By relabelling the pairs, fix p=(1,...,k), and multiply the final sum over q by k!. Shared edges form disjoint intervals of the p-path; on each interval q traverses all edges in one direction, so epsilon is constant throughout that interval.

Expand the product of w_e. For a free r_i, only even powers survive expectation. Every edge term chooses at most one endpoint variable, and every vertex has at most two incident shared edges. Consequently nonconstant surviving terms consist exactly of disjoint pairs of consecutive edges whose chosen endpoint is their common vertex. Think of each such two-edge choice as a domino. Its factor is epsilon^2/4=1/4. No two dominoes can share an edge. They may be adjacent and hence merge into a longer selected-edge interval.

In the equal-terminal case there is one additional possibility: the last edge may select the fixed r_k=1, contributing 1/2. Its sign is positive, because q ends at k and must traverse that terminal shared-edge interval in the same direction as p. This is a single terminal monomer, disjoint from all dominoes. The contribution without a monomer is already counted by the free-variable calculation; the monomer gives the entire correction.

Now sum over q by first fixing a domino configuration on p's k-1 edges. Every maximal selected-edge interval has even length 2s and contracts to a block of 2s+1 vertices. It has two possible orientations in q and weight 4^(-s). Every unused vertex is a singleton block. These blocks occur consecutively in p, and after contraction all m blocks can be permuted freely in q, giving m!. This counts permutations containing the selected edges, allowing any additional shared edges, as required by the product expansion.

The ordinary generating function of a singleton or one nontrivial selected-edge block is

    W(z) = z + sum_(s>=1) 2*4^(-s) z^(2s+1)
         = z R(z).

Therefore the free-variable contribution, summed over q, is

    [z^k] sum_(m>=0) m! W(z)^m.

For the terminal monomer correction, the last selected-edge interval has odd length 2s+1, consists of s dominoes followed by the monomer, and contracts to a terminal block of 2s+2 vertices. Its orientation is forced, it must be last in q, and its weight is (1/2)4^(-s). Its generating function is

    V(z) = sum_(s>=0) (1/2)4^(-s) z^(2s+2)
         = (z^2/2)/(1-z^2/4).

The other m blocks are as before and may be freely permuted. Thus the correction is

    [z^k] V(z) sum_(m>=0) m! W(z)^m.

Since 1+V(z)=R(z), adding and multiplying by k! proves the formal-series identity. Extracting even powers from R yields the stated one-sum coefficient formula.

## Independent exact checks

`verify_second_moment.py` computes the overlap factors directly for every q in S_k and every d, without importing the producer scripts or using the domino proof. It compares that result with the coefficient formula for k=1,...,7. It also compares the formula with the moments in the complete frame census for k=1,...,5. All arithmetic is rational/integer.

## Independent review and reproduction

The independent mathematical review checked the joint sign normalization, constant direction sign on shared intervals, disjoint-domino expansion, terminal-monomer correction, and block-permutation factorials. It independently reran the overlap verifier for k=1,...,7 and compared with the complete external-frame moments for k=1,...,5. The proof and checks were approved on 2026-10-20. See `SECOND_MOMENT_REVIEW.md`.

From the published `second_moment` directory, run

    python3 verify_second_moment.py --census-dir ..

To check the formula independently of the finite census files, omit the option:

    python3 verify_second_moment.py

Both modes use only the Python standard library. All validation failures are explicit errors and remain active under `python3 -O`. The output reports whether each formula value was also compared against a complete external-frame census; k=6,7 never receive that label.
