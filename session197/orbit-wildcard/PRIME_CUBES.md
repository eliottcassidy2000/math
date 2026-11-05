# Prime-cube rigidity as a corollary of the strong-order bound

Separate Session197 corollary. **If H(T)=p^3 for a prime p>=5, then T has trivial automorphism group.** This extends the prime and prime-square conclusions stated in the named S196 symmetry packet. The underlying free-action principle, strong-order estimate, and small-orbit method are inherited. No literature-priority claim is made, and no existence theorem for any p^3 value is asserted.

## 1. Inputs and scope

The [S196 symmetry proof](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/symmetry/PROOF.md) proves:

- |Aut(T)| divides H(T), by its free action on Hamilton words.
- An automorphism of prime order p has vertex orbits of size1 orp. Each nontrivial orbit induces a prime-order circulant and is regular and strong.
- Disjoint induced sets give H(T)>=product H(T[S_i]). This also has the direct shuffle proof in the named [S196 orbit-spectrum packet](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/orbit-spectrum/PROOF.md).
- Five-circulants all have H=15,F2=90; seven-circulants have H175 or189,F2=1638. Here F2 counts ordered two nonempty path covers.
- If a strong module M is properly contained in a strong component, a directed triangle through its contracted vertex gives an induced cyclic extension M->a->b->M, with H=3H(M)+F2(M).

We also use Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics13(2006), N3, Theorem1, DOI[10.37236/1141](https://doi.org/10.37236/1141). Its residue-sensitive strong-order bound implies, for every strong n-vertex tournament,

    H >= 5^((n-1)/3).                                  (1)

This is the same externally checked quantitative input used in the Session197 entropy packet. Official paper: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf . The theorem is applied only to a strong induced orbit tournament, not to weak tournaments with arbitrary singleton padding.

The64 connection sets at p=13 are exhaustively certified below; the32 at p=11 are also recomputed rather than only inherited.

## 2. Proof

Suppose Aut(T) is nontrivial and H(T)=p^3. Its order divides p^3, so Cauchy's theorem supplies an automorphism g of order p. Select a nontrivial orbit O.

For every integer p>=17,

    5^(p-1)>p^9, hence 5^((p-1)/3)>p^3.                (2)

For an entirely integer check, the base inequality is5^16>17^9. The inequality5*17^9>18^9 and monotonicity of p/(p+1) imply5p^9>(p+1)^9 for all integers p>=17; induction proves (2). Applying (1) to O contradicts H(O)<=H(T). Thus only p=5,7,11,13 remain.

For p=11, the complete circulant bank has minimum H92411>11^3=1331. For p=13, the complete bank has minimum H3669497>13^3=2197. Induced monotonicity excludes both.

For p=5, every orbit has H15, so two nontrivial orbits would give H(T)>=225>125. Thus there is exactly one nontrivial orbit M. All other vertices are fixed by g and see M uniformly, making M a module. Let S be its strong component. If S=M, the ordered SCC product makes15 divide125, impossible. If S properly contains M, its cyclic extension has H=3*15+90=135>125, also impossible.

For p=7, two nontrivial orbits would give H(T)>=175^2>343, so again M is the unique nontrivial orbit and a module. If S=M, one of175 or189 divides343, impossible. Otherwise the cyclic extension has H>=3*175+1638=2163>343. This contradiction finishes the proof.

The argument allows arbitrary ambient order and weak tournaments. A large number of fixed vertices is handled by the module and SCC argument, not by an invalid global order cutoff.

## 3. Exact finite certificates and consequences

The p=13 histogram over all2^6 skew connection sets is

    3669497:12, 3683797:12, 3703011:4,
    3704857:12, 3707483:12, 3711175:12.

The p=11 histogram is

    92411:10, 93027:10, 93467:10, 95095:2.

Each connection set chooses exactly one of d and -d for1<=d< p/2; this exhausts all possible orbit presentations. Two independent endpoint/startpoint subset recurrences check every Hamilton count. The five- and seven-orbit F2 values are independently checked from subset H values. All integer inequalities and divisibility exclusions used above are explicit verifier gates.

The boundary p=3 cannot be included: an ordered join of three directed triangles has H27 and nontrivial automorphism group. The verifier constructs this graph. This does not make an assertion about every H27 tournament.

Consequently every ordinary class with H=p^3, p>=5, has exactly n! labelled realizations and contributes p^3 to the path-pointed A000568 mass. These are direct consequences of rigidity, not formulas for the number of such classes.

Run `python3 verify_prime_cubes.py` and `python3 -O verify_prime_cubes.py`. The original wildcard and all other frozen supplements are unchanged.
