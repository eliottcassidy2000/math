# Exact maps to binary automata, tournament parity, and arithmetic lifts

2026-10-19. Status: PROVED elementary statements; exact finite replays described below. These are explicit transfers of the supplied paper's product-character identity, not transfers of its analytic assumptions or its headline.

## 1. Biased binary strips are Pascal-code weight enumerators

Inheritance: 05-knowledge/results/gilbreath_binary_entry_20261019/BINARY_ENTRY_AND_EXACT_MIXING.md proves the binary system

    x(t+1)=A x(t)+e_F u(t),     A=I+S over F_2,

where S shifts right neighbors left and kills the last coordinate, and u(t) are independent fair bits. It proves exact uniformity after F steps. Here retain a general input bias q, with independent Bernoulli(q) inputs. This does not model a biased initial infinite sea: its boundary trace need not be independent.

Put a=1-2q. For t steps define columns b_j=A^j e_F, 0<=j<t, and for v in F_2^F put

    w_t(v)=#{j<t: v dot b_j=1}.

The exact Fourier transform is

    E (-1)^(v dot x(t))
       =(-1)^(v dot A^t x(0)) a^(w_t(v)).               (A1)

This is the paper's Bernoulli subset-product identity on the additive group F_2^F, with character z=(-1)^(v dot b_j). Its energy is exactly 2 w_t(v). No limiting approximation is involved.

Let P be the least power of two >=F. Since S^P=0 and (I+S)^P=I+S^P, A^P=I. Define the length-P binary linear code

    C_F^P={ (v dot A^j e_F)_(0<=j<P): v in F_2^F },
    W_F^P(z)=sum_(v in F_2^F) z^(w_P(v)).

The first F columns form the invertible Pascal matrix, so this code has dimension F. At t=mP,

    P(x(t)=x(0)) = 2^(-F) W_F^P(a^m),                 (A2)
    ||Law(x(t))-Uniform||_TV
       <= (1/2) sqrt(W_F^P(a^(2m))-1).                 (A3)

For (A3), Fourier Parseval gives the squared L2 density distance W_F^P(a^(2m))-1 and Cauchy–Schwarz bounds total variation. The same bound applies to every deterministic initial state; its affine translation changes Fourier signs only.

### Exact recursion and the non-dyadic distance gain

Include F=0 with W_0^P=1. If P=2L then

    W_F^P(z)=W_F^L(z^2),                      0<=F<=L;          (A4)

if F=L+r with 1<=r<=L, then

    W_F^P(z)=sum_(v in C_r^L)
                   (2z)^(wt(v))(1+z^2)^(L-wt(v)).             (A5)

Proof. The Pascal evaluation vector of row index h has jth entry binom(j,h) mod2. For h<L, its two half-blocks are equal. For h=L+h', the first block is zero and the second is the h'th row of length L. Thus for F<=L the words are (u,u); for F=L+r they are precisely (u,u+v), where u is arbitrary in F_2^L and v belongs to C_r^L. A coordinate with v_i=0 contributes 1+z^2 after summing u_i, and v_i=1 contributes 2z. QED.

For P/2<F<P, the minimum nonzero weight is exactly 2. Every Pascal row index h<P-1 has even weight, while h=P/2-1 gives a row of weight 2 and is available because F>P/2. For F=P the code is the full binary space and its minimum weight is 1. Consequently, at complete periods,

    |nontrivial Fourier coefficient| <= |1-2q|^(2m), F<P;
    |nontrivial Fourier coefficient| <= |1-2q|^m,    F=P.       (A6)

This is a two-versus-one exponent in the number of complete periods, not necessarily in physical time: P changes with F. At q=1/2 the inherited exact F-step mixing is stronger. At q different from 1/2 and 0<q<1, no finite t gives an exactly uniform law, because every nontrivial Fourier factor in (A1) is nonzero.

Examples computed by (A4),(A5):

    F=3,P=4: W=1+6z^2+z^4.
    F=5,P=8: W=1+4z^2+22z^4+4z^6+z^8.
    F=7,P=8: W=((1+z)^8+(1-z)^8)/2.

### The complementary-dimension dual is another Pascal code

The coordinating research lane identified this explicit duality. Allow any 0<=F<=P, with P=2^m, and retain the first F Pascal rows. Identifying j with its m binary digits, the row indexed by h is the Boolean monomial

    j -> product_(r: h_r=1) j_r,

by the binomial identity (1+X)^j=product_(r:j_r=1)(1+X^(2^r)) over F_2. The inner product of rows h and k is one exactly when h OR k=P-1: the sum of their product over the Boolean cube is 2^(m-popcount(h OR k)) modulo two.

If h<F and k<P-F then h+k<=P-2, hence h OR k cannot equal P-1. Dimensions now give

    (C_F^P)^perp = C_(P-F)^P.                          (A9)

For the homogeneous enumerator W_C(u,v)=sum_(c in C)u^(P-wt(c))v^wt(c), character orthogonality yields

    W_F^P(u,v)=2^(F-P) W_(P-F)^P(u+v,u-v).             (A10)

Indeed expand each product over coordinates; summing signs over the dual code kills every word outside C and gives |C^perp| for a word in C. This proves the identity without an unverified coding-theory import. It is the classical MacWilliams transform in an explicitly identified control code, not a claim of historical novelty.

The same base-two digits are the Hasse-digit carriers in the repository's tournament-deletion packet, but here they index Boolean monomial evaluation vectors. The common map is the exact binomial/Frobenius digit expansion; the dynamical state and the deletion operator remain distinct objects. This gives a native dual code for the squarefree/cycle lane's character-code transform.

### Fixed-weight forcing exposes the exact missing coset

If a uniform input word has exactly k ones, the state increment is a uniform k-subset sum of the columns b_j. Its Fourier transform is

    K_k(w_t(v);t)/binom(t,k),                       (A7)
    K_k(w;t)=sum_j (-1)^j binom(w,j)binom(t-w,k-j).

The last coordinate of every b_j is one, so

    x_F(t)=x_F(0)+k mod2.                           (A8)

Thus fixed-weight forcing cannot mix over all 2^F states. For t>=F the affine span of the columns is exactly the last-coordinate-one hyperplane, because the first F columns are affinely independent: their pair differences have rank F-1. This is the affine sidecar in FIXED_CARDINALITY_TRANSFER.md in a native dynamical carrier. It identifies the only *affine-linear* obstruction, not an assertion that every compatible small k reaches all states.

Source -> target: input bits -> terminal binary state by the displayed control matrix.
Preserved: every terminal coordinate and its probability.
Lost: the input word, nonlinear pre-entry trajectory, and any prime-gap interpretation.
Sidecar: full input weight/parity when conditioning, plus bounded integer states before entry.
Hostile: F=2,k=1 gives only the last-coordinate-one states.

## 2. Tournament score parity is the cut-character transform

Take a labelled tournament T_0 on v vertices. Independently reverse every edge of its underlying complete graph with probability q. Let epsilon_e be the reversal bits and d(T) the vector of outdegrees modulo two. Each reversal changes the two incident outdegree parities, so

    d(T)=d(T_0)+partial epsilon,
    partial(e_(i,j))=e_i+e_j in F_2^v.              (T1)

For v>=2 the image is the even-weight subspace E of dimension v-1. Hence the possible parity vectors form the affine hyperplane

    sum_i d_i = binom(v,2) mod2.                    (T2)

The characters of E are indexed by vertex subsets S modulo S~S^c. For any representative S,

    E (-1)^(sum_(i in S)(d_i(T)-d_i(T_0)))
       =(1-2q)^(|S|(v-|S|)).                       (T3)

Indeed an edge column pairs nontrivially with S exactly if it crosses the cut. Character inversion gives a completely explicit probability for each score-parity vector:

    P(d(T)=d(T_0)+z)
       =2^(-(v-1)) sum_[S] (-1)^(sum_(i in S)z_i)
                            (1-2q)^(|S|(v-|S|)), z in E.     (T4)

At q=1/2 every feasible parity vector occurs in exactly

    2^(binom(v,2)-v+1)

labelled tournaments. This counts orientations, not tournament isomorphism classes. It is not an assertion about the full score sequence or Hamiltonian-path count.

If exactly k edges are reversed, replace the factor in (T3) by

    K_k(|S|(v-|S|);binom(v,2))/binom(binom(v,2),k).           (T5)

This is another literal instance of (F1). The old modulus energy is 2|S|(v-|S|), while the fixed-size dispersion is

    binom(v,2)-|binom(v,2)-2|S|(v-|S|)|.                    (T6)

The complement inside the absolute value matters: a cut containing every allowed edge creates an exact cardinality obstruction.

### General edge-reversal host: bipartiteness is the affine obstruction

Allow reversals only on a connected simple graph Gamma with v>=2 vertices and at least one edge. The ordinary syndrome span is still E. Its affine-difference span has dimension

    v-2 if Gamma is bipartite,
    v-1 if Gamma is not bipartite.                         (T7)

Proof. A linear functional c on the ambient F_2^v is constant on all edge columns exactly when c_i+c_j is the same value ell on every edge. If ell=0, connectedness forces c constant, which is the trivial functional on E. If ell=1, c is precisely a bipartition coloring; such a coloring exists iff Gamma is bipartite, and is unique modulo the constant vector. Thus there is exactly one extra affine equation in the bipartite case and none otherwise. Duality gives the stated dimensions. QED.

For a bipartition V=L union R, exact-k reversal obeys

    sum_(i in L)(d_i(T)-d_i(T_0))=k mod2.                  (T8)

The triangle is the smallest graph that destroys this additional parity equation. This is a sharp structural control, not an asymptotic mixing theorem.

Connection to existing tournament-statistic work: the score-parity indicator is an honest integer-valued function, so the integral source-free basis in 05-knowledge/results/tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md expands it exactly. But the edge-reversal operation here is not that paper's vertex-deletion operator. The common target is a statistic space; neither operation is substituted for the other. No synthetic tournament on tied phase values is introduced.

## 3. Sparse binary correction digits for arithmetic lifts

Let M be an integer b-by-a matrix, p prime, and vbar in ker(M mod p). Choose an integer lift vtilde and put

    d=M vtilde/p mod p.

An unrestricted first lift exists exactly when

    -d belongs to im(M mod p),                            (L1)

as recorded in the repo's first-obstruction theorem. Now impose the genuinely additional restriction that correction digits w_i are either 0 or 1:

    M(vtilde+p w)=0 mod p^2,
    w in {0,1}^a,                                        (L2)

which is exactly the labelled subset-sum problem

    sum_i w_i c_i=-d in G=im(M mod p),                    (L3)

where c_i is column i of M mod p. If (L1) fails, no subset-product argument can repair it.

If each w_i is independent Bernoulli(rho), the exact transform is

    E chi(sum_i w_i c_i)=product_i(1-rho+rho chi(c_i)).    (L4)

If a correction with exactly k nonzero binary digits is required, (F1)–(F4) apply directly to all a columns when their dispersion condition holds. In that case every target in G has at least binom(a,k)/(2|G|) binary corrections. Alternatively, affine pruning may retain only n<=a coordinates in a coset c+H. Set every removed correction digit to zero. For targets in kc+H, the corresponding mixing condition guarantees at least binom(n,k)/(2|H|) corrections. One must not put the original a back into this pruned binomial coefficient. This is a multiplicity theorem for admissible corrections after the first-lift obstruction has vanished.

An independent review supplied a sharp scope hostile: nine zero columns and one unit column over F_2, with k=8 and target zero. There are exactly binom(9,8)=9 corrections, fewer than binom(10,8)/(2|F_2|)=45/4. Pruning to the nine zero columns makes H trivial and correctly retains the count 9. It does not license the larger unpruned entropy bound.

At p=2, unrestricted digits are already binary, and fair independent correction digits give exactly 2^(a-rank M) solutions for every reachable target. At odd p, binary restriction is substantive; even a reachable unrestricted target need not be a binary subset sum.

Minimal hostile: p=3, M=[1], vtilde=3. Then vbar=0 is in the mod-3 kernel and d=1. The required correction is w=2 mod3, so the unrestricted lift vtilde+3w=9 exists, while neither binary digit 0 nor 1 works. The obstruction class vanishes in the cokernel, yet the restricted digit supply is empty.

Map: selected correction-coordinate labels -> obstruction-killing syndrome -> actual mod-p^2 lift.
Preserved: the full matrix congruence, target, and binary digit support.
Lost if passing to cokernel alone: the correction witness and its digit constraints.
Sidecar: actual columns and the chosen integral defect Mvtilde/p.
Cheapest test: direct enumeration of binary digit vectors versus modular syndrome dynamic programming.

## 4. Precise next questions, not conclusions

1. For binary strips, classify the weight enumerator at arbitrary t=mP+r directly from the truncated Pascal code, and optimize fixed-k mixing on the parity-compatible affine hyperplane. Equations (A1),(A7) already give finite exact algorithms.
2. For sparse edge-reversal hosts, use the graph's cut-size spectrum to give quantitative fixed-k score-parity mixing. Bipartiteness is only the affine obstruction; narrow cuts can cause slow mixing.
3. For arithmetic lifts, determine a natural family of correction-column reservoirs with large phase dispersion and a small affine quotient. The general first-lift criterion alone supplies no such reservoir.
4. For LRC, Fourier kernels can be transported only after identifying a positive ambient phase measure and preserving the target event. No such ambient estimate for the physical owner-labelled frontier is established here; the use of Fejer polynomials alone supplies no LRC(14) result.

These questions reuse exact objects already present in the repository. None asserts randomness of prime gaps, an unrestricted Collatz carry estimate, a Wilson/Wieferich finiteness result, or a new major-conjecture solution.

## 5. Identity-positive reservoirs bypass energy estimates

This additional consequence was independently identified by the squarefree/cycle research lane. Let H be a finite elementary abelian 2-group, let g_1,...,g_n generate H, and select labels independently with probabilities rho_i<=1/2. Every Fourier coefficient

    product_i(1-rho_i+rho_i chi(g_i))

is a product of factors 1 or 1-2rho_i, hence is nonnegative. Therefore

    P(product=1)>=1/|H|,                            (P1)

and the identity is a most likely target. Character inversion proves both statements: every other target inserts signs, while the identity inserts only +1. No pruning or lower energy bound is needed. The same conclusion holds if H is taken as the generated subgroup rather than an unnecessarily larger ambient group.

For equal rho=t/n<=1/2, a Chernoff size-window tail smaller than 1/(2|H|) leaves at least 1/(2|H|) probability in identity subsets of typical size. Thus t much larger than log|H| suffices for this size-window conclusion; an affine obstruction still prevents replacing it blindly by a chosen exact size.

More generally, in any finite abelian group, suppose labels are partitioned into inverse pairs (g,g^(-1)), and the two selection probabilities in each pair agree. Each paired Fourier factor is

    |1-rho+rho chi(g)|^2>=0,

so (P1) and modal identity hold for every rho in [0,1]. This is positivity of an autocorrelation measure. It is a useful independent sufficient structure when an inverse-paired reservoir is actually present.

Minimal failure beyond these hypotheses: in C_3 with two labels both equal to a generator and rho=1/2, identity occurs only for the empty subset, with probability 1/4<1/3. The formal product transform remains valid; positivity has failed. Selected prime residues in the Carmichael paper are not automatically inverse-paired, so this does not replace its high-order Fejer argument.
