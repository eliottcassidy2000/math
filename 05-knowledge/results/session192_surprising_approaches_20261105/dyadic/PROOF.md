# Fixed-exterior residue programming and simultaneous square-root nets

**PROVED elementary statements, with finite exact controls; independently audited.** Session 192, 2026-11-05. The objects are specified one-apex presentations, not uniformly sampled unrooted tournaments. This does not prove exact integer coverage or the global `{7,21}` conjecture.

## 1. Object and inherited product coordinates

For a binary word w=b_1...b_n, let T(w) have vertices x,1,...,n, with i->j for i<j and x->i iff b_i=1. Write H(w) for its Hamiltonian-path count. **n is chain length; the tournament order is n+1.** T(w) is strong exactly when b_1=1,b_n=0 (n>=2).

**Prior result recovered during cross-lane review:** S26 already proves the binary-product approximation, product-of-uniforms limit, variance, exact finite mean, and its quantitative coupling in [the S26 classification](../../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md), sections 3–4. Its sections 5–7 also prove a positive proportion of globally missing one-apex values. Sections 1 and 3 below rederive that inherited mechanism for a self-contained proof; no novelty is claimed for them. The additions are the explicit selector uniform over arbitrary middle words (section 2), insulated fixed-length residue programming (section 4), and their simultaneous conclusion (section 5).

The direct path-slot decomposition, proved in S25 and restated in S191 arithmetic, is

    H(w)=b_1+1-b_n+N_01(w)
         + sum_(j<i) b_j(1-b_i) 2^(i-j-1),              (1)

where N_01 counts adjacent 01 pairs. Briefly: a path has x first, last, or internal. In the internal case, if its predecessor i is below its successor j, the pair must be adjacent; otherwise each of the i-j-1 intervening chain names chooses independently which side of x it occupies.

Introduce the *two* dyadic coordinates

    X(w)=sum_(j=1)^n b_j 2^(-j),
    Y(w)=sum_(i=1)^n (1-b_i) 2^(-(n+1-i)),
    R(w)=sum_(i<j) (1-b_i)b_j 2^(i-j-1).                (2)

Thus X reads from the left, while Y reads the complemented word from the right. Expanding their product and separating j<i, j=i and j>i gives the exact identity

    H(w)=2^n X(w)Y(w)+b_1+1-b_n+N_01(w)-R(w).          (3)

The diagonal term is zero. This is an exact identity over the dyadic rationals, not an asymptotic fit.

For every 01 pair i<j, at least one adjacent 01 boundary t lies between i and j. For any fixed such boundary,

    sum_(i<=t<j) 2^(i-j-1) < 1.

Counting a pair once for every intervening boundary can only increase its weight. Hence 0<=R(w)<N_01(w) when N_01>0; if N_01=0 then R=0. For strong words, therefore,

    2 <= H(w)-2^n X(w)Y(w) <= B_n,
    B_n = floor(n/2)+1.                                (4)

Indeed N_01<=floor((n-2)/2). The lower equality holds exactly when N_01=0, namely w=1^o0^z with o,z>=1. The error is nonnegative. Keeping this sign avoids an unnecessary absolute-error estimate.

## 2. A uniform explicit selector at every fixed order

Let L>=0, a,b>=1, and n=a+b+L. The middle word of length L may be **arbitrary**, including chosen after the two ends. Given any real target h in [2^(n-2),2^n], put t=sqrt(h/2^n), so 1/2<=t<=1. Define

    A=min(2^a-1, floor(2^a t)),
    B=min(2^b-1, floor(2^b t)).                         (5)

Use the length-a binary word for A as prefix, and the reverse complement of the length-b binary word for B as suffix. Both start digits in their respective readings are 1, so the resulting word is strong. For **every** choice of the L middle bits,

    |X-t| <= 2^(-a),    |Y-t| <= 2^(-b).

Since X,Y,t lie in [0,1], |XY-t^2|<=2^(-a)+2^(-b). Thus the output satisfies

    |H(w)-h| <= E(n;a,b),
    E(n;a,b)=2^(n-a)+2^(n-b)+B_n.                     (6)

No search through tournaments is required. For integer h, the floors in (5) are computable with integer square roots: floor(sqrt(h*2^(2a-n))), clearing a power-of-two denominator before applying isqrt.

Taking L=0 and balancing a,b proves that the strong one-apex counts at chain length n form an E_n-net throughout [2^(n-2),2^n], with

    E_n=2^floor(n/2)+2^ceil(n/2)+floor(n/2)+1.          (7)

In particular every gap between consecutive distinct counts in this fixed-order class has length at most 2E_n: apply (6) at its midpoint. S191's independently established first distinct lower gap is 2^(m-2) at n=2m>=6, and 3*2^(m-2) at n=2m+1>=7. Together these give

    maximum fixed-order gap = Theta(2^(n/2)).          (8)

The constants here are not asserted sharp. This is a mesoscopic size statement, not interval filling: the gaps themselves grow. If one elects not to import S191's exact boundary-gap theorem, (7) and the O(2^(n/2)) upper bound are independent of it.

## 3. The inherited S26 product distribution (reproduced, not new)

Choose w uniformly from all 2^(n-2) strong binary presentations, n>=2. Let a=floor(n/2), b=ceil(n/2). The first a bits of w and the last b bits in reverse-complement reading are disjoint independent blocks. Their leading digits are fixed at 1; the remaining digits are independent fair bits. The corresponding truncations X_a,Y_b are independent and uniform on the dyadic grids

    {1/2+j/2^a : 0<=j<2^(a-1)},
    {1/2+j/2^b : 0<=j<2^(b-1)}.

Add independent uniform tails of widths 2^-a and 2^-b to obtain independent U,V, each Uniform[1/2,1]. Since X and U lie in the same length-2^-a cell, and Y and V in the same length-2^-b cell, (4) gives the coupling bound

    |H(w)/2^n - UV| <= delta_n,
    delta_n=2^-a+2^-b+B_n/2^n.                        (9)

Consequently H(w)/2^n converges in distribution, and in every fixed positive integer moment, to Z=UV. Its density is

    f_Z(z)=4 log(4z)  for 1/4<z<1/2,
           -4 log z  for 1/2<=z<1,
           0 otherwise.                             (10)

This follows by integrating the density 4 over the square [1/2,1]^2 along xy=z: f(z)=4 integral dx/x over max(1/2,z)<=x<=min(1,2z). The density is bounded by 4 log 2, so the Kolmogorov distance is at most (4 log 2)*delta_n. In particular

    lim E[(H/2^n)^r] = [2(1-2^(-r-1))/(r+1)]^2,
    E[Z]=9/16,       Var(Z)=55/2304.                  (11)

For n>=3 the finite mean is exactly

    E[H(w)] = 9*2^(n-4)+1/2.                         (12)

To verify (12) directly in (1): the endpoint-endpoint long pair contributes 2^(n-2); the two one-endpoint sums contribute 2*(2^(n-3)-1/2); the wholly internal pairs contribute 2^(n-4)-(n-1)/4. The boundary term has expectation 2+(n-3)/4. Their sum is (12). The n=2 exception is H=3.

This is a distribution of **presentations with multiplicity**. It supplies no limiting density of the distinct integer support, and says nothing about uniformly sampled unlabelled tournament classes. Smooth limiting density is compatible with growing arithmetic gaps.

## 4. Exact local residues and insulated binary switches

Modulo 2^k, (1) retains only pairs with 1<=i-j<=k; farther pairs vanish. The count is therefore a finite-range interaction of the attachment bits, together with the endpoint and adjacent-01 terms.

For arbitrary words u,v, let A be the integer represented by u and p its last bit (p=0 if u is empty). Put

    C(v)=sum_(j:v_j=0) 2^(j-1),
    epsilon(v)=[v is nonempty and v_1=0].

The exact recurrence H(u0)=H(u)+A+p, H(u1)=H(u), with the integer state updating A->2A+b, proves the adjacent-switch law

    H(u10v)-H(u01v)=A+2-p+C(v)-epsilon(v).             (13)

Fix k>=2. For j=1,...,k-1, set A_j=2^j-2 and let P_j be its binary expansion padded on the left to exactly k bits. Define two length-(2k+2) gadgets

    G_j(0)=P_j 01 1^k,
    G_j(1)=P_j 10 1^k.                               (14)

In every exterior, the switch has preceding k bits P_j and following k bits all 1. Hence A=A_j mod 2^k, p=0, C(v)=0 mod 2^k and epsilon(v)=0. Therefore changing G_j(0) to G_j(1) changes H by exactly 2^j modulo 2^k, **independently of all other bits**.

Let

    M(e)=G_1(e_1)...G_(k-1)(e_(k-1)),
    L_k=(k-1)(2k+2)=2k^2-2.                          (15)

For any fixed prefix u and suffix v, all 2^(k-1) choices e obey

    H(u M(e) v)=H(u M(0) v)+sum_(j=1)^(k-1) e_j 2^j
                     (mod 2^k).                     (16)

All H are odd (also immediate from the recurrence), so (16) is a bijection onto the odd residues. This is a **fixed-length, fixed-exterior** residue theorem, stronger in these respects than the inherited existence of some arbitrarily large modular witness. The gadgets may change exact H enormously; only the residue increment is local.

For k=1 take L_1=0 and the empty middle bank. The only odd residue modulo 2 is automatic.

## 5. Simultaneous residue and magnitude, with explicit quantifiers

For every k>=1, every n>=L_k+2, every odd residue r modulo 2^k, and every real h in [2^(n-2),2^n], there is a strong one-apex word w of length n such that

    H(w)=r mod 2^k,
    |H(w)-h| <= E_(n,k),

    a=floor((n-L_k)/2), b=ceil((n-L_k)/2),
    E_(n,k)=2^(n-a)+2^(n-b)+floor(n/2)+1.              (17)

Choose the prefix and suffix by (5), put the gadget bank between them, then select e by the binary expansion of r-H(uM(0)v) modulo 2^k. Equations (6) and (16) apply simultaneously to this same word. For each fixed k, E_(n,k)=O_k(2^(n/2)); the hidden constant grows roughly as 2^(k^2).

This does not silently interchange quantifiers: k is fixed for the square-root asymptotic. It supplies a residue-specific square-root net, not exact coverage. An exact target h would follow from H=h mod 2^k and |H-h|<2^k. The bound in (17) never supplies this inequality: for k>=2, L_k>=k and a,b>=1 force E_(n,k)>=2^(L_k+1)>2^k; k=1 is also insufficient. A sharper coupled correction mechanism, not separate congruence and size claims, is needed to cross that boundary.

## 5.1 The fixed-modulus nets coexist with positive missing density

S26, `APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md` sections 5–7, proves an independent, stronger hostile to any exact-coverage interpretation. If A(X) counts distinct values H<=X attained by any one-apex tournament, including all source/sink padding, then

    A(X) <= (log 2) X/2 + O(sqrt X).

Consequently at least

    (1-log 2) X/2 - O(sqrt X)

positive odd integers through X are absent from this class. Equivalently, the lower proportion missing among odd integers is at least 1-log 2, approximately 0.30685. This inherited theorem counts **distinct values**, unlike the presentation distribution in section 3.

There is no conflict with (17): for each fixed k, every square-root-scale window of the stated sufficient width inside the fixed-order size range contains a value in each specified odd residue, yet a positive proportion of individual odd integers remain missing globally. Even the fixed-order maximum gaps grow at square-root scale by (8). Thus “microscopic obstruction” means a loss below the resolution of a fixed-modulus/mesoscopic observer; it does not mean all holes have bounded length. The modulus-dependent constants and the order of quantifiers matter. Neither the simultaneous nets nor a smooth presentation limit implies asymptotic exact coverage, density one of attained integers, or permission to let k grow without paying its explicit L_k and E_(n,k) costs.

## 6. Cross-thread source and limitations

The source inspiration is the distinction between retained old digits and fresh future digits in THM-2555 / `natural-extension-sheet-charge-and-future-digit-boundary`, together with the exact rather than heuristic digit laws in THM-4473 / `collatz-digit-chains-rotation-repunits`. Their LRC and Collatz conclusions are not used as tournament theorems. The actual map here is the explicit two-ended dyadic map (2); the preserved quantity is H up to the exact correction (3). The word, and thus H, can be recovered from X at finite n, but the continuum limit forgets its microscopic arithmetic. The product law was already proved in S26 and is reused here. The fixed-exterior residue switches and simultaneously controlled selector are the additions; their point is to exploit that earlier product mechanism rather than rebrand its limit theorem.

The direct slot identity and swap identity are inherited from S25 and S191 and are reproved above. The sharp lower boundary gap in (8) is explicitly inherited. No literature-priority claim is made. The S26 all-order exclusion of H=13 from this class, also recovered in S191, remains a decisive hostile against exact-coverage claims: even this simultaneous local/magnitude theorem does not remove it.

Run `python3 verify.py` and `python3 -O verify.py`. The verifier checks the product identity with integer arithmetic, independently counts paths by vertex-subset DP and permutations on small universes, tests the deterministic selector, exhausts gadget residue banks in stated finite ranges, and records exact presentation moments and support gaps. Numerical CDF diagnostics are only diagnostics; (9)-(11) have the elementary proof above.
