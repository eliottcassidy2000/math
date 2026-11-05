# Exact cyclic-context content and a physical three-level divisibility hierarchy

Session197. **PROVED for every nonempty physical tournament.** Fix V and define

    C(V)=gcd_(nonempty tournaments T,U) H(C3[T,U,V]).

Let g(V) be ordinary bin-value content and G(V) the universal intact-module content. If B_V=sum_k u_k q_k for odd order, or sum_k u_k r_k for even order, with the central-factorial basis and gamma_k below, then

    g(V)=gcd_k gamma_k u_k,
    G(V)=gcd_k k! gamma_k u_k,
    C(V)=gcd_k (k!)^2 gamma_k u_k.                       (1)

The first two identities are the separately frozen UNIVERSAL_MODULE_CONTENT theorem. This proof establishes the third from the actual cyclic Gram kernel. It gives

    g(V) | G(V) | C(V) | H(V),
    g(V) C(V) | G(V)^2, G(V)^2 | g(V)^3,
    C(V) | g(V)^2, C(V)^3 | G(V)^4.                     (2)

An actual ten-vertex V has g=9, G=27 and C=81, so the cyclic restriction can strictly overestimate the universal module divisor. It also attains equality in all the displayed exponent bounds. Gcds are not claims that their values are attained by the construction.

## 1. Physical kernel and a finite test bank

Put n=|V| and f_r=r!c_r(V), with f_0=f_(n+1)=f_(n+2)=0. The reviewed CYCLIC_GRAM_FORESTS theorem gives the exact matrix K on unordered covers of T,U, indexed1 through n+1:

    K_(r,r)=(r!)^2 (f_(r-1)+3f_r+f_(r+1)),
    K_(r,r+1)=K_(r+1,r)=r!(r+1)!(f_r+f_(r+1)),

with other entries zero. In particular C(V) is the gcd of the entries of K. One direction follows because every actual count is an integer bilinear form in K. For the reverse direction, use actual inputs Tr1,...,Tr_(n+1). Their unordered-cover matrix is the lower unit-triangular Stirling matrix S(i,r). The physical count matrix is S K S^T, and S and its inverse are integer matrices. Thus the entries of these two matrices have the same gcd. No virtual input tournament is used in this finite test bank.

Every actual Hamilton count is odd, by the tournament odd-cycle collection identity. Therefore C(V) is odd. This parity is essential in the next arithmetic step.

## 2. Central coordinates and the determinant-four boundary

Write d=ceil(n/2), and

    q_k=x product_(j=1)^k(x^2-j^2), gamma_k=(2k+1)!,
    r_k=x^2 product_(j=1)^k(x^2-j^2), gamma_k=(k+1)(2k+1)!.

All u_k are integers. The normalized binomial identities from UNIVERSAL_MODULE_CONTENT give

    f_r=sum_k N_(r,k) gamma_k u_k,

where every N_(r,k) is an integer, is zero for r<k+1, and equals1 at r=k+1. Set a_k=(k!)^2 gamma_k u_k.

Every kernel entry is an integer combination of the a_k. In diagonal row r only k<=r can occur; its coefficient on a_k is ((r!)^2/(k!)^2)(N_(r-1,k)+3N_(r,k)+N_(r+1,k)). The analogous off-diagonal coefficient is r!(r+1)!/(k!)^2 times N_(r,k)+N_(r+1,k). The factorial ratios are integers whenever the relevant coefficient is nonzero. For r=k>=2 the coefficient of a_k in K_(k,k) is exactly1.

The first two entries, including their exact small-index normalization, are

    (K_11,K_12)=(3a_0+a_1, 2a_0+2a_1) for odd n,
    (K_11,K_12)=(5a_0+a_1, 6a_0+2a_1) for even n.       (3)

When d=1 put a_1=0. Both displayed2-by-2 coefficient matrices have determinant4. Hence C divides4a_0 and4a_1. Induction on k using the diagonal K_(k,k) gives C|4a_k for every k. Conversely gcd(a_k) divides every kernel entry. Since C is odd, these two divisibilities prove

    C=gcd_k a_k.

Also a_0=u_0=B_V(1)=H(V), because every higher central basis polynomial vanishes at1. This proves C|H and the chain in (2).

Only the entries with indices at most L=max(2,d-1) were used to recover all a_k. The same Stirling change of coordinates shows that the physical square bank T,U in {Tr1,...,Tr_L} already has gcd C. Its graphs have at most n+2L vertices. This is a sufficient bank, not a minimal-size assertion.

The odd physical hypothesis cannot simply be discarded. The formal positive odd polynomial P=x^3+x has a_0=2,a_1=6, whose gcd is2, but the corresponding kernel has content4: K_11=12, K_12=16 and all its other entries are divisible by4. It is not the bin polynomial of a tournament, because P(1)=2 is even. The determinant-four step explains the exact obstruction to blindly extending (1) to arbitrary formal profiles.

## 3. Sharp valuation bounds and target scope

Fix a prime p and let a=v_p(g), b=v_p(G), c=v_p(C). Choose a coordinate attaining the minimum for g. Since (k!)^2 divides gamma_k and u_k is an integer, v_p(k!)<=floor(a/2). The cyclic formula therefore gives

    a <= b <= c <= a+2 floor(a/2) <= 2a.

Choosing instead a coordinate attaining b, its factorial valuation is at most floor(b/3), because v_p(gamma_k u_k)+v_p(k!)>=3v_p(k!). Hence

    c <= b+floor(b/3).

There is also a stronger mixed divisibility. For every k, write A_k=v_p(gamma_k u_k), K_k=v_p(k!), omitting zero coordinates. We have a<=A_k and c<=A_k+2K_k. Choosing k attaining b=A_k+K_k gives a+c<=2b. Thus g C divides G^2. This is concavity of the minimum of the affine functions A_k+tK_k at t=0,1,2, with no probabilistic or analytic approximation. Together with the inherited G^2|g^3, it also implies the last two divisibilities in (2). The direct bounds above independently verify them. No physical interpretation is asserted for further formal factorial powers. In particular a=0 forces c=0 and a=1 forces c=1. Thus g,G,C have the same prime support, and divisibility by p^2 is equivalent for all three.

This sharper content is not a new prime/7p obstruction beyond ordinary content and induced-count monotonicity. If N is prime or7p and g|N, H(V)<=N, then C|N. For squarefree N, use the preceding valuation-one rigidity. The only remaining case N=49 uses g in{1,7,49}; if g=49, then49|C|H<=49 forces C=49. This limitation does not negate the exact higher-prime-power hierarchy.

## 4. Literal physical separation:9,27,81

Take

    V=[382,420,88,706,810,972,18,597,204,327].

Its unordered path-cover counts c_1,...,c_10 are

    7857,70551,185274,211989,124725,40851,7672,814,45,1.

Its bin polynomial and even central coordinates are

    B_V=x^10+64x^8+1000x^6+4020x^4+2772x^2,
    (u_0,...,u_4)=(7857,10449,2043,94,1).

The three gcd formulas give g=9,G=27,C=81. These are independently checked from the full ordered-cover vector, its factorial weighting, and every entry of the literal cyclic kernel.

Two small actual cyclic graphs already certify that C cannot exceed81:

    H(C3[Tr1,Tr1,V])=164673,
    H(C3[Tr1,Tr2,V])=462591,
    gcd(164673,462591)=81.

Both counts are computed directly on their12- and13-vertex graphs. The kernel formula proves divisibility by81 for every other cyclic pair.

For contrast, substitute V into the four-run matching context Q_4 of UNIVERSAL_MODULE_CONTENT. Its exterior run coefficients are(91,101,44,6), giving the literal16-vertex graph count94405041, which is27 modulo81. This graph is a whole-module context but is not of the two-block cyclic form used to define C. Together with H(V)=7857 it has gcd27, certifying that the universal divisor really is smaller. The verifier directly counts this16-vertex graph as well.

## 5. Reproduction and named dependencies

Run `python3 verify_cyclic_module_content.py` and its `-O` counterpart. It checks both kernel-content formulas and finite physical banks on all labelled tournaments through order5, exact central transformations on signed formal profiles with odd evaluation at1, the even-evaluation hostile, the sharp valuation bounds and all three literal larger witness graphs. The source uses exact integers throughout.

Named dependencies, all in the same published path-polynomial directory:

- CYCLIC_GRAM_FORESTS.md SHA2562fe61a1fcbd4133fbd9bddcdd8a337dfe521d10ca41514f69ecb25f358327c83.
- verify_cyclic_gram_forests.py SHA256f56c9feed4d93b757a35bf7e847fc03a38d623d25214bae59fc9fefa5c34f6f0.
- UNIVERSAL_MODULE_CONTENT.md SHA2564f86f7ada5eff79e61bb2ef4b8139b76df8efedf6996a71dae327bc66b24d9cb.
- verify_universal_module_content.py SHA256c9779056d1a4269d8fa6146a4dc37735c69f6f21b8e379b884794193921e5b78.
- verify.py SHA25610c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb.

The universal-content dependency is separately under independent review. No claim from an exploratory search or unproved profile-realizability assertion is used.
