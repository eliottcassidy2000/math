# Integral cancellation and typical norm behavior along the Cayley--Dickson tower

**PROVED, scoped; same-agent proof review.** These are explicit statements about the standard integral basis and the stated multiplication. They are not a classification of arbitrary zero divisors, a literature-priority claim, or a solution of LRC(14).

## C1. Convention and the integral two-support theorem

Let A_d have dimension N=2^d, basis e_i for 0<=i<N, and Baez's doubling convention

    (a,b)(c,d)=(ac-d conjugate(b), conjugate(a)d+cb),
    conjugate(a,b)=(conjugate(a),-b).

All structure constants are integral, and

    e_i e_j=S_d(i,j)e_(i XOR j),   S_d(i,j) in {+1,-1}.

This is the convention already used in S13 of the repo's [Paley/octonion construction](../collatz_mod6_20260922_paley_fano_octonion_design.md). Write L_i for left multiplication by e_i. It is a signed permutation matrix, hence integral orthogonal. Fix DISTINCT i,j and integers a,b, not both zero. Put x=a e_i+b e_j, t=i XOR j, and

    P=L_i^T L_j.

The unsigned permutation of P is translation by t. Its cycles are exactly the N/2 pairs {v,v XOR t}. In each such pair P has matrix [[0,s],[u,0]], with s,u signs. Thus multiplication by x, after an explicitly known signed permutation of output and a common pairing of input and output, has blocks

    B_epsilon = [[a,b s],[b u,a]],  epsilon=s u in {+1,-1}.       (1)

Let k_d(i,j) be the number of epsilon=+1 blocks. Then its integral cokernel is

    [Z/g + Z/((a^2-b^2)/g)]^k
      + [Z/g + Z/((a^2+b^2)/g)]^(N/2-k),                       (2)

where g=gcd(a,b)>0, Z/0 means Z, and cyclic orders are taken in absolute value. The plus signs in (2) mean direct sums. This is a cyclic decomposition; to obtain the globally ordered Smith invariant factors, sort valuations separately at each prime (or use gcd/lcm operations). Simply sorting the displayed integers is not generally Smith form.

**Proof.** On the pair indexed by v,w=v XOR t,

    s=S_d(j,w)S_d(i,v),  u=S_d(j,v)S_d(i,w).

No other entries occur. The first Smith factor of a nonzero 2-by-2 block is the gcd of its entries, namely g, and the second is its determinant divided by g. All coordinate changes used before these Smith operations are integral unimodular. This proves (1)--(2), including zero coefficients and singular blocks.

For every modulus m>=1 the exact number of solutions to xy=0 is consequently

    [gcd(m,g) gcd(m,(a^2-b^2)/g)]^k
       [gcd(m,g) gcd(m,(a^2+b^2)/g)]^(N/2-k).                 (3)

This includes composite moduli and prime powers. To solve xy=z for a NAMED z, retain the actual paired coordinates and signs, or full Smith transformations. The integer k and cyclic factors alone do not transport z or individual kernel vectors.

Over R, the squared singular values are a^2+b^2 twice per negative block, and (a+b)^2,(a-b)^2 per positive block. For ab!=0, multiplication by x scales every squared norm by a^2+b^2 iff k=0. It has a real kernel iff k>0 and a=+/-b; that kernel has dimension k. These spectral facts are consistent with classical Cayley--Dickson eigentheory; our arithmetic statements use the integral blocks rather than inferring Smith forms from real eigenvalues.

## C2. An O(d) exact recursion for every named pair

Set k_d(0,j)=k_d(i,0)=0. For i!=j and h=2^(d-1), put u=i mod h,v=j mod h. Then

    same half: k_d(i,j)=2 k_(d-1)(u,v);

    opposite halves, uv=0 or u=v: k_d(i,j)=0;

    opposite halves, u,v nonzero distinct:
        k_d(i,j)=h-4-2 k_(d-1)(u,v).                       (4)

The computation follows one shorter pair at each step, using O(d) arithmetic operations; this is not a bit-complexity claim with unit-cost huge integers. It supplies (2)--(3) at arbitrary dimension without building N-by-N matrices. A compressed answer records factors with multiplicities.

**Proof.** Let C be conjugation, and let c_u=1 for u=0 and -1 otherwise. At the next doubling,

    L_u = diag(L_u, c_u L_u),
    L_(h+u) = [[0,-c_u R_u],[R_u,0]],

where R_u is right multiplication in the previous algebra. Conjugation gives R_u=c_u C L_u C, so left and right pairs have the same cycle-square signs. This proves the same-half case, including upper/upper pairs.

For a lower/upper pair (u,h+v), put Q=L_u^T R_v. The two diagonal halves of P^2 are -c_u c_v Q^2. A pair containing zero, or with u=v, yields -I. For distinct nonzero u,v, c_u=c_v=-1. Compare Q to P_old=L_u^T L_v. Anticommutation of distinct imaginary basis units says that right and left multiplication by e_v differ by the column sign

    chi_v(b)=+1 for b in {0,v}, and -1 otherwise.

Consequently Q^2 and P_old^2 differ exactly at the four columns

    H={0,u,v,u XOR v}.

The basis table on H is quaternionic; P_old^2=-1 on those four columns. Away from H, Q^2=P_old^2. Exactly h-2 k_old-4 columns have Q^2=-1. Each lower column corresponds to one of the cross-half two-cycles, whose positive square condition is -Q^2=+1. This is (4).

The elementary identities used here hold at every stage: distinct imaginary basis units anticommute, conjugation reverses products, and L_u^2=R_u^2=-I for nonzero basis u. They follow by induction from the displayed block formulas (use R_u=c_u C L_u C for the simultaneous right assertion). Basis alternativity and anticommutation also force the cyclic quaternion table on H. No assertion of alternativity for a general sum of basis vectors was used.

## C3. Full support census, exact extrema, and its first moment

Let C_d(z)=sum_(0<=i<j<2^d) z^(k_d(i,j)). Start with C_1=1,C_2=6. For d>=3 and h=2^(d-1),

    C_d(z)=2 C_(d-1)(z^2)
           +2 z^(h-4)[C_(d-1)(z^-2)-(h-1)] +3h-2.          (5)

Despite the notation, the right side is a polynomial with nonnegative integer coefficients: it is the disjoint pair classification (4). Within each half there are two copies of the old unordered pairs. Cross-half pairs with distinct nonzero lower indices give two copies of the old unordered pairs after removing the h-1 pairs containing zero. The other cross-half pairs number 3h-2.

For d>=3 this proves

    [z^0]C_d=(6d-5)2^(d-2)+2,                              (6)
    support(C_d)={0,4,8,...,2^(d-1)-4},                    (7)
    C_d'(1)=(h-4)(h-1)(h-2).                              (8)

For (7), at d=3 only zero occurs. At each later stage the same-half transformation supplies multiples of 8, and the cross-half reflection supplies the intervening values congruent to 4 modulo 8; its zero-class coefficient remains positive after removing h-1. For d>=4 the reflected old exponent could contribute to the new constant term only at k=(h-4)/2, which is 2 modulo 4 and does not occur. Thus the constant obeys A_d=2A_(d-1)+3h-2, with A_3=28, giving (6). Differentiate (5): the two old first moments cancel, leaving (8).

| Dimension | Pair census by real kernel dimension for e_i+e_j |
|---|---|
| 8 | 0:28 |
| 16 | 0:78; 4:42 |
| 32 | 0:202; 4:84; 8:84; 12:126 |
| 64 | 0:498; 4:252; 8:168; 12:168; 16:168; 20:168; 24:252; 28:342 |

In particular, this entire TWO-SUPPORT family has maximal annihilator dimension N/2-4, with every intermediate multiple of four attained. This is different from the classical maximum N-4d+4 for unrestricted elements. They agree at dimension16, then diverge. The proportion of pair supports admitting a real zero divisor tends to one, since the exception count (6) is O(dN) among N(N-1)/2 supports.

## C4. Every odd prime becomes an integral obstruction at dimension sixteen

For a primitive coefficient pair gcd(a,b)=1, an odd prime dividing a^2+b^2 must equal1 modulo4: neither coefficient vanishes modulo that prime, and their ratio squares to -1. This is the elementary prime-divisor direction behind the two-square theorem, not its full converse.

At dimension16, k_4(1,10)=4. At every higher dimension the same embedded pair has k_d(1,10)=N/4. Hence its cokernel has N/4 factors of each of a^2-b^2 and a^2+b^2, with the extra unit factors understood.

For ANY odd prime ell take a=1,b=ell-1. Then

    a^2-b^2=ell(2-ell),    a^2+b^2=2 modulo ell.

The multiplication is invertible over R, but its kernel modulo ell has dimension N/4. Thus every odd prime, including3 modulo4, occurs as a modular cancellation obstruction in this fixed support family. No prime search or unproved distribution assertion is needed.

At a=1,b=2 in dimension16 the cokernel is (Z/3)^4 + (Z/5)^4, equivalently (Z/15)^4; the ranks modulo3 and modulo5 are both12. A primitive two-support octonion multiplication with these coefficients instead has cokernel (Z/5)^4. Equal coefficient norm5 does not predict the exceptional primes once positive blocks occur.

At a=b=1 the same sixteen-dimensional example has cokernel Z^4+(Z/2)^4. Nevertheless x conjugate(x)=conjugate(x)x=2. A two-sided element inverse exists over R while left multiplication is singular. Reassociation is the invalid step in trying to use that inverse for cancellation. The checker retains an explicit nonzero two-support y with xy=0.

## C5. Typical dense norm behavior improves while sparse failures proliferate

Let X,Y be independent random vectors in A_d whose coordinates are independent uniform signs divided by sqrt(N). Both have squared norm1. For d>=3,

    E ||XY||^2 = 1,
    E (||XY||^2-1)^2
        = 16(h-4)(h-1)(h-2)/N^4,  h=N/2.                 (9)

Thus ||XY||^2 tends to1 in probability as d grows, with exact variance asymptotic to2/N. For example, dimensions16,32,64 give respectively21/512,315/8192,3255/131072. The variance need not decrease from the octonions onward: it is zero at8 and becomes positive at16.

**Proof.** Use the inherited rectangle coefficient calculation from [NORM_SIGNS, section1](../session112_counting_interfaces_20261030/norms/NORM_SIGNS.md) and its mixed-monomial orthogonality argument from [NORM_DEFECT, D2](../session112_counting_interfaces_20261030/norms/NORM_DEFECT.md). A violated rectangle is exactly a positive two-cycle in C1. Every mixed coefficient of the unnormalized norm defect is either0 or+/-4; different mixed monomials use distinct unordered coordinate pairs and are orthogonal in expectation. Their total number is sum_(i<j) k_d(i,j)=C_d'(1). The raw second moment is16 C_d'(1). Normalizing both vectors divides a squared biquadratic defect by N^4. Its mean is zero. Insert (8).

The proof actually needs only centered, mutually independent, variance-one raw coordinates, and independence of the two vectors. Fourth moments are not needed because each mixed monomial uses each coordinate only once. The unit-length and displayed probability interpretation use the random-sign choice. Chebyshev gives the explicit upper bound variance/epsilon^2 for a squared-norm error greater than epsilon.

This is a two-random-input statement. It does not say that a random multiplier approximately preserves EVERY vector, nor that these algebras are quasirandom tournaments. Equations (6) and (9) show a precise difference between typical mixed observations and worst-case algebraic identities.

## C6. A compiler for every finite support, with its phase retained

For any nonempty finite support S in F_2^d, choose i0 in S and let

    U=span_F2 {i XOR i0 : i in S},   r=dim U.

For arbitrary coefficients a_i over ANY commutative ring, left multiplication by x=sum_(i in S) a_i e_i maps the coordinate subspace on b+U into the coordinate subspace on b+i0+U. These disjoint input and output cosets give exactly2^(d-r) square blocks of order2^r. On columns b+u, rows b+i0+w, the contribution of i=i0+v is

    a_i S_d(i,b+u) when w=u+v.                            (10)

This is an exact decomposition, not a graph connectivity surrogate. It proves that determinant, Smith, rank and linear-system computations reduce to these blocks, retaining the coset, row/column orders and signs. For s support elements, r<=s-1. Computing all sign entries directly takes O(d N s) elementary sign-recursion operations; dense cubic elimination per block would use O(N 4^r) arithmetic operations over a field. This is not a uniform bit bound for integer Smith algorithms. A single supplied block is smaller still.

The proof is simply the support law i+(b+u)=b+i0+(i+i0+u). This works in every signed-XOR algebra, even ones that fail composition. It explains why low-support integral responses remain tractable far up the tower. It does not imply that blocks on different cosets have identical signs or cokernels.

## C7. Characteristic two restores the laws but destroys orientation transport

Over F_2 every sign disappears. There is an explicit algebra isomorphism

    A_d tensor F_2 = F_2[(Z/2)^d]
       = F_2[t_1,...,t_d]/(t_1^2,...,t_d^2),               (11)

given by e_(2^j) <->1+t_(j+1). This is associative and commutative at EVERY dimension, with nilpotents. Frobenius sends every x to (sum_i x_i)e_0. In particular it kills the entire augmentation ideal, rather than acting as an invertible field clock. Coherent real orientation signs cannot be recovered from this quotient.

For x=e_i+e_j, i!=j, multiplication modulo2 has rank N/2 whatever k_d(i,j) is. Its Gram matrix L_x^T L_x is zero modulo2. Thus a Gram statistic loses even the surviving multiplication information. The integer block presentation, with its signs and valuation layers, avoids this loss.

For the language consisting of all named basis-left-multiplications and scalar-coordinate observation, any exact linear carrier over F_2 needs dimension N: observing the scalar coordinate after multiplication by e_i returns exactly x_i. This family of observations separates every pair of states. There are2^N distinguishable states even for a nonlinear finite carrier. This minimality statement concerns arbitrary binary states and this operation language, not a single orbit or a restricted support promise.

## C8. Every fixed support shape has an all-dimension joint transfer recurrence

Fix r and a subset D of F_2^r containing0 and spanning F_2^r. Retain an injective linear frame G:F_2^r -> F_2^d, a support origin i0, and a column address b. The full joint sign response of its block is

    S_d(i0+G s, b+G u),  s in D, u in F_2^r.              (12)

For EVERY such finite D, the counts of these framed joint sign patterns, as d increases, have an exact finite transfer matrix over Z and hence an integer linear recurrence and a rational generating function. This statement counts framed data with multiplicity; it is not an unnamed-orbit count.

Here is a constructive finite state set. For each query (s,u), retain four bits: the forward and reverse multiplication exponents f,g, and the indicators z_a,z_b that its two current indices are zero. Also retain the row span W of the current matrix G in F_2^r. Start at d=0 with f=g=0,z_a=z_b=1,W={0}. To append a highest coordinate, choose bits o,beta and a row vector v in F_2^r. For each query set alpha=o+s dot v and gamma=beta+u dot v. Its new bits are

    f' = f+gamma(1-z_a) if alpha=0, else g+gamma z_a;
    g' = g+alpha(1-z_b) if gamma=0, else f+alpha z_b;
    z_a'=z_a(1-alpha),   z_b'=z_b(1-gamma),                (13)

with sign equations modulo2, and W'=span(W,v). These are exactly the Cayley--Dickson sign recursion, retaining the zero information needed for conjugation. Accept only W=F_2^r, and read the forward signs. The state count is at most16^(|D|2^r) times the number of subspaces of F_2^r; reachable states can be much fewer. The input alphabet has2^(r+2) letters. The zero letter is an identity transition.

This proves a constant finite transfer for each fixed shape, without a bound on d. The total accepted mass is independently known:

    2^(2d) product_(j=0..r-1)(2^d-2^j).                  (14)

It counts ordered independent frame columns, and both named origins. Changing b within its coset changes the named frame and is deliberately counted. For a property invariant under that change, each actual coset has2^r such addresses; for phase-sensitive queries one must retain them. A specified block determinant, modulus, or Smith-type predicate is a terminal observation on this same transfer. It does not require independent transports for each coefficient.

The Cayley--Hamilton recurrence from this construction is sufficient, not claimed minimal. It is a recursion in algebra dimension d; its operators are bit append and terminal observation. Arbitrary metagraph complements and Collatz iteration are not among its operations.

## C9. Complete three-support sign realizability, certified by finite closure

**FINITE-EXACT with an all-size reduction supplied by C8.** For D={0,1,2} in F_2^2, denote the twelve forward sign bits by f_s(u), with s=0,1,2 and u=0,1,2,3. Exactly256 of the4096 candidate arrays are impossible at every dimension. They are precisely the arrays satisfying ALL FOUR equations

    sum_u f_0(u)=1,   sum_u f_1(u)=1,   sum_u f_2(u)=1,
    f_0(1)+f_0(2)+f_1(0)+f_1(2)+f_2(0)+f_2(1)=0.          (15)

Every other array occurs by dimension64 and at every higher dimension. The four displayed affine parity equations are independent; their simultaneous solution set has size2^8. Their complement is the realizable set, not an affine code itself.

The certificate enumerates reachable C8 states by exact integer transitions. Their numbers at d=1,...,7 are

    16,86,206,1422,4286,4534,4534.

The checker compares the FULL sets at the final step. Because the zero letter is identity, this equality establishes closure under every letter and excludes later new states. The full-rank accepted forward-sign arrays number32,64,1076,3592,3840 at d=2,...,6. Exhaustive comparison with (15) tests all4096 candidates. A second implementation directly constructs every frame and both origins through dimension16 and matches the complete sign histogram, including multiplicities.

This is a complete uniform realization guard for the stated three-support block interface. It does not classify all integer coefficient triples up to Smith equivalence, nor arbitrary four-support blocks or unnamed tournament attachments. The source frame, full joint signs and zero guards are what make the finite closure valid.

## C10. Minimal closed carriers for the named bit-append language

**FINITE-EXACT via deterministic partition refinement on the complete reachable graphs.** Give C8's machine a terminal output equal to the full sign array when the frame has full rank, and a single rejection symbol otherwise. Keep every bit-append letter named. Then its smallest deterministic carrier has

    38 states for r=1, D={0,1};
    4534 states for r=2, D={0,1,2}.                        (16)

Both complete reachable machines are already minimal. Begin by identifying states with the same current terminal output, then refine by the output and all named one-letter successor classes. The partition sizes are

    17 ->37 ->38 ->38;
    3841 ->4529 ->4534 ->4534.

Every strict refinement has a distinguishing word of at most the new depth, by induction. At the fixed point all futures agree within each class, again by induction on word length. Thus these are exact minimalities, and every pair of distinct reachable states is distinguished by some continuation of at most TWO letters. The checker constructs the complete graph by independent breadth-first discovery, verifies transition consistency, and saves a concrete pair with equal accepted current reading but different readings after ONE legal letter.

The minima are for arbitrary reachable prefixes, including currently rank-deficient prefixes which may later become valid frames. They are not minimalities for already-full-rank sources alone, a single trajectory, unnamed basis relabellings, or tournament complement/re-marking. They show in this precise interface why current joint signs alone do not form a closed carrier: the retained conjugation zero tests and reverse signs encode distinctions a later coordinate can expose.

## C11. A literal arithmetic bridge: odd-prime dilation is the binary tower

This is a corollary of C7 and the already proved [prime-dilation D5--D6](../prime_dilation_response_20261030/PROOFS.md), credited there. Fix distinct ODD primes p_1,...,p_d and a positive source n=u product p_i^a_i. Let C=sigma(u). On its binary aliquot-response cube let T_i be multiplication of the source by p_i. The inherited recurrence gives

    (T_i-1)(T_i-p_i)=0, hence T_i^2=1 over F_2.

There is therefore the explicit homomorphism

    A_d tensor F_2 -> binary prime-dilation operators,
    e_v -> product_i T_i^(v_i).                            (17)

If C is odd, the inherited response pairing is nonsingular and this map is a faithful isomorphism onto the full2^d-dimensional operator algebra. It identifies the radical generator e_(2^i)+1 with the finite-difference operator T_i+1. The observable is retained as lambda_n(e_v)=s(n product p_i^(v_i)) mod2; it is NOT the Cayley--Dickson norm.

For P={7,11,13}, the sources1001 and343 have u=1 and give faithful dimension8 realizations. At57 the P-free core is57, with sigma(57)=80 even; the observable quotient has rank1. The checker independently computes the actual integer aliquot responses and verifies the shift intertwiners, recovering the inherited ranks8,8,1. The source-side change, not the ambient dimension, explains that collapse.

This bridge preserves the composition of named prime-dilation histories, their binary linear combinations and the source-tagged scalar response. It discards real multiplication signs and higher integral precision. Applying an aliquot step changes the source functional, so it is not an additional closed operation of (17). In particular it does not identify the orbit1001 ->343 ->57 with a Cayley--Dickson doubling, or transfer a Collatz convergence statement.

## Questions left open by this packet

1. Given the complete sign-realization list C9, classify its four-by-four integer blocks for all coefficient triples up to actual integral transport. Find minimal counting recurrences for their Smith strata; C8 already supplies a sufficient finite recurrence.
2. Simplify the realization guards for other fixed affine support shapes and preserve simultaneous phases when several operators are composed. Individual Smith forms discard the common coordinate alignment.
3. Determine stronger concentration or limiting distributions for the centered norm in C5. The exact second moment alone does not supply a Gaussian limit or an operator-norm bound.
4. Identify an actual guarded map from these sparse integral modules to strong-atom attachment modules. The common use of Smith form is a transferable method, not such a map. The 81/27 transporter pair remains a mandatory hostile if one is proposed.

LRC(14), Collatz convergence, unrestricted metagraph complement closure and a new A000568 value remain open here.
