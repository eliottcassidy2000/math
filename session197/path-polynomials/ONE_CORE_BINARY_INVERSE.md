# One stationary core between two binary histories

Session 197. This is an exact physical two-integer normal form and finite inverse for a specified tournament language. The bounded certificate below is not an all-order spectrum theorem.

## 1. Sources, conventions and the physical language

For a nonempty rooted tournament (A,a), write h=H(A), e for its Hamilton paths ending at a, d=H(A-a), and R for its labelled two-bin path configurations in which the first bin is nonempty and its path starts at a. The other bin may be empty. H(empty)=1. A two-bin configuration partitions all vertices, and each occupied bin carries one directed Hamilton path of its induced vertices.

The stationary singleton moves used here are inherited from S195 combined-growth:

https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/combined-growth/PROOF.md

Call them O and Z here, to keep the binary encoding visible. O appends a universal sink. Z appends a vertex beating only the retained root a. Their exact state laws are

    O(h,e,R,d) = (h,0,2R,d),
    Z(h,e,R,d) = (h+R-e,d,2R-d,d).                 (1)

The positive four-state interface and its rank certificate are inherited; no novelty is asserted for them. The arbitrary-block stationary sewing law, as opposed to moving the root into the new block, is proved in RETAINED_ROOT_SEWING.md. Its attribution is clarified in RETAINED_SEWING_SOURCE_SCOPE.md. Both are dependencies, not new claims here.

The scalar binary recurrence used below is also inherited. It is exactly f(q)=g(q+1) in the reviewed modular-gadget theorem:

https://github.com/eliottcassidy2000/math/blob/48e10c6a4773de23f9ca19ef120b17012a0d7812/session197/modular-gadget-universality/PROOF.md

That source uses its own canonical start-one physical compiler. We do not identify that rooted state with the stationary prefix state below merely because their scalar counts agree.

Fix one actual nontransitive rooted core (B,b). Starting from a rooted singleton a, apply any finite O/Z prefix, sew B on the right by orienting all old vertices toward B except reversing a->b to b->a, retain a, and apply any finite O/Z suffix. Exactly one nontransitive block is attached. Different words or cores can give isomorphic graphs or the same count; this is a construction language, not a unique tournament decomposition. We also allow no core as a separate language.

A stationary attachment of a transitive right block of order m, marked in position j from0, is literally O^j Z O^(m-1-j), with vertices labelled in creation order. Conversely every prefix can be implemented by transitive blocks of unrestricted order, allowing harmless initial Z padding only at the singleton seed. A suffix containing a Z partitions literally into such transitive blocks. A suffix consisting only of O steps is scalar-neutral only because it is final, and is removed as in Section3. Inserting an extra Z after the core is not generally neutral. Thus the language allows arbitrary binary histories; it is not limited to the order-at-most5 transitive letters of the earlier finite bank. Internal O or Z steps are never discarded.

## 2. The prefix integer and its physical inversion formula

Define

    F(0)=1,
    F(2q+1)=F(q),
    F(2q)=F(q)+q+(q mod2), for q>=1.               (2)

Read the binary digits of q from most significant to least, with O=1 and Z=0. The empty word represents q=0. Applying (1) inductively from (1,1,1,1) gives the full state

    (h,e,R,d) = (F(q), (q+1) mod2, q+1, 1).        (3)

Leading zero digits are full-state-neutral Z steps at this seed. They may be removed before any later core, because (3) preserves all four coordinates required by stationary continuation. This is not deletion of arbitrary padding from inside a history.

There is a direct cycle interpretation. Write b_i for bit i of q. All vertices outside a form a transitive chain in creation order. Every directed cycle contains a. Its first chain vertex is an earlier O, its last is a later Z, and the intervening chain vertices can be any subset of even cardinality for an odd cycle. Consequently, by the tournament odd-cycle collection identity,

    F(q)=1+2 sum_(i>j) b_i(1-b_j) w_(i-j),         (4)
    w_1=1,  w_l=2^(l-2) for l>=2.

There are no disjoint pairs of these cycles because they all contain a. Formula (4) in particular proves positivity and oddness, and independently reproduces (2).

Taking the converse reverses the chain order and complements the incidences with a. For any fixed width l with 0<=q<2^l, including leading zero padding,

    F(q)=F(2^l-1-rev_l(q)),                       (5)

where rev_l reverses exactly l bits. This asserts equality of scalar counts under actual converse graphs, not equality of their rooted endpoint states.

## 3. Every suffix has one scalar integer

For an arbitrary rooted input (h,e,R,d), the scalar output of any O/Z suffix is

    h+t(R-d)+d(F(t)-1)+(t mod2)(d-e),              (6)

where t is the integer whose bit i is1 exactly when suffix step i, starting at0, is Z. Thus the suffix reads the bits of t from least significant to most, with Z=1 and O=0. Any trailing O steps are scalar-neutral and may be removed only because the construction ends there. For t=0 choose the empty suffix.

Proof: by (1), the coefficient of R in the final scalar is sum_(Z steps i)2^i=t. The coefficient of the initial e is minus1 precisely when the first step is Z, hence minus(t mod2). The coefficient of h is1. To determine the remaining coefficient of d, start with (h,e,R,d)=(1,1,1,1). The chronological O=1/Z=0 word is the reverse complement of the bits of t. Equations (3) and (5) show that its count is F(t), which determines the coefficient of d and gives (6). All state transfers are linear, so this determines the formula for every input, without assuming that arbitrary coordinate quadruples are physical.

## 4. The exact one-core formula

For the core define

    H=H(B), E=end_b(B), S=start_b(B),
    Q=R_b(B), L=L_b(B), C=B_B(2).

L is the analogous two-bin response with one designated nonempty bin ending at b, and C is the complete labelled two-bin count. Set epsilon_q=(q+1) mod2. The stationary sewing law and (3) give the middle state

    h_q=H F(q)-S epsilon_q+L(q+1),
    e_q=E, R_q=C(q+1)-Q, d_q=H.                  (7)

Combining (6) and (7) proves

    N_B(q,t)=H[F(q)+F(t)-1]+(L+Ct)(q+1)
             -S epsilon_q-t(Q+H)+(t mod2)(H-E).  (8)

Every q,t>=0 has the literal construction just given, and every word with this core has some such pair with the same count. This is a scalar normalization; omitted leading prefix Z or terminal suffix O vertices need not preserve isomorphism or graph order. The canonical witness uses

    1+|B|+bit_length(q)+bit_length(t)              (9)

vertices, with bit_length(0)=0. Deleting the retained root leaves an ordered join of two transitive chains and B, so its count is H. This language is not claimed to lie in the earlier C2 class: its outputs include67.

## 5. A complete target inverse and explicit bounds

The elementary physical inequalities are

    H>=S,E,  L>=H,
    C>=(|B|+1)H,  Q<=C/2.                       (10)

For the cut bound, split any Hamilton word into two consecutive labelled pieces, allowing empties; the word and cut are recovered by concatenating the bins. For Q<=C/2, forget the bin labels on any configuration counted by Q; its unique path containing b is prescribed to start at b, and the two choices of bin labels are distinct full configurations. The remaining inequalities follow from extending or isolating an endpoint in the labelled-bin interpretation and the inherited sewing proof.

Since a nontransitive tournament has at least3 vertices,

    a=C-Q-H>0,
    A_q=C(q+1)-Q-H >= a(q+1).                    (11)

No zero denominator occurs for this language. Singleton or transitive cores are not needed: the no-core language is treated separately. From (6), F(t)>=1 and H-E>=0, while (7) gives h_q>=L(q+1). Therefore

    N_B(q,t)>=h_q+t A_q
             >=(q+1)(L+a t).                    (12)

For every positive bound M, all outputs N<=M are found by the finite exact test

    0<=q<floor(M/L),
    retain q only if h_q<=M,
    0<=t<=floor((M-h_q)/A_q),
    evaluate (8) and retain N_B(q,t)<=M.          (13)

Necessity follows from (12); sufficiency follows from the actual graph constructor. No recurrence inverse or guessed physical state is used. For an exact target, keep only equality in the final test. For K=floor(M/L), at most

    K+(M/a)(1+log K)

pairs are tested when K>=1, since sum_(j=1)^K 1/j<=1+log K. This is an arithmetic-operation count apart from the O(log M) binary evaluations of F, not a bit-complexity assertion. A table of F values removes those repeated evaluations.

For q,t>=1, (12) implies qt<=M, and hence bit_length(q)+bit_length(t)<=floor(log2 M)+2. When one index is zero the sum is at most floor(log2 M)+1. Thus (9) is at most

    |B|+3+floor(log2 M).                          (14)

All physical inputs are retained in the construction; only canonical leading or final neutral steps identified above are replaced.

For no core, strip terminal O steps using F(2q+1)=F(q). The remaining integer is0 or even, q=2v, and F(q)=F(v)+v+(v mod2)>=v+1. Therefore all no-core values <=M are exactly obtained from q=0 and the even integers 2<=q<=2M-2. This separate finite inverse is necessary because unrestricted q has infinitely many neutral representations of1.

## 6. Complete finite bank certificate and limits

The checker independently inventories all labelled tournaments of orders1 through5 and all their roots, using literal Hamilton words and labelled cuts for the six core statistics. There are50 distinct nontransitive signatures, inherited as a subset of the65-signature stationary bank. One representative graph per signature suffices for (8), since all q,t depend only on these six numbers.

Using all50 cores and the no-core language, the exact inverse through M=1000000 has499,997 odd outputs. Its only missing odd values are7,21 and462661. It evaluates14,197,377 core-index pairs; the no-core part contributes249,077 distinct values. Thus462661 is the first additional gap in this language, not merely an unsuccessful search at one chosen pair. The output receipt includes the full core inventory, sorted support digest, per-core pair counts, and literal adjacency witnesses for selected difficult values. Separate subset path counting checks those actual witnesses and all small prefix/core/suffix graphs.

The gap is repaired by two actual nontransitive cores. Use prefix q=2, then B1 with rows[10,0,3,6], mark0, then B2 with rows[8,1,3,6,15], mark4, always retaining the original root, then suffix t=112. The signatures are respectively (3,0,1,6,3,24) and (5,0,5,32,5,64). The prefix state is (3,1,3,1); the first core gives (17,0,66,3), and the second gives (415,0,4128,15). Since F(112)=107, formula (6) gives

    415+112(4128-15)+15(107-1)=462661.

The literal graph has19 vertices and is independently counted by a full subset Hamilton-path DP in the checker. Therefore this is a strict difference between physical construction languages, not a forbidden Hamilton count. The missing value factors as127*3643 and is not divisible by7; consequently the one-core bank contains every odd prime and every seven-times-odd-prime value through999999 other than7 and21. No primality distribution or infinite coverage is asserted.

These finite claims do not assert that every larger odd integer except7 and21 is attained, that any fixed core is sufficient, or that either parameter is an affine digit operation on the final count. They also do not identify this stationary language with the start-one source oracle. The exact inverse makes a concrete all-order coverage problem in two nonnegative integers and a finite physical core bank; it does not solve that problem.

Run `python3 verify_one_core_binary_inverse.py` and its `-O` form. The proof and verifier depend on the frozen retained-root packet and helper, whose hashes are recorded in the receipt. The new verifier does not import any exploratory census.
