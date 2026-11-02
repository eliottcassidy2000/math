# A minimal one-edge arithmetic reflector outside the induced-extension cone

**Status: PROVED (universal graph and polynomial identities); FINITE-EXACT
(minimal exterior budget and unique minimal pair up to reversal/relabeling).**
The full odd Hamiltonian-path spectrum and bounded exact spectral width remain
OPEN. No priority claim is made.

## 0. Inheritance, board, and scope

The inherited Session137 constructors `THEORY.md` proves universal fresh-block
source maps are nonnegative integer maps in actual source coordinates, and
that sufficiently large original-port cover cost seals a right-local history
into multiples of its original H. Its `BOUNDED_INTERIOR.md` gives dyadic normal
forms at bounded parent height and interior budget. Session134
`descent/THEORY.md` proves that every fixed positive **induced** increment has
bounded endpoints, hence a finite union of nested induced rays has density
zero. Session137 `coverage/PHASE2.md` shows full-star unit connectivity can fail.

Our anchor is actual exact realization, not further modular coverage. The
live board is: actual source ancestry; old-edge changes; two-path-cover cost;
profile cancellation; exact versus local width. The hostile near miss is
calling a small scalar difference an induced insertion, or a scalar identity
an identity of future profiles. The least-used sidecar here is the ordered
pair `(H(U), E(U))`, where E counts ordered two-path covers allowing empty
paths. We change **one old exterior arc**. No addition/deletion of old vertices
or re-marking is presented as the same operation.

For nonempty U let p_U(c) count covers by c unordered nonempty directed paths,
and put

    F_U(z)=sum_(c>=1) p_U(c)(z)_c,
    h=F_U(1)=H(U), E=F_U(2)=2(H(U)+p_U(2)).

The inherited homogeneous substitution kernel is linear in this polynomial.
The polynomial odd-cycle formula, proved explicitly in Session137
`width_phase2/PROOF.md`, is

    F_T(z)=sum_C 2^|C| z^(|T|-sum_(gamma in C)(|gamma|-1)),

where C ranges over vertex-disjoint collections of directed odd cycles.
We use both identities with their full graphs retained.

## 1. Two literal contexts differing in one arc

Rows below are outgoing-neighbor bitmasks on vertices 0,1,2,3,4. Distinguish
vertex 0 as the substitution port. Define

    Q_minus = (16,17,3,7,12),
    Q_plus  = (16,25,3,5,12).

They differ only by reversing `3->1` to `1->3`. The port is dominated by
1,2,3 and dominates 4. The unchanged exterior arcs are

    2->1, 1->4, 3->2, 4->2, 4->3.

For any nonempty tournament U, replace port0 homogeneously by U. Call the
resulting literal tournaments R_minus(U),R_plus(U). They have exactly the
same vertices and edges except for the single old exterior edge `{1,3}`.
Both quotients are strong, so both substituted tournaments are strong even
when U itself is reducible.

### The universal arithmetic reflector

For every nonempty tournament U,

    H(R_minus(U)) = 11h+12p_U(2) = 6E-h,
    H(R_plus(U))  = 13h+12p_U(2) = 6E+h.              (1)

In particular that one edge reversal has exact increment `2H(U)`.

**Elementary coefficient proof.** Maximal U-run splitting expresses each
response as `sum c! W_c p_U(c)`. Every return to U must pass through its sole
exterior outneighbor4, so no quotient word can use U more than twice. The
single-run counts are H(Q_minus)=11 and H(Q_plus)=13. The six two-run words
for Q_minus are

    043210, 104320, 204310, 210430, 304210, 310420;

those for Q_plus are

    042130, 043210, 104320, 130420, 210430, 304210.

Multiplying by2! gives the coefficient12. This is an exact count of a fixed
finite context, not an extrapolation from sampled U. Equation(1) follows.

The minus context is already a width-three construction:
`C3[L2, K1->U, K1]` with cyclic block order as displayed. The plus context is
an actual order-five quotient. At U=K1 the counts11 and13 recover the inherited
exact width contrast: 13 has width5, while11 has width3. No claim that all
other plus outputs have exact width5 is made (25 already has a width3 join
witness).

## 2. Full polynomial law and why cancellation is genuine

Write F=F_U and F_+=F(z+1),F_-=F(z-1). The actual full-profile operators are

    F_Rminus(z) = (z^4-2z^2)F
                +z(z+1)(2z+1)F_+
                -z(z-1)(2z-1)F_-,                    (2)

    F_Rplus(z)  = z^4 F
                +2z(z^2+z+1)F_+
                -2z(z^2-z+1)F_-.                    (3)

Their difference is the particularly short exact operator

    Delta F(z)
      = z[F(z+1)-F(z-1)]
        -z^2[F(z+1)+F(z-1)-2F(z)].                  (4)

The signed coefficients are **evaluation formulas for actual graph
operations**, not permission to subtract tournaments or admit nonactual
source states. At z=1, F(0)=0, and (4) gives precisely2h.

**Proof on a polynomial basis.** A transitive module L_k has F(z)=z^k,
k>=1. Every cycle meeting that module must contain exterior vertex4, and
its U part is one increasing nonempty subset. No two odd cycles can be
disjoint: any U-meeting cycle meets4 and also at least one of1,2,3; the only
exterior triangle avoiding4 is123. The induced exterior tournament has
exactly two triangles in each context.

After a chosen increasing module subset, a cycle goes to4 and follows a
nonempty directed path among1,2,3 starting at2 or3. The numbers of such
paths of lengths t=1,2,3 are respectively

    Q_minus: (2,3,1),       Q_plus: (2,2,2).

For a module subset of size s, the cycle has s+1+t vertices and is odd
exactly when s+t is even. Consequently, letting

    O=sum_(s odd) binom(k,s) z^(k-s)
      =[(z+1)^k-(z-1)^k]/2,
    V=sum_(s even,s>=2) binom(k,s) z^(k-s)
      =[(z+1)^k+(z-1)^k-2z^k]/2,

the minus polynomial is

    z^(k+4)+4z^(k+2)+2[(2z^3+z)O+3z^2 V].

This simplifies to(2). The changed path counts give
`Delta F=2(zO-z^2 V)`, which is(4), and adding gives(3).
The substitution kernel is linear in p_U, hence a linear operator on F_U.
The transitive polynomials z^k, k>=1, form a basis of polynomials with zero
constant term. Therefore the identities hold for every actual U, not merely
for transitive modules. QED.

### A hostile to profile monotonicity

For U=L_k, k>=1, formula(4) at z=2 gives

    Delta F(2)=8*2^k-2*3^k-6,
    Delta p(2)=4*2^k-3^k-5.                            (5)

At k=4 this changes H by +2 but p(2) by -22. Higher context responses are
not monotone; `(H,E)` is sufficient for the scalar outputs(1), not for
arbitrary subsequent substitution. This is precisely how an old-edge edit
escapes positivity statements whose hypothesis is universal induced
fresh-block extension.

## 3. Unbounded exact fixed-gap strong pairs

For ordered joins, both H and E multiply: restrict the one or two ordered
paths to each block. A singleton has H=1,E=2. For a fixed actual nonempty U
and r>=0, put `U_r=U->L_r`, allowing L_0 empty only as padding notation. Then

    H(U_r)=h,              E(U_r)=2^r E(U).

The two actual strong tournaments R_minus(U_r),R_plus(U_r), differing in one
exterior arc, have counts

    6E(U)2^r-h,            6E(U)2^r+h.                 (6)

Their gap is the fixed positive integer2h while both counts tend to infinity.
Thus every **attained half-gap** h has unbounded strong one-edge-flip
realizations. The condition that h is attained is explicit; this does not
classify all possible fixed one-edge-flip gaps.

In particular, for U=L_k, k>=1,

    H(R_minus(L_k))=6*2^k-1,
    H(R_plus(L_k)) =6*2^k+1.                           (7)

This is an infinite exact family of neighboring odd counts connected by one
arc reversal. It is not a global interval theorem: the family is dyadic and
has only O(log X) outputs below X. At a fixed k the two tournaments have the
same order, so neither is a proper induced extension of the other. Across k,
the natural nested padding paths have exponentially growing gaps. There is
no conflict with the bounded-induced-gap or nested-ray barriers.

## 4. A source-certified exact inverse, and its resource cost

For any core C with no singleton SCCs, let `(h_C,e_C)=(H(C),E(C))`.
Allow the empty core with `(h_C,e_C)=(1,1)`. Removing all singleton SCCs
from U preserves H, and E multiplies over the ordered SCCs. Consequently

    (H(U),E(U)) = (h_C, e_C*2^r),                       (8)

where r is the number of removed singleton SCCs. Conversely C followed by r
singletons realizes the pair. If C is empty, require r>=1 because U is
nonempty; otherwise r>=0. Order and positions of the deleted singleton
components do not affect this pair, although they affect the original graph.

The inherited Session133 `core_reduction/THEORY.md` proves the exact budget

    K(B)=max {3a+4b : 3^a*5^b<=B}
        =max_(b=0,1,2;5^b<=B) [4b+3 floor(log_3(B/5^b))].

Every actual core with H<=B has at most K(B) vertices. This is an inherited
bound, not a claim proved by the reflector construction.

Therefore, at a prescribed input height h, enumerate actual cores of order
at most K(h) with H=h, retaining their adjacency matrices and e_C. A target N
is a reflector output from that height **iff**, for some retained core, sign
sigma in {-1,+1}, and allowed r,

    (N-sigma*h)/(6e_C) = 2^r.                           (9)

Acceptance returns the literal graph C->L_r followed by the chosen context.
No formal `(h,E)` pair is accepted without an actual core witness. Enumerating
h<=floor(N/11) is already a complete variable-budget inverse, because (1)
and E>=2h imply R_minus>=11h and R_plus>=13h. This is an exact decision for
this constructor image, not a positive inverse for every odd target.

A useful sharper small-target bound follows from `E(U)>=(|U|+1)H(U)`:
each Hamiltonian path together with one of its |U|+1 cuts is a distinct
ordered two-path cover. If h>1, U is nontransitive and |U|>=3, hence every
minus output is at least23h and every plus output at least25h. Thus all
reflector outputs through200 use either h=1 or h<=8. Here K(8)=4, and complete
core enumeration gives just `(h_C,e_C)=(1,1),(3,12),(5,32)`. The exact
all-order image through200 is

    11,13,23,25,47,49,69,75,95,97,141,147,187,191,193,197.

The verifier retains actual core-plus-padding witnesses for all16 values.
The omission of other targets is only from this reflector image.

### Quantitative budget obstruction

Let A(X,B) be the number of distinct reflector outputs <=X whose input U has
H(U)<=B. Counting every labelled core (an intentional overcount), equation(8)
gives the explicit bound

    A(X,B) <= 2[1+floor(log_2(X+B))]
              sum_(n=0)^K(B) 2^(n(n-1)/2).             (10)

Every sign/base pair contributes at most the displayed number of padding
exponents. The empty-core nonemptiness condition can only reduce this bound.
Since K(B)<=3 log_3 B, if

    log B(X) = o(sqrt(log X)),

then A(X,B(X))=X^o(1). More explicitly, with base-two logarithms, when
`B(X)<=2^(c sqrt(log_2 X))`, the upper exponent is

    A(X,B(X)) <= X^(9c^2/[2(log_2 3)^2]+o(1)).          (11)

For every fixed `c < sqrt(2)*log_2(3)/3`, this is a zero-density target set.
Thus the old-edge operation genuinely bypasses the fixed-induced-gap barrier,
but a reflector compiler covering a positive-density target population
cannot keep its module height below this stretched-exponential-in-bit-length
threshold. These are necessary resource conditions, not existence claims.

## 5. Minimal nontrivial cancellation context: exact finite classification

Call a pair of one-port contexts a **nontrivial unit isolator** if its entire
H-response coefficient rows differ by `(2,0,0,...)` and at least one common
higher-run coefficient is positive. The latter condition rules out trivial
ordinal multiplication contexts that have no growth at transitive modules.
The contexts must have the same number of exterior vertices, but this
minimality statement does not initially restrict how their exterior edges
differ.

**FINITE-EXACT theorem.** Four exterior vertices are necessary and sufficient.
Among contexts on at most four total vertices, the possible coefficient rows
are

    order1: (1),
    order2: (1),
    order3: (1,0), (3,2),
    order4: (1,0), (3,0), (3,2), (5,4), (5,6).

No two rows with a positive common higher coordinate differ by2 only in the
first coordinate. The lists exhaust every labelled orientation with port0;
fixing the port loses no context because labels may be permuted.

At order5, the displayed Q_minus,Q_plus supply the isolator. Among all1024
labelled order-five contexts, and all six possible exterior edge flips,
exactly48 **oriented low-to-high** moves are nontrivial unit isolators. Every
one has response rows `(11,12,0)` and `(13,12,0)`. They form two port-preserving
relabeling orbits of24 moves each, interchanged by global edge reversal.
Thus the pair is unique up to relabeling and global reversal in the minimal
one-edge setting. These classifications are exact finite facts supported by
complete independent enumeration, not an all-order uniqueness assertion.

## 6. Verification, loss ledger, and next boundary

`verify.py` is self-contained standard-library code. It independently uses:

* literal distinct quotient words for the response rows;
* subset Hamiltonian DP and set-partition cover DP on actual substitutions;
* integer polynomial shifts for(2)-(4);
* exhaustive labelled context/edge universes for minimality;
* explicit adjacency witnesses for both sides of(7), with strongness checks;
* mutation controls that reject the wrong edge, wrong coefficient, and the
  false profile-monotonicity extension.

Run `python3 -B verify.py` and `python3 -B -O verify.py`; receipts must agree.
No source code from inherited packets is imported by the verifier.

Connection contract: actual U -> two actual same-vertex strong contexts,
with sidecars h,E and full polynomial F. The graph map preserves U exactly,
changes only one exterior arc between its two outputs, and has controlled
scalar gap2h. Scalar output loses the full profile, source ancestry and
spectral-width minimality; equations(2)-(4) retain the profile when needed.
The target inverse still asks whether an actual h,E can land at a desired
integer N=6E+/-h. No arbitrary integer pair is admitted as an actual module,
and this construction does not prove all permitted odd integers attainable.
