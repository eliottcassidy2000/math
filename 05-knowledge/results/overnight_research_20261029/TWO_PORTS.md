# Two port groups and the cost of preserving a clock

Status: **PROVED elementary** for P1, P2 and P4 below; **FINITE-EXACT**
for the complete fixed-block tournament census P3. No arithmetic
realization as a Tate-Shafarevich group is asserted.

Write `J=[[0,1],[-1,0]]` and `H_m=(Z/m)^2` with its standard perfect
alternating pairing. The case `H_1` is trivial. A nonsingular integral
alternating matrix K gives the discriminant pairing
`b([x],[y])=x^T K^-1 y mod Z` on coker K. Integral alternating reduction
and hyperbolic peeling are inherited from
[Session96](../session96_sha_pairing_growth_20261029/tournaments/PAIRING_GROWTH_MODES.md).

## P1 Complete classification for two active hyperbolic blocks

Fix positive integers d,e and let B range over **all integer 2 by 2
matrices**. For

    K = [dJ, B; -B^T, eJ],       t = de - det B != 0,

put `s=gcd(d,e,all entries of B)`. Then, as an unmarked paired group,

    coker K = H_s orthogonal-sum H_(|t|/s).                 (P1)

Consequently the attainable groups, with d,e fixed, are exactly

    { H_s orthogonal-sum H_(s n) : s | gcd(d,e), n >= 1 }.  (P2)

**Proof.** The Pfaffian is t. The first alternating invariant factor is
the gcd of all entries, namely s; the product of the two invariant
factors is |t|. Every entry is divisible by s, so s^2 divides t. These
facts give (P1), including the divisibility chain. Conversely choose

    B = s diag(1, (d/s)(e/s)-n).

Its content together with d,e is exactly s, and its Pfaffian is s^2 n.
This realizes every group in (P2). Hyperbolic reduction classifies the
pairing as well as the underlying group. The basis chosen in reduction
does not canonically preserve the input ports. QED.

In particular s=1 gives **every H_n and hence every prime**, even with
d=e=3 fixed. Orthogonal sums of these constructions realize every finite
perfect alternating paired group. A single four-dimensional attachment
has at most two hyperbolic factors; it does not realize arbitrary numbers
of generators in one step.

| target | B for d=e=3 | Pfaffian | paired group |
|---|---|---:|---|
| 25 | diag(1,-16) | 25 | H_25 |
| 27 | diag(1,-18) | 27 | H_27 |
| 57 | diag(1,-48) | 57 | H_57 = H_3 orthogonal-sum H_19 |

These are weighted lattice attachments. Their cross entries are not
tournament orientations. The decomposition at 57 uses coprime primary
parts; 25 and 27 instead increase a prime-power exponent. This makes the
three numbers different arithmetic tests despite the common construction.

**Hostile to keeping only the determinant.** With d=e=3, B=0 gives
H_3+H_3 and B=diag(1,0) gives H_9. Both have order 81. Named ports and
the integer lattice, not just the order, remain necessary.

## P2 A recurrence that retains every named joint port

Let M be any nonsingular integral alternating matrix, and let the columns
of P be the named integer ports. Retain `G=P^T M^-1 P`. For an integral
alternating coupling matrix E, set `K=M+P E P^T`. Provided K is
nonsingular,

    G_new = (I+G E)^-1 G,                                 (P3)
    det K = det M det(I+E G).                             (P4)

To append two objects use block-diagonal M and P and
`E=[0,H;-H^T,0]`; a two by two H couples two ports on each side.
To retain or re-mark ports by an integer matrix R, update
`G_new -> R^T G_new R`. Thus (P3) is closed for subsequent updates whose
ports remain in the retained span. Proof: the resolvent identity
`K^-1 P=M^-1 P(I+E G)^-1`, followed by multiplication by P^T; the
determinant identity is Sylvester's identity. Neither requires G or E
to be invertible.

The recurrence saves work on future rational pairing queries. It does
not reconstruct the full integral cokernel: lattice indices, named
vectors, and their embedding remain a sidecar for that task. With a
clock T on the lattice and dual port action C satisfying `T^T P=PC` and
`T^T M T=M`, an attachment preserves that clock precisely when
`T^T P E P^T T=P E P^T`. It is safest to test this equality on the
actual matrices; covariance of port vectors must not be confused with
covariance of linear functionals under transpose.

For a direct join `[M,B;-B^T,N]` with individual isometries T_M,T_N, the
exact clock condition is simply `T_M^T B T_N=B`. Losing it loses the
clock even when the paired groups happen to be isomorphic.

## P3 Every legal cut between the fixed four vertex blocks

Use the prior block

    M = [ 0  1  1  1 ]
        [-1  0 -1  1 ]
        [-1  1  0 -1 ]
        [-1 -1  1  0 ].

In this displayed order **Pf(M)=-3**, not +3. The previous note's sign
is corrected in place; its paired group H_3 and the product Pf(M)^2=9
were correct. Every B in {+1,-1}^{4 by 4} gives a labelled eight-vertex
tournament K. Exactly 128 cuts have rank one and 6,272 have rank two.
The other ranks have 36,864 and 22,272 cuts respectively.

For rank two the full paired groups and multiplicities are:

| group H_m | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 21 | 27 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| labelled cuts | 576 | 768 | 1152 | 576 | 64 | 576 | 1152 | 384 | 576 | 384 | 64 |

All rank-one cuts give H_9. Thus the attainable prime support here is
`{3,5,7,11,13,17}`. H_25 and H_57 are absent from this fixed-block,
rank-at-most-two universe; H_27 is present. This is not an all-order
obstruction to those absent groups in other tournaments.

The producer exhausts 65,536 bit patterns. An alternative audit generates
rank-two cuts from pairs of normalized row directions: there are
`choose(8,2)*(4^4-2*2^4)=6272`. This construction is complete because a
sign vector in the span of two independent sign vectors is one of their
four signed copies. Direct Pfaffians and Pfaffian minors independently
check the compressed Schur computation and paired factors.

## P4 The simultaneous three cycle removes all new primes

Let P fix vertex 0 and rotate 1,2,3. Among the above rank-at-most-two
cuts, impose `P^T B P=B`. There are sixteen:

| cut rank | group | number |
|---:|---|---:|
| 1 | H_9 | 8 |
| 2 | H_9 | 4 |
| 2 | H_27 | 4 |

Here the exclusion of primes other than 3 also has a short proof. A
commuting B has a constant first row/column off the corner and a
three by three circulant lower block. A nonconstant sign circulant has
rank two on the nontrivial rational C3 representation. Its nonzero
invariant part then forces total rank at least three. Hence rank at most
two forces a block-constant cut

    B = [a, b b b; c, d d d; c, d d d; c, d d d],
    a,b,c,d in {+1,-1}.

For u=(0,1,1,1), `e_0^T M^-1 u=-1`. The two-port Pfaffian formula
from [J3](../sha_marked_join_f1_20261029/JOIN_RECURRENCES.md) gives

    Pf K = 9(1-ad+bc).                                    (P5)

The determinant `ad-bc` is 0,2 or -2, with multiplicities 8,4,4.
The finite minor check supplies the full group H_9 or H_27, rather than
merely its order. Preserving this named clock is the operation that
eliminates the 5 and 11 constructions from the preceding session.

## P5 An all weight counting recurrence and a fast exact census

Restrict every entry of B to [-L,L]. Define the product multiplicity
`f_L(k)=#{(a,d) in [-L,L]^2 : ad=k}` and

    D_L(t) = sum_k f_L(k) f_L(k-t).

This counts the labelled two by two matrices of determinant t. Since
f_L is even, D_L is a polynomial square: if

    F_L(z) = sum_(k=-L^2)^(L^2) f_L(k) z^(k+L^2),

then `D_L(t)=[z^(t+2L^2)] F_L(z)^2`. In particular
`sum_t D_L(t)=(2L+1)^4`.

For fixed d,e and s dividing gcd(d,e), the exact number of nonsingular
cuts producing H_s+H_(s n), n>=1, is

    sum_(epsilon=+1,-1) sum_(h | gcd(d/s,e/s), h^2 | t_epsilon)
         mu(h) D_floor(L/(s h))(t_epsilon/h^2),            (P6)
    t_epsilon = (d/s)(e/s) - epsilon n.

**Proof.** Write B=sC. Its Pfaffian has sign epsilon and absolute value
s^2 n exactly when det C=t_epsilon. The condition that its first
alternating factor is exactly s is
`gcd(d/s,e/s,all entries C)=1`. Divisor Mobius inversion enforces this
primitive-content condition; dividing C by h divides its determinant by
h^2 and changes the box to floor(L/(s h)). This proves (P6).

This gives a literal role for the Mobius function: remove common content
before classifying the paired group. Liouville's total factor parity
would not perform that primitive-content inversion. The square in h^2
comes from the determinant's degree two, not from an unproved arithmetic
duality.

The implementation uses carry-free Kronecker encoding. Choose an integer
base larger than (2L+1)^4, pack the coefficients of F_L, square that one
integer, and unpack. It constructs the entire determinant distribution
with O(L^2) coefficient work plus one integer square on O(L^2 log L)
bits, instead of enumerating (2L+1)^4 matrices. This bound is for the
determinant distribution; divisor queries in (P6) have their additional
explicit cost. Full labelled cuts and clocks cannot be reconstructed
from this histogram alone.

`weighted_counts.py` agrees with direct enumeration on 109,080
nonsingular attachments, including zero-count classes, and accounts for
all 104,060,401 matrices in the L=50 box by one packed square. From
H_3,H_3 the least L producing H_25,H_27,H_57 is respectively 4,4,6.
For example determinant -18 occurs already at L=3, but every such cut
has content 3 and yields H_3+H_9 rather than H_27. The Mobius correction
detects this precise distinction.

## What this contributes to Selmer arithmetic

P1 removes a **representation** obstruction: weighted two-port lattices
can model every finite alternating group. P3 and P4 quantify how native
entries and clocks restrict that freedom. Actual Sha finiteness still
requires arithmetic input. For a fixed curve and prime ell, the exact
sequences

    0 -> E(Q)/ell^k E(Q) -> Sel_(ell^k)(E) -> Sha(E)[ell^k] -> 0

show what must be certified: local conditions, compatible transition
maps, and a bound on the quotient after removing the Mordell-Weil part.
A uniform finite bound on |Sha[ell^k]| forces that primary group to be
finite; full Sha finiteness additionally needs only finitely many
nonzero primary parts. These are arithmetic obligations, not consequences
of the formal recurrence. In particular an F1 specialization of counting
polynomials supplies neither obligation.

There is a more quantitative primary criterion. If
`d_ell=dim_Fell Sha(E)[ell]` and every element of Sha(E)[ell^infinity]
is killed by ell^a, then its order is at most ell^(a d_ell). Prove this
by induction using the multiplication-by-ell kernel, or by applying that
argument to each finite generated subgroup; their ell-torsion ranks are
at most d_ell. Actual d_ell is finite by finite Selmer theory. Thus a
source-preserving, all-level exponent certificate would suffice at one
prime. The preceding one-port recurrence does not supply it: its copy
parameter grows the mod-ell dimension and is not the Selmer tower of a
fixed curve. P1 also shows that unrestricted two-port gluing gives no
automatic finite-prime-support bound.
