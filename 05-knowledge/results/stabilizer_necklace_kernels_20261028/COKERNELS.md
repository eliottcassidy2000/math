# The complete integer cokernel is a tower of golden residue rings

**PROVED** C1-C3; local labels are qualified by this filename. This concerns
the specific recursive sign word **in its specified cyclic coordinate**.
It is not an invariant of the unmarked tournament or of the colouring orbit.
The clock-preservation hostile is in C3.

## C1. Exact decomposition, not just an equal determinant

Put N=3^k, P(z)=1+z-z^2, D_k(z)=product_(j=0..k-1)P(z^(3^j)), and
R_N=Z[z]/(z^N-1). For the integer circulant M_k from [KERNELS](KERNELS.md),
identify coordinate i with z^(-i); multiplication by D_k then has exactly
the matrix entries d(j-i). Let O=Z[phi], phi^2=phi+1. There is an additive
group isomorphism

    coker_Z(M_k)
       = direct_sum_(j=0..k-1) [O/(phi^(3^(k-j))-1)]^(3^j).       (C1)

At k=0 both sides are the trivial group. For k>=1 this is a finite group.

Proof of the CRT step: the ideals(P(z^(3^i))) in R_N are pairwise comaximal.
Suppose two with i<j belonged to one maximal ideal. The residue field has
prime characteristic p, and the image zeta of z satisfies zeta^N=1.
For p=3, zeta=1, but P(1)=1. For p!=3, zeta has3-power order.
Set a=zeta^(3^i), b=a^(3^(j-i)); both are roots of P.
If a=b, then a has order dividing3^(j-i)-1 as well as a power of3,
forcing a=1, impossible. If a,b are different, their product is-1. For
odd p this is impossible because the product of odd-order elements has
odd order. For p=2 the roots of P are the primitive cube roots, so
a^(3^(j-i))=1, again impossible. This covers repeated-root characteristic5
as well through the a=b case. Every maximal ideal has prime residue
characteristic because R_N is finite integral over Z. Thus the ideals are
comaximal, proving the CRT decomposition of R_N/(D_k).

For one factor put d=3^j and m=N/d. Writing w=z^d makes R_N free of rank d
over Z[w]/(w^m-1), with basis1,z,...,z^(d-1). Imposing P(w)=0 identifies
w with phi, since P is minus its monic minimal polynomial. Its additive
quotient is therefore[O/(phi^m-1)]^d, proving(C1).

The cyclic clock is retained, not reset independently on each copy. On the
j-th block, multiplication by z acts as

    (u_0,...,u_(d-1)) -> (phi*u_(d-1),u_0,...,u_(d-2)).          (C2)

After d steps it multiplies every entry by phi; after N steps it is identity.
In the original index convention this is the shift i->i-1. Independent
scalar residue lists without this twisted shift would discard the clock.

## C2. Two new closed recurrences and all Smith factors

For odd m, |O/(phi^m-1)|=L_m. Hence

    Delta_0=1,
    Delta_k=Delta_(k-1)^3 L_(3^k)
           =product_(r=1..k) L_(3^r)^(3^(k-r)).                 (C3)

These determinants are positive: the constant eigenvalue is1 and the
remaining real-circulant eigenvalues form nonzero conjugate pairs. The
first values are1,4,4864,50541367159422976.

The inherited [odd-period ring theorem](../weighted_crt_clocks_20261028/CRT.md), C1,
gives, for every m=3^r,r>=1,

    O/(phi^m-1) = F4 x Z/(L_m/4).

Write q_r=L_(3^r)/4. The Lucas triple-angle identity for odd m,
L_(3m)=L_m^3+3L_m, gives

    q_1=1,   q_(r+1)=16 q_r^3+3q_r
                       =q_r(16q_r^2+3).                        (C4)

In particular q_r divides q_(r+1), all q_r are odd and1 mod3, and the
new factor16q_r^2+3 is coprime to q_r. We do not assert these factors
are prime or squarefree. The first q_r are1,19,109801,... .

Equation(C1) now determines the **entire Smith normal form**. For k>=1 put
b_k=(3^(k-1)-1)/2. The N diagonal entries, in divisibility order, are

    1;
    2 repeated N-1-b_k times;
    2q_r repeated 3^(k-r) times, for r=2,...,k.                  (C5)

Indeed the F4 summands contribute N-1 copies of Z/2, and the nontrivial
odd cyclic summands have orders q_r in a divisibility chain. Pair one
2-factor with each odd cyclic factor; all remaining2-factors stay alone.

Examples:

    N=9:  1, 2 (seven copies), 38;
    N=27: 1, 2 (22 copies), 38 (three copies), 219602.

An independent Euclidean integer Smith computation verifies(C5) through27;
Bareiss determinants independently verify(C3) through81. The attempted
generic Smith computation at81 was stopped because of coefficient growth;
no receipt claims that larger literal calculation completed. See
[check_cokernels.py](check_cokernels.py) and [its JSON](check_cokernels.json).

Since19 first appears in q_2 and the subsequent cofactors are coprime to19,

    nullity_F19(M_k)=b_k.                                    (C6)

This equals the minimum regular-subset size for C_k in(S6). Both equal the
same geometric depth sum1+3+...+3^(k-2). On the linear side these are the
characters over an extension field with z^(3^j)=5, where5 is the order-nine
root of P in F19. On the symmetry side they are the minimum contributions
of independent paired levels. These source objects and their symmetries
differ; an equivariant identification between them is not claimed.

## C3. The cyclic coordinate is a necessary sidecar

For the depth-two sign word, the recursive colouring stabilizer is trivial,
but it has81 isomorphic A_2-translates. Keep the numerical cyclic order
fixed while translating the colouring. Their determinants have the census

    4864:30,   1792:18,   1024:18,   256:15.

Thus(C3) is not a colouring-orbit invariant. Even the three translates
under C_2 give different observations:4864,1024,1024. Their complete
three-entry future profiles under C_2 are all distinct. Within this exact
three-state source orbit, the smallest carrier retaining the fixed-clock
determinant under re-marking therefore has3 states, although the unmarked
colouring quotient has1. The hostiles are checked as actual permutations.

The repair is to retain and transport the cyclic clock and its origin along
with the colouring. Simultaneous relabelling then conjugates the resulting
matrix and preserves its Smith form. Transporting only the colouring while
silently keeping the old clock is a different operation. This makes the
new recurrence usable as a typed instrument without confusing it with a
solution to arbitrary tournament-composition closure.
