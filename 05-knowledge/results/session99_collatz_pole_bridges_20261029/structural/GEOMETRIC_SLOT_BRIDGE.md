# A geometric-slot bridge, not a Collatz conjugacy

Status: PROVED elementary identities and a scoped equivariance obstruction.
Source pin: `d9d7e33349653b6b65dbad477841c9ed518c081a`.
Inherited inputs: Session98 `limits/ALL_POLES_CORE_EXPANSION.md` and Session94
`dynamics/VALUATION_TREE_CLOCKS.md`. The former is used from the verified-final
Session98 packet. No originality claim for stars-and-bars or branching formulas.

## 1. The exact common object

For k>=0 let M=3k and N=M+1. The top core at the tournament pole 4^k is
mu=3^k, and its contribution is

    T_k(z)=c_k z^M/(1-z/4^k)^N,
    c_k=(3k)!/(3^k k!) 2^(-3k^2+k).

The empty core has c_0=1. Fix a permutation tau of type 3^k on M ordered
labels. Extending it by f fixed labels means choosing its M-label support
from M+f labels and transporting tau increasingly to that support. The
chosen support is equivalent to a weak composition

    (g_0,...,g_M), sum g_i=f,

where g_i counts fixed labels before, between, or after successive support
labels. Thus the pole order N is the number of insertion slots, not the
number of cycles. Each fixed label costs 2k independent orientation bits,
so every slot has geometric weight t=z/4^k. Core orientation weight and
core multiplicity are entirely in c_k. This gives the denominator (1-t)^N.

For the inherited odd Collatz map S(x)=(3x+1)/2^a with a=v2(3x+1), odd-Haar
valuation digits are iid with mass 2^-a on the punctured coding domain.
Write a=1+b+2G: b has P(b=1)=1/3; independently
P(G=g)=(3/4)4^-g. An even digit creates three children, an odd digit none.
A stopped full ordered ternary tree with k internal vertices has M=3k
edges and N=3k+1 vertices by the rooted-tree identity vertices=edges+1.
Its minimum division clock is 4k+1. Conditional on a fixed shape, the N
vertex marks G are independent, so their sum R has

    P(R=f | k,shape)=binom(3k+f,f) (3/4)^N 4^-f.       (1)

This is EXACTLY the same weak-composition count. Fix any one such tree
shape and enumerate its vertices breadth-first. Sending the insertion-slot
vector to the vertex-mark vector is a degree-preserving bijection. This is
an explicit bridge after fixing a core tau, its ambient label order, and a
tree shape; it is not a bijection of all cores with all shapes.

For a precise normalized GF identity, set

    Q_N(t)=(1-t)^(-N),
    J_k(z)=T_k(z)/(c_k z^(3k))=(1-z/4^k)^(-N).

(The quotient is interpreted by its removable value at zero.) For 0<=q<1,

    J_k(4^k q s)/J_k(4^k q)=((1-q)/(1-qs))^N.          (2)

At q=1/4, (2) is the PGF of the Collatz residual R in (1). This is a
normalized *single top-core tail*, not the entire B, and the normalization
point is z=4^(k-1). The full tree clock PGF conditional on k is

    w^(4k+1) ((3/4)/(1-w^2/4))^N.

Its poles as a function of w are +/-2; these are NOT the tournament pole
locations 4^k. The matching variable is t=w^2/4=z/4^k, after removing
prefactors. In particular a common 4 arises from two binary precision bits
per residual mark, whereas tournament 4^k comes from 2k lost edge bits per
fixed label. The common 3k+1 is a slots/Euler-count identity, not the
arithmetic affine constant in 3x+1.

## 2. Recurrences that transfer, and ones that do not

Let q_(k,f)=binom(3k+f,f), zero for f<0. Then

    (f+1)q_(k,f+1)=(3k+f+1)q_(k,f),
    q_(k+1,f)=sum_(j=0)^f binom(j+2,2) q_(k,f-j).       (3)

The second recurrence adds three geometric slots. For probability weights
p_(k,f)=(3/4)^(3k+1)4^-f q_(k,f), the same recurrence is convolution with
an independent NegativeBinomial(3,3/4) variable. Thus

    E[R|k]=k+1/3, Var(R|k)=4(3k+1)/9.

Both the deterministic core-gap class and a fixed tree-shape mark class
obey this recurrence. It does not say that Collatz evolution increments k.

The leading tournament amplitude A_k=c_k 4^(3k^2) satisfies

    A_0=1,
    A_(k+1)=(3k+1)(3k+2)2^(6k+4) A_k.                (4)

Hence in the shared t coordinate, U_k(t)=T_k(4^k t),

    U_(k+1)(t)=(3k+1)(3k+2)2^(6k+4)
                  * t^3/(1-t)^3 * U_k(t).           (5)

By contrast tree-shape counts F_k=binom(3k+1,k)/(3k+1) satisfy
F_0=1 and F_(k+1)=sum_(i+j+l=k)F_i F_j F_l. The cubic tree recursion does
not enumerate the tournament top cores or their orientation weights.
For k=2 there are 3 ordered full ternary shapes but 40 permutations of
type 3^2 on six labels. Non-top cores also occur in the same pole bank:
mu=(5) has order 6 at pole16, versus order7 for mu=(3,3). The bridge to
N=3k+1 concerns the unique maximal-support core, not every bank summand.

## 3. A precise equivariance obstruction

For k>=1 let X_k be the conjugacy class of permutations of type 3^k in
S_(3k). Let Y_k be ordered rooted full ternary trees with k internal
vertices and their 3k edges bijectively labelled by [3k]. S_(3k) acts on
X_k by conjugation and on Y_k by relabelling edges. There is NO equivariant
map X_k -> Y_k.

Proof: the stabilizer of tau in X_k contains tau itself, a nonidentity
order-three permutation. The action on Y_k is free: an automorphism of
an ordered rooted tree fixes the root and each first/second/third child
inductively, hence every edge; a label relabelling can fix it only if the
relabelling is identity. Equivariance forces Stab(tau) to be contained in
Stab(F(tau)), contradiction. This also prohibits an equivariant bijection.
It does NOT prohibit non-equivariant coding using extra order/root data,
or constant maps to unlabelled statistics. The insertion-slot bijection
above deliberately chooses such extra order and shape data.

There is likewise no arithmetic conclusion from the normalized law:
positive source1 has constant valuation2 and an infinite ternary queue,
although its Collatz orbit is already at1. The terminating word(2,1,1,1)
has source41 mod64 and sends41 to107, so tree extinction is not descent.
Retaining the sum R loses the ordered marks and therefore loses affine
carry/source-cylinder information. The bridge transports the residual
clock law only, not integer convergence, first descent, ROOT forcing,
source guards, or a tournament group action.

## 4. Reproduction

`python3 check_geometric_slots.py` and `python3 -O check_geometric_slots.py`
write identical JSON. Exact finite checks independently enumerate support
subsets/gaps, ternary tree words/marks, conjugation stabilizers at k<=2,
and affine source controls, then verify recurrences with Fraction arithmetic.
These finite checks support the displayed identities; the proofs above
supply the all-k and impossibility statements.
