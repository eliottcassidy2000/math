# The 1–2–4 reflection refinement: counts, carries, and the next bit

**Status: PROVED formal consequences of the cited OCF and the inherited reflection lemma; FINITE-EXACT controls are separate. No literature-priority claim.**

## 1. Inheritance and the three different meanings of 1–2–4

The primary odd-cycle collection formula is due to Grinberg–Stanley,
[*The Redei–Berge symmetric function of a directed graph*](https://arxiv.org/abs/2307.05569),
and is restated in Irving–Omar,
[*Revisiting The Rédei-Berge Symmetric Functions via Matrix Algebra*, v3, Corollary 20](https://arxiv.org/html/2412.10572v3).
We use the tournament specialization

    H(T) = sum_(k>=0) 2^k alpha_k,    alpha_0=1,                 (1)

where H counts ordered directed Hamiltonian vertex words and alpha_k counts
sets of k vertex-disjoint literal directed odd cycles of length at least
three. Cycles are identified by rotation only, not reversal or support.
The weights 1,2,4 are the first three *packing-size weights*, not a claim
that the coefficients alpha_k are binary digits.

Three precise antecedents were recovered by filename enumeration followed
by targeted `git show`; source snapshots and blob pins accompany this note:

* `aliquot_power_124_20261026/tournaments/TOURNAMENTS_124.md`, section 3:
  7=1+2+4 is the proper-divisor sum of 8, with monomial profile 1+X+X^2.
  This is not a graph independence polynomial: its one size-two face would
  require two distinct singleton faces. The valid graph profile 1+3X also
  evaluates to 7; realizing that graph as the *whole* tournament cycle
  conflict graph is a separate gate. We do not identify the two profiles.
* THM-466 already proves the weighted digit tower modulo powers of two.
  Its 2026-10-26 correction explicitly denies recovering alpha from H.
  The carry recursion below makes that inherited forward map explicit.
* Session 130 proves the reflection refinement (2) below. The present note
  supplies actual binary digits, a modulo-16 lift, and hostile controls.

An early `zeckendorf_tournament.py` exploration contains two superseded
claims: that H determines its OCF profile, and an incorrectly shifted path
formula. We do not propagate either. The coefficient-nonuniqueness
correction already exists in THM-466; our actual-tournament check is a replay,
not a new discovery.

## 2. The actual second and third binary digits

Let rho be an involutory anti-automorphism of T. Its action on cycles is
rho_*(C)=reverse(rho(C)). Let q be its number of fixed literal odd cycles,
i its number of nonfixed mirror pairs whose supports intersect, and j its
number of nonfixed mirror pairs whose supports are disjoint. These are
nonnegative **counts**, not bits.

Session 130 gives

    H = 1 + 2q + 4i (mod 8).                                  (2)

For b_r=floor(H/2^r) mod2, this means precisely

    b_0=1,
    b_1=q mod2,
    b_2=(floor(q/2)+i) mod2.                                  (3)

**Proof.** Write q=2u+epsilon, epsilon in {0,1}. Then
1+2q+4i=1+2epsilon+4(u+i). Neither the constant 1 nor 2epsilon
creates a lower-position carry. Hence (3). QED.

Thus the exact low-bit condition is

    H=7 (mod8)
    iff q is odd and floor(q/2)+i is odd
    iff (q=1 mod4 and i odd) or (q=3 mod4 and i even).           (4)

This is a residue criterion, not an assertion H=7. Order-five examples
from Session 130 already give

| adjacency mask | H | q | i | low binary word |
|---:|---:|---:|---:|:---|
| 8 | 9 | 2 | 1 | 001 |
| 40 | 15 | 3 | 2 | 111 |

Bits in the displayed word run high-to-low (b2,b1,b0). The mask uses bit
k for the upper-triangle arc u->v, with pairs (u,v) in lexicographic order.
The first row refutes replacing q and i by their parities before adding;
the second refutes identifying b2 with i mod2 and refutes forbidding residue7.

At even order no fixed odd cycle exists, so q=0 and (3) reduces to
b1=0, b2=i mod2. That special case must not be promoted to odd order.

## 3. General carry theorem and the realizability gate at exact seven

For every finite nonnegative coefficient sequence alpha with alpha0=1,
set c0=0 and recursively

    b_r=(alpha_r+c_r) mod2,
    c_(r+1)=floor((alpha_r+c_r)/2).                             (5)

Extend alpha by zeros. Then these b_r are exactly the binary digits of
sum 2^r alpha_r. In particular,

    b1=alpha1 mod2,
    b2=(alpha2+floor(alpha1/2)) mod2,
    b3=(alpha3+floor((alpha2+floor(alpha1/2))/2)) mod2.           (6)

**Proof.** At each r the invariant is
sum_(k<r)2^k alpha_k = sum_(k<r)2^k b_k + 2^r c_r.
Adding 2^r alpha_r and applying Euclidean division by two gives the
next invariant. Eventual zero coefficients make the carry terminate. QED.

Equation (5) is the correct relation between packing coefficients and
ordinary binary expansion. Binary 7=1+2+4 does not license alpha=(1,1,1).
If H is *exactly* seven, positivity in (1) gives alpha1+2alpha2=3 and
alpha_(k>=3)=0. The two integer candidates are (alpha1,alpha2)=(3,0)
and (1,1). A disjoint cycle pair includes two distinct cycles, excluding
(1,1). Therefore exact H=7 would require the *whole* conflict graph K3.

In the reflection setting it would further require j=0 and

    (q,i)=(3,0) or (1,1).                                    (7)

These are necessary orbit profiles, not constructed tournaments. The
separate inherited all-order H7 exclusion has an additional tournament
realizability proof (THM-343 or THM-1370); it does not follow merely from
(2), and no all-order 7/21 exclusion is newly asserted here.

The map alpha -> H is not injective even among actual tournaments:
H=9 occurs with alpha=(1,4) at order five and alpha=(1,2,1) for an ordered
stack of two directed triangles at order six. This is exactly the
THM-466 corrected hostile. Numerical equality does not identify supports,
coefficients, graph carriers, or binary digits before carrying.

## 4. A modulo-16 packing-orbit theorem

Let p be the number of **nonfixed two-element orbits of disjoint two-cycle
packings** under rho_*. These are orbits of packings, not orbits of cycles.
Let t be the number of fixed three-cycle packings. Then

    H = 1+2q+4i+8(j+p+t) (mod16),                             (8)

and the actual fourth binary digit is

    b3=(floor((q+2i)/4)+j+p+t) mod2.                           (9)

**Proof.** An anti-automorphism fixes at most one vertex: two fixed vertices
would force their mutual arc both ways. Every rho-invariant odd support
contains a fixed vertex. Thus all fixed odd cycles meet; in even order
there are none. A fixed two-cycle packing consequently cannot consist of
two separately fixed cycles; it is exactly one disjoint mirror pair.
Therefore

    alpha1=q+2i+2j,       alpha2=j+2p.

The involution on three-cycle packings gives alpha3=t mod2. Truncating
(1) modulo16 and substituting gives

    1+2(q+2i+2j)+4(j+2p)+8t
      =1+2q+4i+8(j+p+t) (mod16).

Finally floor((1+2q+4i)/8)=floor((q+2i)/4), since q+2i is an integer.
This proves (9). QED.

A fixed three-packing consists necessarily of one fixed odd cycle and one
disjoint mirrored pair, all three supports mutually disjoint. Conversely
any such choice is fixed. This provides a direct combinatorial formula for t.
In even order q=t=0, so

    H=1+4i+8(j+p) (mod16),
    b3=(floor(i/2)+j+p) mod2.                                 (10)

The new 8-weight is not a fourth unconstrained bit either. Counts j,p,t
must again be reduced only after the lower carry in (9) is added.
In fact the first orbit decomposition gives an exact nonnegative lift
at every order:

    H=1+2q+4i+8R,
    R=j+p+sum_(k>=3)2^(k-3) alpha_k.                         (10a)

Thus the disjoint mirror pairs do not disappear combinatorially. Their
4j contribution from single cycles and their 4j contribution from fixed
two-packings combine to 8j: a genuine upward carry, zero only after
reduction modulo8. The modulo16 theorem is R=j+p+t mod2.
For n<=8, t=0 and all higher packings vanish, so the identity strengthens to

    H=1+2q+4i+8(j+p)               (exact for n<=8).            (11)

Indeed three odd cycles would use at least nine vertices. Thus at these
orders j+p is exactly the higher integer quotient left after subtracting
1+2q+4i; it is not merely an observed modular correction.

The checker finds actual same-order, same-(q,i) tournaments with different
H modulo16. Hence the new packing-orbit information cannot generally be
recovered from (q,i). Separate hostile fixtures show that each of j,p,t can
change the fourth digit when omitted from (8).

## 5. The central cut as a binary carry source

In Session 130's centered odd-order setting, delete the fixed center z to
obtain B. Write r_z for the nonfixed mirror cycle pairs containing z, and
K for the exact nonnegative weighted cut charge, so

    H(T)=H(B)+2K,       K=q+2r_z (mod4),
    i(T)=i(B)+r_z.

Since B has even order, H(B)=1+4i(B) mod8. Consequently

    b1(T)=K mod2,
    b2(T)=(i(B)+floor(K/2)) mod2.                              (12)

The apparent separate cut and cycle descriptions agree because
floor(K/2)=floor(q/2)+r_z mod2. This identifies the actual bit, rather than
mistaking the cut count K or the fixed-cycle count q for a binary digit.

## 6. What does and does not transfer from Zeckendorf

A Zeckendorf expansion uses nonconsecutive Fibonacci *indices*: its
selection indicator is an independent set of a path. A cycle packing is
an independent set of Omega(T). This is a common independence-language
interface, not an identification of numeration systems.

For the path P_m, deletion of an endpoint gives

    I(P_m,x)=I(P_(m-1),x)+x I(P_(m-2),x),
    I(P_0,x)=1, I(P_1,x)=1+x.

At x=1 this counts nonconsecutive selections with Fibonacci growth. At x=2
it counts the same selections with two decorations for every chosen vertex:

    I(P_m,2)=(2^(m+2)+(-1)^(m+1))/3.                           (13)

**Proof.** The recurrence has characteristic roots 2 and -1; the stated
formula has the two required initial values 1,3 and satisfies the
recurrence. QED. The values are 1,3,5,11,21,43,...; the early exploration's
formula (2^(m+1)+(-1)^m)/3 has the wrong index and already fails at m=1.

No conclusion here says Omega(T) can be an arbitrary path. Nor do
Zeckendorf's nonadjacent selection bits force nonadjacent binary digits of
H: a directed triangle gives H=3=(11)_2, and a verified five-vertex
self-converse tournament gives H=15=(1111)_2. The decisive data remain the
realizable cycle complex, its reflection action, and integer carries.

## 7. Verification scope

`check_binary_carries.py` is self-contained, uses explicit RuntimeError
checks rather than asserts, and imports no producer or repository code.
Hamiltonian paths are counted by subset DP; literal cycles are enumerated
by rooted DFS; every disjoint packing and its reflected image are checked.
It exhausts all labelled tournaments compatible with the standard reversal
anti-involution at orders1..7, not all labelled tournaments, then checks
seeded samples at orders8,9 and ordered triangle stacks through order12.
Any involutory anti-automorphism of a tournament is conjugate to this
standard permutation by its fixed-point count, so the small-order census
also covers every isomorphism class admitting such an anti-involution.
Repeated labelled presentations remain distinct in this census.

The receipts record exact universes, counts and adjacency-row witnesses.
Theorems (3), (5), (8), (9), (11), (12), (13) rest on their proofs, not
extrapolation of the finite tests. No H-spectrum completeness, general
nonrigid path-orbit spectrum result, or Collatz conclusion is claimed.
