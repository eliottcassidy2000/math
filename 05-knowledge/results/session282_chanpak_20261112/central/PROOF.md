# Central fusion banks: prime ownership and the exact absolute carry

**Session282. PROVED, independently reviewed and freshly replayed after restoration.** Chan–Pak remains
OPEN. The arithmetic below is elementary bounded coin/CRT algebra; novelty is
not claimed for that algebra. The contribution is its exact physical carrier
interface, a target decoder retaining the absolute carry, and the restriction
it puts on repeated central gadget designs.

## Inheritance and connection

Use the physical rooted-fusion matrix law from
[Session264](../../session264_constructive_20261111/fusion/PROOF.md) and the
complete central-gate criterion from
[Session281](../../session281_chanpak_20261112/wildcard/PROOF.md), the latter
already published at6301e35a before this session. A nontrivial rooted
matrix has diagonal `(h,Q,0,h)` and its `(1,4)` entry is `H-h`.
A central bank is a finite list of actual rooted tournaments whose seven
coordinates `(H,h,Q,U,V,d,e)` agree except for H. The singleton boundary
is harmless for the algebra; the physical banks below are nontrivial.

The distant source is
[THM-5103, faithful Fermat arithmetic](../../../../01-canon/theorems/THM-5103-faithful-fermat-path-count-arithmetic-and-sharp-coprime-bank-cost.md).
There, private prime ownership makes a multiplicative encoding faithful.
Here the map sends a fusion presentation to its normalized count difference:
prime ownership instead determines whether the available additive mesh is
primitive. This map loses the baseline and the bounded digit sets. Restoring
both gives the exact carry test below. No Fermat-prime or external conjecture
is assumed.

## 1. Exact central convolution

Let bank i have common deletion count h_i and counts `b_i+2d`, for d in a
finite integer set E_i containing0. Fix an order of fusion and arbitrary fixed
physical contexts between the banks. Let H0 be the count with every d=0,
and D the product of all deletion counts, including the fixed contexts.
Then every choice has

    H = H0 + 2 sum_i (D/h_i) d_i.                    (1)

Every other carrier coordinate is independent of the choices. Indeed, put
E=E14. For every physical matrix M, `ME=h_M E`, `EM=h_M E`, and `E²=0`.
Expand the ordered product after replacing each selected factor by
`M_i+2d_i E`. All terms using two E factors vanish, including those with
intervening matrices. A single E term has coefficient D/h_i. This proves
(1) without commuting the physical factors.

Consequently (1) is a **renewable full-carrier bank**, not just an output
interval. Every later common two-sided fusion multiplies its differences
by the deletion counts of the added contexts.

Let g_i=gcd{|d|:d in E_i}; ignore a singleton bank E_i={0}. The gcd of all
normalized output differences is exactly

    G = gcd_i ((D/h_i) g_i).                          (2)

The differences generate the ideal in(2), and selecting only one nonzero digit
shows the reverse inclusion. A gcd of1 does not imply the finite support is
an interval.

If every g_i=1 and there are no fixed nonunit contexts, then

    G = D/lcm(h_1,...,h_r).                           (3)

For a prime p the valuation is `sum_i v_p(h_i)-max_i v_p(h_i)`.
Thus G=1 **iff no prime divides two different h_i**, equivalently the
h_i are pairwise coprime. Repeating a nonunit deletion type destroys this
primitive mesh, however many copies are used. Fixed contexts contribute
their deletion-count product as an additional unavoidable factor.
This is a per-presentation assertion, not a sparsity theorem for unions of
different baselines or arbitrary tournament operations.

### Physical repeated-gate control

The inherited literal pair

    A=(56,13,25,48,34,6), B=(52,13,24,49,34,6)

has H=41,43 and agrees in every other carrier coordinate, with h=15.
Fusing k independent choices gives exactly k+1 distinct counts

    H(A*...*A) + 2*15^(k-1) j,       0<=j<=k.          (4)

This is an arbitrarily long renewable arithmetic progression, but its mesh
is15^(k-1), not1. Merely adding more interchangeable copies cannot make it
an odd interval. The controls literally count the order11 two-factor case.

## 2. A complete absolute-target test when deletion types are coprime

Now omit nonunit fixed contexts and assume pairwise coprime h_i>1. Write
D=product h_i and w_i=D/h_i. For an odd prescribed target N put
`z=(N-H0)/2`. Compute its unique residue digit

    r_i = z*w_i^(-1) mod h_i,      0<=r_i<h_i.         (5)

Define the available carry sets and required absolute carry

    K_i={ (d-r_i)/h_i : d in E_i, d≡r_i mod h_i },
    c=(z-sum_i w_i r_i)/D.                            (6)

CRT guarantees c is an integer. Equation(1) realizes N **iff**

    c belongs to K_1+...+K_r.                         (7)

Proof: reduce(1) modulo each h_i to obtain(5). Every allowed digit is
uniquely r_i+h_i k_i; substitute and divide by D. Conversely such a choice
of k_i reconstructs actual bank members, hence an actual fusion tournament.

This test includes the absolute baseline H0; passing every residue test
only asserts each K_i is nonempty and does not prove(7). For E_i containing
exactly one representative of each residue, there is exactly one carry
value for each target residue. The many other integers in that residue
class are not represented by this fixed bank.

When E_i=[0,m_i] with m_i>=h_i-1, each K_i is the complete interval
`[0,floor((m_i-r_i)/h_i)]`. Hence(7) reduces to the two integer inequalities

    0 <= c <= sum_i floor((m_i-r_i)/h_i).              (8)

A greedy bounded allocation of c gives the digits and the explicit physical
construction. Unlike modular surjectivity, (8) is an exact target selector.
It is conditional on availability of the specified physical central banks.
No arbitrary integer carrier tuple is assumed realizable.

## 3. Guaranteed interior and a quantitative constructive supply target

For the complete digit banks above let

    S=sum_i w_i m_i,
    C=(r-1)D-sum_i w_i+1.                             (9)

Every integer z in `[C,S-C]` is represented. If this interval is empty,
there is no assertion. Thus a physical central odd interval of
`S-2C+1` values is supplied whenever this number is positive.

Proof: R=sum_i w_i r_i <=rD-sum_i w_i=C+D-1. If z>=C,
R-z<=D-1. Since R-z is divisible by D, it cannot be positive, so c>=0.
Apply the same argument to S-z and complementary digits m_i-d_i: its
nonnegative carry condition is exactly the upper inequality in(8).
Therefore both conditions hold. For r=1, C=0 and the result is the original
bank. For r=2 this is the familiar bounded two-coin conductor
`(h_1-1)(h_2-1)`, expressed in physical deletion units.

This separates three obligations that should not be merged:

1. **Primitive mesh:** pairwise disjoint prime support and unit digit gcd.
2. **Digit supply:** enough actual central realizers in each residue.
3. **Absolute ownership:** the carry inequalities, including H0.

The current known central41/43 gate has only two digits and h=15, far below
the complete-residue condition. Neither(9) nor existing modular completeness
supplies the missing physical banks. The theorem is a constructive interface
specification, not a Chan–Pak proof.

## 4. Small hostile to residue-only reasoning

Take arithmetic bank data h=(3,5), E_1=[0,2], E_2=[0,4]. The sum
`5E_1+3E_2` has one representative of every residue modulo15, but misses1
(and many other positive integers). At z=1, the canonical residues are
r=(2,2), and c=(1-10-6)/15=−1. The exact selector rejects it immediately.
These are arithmetic controls, NOT claimed physical central banks of those
deletions. An arbitrary bank of matrices is not a tournament.

## 5. Physical cross-lane application: primitive mesh still needs ownership

The new [carry-loop pair](../wildcard/PROOF.md) has counts2747/2749 and common
full carrier deletion169. Compose it on the right of the inherited41/43
central pair. The four actual sixteen-vertex tournaments have counts

    474053,474083,474391,474421.

Their shared deletion is2535=15·169; normalized choice weights169 and15 have
gcd1. The exact target test selects precisely these four outputs with
E_1=E_2={0,1}, and rejects474055. Independent path recurrence checks all four
literal counts. This physical example separates primitive mesh from interval
coverage; it does not realize the complete digit banks hypothesized in§3.

This section was reconstructed after environment rollback and freshly replayed;
the final session audit binds the exact restored files.
