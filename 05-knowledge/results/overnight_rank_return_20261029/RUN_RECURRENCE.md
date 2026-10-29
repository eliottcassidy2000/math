# Strong targets and a joint-port recurrence recovered from the March family

2026-10-29, 11:00 UTC continuation. **PROVED scoped all-size formulas;
FINITE-EXACT strong-target census through seven.** Same-agent review.
The interval family is not closed under actual path complement. The
all-size integral selected-path conjecture and general attachment
realization remain open; none of the following sequences is A000568.

## R1. A successful complement can stay strong through seven

Use the three-run operator O_(T,P)=(M_T;M_(T^op xor E(P))) and arbitrary
port weights from [COMPLEMENT_KERNEL K1--K4](COMPLEMENT_KERNEL.md).

**FINITE-EXACT.** Every one of the396 strong classes with3<=n<=7 admits
recovery at every modulus with at most one actual complement, while every
performed complement has a strong target. If M_T is already unimodular,
no move is needed. Otherwise strong_target.json supplies a Hamilton path
with a strong target and coprime maximal minors of O_(T,P).

| Order | Already unimodular | One strong-target complement |
|---:|---:|---:|
|3|1|0|
|4|1|0|
|5|3|3|
|6|7|28|
|7|16|337|

The source bank and its independent labelled-mass exhaustiveness proof are
the preceding checkpoint's. The new verifier replays the actual source and
target matrices, their integer minors and strong connectivity. It does not
infer strong-target success from a class count or from rational rank.

Zero moves are essential to the statement: all three Hamilton complements
of C3 are transitive, so a requirement to perform a positive-length move
while staying strong already fails at three. Its three-singleton response
is invertible before moving. The meaningful all-size strong-target question
is therefore about sources that still need an additional observation.

## R2. Inheritance and the actual interval source

Let S_n have vertices0,...,n-1 and arcs

    i -> i+1;       j -> i whenever j>=i+2.               (R1)

The March5 theorem, recorded in
[tiling-class-structure F3](../../../03-artifacts/drafts/tiling-class-structure.md)
and [verify_tribonacci_proof.py](../../../04-computation/verify_tribonacci_proof.py),
already proves its Hamilton count h_n, with

    H(z)=sum_(n>=0) h_n z^n=(1-z^2)/(1-z-z^2-z^3).        (R2)

The correction lineage is MISTAKE-645 and
[THM-4731, statement6](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md).
The old guess2^(n-2)+1 is false from seven. No new priority is claimed for
tribonacci counting. The October20
[future-boundary propositions3--4](../tournament_future_boundary_20261020/README.md)
already retain the run endpoints and the actual complemented relation.

The inherited bijection is particularly useful here: a Hamilton path of
S_n consists of consecutive label intervals in descending interval order,
each traversed increasingly. Its interval lengths form a composition of n
with no two adjacent parts both1. If the top and bottom lengths are a,b,
its start and end are(n-a,b-1).

Closest mechanism: that run/port bijection and the previous integer recovery
criterion. Canonical hostile:81/27, supplemented by the native six-vertex
exit in R6 below. Corrected near miss: a scalar recurrence does not retain
joint endpoints or the target relation. Least-used sidecar: the first and
last interval of EACH colour. The board is interval source / named colour /
joint endpoints / integral clock / strong target / literal complement.
Anchor: strong-target recovery; niche: March run arithmetic; wildcard:
evaluate one exact response without constructing its large matrix.

## R3. Complete endpoint kernels and forced unit matching

Put D(z)=1-z-z^2-z^3 and define three coefficient sequences

    h(z)=(1-z^2)/D(z),
    a(z)=(1-z)/D(z),
    b(z)=z^2/D(z).                                      (R3)

Let E_n(s,t) count Hamilton paths of S_n from s to t. The following is a
direct endpoint refinement of the inherited composition bijection:

* E_n(0,n-1)=1.
* For s>t, set m=s-t-1. Use h_m if t!=0 and s!=n-1;
  use a_m if exactly one of t=0,s=n-1 holds; use b_m if both hold.
* Every other entry is zero.

The first rule includes n=1. For n=2 only E_2(0,1) is nonzero.
To prove the second rule, remove the prescribed top and bottom parts.
The remaining composition has total m. A boundary part of length1 forbids
an adjacent singleton. With L=z^2/(1-z), the two-state transition matrix

    W(z)=[[L,z],[L,0]]

records whether the preceding part is a singleton. The four boundary
series e_alpha^T(I-W)^(-1)(1,1-beta)^T are h,a,a,b. For m=0 the two
singleton boundaries give b_0=0, so this formula also handles that guard.

**Corollary (PROVED).** For every n>=3, E_n has determinant(-1)^(n-1) and
permanent1. Row0 forces column n-1; then row1 forces column0, row2 forces
column1, and so on. Every selected entry is1. Rotating the last column to
the front makes the matrix unit lower triangular. Thus its single-run
weight response is invertible over every coefficient ring.

This uses the same forced-matching mechanism as the independently audited
[balanced endpoint extension](../balanced_fibre_endpoint_extensions_20261020/AUDIT.md#3-full-endpoint-unimodularity),
on a different explicit source. It does not make its three-run matrix
unimodular: S6 has three-run rational rank100 of120.

The marginals needed below also have constant-time formulas after prefix
sums of h and a are precomputed. The start count f_n(0)=1. For1<=s<n-1,
f_n(s)=a_(s-1)+sum_(m=0)^(s-2)h_m; for s=n-1 it is
b_(s-1)+sum_(m=0)^(s-2)a_m. The terminal count at t is f_n(n-1-t).

## R4. All ordered cover counts have an exact minimal recurrence

Let C_(n,k) count spanning covers of S_n by k ordered, nonempty directed
paths. Assign each vertex the colour of its path. A monochromatic subset
splits into maximal consecutive intervals. Between different intervals of
one colour, every arc points from the higher interval to the lower one.
Therefore the path on that subset visits those intervals in that order,
and its count is the product of their h_lengths.

First allow j colours, including unused colours, and write F_j(z) for the
resulting weighted word series. Each maximal colour run contributes
R(z)=H(z)-1. There are j choices for its first colour and j-1 thereafter:

    F_j(z)=1+jR/(1-(j-1)R)
          =(1-z^2)/(1-jz-z^2-jz^3).                     (R4)

For j=0 this reduces to1. Excluding unused colours gives the all-size law

    sum_(n>=0) C_(n,k) z^n
      =sum_(j=0)^k (-1)^(k-j) binom(k,j) F_j(z).         (R5)

For fixed k>=1 the minimal eventual homogeneous constant-coefficient
recurrence over Q has order exactly3k, with characteristic polynomial

    product_(j=1)^k (X^3-jX^2-X-j).                     (R6)

Proof of minimality: the k denominator cubics are pairwise coprime. A
common root of two would satisfy z(1+z^2)=0, whereas their values there
are1 or2. Neither denominator vanishes at z=+1 or-1, so its numerator
1-z^2 cannot cancel. The nonzero inclusion-exclusion coefficients leave
every pole. The j=0 constant affects initial terms but not the denominator.

For k=3, the initial nonzero values at n=3,... are

    6,36,186,828,3414,13428,51162,190764,700590,2545236,...

This is a degree-nine recurrence for the actual three-run mass of S_n,
not an extrapolation from a fit. It evaluates totals using O(kn) integer
arithmetic operations through size n, or companion powering for isolated
large n. This statement counts arithmetic operations, not unit-cost bit
operations or all tournament isomorphism classes.

## R5. A joint-port compiler with three local coordinates per colour

Fix distinct ordered starts(s_1,...,s_k) and ends(t_1,...,t_k). The same
interval decomposition evaluates the exact native entry M_k[s,t], without
constructing the (n)_k-by-(n)_k matrix.

For colour i retain three coordinates: e for no interval yet, G for a
terminal-only response, and F for the joint start/terminal response.
On a new interval I processed from LOW labels to HIGH labels, define

    H_I=h_|I|,
    g_I=number of interval paths ending at t_i,
    f_I=number of interval paths starting at s_i,
    E_I=number starting at s_i and ending at t_i.

An endpoint outside I gives zero. The exact update is

    (e,G,F) -> (0, g_I e+H_I G, E_I e+f_I G).            (R7)

The new higher interval becomes the first interval of the path. Hence it
replaces the source endpoint, while retaining the old terminal. In
particular the update does not multiply unrelated start/end marginals:
the first interval uses the actual joint E_I. Subsequent intervals have
uniform inter-interval arcs, which is precisely why their product is legal.

Tensor these three-coordinate updates across the k named colours. Sum
over interval partitions, keeping the last colour to forbid equal adjacent
colours, and extract the all-F coordinate at position n. Start with all-e
at position zero. This constructs every colour word once and every native
cover with its correct endpoint multiplicity, proving the algorithm.

R3 makes every local coefficient available in O(1) arithmetic operations
after O(n) preprocessing. The direct dynamic program uses
O(k^2 3^k n^2) arithmetic operations and O(k3^k n) stored integers. The
3^k local coordinates are a sufficient representation, not a proved minimum.
They describe an evaluator for one fixed query, not a faithful compressed
state containing the entire arbitrary joint-port tensor.

Exact long queries with starts(n-3,n-2,n-1), ends(0,1,2) give

| n | Native three-run entry |
|---:|---:|
|20|20970288|
|40|1840936361018199878|
|100|1245494853107606036773627246585960940439542446873568|

These are evaluations of the proved recurrence. Direct native enumeration
independently checks every three-run entry through five,96 fixed-seed entry
pairs at each size six through nine, and32 two-run pairs at each size two
through seven. All cover arities are independently checked by permutation
cuts through eight. No direct size100 enumeration is claimed.

## R6. The actual complement relation still exits this family

At n=6, S6 has rows[2,4,9,19,39,15]. The legal path
P=(2,3,4,5,0,1) gives the strong target

    U=[30,56,58,48,32,1].                               (R8)

This is the successful integral strong-target witness in R1. Both sources
have H=17 and the same score multiset, but they are not isomorphic. The
multisets of(degree, sum of out-neighbour degrees) are respectively

    S6: (1,1),(1,2),(2,4),(3,6),(4,7),(4,8);
     U: (1,1),(1,4),(2,2),(3,4),(4,7),(4,10).

Thus the scalar recurrence cannot substitute for the actual relation after
this very operation. The inherited exact target formula is

    tau_P(S_n)=T_up xor (I xor E(P)),                    (R9)

where I is the set of adjacent-label pairs. It retains every cut edge and
every long hinge of the chosen run composition. R5 sums endpoint responses;
it does not retain each internal path witness. General quotient-walk
assembly still needs the full complemented relations and their alignment.

## R7. The recurrence clock has one exact arithmetic loss direction

The homogeneous recurrence in R4 uses the integer matrix

    A_j=[[j,1,j],[1,0,0],[0,1,0]],       j>=1.           (R10)

**PROVED.** For every r>=0 its Smith form is

    Smith(A_j^r)=(1,1,j^r).                             (R11)

For r=0 this is the identity. For r>0, det(A_j^r)=j^r. At every prime
dividing j, A_j reduces to A_0, whose powers have rank two: the leading
two-coordinate block is a coordinate swap and remains invertible. At all
other primes the matrix is invertible. Therefore no prime divides either
of the first two invariant factors; the determinant supplies the third.

For an arbitrary initial three-coordinate boundary state modulo m, the
exact image size is m^3/gcd(m,j^r). At m=p^a the image needs
2a+max(a-r*v_p(j),0) base-p digits. One step has the exact integer image
guard y_0-y_2=0 modulo j and inverse

    (u,v,w)=(y_1,y_2,(y_0-jy_1-y_2)/j).                 (R12)

Thus loss depends on both precision and elapsed clock r. This is a literal
integer recurrence operator, with compatible reduction maps. It is not
Frobenius, which is invertible on a finite field, and not the C3 group
action on run names. The full source state is arbitrary here; the one
enumerative sequence F_j has a particular prescribed initial state.

Incoming [primary clock completion](../overnight_research_20261029/PRIMARY_CLOCK_COMPLETION.md)
separates adjoining roots from preserving the original source under quotient
maps. R12 is a small explicit companion example of why division needs a
source/image guard. No trajectory map from this recurrence to Collatz is
supplied; the five-cylinder source-preserving return problem remains open.

## Reproduction and frontier

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/strong_target.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/strong_target.py
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/run_recurrence.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/run_recurrence.py
```

The recurrence checker compares independent endpoint DP through eleven,
all-arity masses through eight, the joint-entry tests above, and the clock
Smith powers for1<=j<=12 and0<=r<=12. A separate full finite image census
checks1<=j<=6,1<=r<=4,2<=m<=8. Normal and optimized checks use explicit
exceptions; they do not depend on assertions disabled by Python -O.

The next structural questions are whether the all-size strong-target
criterion holds whenever an extra probe is needed, and whether the R5
compiler can be extended to the actual sparse-defect targets in R9 while
retaining all joint ports and global path connectivity. The native81/27
control is replayed unchanged. No extension of A000568 or complete
size-eight class census is claimed. Meta-patterns used: recover the original
source before composing; retain the first/last boundary; test the actual
integral consumer instead of only its scalar mass.

**12:00 UTC continuation:** [DEFECT_TRANSFER D2--D4](DEFECT_TRANSFER.md)
answers the single-query compiler question for every first actual interval
complement: its backward-arc cutwidth is at most two, giving O_k(n)
arithmetic evaluation with paired components. The original interval source
has width one in reversed order. Literal target distributions and bounded
width under repeated complements remain separate obligations.
