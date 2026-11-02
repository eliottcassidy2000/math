# Collision partitions and marked-port recurrences

**PROVED elementary identities; FINITE-EXACT controls; no literature-priority
claim.** Classical Burnside enumeration and collision-free sieve recurrences
are recovered and credited below. The extension here is an explicit
collision-block transition valid before the usual collision-free threshold,
and its comparison with marked tournament counting.

## 1. Exact wheel transitions at every new prime

Let Q>=2, p a prime not dividing Q, and let w=(g_1,...,g_j) be an occurrence
of consecutive gaps in the cyclic reduced-residue system modulo Q. Fix its
two endpoint ports; write s_0=0,s_i=g_1+...+g_i. For a fixed total g=s_j,
allow cyclic traversal as necessary. Partition the indices0..j by equality
of s_i modulo p; call this partition Pi_p(w). For an internal block B,
w/B is the word obtained by deleting all boundaries indexed by B and
merging the adjoining gaps. Then the exact endpoint-preserving operator is

    T_p[w]=(p-|Pi_p(w)|)[w]
              +sum_(B in Pi_p(w), 0,j not in B)[w/B].       (W)

**Proof.** An occurrence starting at a lifts to p copies starting at a+tQ.
The boundary i is removed exactly in copy t=-(a+s_i)Q^(-1) mod p.
Thus two boundaries disappear in the same copy iff they are in the same
block of Pi. There are p-|Pi| untouched copies. A block containing either
endpoint destroys that pair of ports; every other block contracts exactly
as stated. Every new occurrence of total g has a unique pair of old
endpoints and hence a unique old driving word. Summing(W) over all old
occurrences is therefore exact, with no multiplicity ambiguity. QED.

This rule needs no g<2p hypothesis. It preserves ordered gaps and both
ports. Counting forgets the starting phase a; recovering exact locations
requires a mod p and Q^(-1) mod p. A count-level recurrence is not a
location-level controller.

The exceptional primes for one word divide the nonzero integer

    discriminant(w)=product_(0<=i<j<=length(w))(s_j-s_i).

Outside these primes all blocks are singletons. The full partition, rather
than merely the number of its blocks, specifies the contractions. Common
endpoint residues also matter: a block containing both endpoints loses one
copy, not two.

## 2. Generating functions after collisions stop

Suppose Q is even and g<2p. All differences of offsets are even, so none
can be a nonzero multiple of the odd prime p within this range. Thus Pi
is discrete. If n_j counts driving words of length j and total g, (W) gives

    n'_j=(p-j-1)n_j+j n_(j+1).                           (L)

For R_g(z)=sum_(j>=1)n_j z^(j-1), equivalently

    R_g,new(z)=(p-2)R_g,old(z)+(1-z)dR_g,old/dz.           (G)

In the basis(z-1)^r, the operator is diagonal with eigenvalue p-2-r.
Consequently, along any sequence of new primes where the collision-free
hypothesis holds,

    R_final(z)=sum_r [R_initial^(r)(1)/r!](z-1)^r
                         product_p(p-2-r).              (D)

These are the standard driving-term mechanism and constant Pascal
eigenvectors of Holt and Rudd, rederived here from(W), not claimed as new.
Their source studies sieve candidates, not a proof of prime-gap conjectures:
[Eratosthenes sieve and the gaps between primes, sections3--4](https://arxiv.org/pdf/1408.6002).

The guard is necessary. At Q30, p7, g14 the old length counts are
n_3=1,n_4=2. Formula(L), applied illegally, predicts(n'_2,n'_3,n'_4)=(2,9,4).
Literal wheel enumeration and(W) give(2,10,6). The endpoints0 and14
collide modulo7. This is an exact minimal-stage hostile in our tested
sequence of primorial stages, not a claim of global minimality over all Q.

Nor does length alone determine a single word's transition before that
threshold. At Q30,p7, the two actual words(2,4,2,4,6) and(4,2,4,6,2)
both have sum18 and length5. Their surviving length distributions are
respectively{4:4,5:1} and{4:3,5:2}. Their residue collision partitions
differ. This is a direct witness against the unguarded length quotient.

## 3. Pair counts, primorial gates, and the champion4

For a squarefree Q and a fixed displacement d, CRT gives the pair count

    A_Q(d)=product_(p|Q)(p-nu_p({0,d})),
    A_(pQ)(d)=(p-1)A_Q(d) if p|d, else(p-2)A_Q(d).       (P)

More generally a finite set of offsets S has exactly
product_(p|Q)(p-|S mod p|) simultaneous unit translates. This counts
pairs or patterns, not necessarily consecutive survivors. The extra
internal boundary information in(W) is precisely what(P) omits.

For Q210, exact pair counts are

|d|2|4|10|30|14|42|54|162|
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|A_Q(d)|15|15|20|40|18|36|30|30|

Thus the pairs5/15 and7/21, after doubling to even displacements, gain a
factor2 when prime3 becomes a divisor. The pair27/81 has the same prime
support and no such change. This is an explicit map and preserved
predicate; it does not identify this factor2 with a tournament automorphism
ratio3. The81/27 tournament control still has different attachment counts.

At every primorial Q divisible by6, the cyclic gap counts for2 and4 are
equal: both start at1 for Q6 and multiply by p-2 at each later prime.
A driving word(2,2) cannot occur, because one of its three offsets is
divisible by3. Nevertheless, direct enumeration of primes through1000 gives
unique champion4 precisely on[131,139) and[491,547), with ties at139 and547.
The interval lengths8 and56 are finite prefix-clock facts; they do not
give a bijection with the Fano56 objects or the56 order-six classes.

At131, gap counts for2,4,6 are10,11,7; at139 they are11,11,8.
At491 they are24,25,24; at547 they are25,26,26. These exact records explain
the ties and demonstrate why cyclic density data cannot recover a prefix
champion. Goldston--Ledoan's eventual primorial theorem assumes specified
Hardy--Littlewood pair and triple conjectures; it is **CONDITIONAL**, not a
resolution here: [primary paper, theorem and Conjectures A/B](https://arxiv.org/pdf/1102.4879).

The old LRC gate has a different but explicit common arithmetic input:
[THM-1286, D-graded gate tower](../../../01-canon/theorems/THM-1286-D-graded-gate-tower-primorial-cascade.md),
completed by [Session114's Farey proof](../session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md).
For C_q={1,...,q-1,q+1}, s(q) is the least prime not dividing q and

    M(C_q union {qm})=max{1/(q+2),m/(qm+s(q))}.

The family is tight iff every prime<=2m divides q (m>=2). Primorials enter
through a least-missing-prime gate. Both this and(P) depend on prime support;
their targets differ. LRC additionally retains the binding owner and its
inverse modulo q, while(W) needs ordered internal boundaries for gaps.

## 4. Tournament cycle interactions and an exact marked generating function

For an odd-part partition lambda, let m_j be its multiplicities, ell(lambda)
its number of parts, z_lambda=product j^(m_j)m_j!, and

    e(lambda)=sum_i(lambda_i-1)/2+sum_(i<j)gcd(lambda_i,lambda_j).

Classical Burnside gives A000568 as a_n=sum_(lambda odd,|lambda|=n)
2^e/z_lambda. Even permutation cycles force an orientation reversal on an
edge orbit, so contribute zero. The gcd term counts edge orbits between
distinct odd cycles; it cannot be discarded as if cycles were independent.
This formula and its packed-response acceleration are inherited from
[the degree-local enumeration packet](../degree_local_tournament_enumeration_20261021/session22-enumeration/README.md).

Keeping cycle markers, define the formal series

    Z(x;y)=sum_(lambda odd) 2^e/z_lambda x^|lambda| product_j y_j^m_j.

Marking and removing one cycle gives the exact functional equation

    dZ/dx=sum_(i odd)2^((i-1)/2)x^(i-1)y_i
                  Z(x;(2^gcd(i,j)y_j)_(j odd)).          (T)

**Proof.** For a partition containing a distinguished i-cycle, removal
reduces e by(i-1)/2+sum_j m_j gcd(i,j), including its interactions with all
remaining cycles. The denominator changes by i times its old multiplicity.
The derivative's factor |lambda| is the sum of i times multiplicity over
all marked choices. The factors cancel exactly to give(T). Constant term
Z(0;y)=1 determines coefficients recursively. This is a coupled functional
recurrence; it is not a scalar Euler product or a bounded-state classifier.

For r ordered, individually fixed port vertices with their internal
tournament prescribed, let b_(m,r) count extensions by m unlabelled body
vertices, up to permutations of the body alone. Then

    b_(m,r)=[x^m]Z(x;(2^r)_(j odd)).                    (B)

Indeed every body cycle has one free orientation bit towards each port,
so a body permutation of type lambda fixes2^(e+r ell) extensions. The
port tournament can be any prescribed orientation: the same count holds.
If its orientations are also allowed to vary, multiply by2^binomial(r,2).
Port permutation is NOT allowed; quotienting it would need another action.

For m=0..6 the one-port counts are1,2,4,12,48,296,3040; the two-port counts
with their edge prescribed are1,4,16,88,704,8912,188288. recurrences.json
contains these counts through body order55 and A000568 through total
order55. The unmerged integer Burnside path agrees with the inherited
packed-response DP for every coefficient0..55.

These are exact extension-orbit counts, not the metagraph's edge adjacency.
The literal81/27 pair still shares67 labelled cyclic triples but has1620
versus1728 regular binary attachments. Formula(T) counts isomorphism
classes; it does not recreate which representative, stabilizer, marked
path, or legal re-marking realizes a transition. Those ports and actions
remain necessary for the full metagraph.
