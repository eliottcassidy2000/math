# Complete classification of ternary-ready roots

**PROVED for all finite strong tournaments.** Let

    h=H(T), s=s_T(r), q=H(T-r), F=E_T,r(2).

E counts ordered spanning path covers in two bins, first bin nonempty and
ending at r; the other bin may be empty. Call a root b-ready when F-h=b.
The inherited two-birth reset has output `h'=h+2F=3h+2b`. This term names
a precise interface, not every conceivable way to realize a ternary value.

## 1 The defect splits the classification into two cases

The [weighted-cut law and rooted defect](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md#3-the-defect-detects-indegree-and-the-legal-inverse)
give, with I the in-neighbours of r and d=|I|,

    delta = F-h-q+s = sum_(a!=b in I) W_ab(T-r),
    delta = 2 sum_{a,b subset I} b_ab >= d(d-1),
    b_ab >= H(T-{r,a,b}) >= 1.                         (7)

W_ab counts two nonempty spanning paths, first ending at a and second
starting at b. Also q>=s by deleting the initial r. Thus a1-ready root
has d=1, delta=0 and q-s=1. A2-ready root has exactly one of:

* d=1, delta=0 and q-s=2;
* d=2, delta=2 and q=s.

For d=1 and q>s the inherited inverse theorem says T-r is strong; its
root is the sole in-neighbour a, with s_(T-r)(a)=q-s. Conversely, adding
a vertex beaten only by such a to a strong parent creates the required
ready root. Therefore the complete1-ready family is the newborn image
of the unique-start family, and the d=1 part of2-ready is the newborn
image of TWO_ENTRY's family.

## 2 The exceptional two-in-neighbour case is a five-word language

Suppose d=2, with the ports oriented a->b. Equality in(7) forces
`H(T-{r,a,b})=1`, hence this residual tournament is transitive. Write its
unique order `u1->...->um`. Encode each ui by the two bits

    (A_i,B_i) = ( [ui->a], [ui->b] ).

The root r beats every ui and loses to a,b. Strongness forces m>=1 and
the last letter not00. Since q=s, every Hamilton path of U=T-r starts
outside {a,b}. In the initial SCC each vertex can start a Hamilton path
of U. Hence this initial SCC avoids a,b; being contained in a transitive
tournament it is a singleton. It must be u1, and u1 dominates a,b. Thus
the first letter is11. Conversely first11 gives q=s, and last not00
gives strongness: r reaches the chain and all chain vertices reach one
of a,b, while r also reaches both ports through u1.

We need exactly `C=W_ab(U)+W_ba(U)=2`. First11 means u1 cannot lie in the
path starting at b or a; it lies in the other path. In a cover counted by
W_ab, assign each chain vertex either before a or after b. The first
assignment is forced to be before a; the first after-b vertex must have
B=0, and the last before-a vertex must have A=1. The reverse cover
interchanges A and B.

Let G1,D1 count mixed assignments with some after-b vertex, classified by
whether the last before-a vertex has A=1 or0. Let G2,D2 be the reversed
counts. They start at zero after the first11. Appending (A,B), with old
last bits (a,b), gives the positive recurrence

    G1' = A(G1+D1)+G1+(1-B)a,
    D1' = (1-A)(G1+D1)+D1+(1-B)(1-a),
    G2' = B(G2+D2)+G2+(1-A)b,
    D2' = (1-B)(G2+D2)+D2+(1-A)(1-b),
    C   = A+B+G1+G2.                                  (8)

Both G coordinates are nondecreasing. Initial11 repetitions leave all
four coordinates zero. A00 makes G1+G2=2 and can never end in an accepted
strong word. After the first10 the state is (1,0,0,0); the only possible
additional letter is01, giving (1,1,0,1), after which every additional
letter makes G1+G2>=2. Interchange the two coordinates for the01 case.
This proves the complete language

    11^k,  11^k 10,  11^k 01,
    11^k 10 01,  11^k 01 10,             k>=1.         (9)

Here11 is a single two-bit letter, and11^k means k repeated letters.
Every word is an actual strong tournament; the rooted isomorphism class
recovers its ordered ports and unique chain, so distinct words are
distinct rooted classes.

## 3 All-size counts and actual Hamilton values

The five tails in(9) have the following values:

|Tail after11^k|Order|H|s=q|F-H|
|---|---:|---:|---:|---:|
|empty|k+3|3*2^k-1|1|2|
|10|k+4|6*2^k+1|3|2|
|01|k+4|5*2^k-1|1|2|
|10,01|k+5|13*2^k+3|5|2|
|01,10|k+5|14*2^k+3|5|2|

For an independent derivation, let W comprise a,b and the final zero,
one or two chain vertices. Prefix W by a transitive k-chain, and add r
losing just to J={a,b}. All paths starting at r number H(W). Deleting r
from any other path gives a cover ending in J, except the covers whose
second path starts in J. Thus

    H = H(W) + 2^k E_J,W(2) - C_J,W(2).                (10)

All prefix vertices must lie in the first path of a C cover, whereas
they can be assigned freely to either bin in E. Here C=2 and the five
pairs `(H(W),E_J,W(2))` are `(1,3),(3,6),(1,5),(5,13),(5,14)`.
These constant-size enumerations give the displayed formulas at every k.

Let U=x^3/D be the inherited one-start rooted series, V the two-start
series, and E the exceptional series. We have the complete ready-class
generating functions

    E = x^4(1+2x+2x^2)/(1-x),
    Ready_1 = xU = x^4/D,
    Ready_2 = xV + E.                                 (11)

All these rooted classes are rigid. This is a classification of the exact
input locus of the reset, not a proof that its H projection covers every
required integer. In particular the exceptional fifth family supplies
`31,59,115,...`; its reset gives `97,181,349,...`.

## 4 Re-marking cannot make the standard reset reusable

**Theorem.** Start with any strong T of order n>=4, rooted at r. Add x
beaten only by r, then y beaten only by x. In the resulting tournament,
no vertex has F-H equal to1 or2. This includes arbitrary re-marking of
every actual old and new vertex.

Proof. An old v!=r acquires two in-neighbours, hence has indegree>=3 and
is excluded by(7). The vertex x has gap zero by the reset identity.
The newest y has gap h-s of the old source. Inserting r after the last
in-neighbour along each path of T-r injects q paths into those not
starting at r, so h-s>=q. If q>=3 this gap is at least three. If q=1,
then s<=1 and the Hamilton cycle in T gives h>=n, so h-s>=3. The parity
of tournament path counts excludes q=2.

Finally r has indegree d_old+1. If d_old>=2 it is excluded. Otherwise
write a for its old sole in-neighbour. At the final r the two in-ports
are y->a. The residual tournament includes x, which dominates all other
residual vertices. If it is not transitive the exceptional criterion
fails. If it is transitive its first letter is x=11, followed by n-2>=2
old vertices, each having first bit0 because y beats it. Language(9)
allows at most one such final letter, namely01. Therefore this case also
fails. QED.

The n>=4 boundary is real: resetting C3 gives an order-five H9 tournament
with two2-ready roots. The theorem concerns this specific reset and
subsequent re-marking. Actual complements, controlled source changes,
new families and other reset schedules remain available. It explains
why the 31 branching obstruction needs a real change of source rather
than a more ingenious selection of a root in the same output.

## 5 The exceptional library is exactly the ready locus stable under R

There is a positive growth law before the ternary sealing step. Let R
denote one indegree-one birth retaining the old root. In every rooted
source its output gap is delta, since the new s equals the new q. Thus
among all1/2-ready strong roots, precisely the five exceptional families
remain1/2-ready after R. In fact they remain2-ready indefinitely.

This has a literal explanation: the newborn beats both in-ports and every
old chain vertex, while the old root beats it. Hence it is a new first
11 in (9). The tail stays unchanged and k increases by one. Their exact
profiles are `(a*2^k+b,b+2,b+2,a*2^k+b+2)`, with the five(a,b) pairs in
the table. This is a complete preparation library stable under R, rather
than a claim that the subsequent ternary output returns to it.

## 6 A singleton core gives a second, smaller two-lane ladder

The short-tail search in TRANSPORT used nontrivial cores of orders3..5.
Returning to the discarded singleton boundary gives another useful
construction. Let E_k be the exceptional word11^k, with physical ports
(r,b), where a->b are its two in-neighbours. It has H=3*2^k-1 and is2-ready
at r. Its converse is1-ready at the **retained** b. To see this literally,
start with a singleton a, add the all-one k-tail funnel rooted at r, then
add a vertex x beaten only by r. Converse, identify b with x, and reverse
the tail order. This gives exactly E_k with the same old names transported.

Consequently a source of count3*2^k-1 can choose both branches

    3H+2 = 9*2^k-1,        3H+4 = 9*2^k+1,           k>=1.

The five-family R growth and this converse/re-marking exchange give an
explicit two-lane preparation ladder at every k. This exchange uses global
converse; it is not being identified with one native off-spine complement.
If native realization is required, the inherited stabilized arc compiler
implements converse by at most n(n-1) native moves, n(n-1)/2 temporary
births and deletions, peak n+1. That generic bound does not assert H or
strongness at intermediate stages, or minimum fixed-order native distance.

The full golden coordinates at the selected ports are

    E_k^op at b: (h, 3^k+3*2^(k-1)-2, 3^k+3*2^(k-1)-1),
    E_k at r:   (h, 2^(k+1)-2,        2^(k+1)-1).

They follow from the singleton all-one funnel polynomials and the actual
birth formulas; E(1)=1 on both sides. Their difference in X and Y is
`3^k-2^(k-1)`. With the F49 convention of TRANSPORT, the exchange is
phase-neutral **iff k=2 mod6**, since5^k=4 mod7 exactly on that clock.
Both phases there are `(alpha,g)=(1,4)`. For k=6m+2, m>=0, this gives

    H=12*64^m-1,
    3H+2=7[36(64^m-1)/7+5],
    3H+4=36*64^m+1.

The first phase-fixed source has H11 and outputs35 and37. The quotient
and second output obey `q_(m+1)=64q_m+9`, `p_(m+1)=64p_m-63`. These
are actual source families, with no assertion of infinitely many prime
outputs. Their ternary children still face the reset obstruction above.
