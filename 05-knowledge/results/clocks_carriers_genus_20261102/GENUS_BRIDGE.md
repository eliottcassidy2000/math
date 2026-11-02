# Fermat rigidity, modular genus and a literal Collatz certificate

Status: CITED inputs and elementary PROVED corollaries, with exact replay
and same-agent proof audit. No independent audit or external novelty claim.

## G0. Repair the objects before comparing them

The pasted proposal confuses two curves. The prime levels with genus-zero
X0(p) are {2,3,5,7,13}. Ogg's fifteen primes instead characterize genus-zero
X0(p)^+=X0(p)/<w_p> and are

    O={2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}.

The quotient matters: X0(17) has genus 1 although X0(17)^+ has genus 0.
The complete intersection of O with Fermat primes is exactly {3,5,17};
no assumption about undiscovered Fermat primes is needed, since the next
Fermat number is 257>71. The intersection of the genus-zero prime levels
of X0 itself with Fermat primes is {3,5}.

CITED primary exposition: [Aricheta, A supersingular observation of Ogg,
Theorem 1](https://celebratio.org/Ogg_AP/article/1157/). This supplies the
fifteen-prime classification and the Fricke quotient, not a Fermat bijection.

The asserted triadic pattern is not a proof of an isomorphism. An involution
on a finite set has f fixed points and t exchanged pairs, with |X|=f+2t.
There is no universal requirement f=1,t=2. The three regular Euclidean
plane tilings have one self-dual member and one exchanged pair; the five
Platonic solids have one self-dual member and two exchanged pairs. These
are different involution-set types. The numbers count objects, not natural
labels for individual primes. Restricting to {3,5,13} changes the proposed
source set, and the displayed assignment13->17 remains arbitrary.

An explicit field obstruction makes one proposed strengthening impossible:
Q(2cos(2pi/13)) and Q(2cos(2pi/17)) have degrees6 and 8. An invertible
rational-coefficient Mobius transform preserves Q(x), so it cannot take
one of these generators to the other. This does not exclude other
correspondences between modular curves and cyclotomic data.

## G1. A genuine common predicate: rotational rigidity

The inherited audited [THM-4644](../../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md)
proves, for odd prime p,

    p is Fermat <=> every circulant tournament of order p has Aut=C_p.

Its multiplier statement is Akin's classical result, not new here. Thus
{3,5,17} is exactly the set of odd Ogg primes whose every circulant
tournament is rotation-rigid. This is an intersection of two independently
meaningful predicates, not a bijection between the entire prime families.
THM-871 already distinguishes the unique rotational classes at 3 and 5 from
the 16 classes at 17. Dropping that multiplicity would erase real structure.

## G2. The modular genus is a certified Collatz trunk number

Let p=F_m=2^(2^m)+1 be prime, m>=2, and g=g(X0(p)). Then

    g=(p-5)/12=(2^(2^m-2)-1)/3,
    3g+1=2^(2^m-2).                                  (G1)

Consequently g is a positive odd integer, its accelerated odd Collatz
successor is 1, and its standard Collatz trajectory reaches1 using one
odd step followed by 2^m-2 halvings. For g=1 this describes one return
around the trivial cycle; its first hitting time is 0.

Proof from the classical genus formula: for prime p>3,

    g=(p+1)/12 - e2/4 - e3/3,
    e2=2 if p=1 mod 4, otherwise 0;
    e3=2 if p=1 mod 3, otherwise 0.

A Fermat prime p>=5 is 5 mod 12, hence e2=2,e3=0. Substitution gives(G1).
Source: [Snowden, Lecture 12, genus formula and prime-level example](https://public.websites.umich.edu/asnowden/teaching/2013/679/L12.html).
This imports the genus formula; the following Collatz/constructibility
consequences are elementary deductions, not statements attributed to that source.

| p | genus g | 3g+1 | odd trunk index k, g=(4^k-1)/3 |
|---:|---:|---:|---:|
|17|1|4|1|
|257|21|64|3|
|65537|5461|16384|7|

The index is always k=2^(m-1)-1. For m>=3 it is odd and exceeds 1.
Therefore **the genus is not a constructible polygon order**, even though
p itself is constructible. This follows from the stronger inherited
[THM-4801](../../../01-canon/theorems/THM-4801-collatz-meets-the-constructible-polygons-trunk-five-initial-segment-principle-collatz-closed-family.md),
which classifies all constructible trunk numbers. Alternatively choose an
odd prime r|k. Any prime q|2^r-1 has order r modulo q, so r|q-1 and q
is not Fermat. Also q!=3. Thus q divides g, obstructing Gauss-Wantzel.
For example5461=43*127. Convergence and constructibility are distinct.

The prime hypothesis is essential. F5=4294967297=641*6700417 is composite.
The formal number(F5-5)/12=357913941 is a Collatz trunk number, but it is
NOT g(X0(F5)). The squarefree-level genus calculation gives

    g(X0(F5))=(642*6700418)/12-2=358472361,
    3g+1=1075417084,

which is not a power of 2. Here the index is 642*6700418, there are 4 cusps,
e2=4 and e3=0; CRT gives the latter counts. This checks the exact boundary
between a numerical substitution and the actual modular curve invariant.

## G3. The larger four-channel family

More generally for any prime p>3, g(X0(p)) is a positive odd Collatz trunk
number (4^k-1)/3, k>=1, iff p belongs to one of the four forms

    p=4^(k+1)+9, +1, +3, or -5,                       (G2)

respectively at residues1,5,7,11 modulo12, with p required prime.
This is simply the four exact genus formulas solved for p. The +1 channel
is the Fermat channel: primality forces k+1 to be a power of 2, by the
factorization of x^r+1 when r>1 is odd. The other three channels need not
be Fermat. For k1 the primes are 17,19,11; for k2 they are 73,67,59;
for k3 they are 257,251. No infinitude in any channel is claimed.

Map/predicate/loss: prime modular curve -> integer genus -> exact odd
Collatz certificate. The map retains the displayed arithmetic relation;
it forgets the curve, its modular interpretation and most prime data.
This proves convergence of this explicit image, not of all positive integers.

## G4. The two/three branching can be checked as an actual permutation map

On P1(F_p), use S(x)=-1/x and R(x)=1-1/x, with the usual infinity
conventions. S has order 2, R order 3, and R after S is x+1. The associated
modular permutation map has d=p+1 edge labels, c(S) white vertices,
c(R) black vertices and c(RS)=2 faces. The action is transitive: translations reach every
finite point and S connects infinity to 0. Thus

    2-2g=c(S)+c(R)+2-(p+1).                           (G3)

Here c(S)=(p+1+e2)/2 and c(R)=(p+1+2e3)/3, recovering the genus formula.
This supplies an independent permutation-cycle computation, with the
modular-curve interpretation inherited from the standard coset action.
At p13 the counts are(8,6,2;14), while at p17 they are(10,6,2;18):
two extra white vertices and four extra edges change genus 0 to genus 1.
The relevant small structure is the actual order 2/order 3 permutation
pair and its product, not an arbitrary assignment of primes to solids.
