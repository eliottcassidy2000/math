# Session265: the converse shadow has geometry, but arithmetic needs chirality

**11 November2026. Final mathematical synthesis. Chan–Pak and Collatz remain OPEN.**

This session merges the user's weighted-path construction with7/21/63,
powers of two, Cayley–Dickson/Fano presentations, Collatz carries, and the
physical source and root-sharing operations from Session264. It concentrates
on new constructions and proofs, accepting old proved results as inputs.
Inherited facts, new proofs, exact finite inputs and remaining questions are
distinguished in the linked lane documents. No external-priority claim or
canon promotion is made.

## 1. The weighted graph recovers every strong component

Let W(u,v) count Hamiltonian paths using u,v consecutively, and let H be the
total path count. Path reversal gives W(T)=W(T*), exactly as suggested.

The new reconstruction theorem is stronger than a visual clustering claim:

    sum_(u in S,v outside S) W(u,v) = H
      iff S is an SCC prefix or suffix.

Every other nontrivial cut has strictly greater weight. Hence the undirected
weighted graph determines the entire SCC partition and its ordered chain up
to reversal. Pairwise minimum-cut queries suffice to recover the partition.
For strong T on n vertices, every cut has the quantitative lower bound
H+(n-2)HC(T), where HC counts directed Hamilton cycles modulo rotation.

Between consecutive SCCs C_i,C_(i+1), the connector is the exact rank-one matrix

    W(u,v)=H·[end_i(u)/H(C_i)]·[start_(i+1)(v)/H(C_(i+1))].

Its TOTAL weight is H. The individual edges need not be uniform or equal to H.
Nonconsecutive SCC connectors vanish. If there are at least two SCCs, boundary
margins and internal weighted degrees recover both normalized endpoint vectors
of every SCC. A single ordinary source observer orients the recovered chain.

This is an actual sufficient interface for ordinal composition, not a claim
that W reconstructs every internal arc. Two explicit strong six-vertex rooted
tournaments have the SAME labelled W and H=45, but the SAME selected-apex
context produces93 versus111. The missing quantity is precisely the existing
rooted response U-d, not a new ad hoc invariant.

[Tight cuts, interfaces and explicit shadow twins](shadows/PROOF.md).

## 2. Seven Fano controls give21 arcs and63 folded character directions

The seven imaginary labels in an eight-dimensional Cayley–Dickson/Fano
presentation organize seven cyclic lines and21 pair locations. Reversing
selected lines gives a literal cube F_2^7 of128 tournaments. Converse adds
the all-one word, so its quotient has64 points. Its function space has exactly
63 nonconstant even Fourier characters. These are precise, DIFFERENT roles
for7,21,63, connected by explicit maps rather than numerical identification.

Hamiltonian counts and every labelled W entry lie in the converse-even sector.
The rooted start-minus-end response lies in the odd sector, which folding
erases. An explicit nonalternative Fano line table A and its converse B have
H=171, identical labelled rooted W, and identical complete unmarked and
root-deleted cover profiles, but

    P_A(w)-P_B(w)=2w²+w⁷.

This is an unusually small chiral response packet: total mass differs by3,
odd-minus-even mass by-1. Standard octonion classification and pending
higher-dimensional claims are not needed for this finite input certificate.

Root-sharing fusion consumes that hidden response. Swapping AB to BA changes
H with the SAME strict sign in every physical left/right context, including
contexts with nonzero parity channels. A fixed inventory {A,B,C3,C3} gives
17-vertex outputs with adjacent-swap gaps5902 and12760, whose gcd is2. Thus
this particular strong-module inventory has no forced congruence beyond
parity. This is not an interval theorem, and it is not the first parity-only
inventory in the unrestricted fusion grammar.

[Fano certificate, Fourier map and positive sorting theorem](doubling/PROOF.md).

## 3. The hidden packet has a literal dyadic amplifier and a step-of-two reader

Appending an m-vertex transitive tail multiplies a rooted path polynomial by
(1+w)^m. The inherited endpoint-preserving wrap turns its coefficients into
actual Hamiltonian paths counted by a gate's position. Applied to A,B, this
builds pairs on9+m vertices with identical complete unmarked and root-deleted
cover profiles, but root-start counts

    260·2^m and257·2^m.

Their difference is coefficientwise positive, with exact gate polynomial
w(2w²+w⁷)(1+w)^m. One tail kills the parity difference while doubling the total
response. This is a physical channel amplifier, not a scalar projector.

A more direct arithmetic reader uses the same two-sided context on A and B:
a root-source two-vertex module on the left and a rooted triangle on the right.
It gives two explicit STRONG ten-vertex tournaments with

    H=1771 and H=1769.

The fusion law proves the exact gap2 by (3·3-1)-(3·2-0)=2. Both literal graphs
and an independent path enumeration are supplied. A one-sided fusion reader
cannot extract a nonzero gap below4 from this packet, so the two-sided
construction genuinely uses both retained response channels.

The values1769 and1771 were already attainable; the new result is their common
context and its invisible input difference. This is one physical unit step,
NOT a reusable +2 operation at an arbitrary count.

[Cross-lane proofs, explicit graphs and controls](BRIDGES.md).

## 4. A single positive operation is a complete binary slope controller

The physical F5 operation acts on the exact source state (r,M,q,E), flips the
unit offset r, and adds four vertices. Its pair iterate G=F5² has characteristic
polynomial

    (t-1)²(t²-745t+64).

One2-adic direction contracts with valuation6. On the stable image, G is
identity modulo2. A small explicit starting calculation and an all-level
squaring argument prove

    v2(M(n)-M(m))=1+v2(n-m)

on the projected orbit. Thus pair-iteration time and the odd slope are
isometric binary coordinates: the projected dynamics is an adding machine.

The projection is only a proof device. After ceil(k/6) ACTUAL pair iterates,
every block of2^(k-1) successive actual iterates realizes every odd slope
modulo2^k exactly once, at either chosen exact offset sign. A bit compiler
finds the required positive run. Its vertex cost is O(2^k), and the integer
counts can be enormous. For example the target7 modulo64 can be reached
at actual slope1,570,759, illustrating the ordinary-height gap.

The same investigation lifts full-state prime control to every p^k for odd
p except5,761. Explicit clocks protecting the other prime factor yield full
control over EVERY3^a7^b. A complete18,522-state modulus21 base plus an
algebraic digit-lifting proof supplies that mixed tower; modulus63 was also
checked over all500,094 source states.

[All-order proofs, bit compiler and exact controls](collatz/PROOF.md).

## 5. Collatz carry is a controlled readout, not a tournament operation

Closing a unit source gives H=r+2^tM. The dyadic controller prescribes the
exact valuation v2(3H+1). Arbitrarily deep halving occurs only on the two
critical faces H=2M-1 and H=1+4M. Every desired halving depth can be realized
inside this unary physical family. Scalar existence of such residues was
already implicit in Session199's broader modular-density result; the new gain
is this one-letter source realization and its explicit dynamics.

The numerical ladder is separated carefully:

    7=2³-1, 21=(2⁶-1)/3, 63=2⁶-1, and3·21+1=64.

Ordinal joining with a triangle implements tripling an EXISTING count.
It cannot start from7 or21, whose exclusions are inherited. The realizable63
does not need those unrealizable intermediate counts. Indeed every H63
tournament has exactly one nontrivial SCC, with count63: every proper odd
factorization would otherwise require7 or21. The weighted-cut theorem locates
that core. Thus63 is irreducible in the realizable-count monoid while composite
as an integer. Its actual realization is inherited, not discovered here.

Applying the Collatz numerical map to one of our counts does not automatically
produce another tournament. Convergence and height-controlled target selection
remain completely separate obligations.

## 6. An electrical-network bridge and the remaining frontier

The same W is the SUM of the Laplacians of all actual Hamiltonian paths. This
gives a literal tournament-to-network map beyond the old spanning-tree analogy:

- lambda_2(L_W)>=2H(1-cos(pi/n)), with equality exactly for transitive T;
- the weighted spanning-tree count satisfies tau(W)>=H^(n-1), with the same
  rigid equality case;
- effective resistance R_W(u,v)<=sum_P distance_P(u,v)/H².

These are exact elementary consequences, not an identification tau=H. A
scaled path matrix passes the basic linear cut/degree tests but fails the
equality rigidity, exposing a concrete inverse-shadow obstruction.
[Network proofs and hostile](geometry/PROOF.md).

The main new constructive objects are now composable: a recoverable SCC
skeleton, a Fano-hidden two-length packet, a positive sorting gate, a physical
unit-step reader, and a unary binary source controller. What is still missing
is a small-height exact target constructor or an unconditioned scalar channel
selector. Primitive gcd, residue surjectivity and one nearby pair do not fill
the odd integers. [Next obligations](OPEN_QUESTIONS.md), [sources](SOURCES.md),
and [reproduction](REPRODUCE.md) preserve that boundary.
