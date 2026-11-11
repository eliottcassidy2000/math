# Chan–Pak through primitive triples, positive compilers, and finite response geometry

**Session256, 11 November 2026. COMPLETE research deliverable.** New
proved-in-note constructions and finite complete certificates, with historical
inputs retained at their stated scope. **Chan–Pak remains OPEN.** No broad
re-audit of earlier claims or literature-priority claim was undertaken.

## 1. A real Berggren-tree lift into tournaments

The three primitive-Pythagorean branches now act on actual decorated graphs.
Write the Euclidean parameters as x=m−n, y=n. Starting at the binary cut
word OI, use

    A: O→OI, I→I;
    B: O→OI, I→OIO;
    C: O→O,  I→OOI.

The word specifies a tournament: its base is transitive in word order,
and a distinguished apex beats the O vertices and loses to the I vertices.
The letter counts are exactly (#O,#I)=(m−n,n), so the decorated tournament
recovers the primitive triple (m²−n²,2mn,m²+n²). These are literal graph
operations, not merely three named arrows on integers. The base order and
apex are essential parts of the state.

A new four-coordinate positive transfer computes its actual Hamiltonian
count. With H the count, P the binary word (O=1,I=0), and e its last bit,

    append O: (H,P,e) → (H,2P+1,1),
    append I: (H,P,e) → (H+P+e,2P,0).

A pair of 4-by-4 matrices then evaluates a whole substitution address without
expanding the exponentially long cut. This gives regular congruence
languages of Berggren addresses for every fixed modulus, with the reversed
composition order stated explicitly.

The lift is not a scalar-H tree. Addresses B,C,AA all produce H=9, but their
A-children have H=107,73,17 respectively. Infinite collisions also occur:
H(C^k)=H(A^(2k)). Thus the needed sidecar is visible in actual small graphs.

[Graph construction, proof, compressed compiler and controls](spinor/PROOF.md).

## 2. Sorting the same spinor loses a whole congruence class

A second, simpler graph lift sorts the cut to O^x I^y. Its count is

    H = 2+(2^x−1)(2^y−1),

so it is exactly the safe-inflation family C3[TT_x,TT_y,K1] from Session251.
In the positive coordinates p=2^x−1, q=2^y−1, all three branches are
polynomials with nonnegative coefficients. Its complete target inverse is
factorization of N−2 into Mersenne factors, with gcd(x,y)=1 and x odd.

But x odd forces p≡1 mod3, so the entire sorted primitive-spinor image has
H≡0 or2 mod3. It never attains class1. The unsorted morphic CA graph has
H=73, whereas sorting that same spinor gives107. Retaining cut order really
changes attainable congruence classes; it is not redundant decoration.

This recovers the old Mersenne arithmetic as a quotient of the new graph
carrier and explains precisely what the quotient loses. Both families
remain restricted and do not cover the Chan–Pak spectrum.

[Sorted lift, exact factorization inverse and congruence boundary](spinor/PROOF.md#8-sorting-the-cut-gives-a-smaller-polynomial-lift-and-a-permanent-obstruction).

## 3. A terminating all-height optimizer, then an exact morphic inverse

Session251's five-core formula depends on two flank addresses A,B and a
finite coloured core signature. The new positive state is

    (g(A),g(B),A,B,AB,A mod2,B mod2,(A mod2)(B mod2),1).

Appending a pair of binary digits acts by one of four nonnegative integer
matrices. For a fixed modulus and any finite common driver, keep the
componentwise-minimal reachable states separately in each residue/controller
label. Dickson's lemma proves this Pareto process terminates at all heights.
Within a label the signed Boolean terms of the count are fixed, so the
retained bank contains a true minimum H, with an actual graph witness.

The process must retain tradeoffs. At modulus6, addresses37 and43 have the
same relevant residues, but one has smaller address and the other smaller
g. Keeping only the shortest representative would lose the optimization
certificate. This is why the output is an antichain rather than one residue
representative.

The same proof works for nonnegative polynomial transitions. Applying it
to the positive matrix-pair substitution carrier from Section1 yields a
terminating exact inverse for the morphic Berggren family. For target N,
track the full matrix state modulo N+1, accept H≡N, and minimize. The minimum
is N if and only if an exact graph witness exists. This requires neither a
depth cutoff nor the assertion that expanded morphic words form a regular
language. No practical worst-case runtime bound is claimed.

**Finite complete outcome:** at modulus16, the minimum H in residues
1,3,5,7,9,11,13,15 is respectively

    17, 3, 5, 408071, 9, 107, 3181, 1615.

The minimizing chronological addresses are AAA, empty, A, BCB, B, BA,
AAAB, BB. Thus this restricted graph lift misses13 and15 at **every depth**.
These are finite complete semigroup certificates, not an extrapolation from
a bounded-depth scan. They do not exclude13 or15 from general tournaments.

The synchronized five-core driver is a different construction: the SAME
ternary word updates a primitive triple and both flank addresses. Under the
specified modulo5 spinor and positive-address restrictions, its exact
modulo4 minima are97 and79, with triples (561,1240,1361) and (341,420,541).
This is a concrete shared-driver result, not independent residue witnesses
mistaken for one common address.

[Positive machine, termination, witnesses and morphic inverse](compiler/PROOF.md).

## 4. Two finite response interfaces replace scalar count labels

Two new exact interfaces give a practical meaning to “keep the sidecar.”

**Five-core attachment interface.** Although Session251 lists nine numerical
core coordinates, its full two-flank response sees only eight combinations:
the two j-coordinates occur through their sum. Eight explicit flank probes
recover all eight combinations. A ninth probe gives a consistency check.
Agreement at the eight probes is equivalent to agreement at every pair of
flank addresses, within that response formula. Minimality is proved for the
formal eight-dimensional response space, not for the set of realizable cores.

**Cyclic substitution interface.** For a fixed c-vertex observer U,

    K_U(S,T)=H(C3[S,T,U])

is a positive Gram kernel of sharp rank c+1. Its coordinates are the first
c+1 ordered path-cover numbers of S and T. The c+1 actual probes
C3[S,TT_a,U], a=1,...,c+1, recover that entire truncated profile. The
Cauchy–Schwarz inequality has exact equality precisely when the profiles
are proportional.

These are different interfaces for different graph operations. Neither can
be silently used as the other, and neither says every formal vector is
realizable. Together they provide finite tests for contextual equivalence
before attempting a count-preserving source switch.

[Eight-probe inverse](history/README.md#2-new-transfer-a-a-nine-probe-response-complete-atlas-for-five-core-constructions),
[positive kernel and sharp observer rank](gram/PROOF.md).

## 5. The cyclic kernel adds controlled packing and the right norm

The Gram formula is a new regrouping of the inherited THM-3121 quotient-walk
kernel. For U=K1 it becomes

    K(S,T)=2H(S)H(T)+(H(S)+f_S(2))(H(T)+f_T(2)),

where f_S(2) is twice the number of spanning two-path covers. The diagonal
is an x²+2y² norm, and the two-source Gram determinant is

    2(H(S)f_T(2)−H(T)f_S(2))².

Thus the correct arithmetic comparison is initially Q(sqrt(−2)), not the
ordinary Gaussian/Pythagorean norm. The quadratic form is another essential
sidecar. An integer square Gram area cannot be manufactured by forgetting
that factor2.

For a fixed observer of order c, the transitive-source count sequence has
the common recurrence with characteristic polynomial

    product_(j=1)^(c+1) (t−j).

The common order c+1 is sharp. If all blocks grow, the observer rank grows,
so this does not contradict the inherited moving-jet non-C-finite theorem.

In particular, U=TT2 gives an explicit all-size formula for
C3[TT_a,TT_b,TT2]. This family has cycle-packing number at most2 and includes
the H45 graph that defeats an illegal intercyclic formula. It implements
Session251's proposed departure from packing number1 while keeping every
extra cycle contribution.

[Factorization, recurrence, packing-two formula and scope](gram/PROOF.md).

## 6. What was recovered from the older repository

The [inheritance atlas](history/README.md) connects22 genuinely distinct
mechanisms, not22 versions of the same tree analogy. It includes:

- Gaussian square/Farey/Lorentz geometry and equal-hypotenuse factor torsors;
- primitive-content cocycles and source-dependent ancestry repair groupoids;
- Pell/Markov selectors, positive two-cube rays and an elliptic-rank carrier;
- sibling determinant circuits, Walsh sectors and Fibonacci owner transport;
- factorial simplex moments and arithmetic decimation;
- Grassmannian carry banks, Rule30 transverse jets and finite observers;
- Stern depth/character packets, harmonic languages and Keller norm−8 renewal.

Each entry names the map, invariant, discarded data, restoring sidecar and
an actual potential consumer. Historical proofs are imported at their
recorded scope. The session's new complete interfaces and graph lift are
linked separately from still-open repair and descent proposals.

## 7. The new incoming main and the remaining constructive gate

Concurrent main `206e445797338638604337b8fc1f581f30be0e00` contributes an
[exact golden content monoid and odd carry lifts](../central_modular_chanpak_20261110/README.md).
It classifies universal primitive content multipliers as5^b9^k21^c in its
specified grammar, and gives nonnegative carry matrices for every fixed odd
multiplier. These are inherited concurrent results, not claims of this session.

Its actual same-chart H611 collision has right addresses37 and39 and carry
difference(-10,2,0,0). That incomparable pair is fully compatible with our
Pareto optimizer: optimization keeps an antichain and never assumes all
same-H sources are comparable. A minimum-value compiler is still not an
exact source-preserving equality-repair compiler.

The incoming bounded-suffix obstruction also matters. A uniform ternary
return implemented only by bounded suffix/cylinder rules forces the deletion
core count to scale by3; a fixed bounded core bank cannot repeat that
forever in that policy class. Our graph lifts and optimizer avoid making
that assertion. A successful five-core coverage proof must instead exploit
unbounded carry conditions, a non-dyadic rebase, more flexible source growth,
or a genuinely different descent.

A second incoming checkpoint, `0f674ffece36a42512114311f6d4d962e86b3df5`,
adds an [endpoint-symmetry extension obstruction](../modular_centers_collatz_pell_20261111/TOURNAMENTS.md).
If every set of at most2k old vertices has a pointwise-fixing subgroup of
order d, every k-vertex extension count is divisible by d, even when the
child is asymmetric. It excludes all one-vertex ternary branches from
nonprime Paley bases and bounded-size extensions of recursive prime-block
bases in its stated range. This is an inherited source-level restriction,
not a scalar-H restriction. A Pythagorean address label or a child whose
symmetry is broken cannot remove it. Our ternary *operation tree* is not
that recursive ternary *substitution tournament* family; the physical
source and named endpoints distinguish the two constructions.

The strongest next questions are in [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md).
The preferred target remains a reusable positive chart switch attaining
3H+2 or3H+4, with source and full response retained, or an all-target proof
for the five-core binary images. Neither is established here.

## 8. Reproduction and publication scope

Run `python3 -B reproduce_all.py` in this directory. It replays the new
engines and the finite morphic certificate. The saved receipts specify their
universes and include direct graph controls. Unbounded claims are supported
by the displayed proofs or a complete finite closure certificate, not by
large numerical coverage scans.

Only this packet and small navigation links are published. The dirty legacy
checkout was not modified. Checkpoints use the current remote tree and a
matching expected main head; concurrent work is preserved, without force
push. See SOURCE_PINS.json, MANIFEST.json and REPLAY.json for provenance and
artifact integrity. No theorem namespace was reserved or promoted to canon.
