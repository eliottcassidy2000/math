# A terminating scalar-fibre optimizer and a synchronized Pythagorean driver

Session256, 11 November 2026. **New proved-in-note mechanism:** a positive
nine-coordinate machine for the actual five-core grammar, and an all-height
terminating Pareto optimizer for its modular/regular fibres. **Explicit bridge:**
a ternary primitive-Pythagorean parameter driver synchronized with graph addresses.
**Not proved:** five-core coverage, a uniform target-height theorem, or Chan–Pak.

## 1. What is inherited, and what changes

Session251's fixed-core theorem is an input. For a core C, an apex cut O,
and addresses A,B >= 0, it constructs the tournament

    TT_left -> C -> TT_right, plus the apex with the encoded cut.

Its signature is (d,k+,k-,i+,i-,o+,o-,j+,j-), and its count is

    d(g(A)+g(B)-1)+(d+k+)AB+(d+k-)(A mod2)(B mod2)
      +i+A-i-(A mod2)+o+B-o-(B mod2)+j++j-.

Here g(0)=g(1)=1, g(2q+1)=g(q), g(2q)=g(q)+q+(q mod2).
The source proof is
`05-knowledge/results/session251_creative_completion_20261111/ternary/PROOF.md`,
GitHub blob `221d323e6ba678471e471ddd393884a3e623fc3c`, fetched this session.
We use the five source row masks literally. No claim is made that their scalar
H-values alone determine attachment behaviour.

Session251 identifies scalar-fibre optimization as the productive alternative
to forcing short representatives for every joint residue pair. Its digit
pullbacks have growing slopes and are not a finite exact-value automaton.
The new theorem solves an optimization problem without bounding slopes or
claiming a finite automaton for all integer targets.

## 2. The positive nine-coordinate digit representation

Define

    v(A,B)=(g(A),g(B),A,B,AB,p,q,pq,1), p=A mod2, q=B mod2.

For s,t in {0,1}, append one bit to BOTH addresses. The exact update is

    gA' = gA+(1-s)(A+p),    gB' = gB+(1-t)(B+q),
    A' = 2A+s,             B' = 2B+t,
    AB' = 4AB+2t A+2s B+st,
    p'=s, q'=t, (pq)'=st, 1'=1.                         (1)

Consequently four fixed nonnegative integral 9-by-9 matrices act on v.
The seed is v(0,0)=(1,1,0,0,0,0,0,0,1). Equation (1) follows directly
from the g recurrence and multiplication. It holds at zero as well.

Every pair (A,B) is represented by a synchronized word, by padding the
shorter binary expansion with zeroes. The graph count is a fixed linear
functional of v, with possibly signed coefficients only in the constant
and Boolean coordinates. The engine stores each of the four matrices.

This is an actual graph-address operation, not a claim of a one-vertex
operation on the already finished child. Appending bits changes the old
flanks and their apex cuts; it is legal because the grammar reconstructs
the whole tournament from C,O,A,B. Core rows and cut remain sidecars.

## 3. All-height Pareto optimization theorem

Fix a positive even modulus M and a finite deterministic controller Q on
the four digit pairs; transitions may be forbidden. The controller may
track regular address restrictions and any other finite synchronized
register. Label a reachable state by

    (controller state, v mod M).

Within each label retain just the componentwise-minimal reachable integer
vectors. If a new vector is dominated by a retained vector, discard it.
Otherwise retain it, delete vectors it dominates, and queue its four legal
successors. Initialize with the seed and the controller's initial state.
Keep a word with every retained vector.

**Theorem.** This algorithm terminates, without any input bound on address
length or graph order. At termination, for each label every reachable
vector is dominated by a retained reachable vector with that label.
For any fixed core/cut, its exact minimum Hamiltonian-path count in any
union of accepting labels is attained by one of the retained witnesses.
The same is true for the minimum across any finite core/cut bank.

**Proof of termination.** There are finitely many labels. Whenever a
vector is added to one label, no earlier retained vector dominates it.
Nor can any previously deleted vector dominate it: deletion always
replaces that vector by a smaller one, so a currently retained vector
would then also dominate the proposed addition. Thus an infinite run
would give, at some label, an infinite sequence with no earlier term <=
a later term. Dickson's lemma for N^9 forbids this. There are finitely many
queued additions, each with finitely many successors, so the queue empties.

**Proof of domination completeness.** At termination the retained union
is downward representative of the seed and closed under each legal
transition up to domination. Matrices (1) preserve componentwise order.
States with the same label have the same controller successor and the
same successor residues. Induction on the length of an arbitrary word
therefore produces a retained representative below its endpoint. Every
retained endpoint was actually reached, with a stored word.

**Proof of exact minimization.** In a fixed label the Boolean coordinates
p,q,pq are equal as INTEGERS, not merely congruent: each is zero or one
and M>=2. The constant coordinate is also equal. The nonconstant remaining
coefficients in H are d,d,d+k+,i+,o+, all nonnegative. Hence domination
within a label implies no increase in H. Every feasible graph is dominated
by a feasible retained graph; the retained graphs are themselves feasible.
Their minimum is therefore the true minimum. QED.

Only four residues need actually be stored: (gA,gB,A,B) mod M determine
all nine for even M. Within a label the Boolean coordinates are fixed,
and A'<=A,B'<=B imply A'B'<=AB. Thus the Dickson argument can equivalently
be made in dimension FOUR. Since g is always odd, there are at most
|Q| M^4/4 populated labels. This bounds the number of residue labels,
NOT the number of Pareto representatives or iterations.

For an odd requested modulus m, run the theorem at 2m and forget the extra
parity information only at acceptance. It is not valid to ignore parity
before the monotonicity argument.

### Why one representative per residue is insufficient

At M=6, the unrestricted engine retains both

    A=0,B=37,gB=15 and A=0,B=43,gB=9

in one residue label. Their words are 100101 and 101011, using pair-digit
2s+t. The first has smaller address, the second smaller g. Which is better
can depend on the coloured core's coefficient of B. A minimum-length or
minimum-address residue representative is not in general an H optimizer.
The full antichain pays exactly this missing tradeoff.

### Exact-value criterion, and the remaining obstacle

For target N, choose an even M>N and accept H == N (mod M), with any
specified finite controller. Counts are positive, so the minimum in this
fibre is N if and only if N is attained. Thus this is also an all-height
membership procedure, although Session251 already supplied a terminating
inverse by another route. The new contribution is reusable simultaneous
optimization over regular and arithmetic fibres, with Pareto compression.

A proof that the least value in each required fibre is N, uniformly in N,
would settle five-core coverage. Neither Dickson termination nor modular
surjectivity proves that inequality. There is no asserted polynomial
runtime, uniform carry bound, or finite-state language of exact H values.

## 4. A literal primitive-Pythagorean synchronization

Use Euclid parameters m>n>0, gcd(m,n)=1, m-n odd, initially (2,1). The
three updates are

    T0(m,n)=(2m-n,m),
    T1(m,n)=(2m+n,m),
    T2(m,n)=(m+2n,n).                                  (2)

Each preserves positivity/order, coprimality, and opposite parity. The
triple (m^2-n^2,2mn,m^2+n^2) therefore remains primitive. These familiar
ternary maps are used here for this preservation fact; no completeness
or unique-parent theorem for the triple tree is needed.

Synchronize T0,T1,T2 with digit pairs 00,01,10 respectively, forbidding 11.
The SAME word now carries both an actual primitive triple and an actual
tournament. Reduce the Euclid parameters modulo any ell to obtain a
finite controller of at most ell^2 states; include any desired finite
word-condition flags. The negative entry in T0 is harmless: the Euclid
register is finite, while domination is applied only to the positive
compiler vector. It would be wrong to extend the positive-order proof
to the full unreduced Euclid state without additional work.

**Corollary.** For any congruence predicate on Euclid parameters or on
their primitive triple, any regular restriction on the shared ternary
word, and any Hamiltonian-count congruence, the least H in this synchronized
five-core family is effectively attained and constructible, or the family
is detected to be empty. No address-height cutoff is needed.

This is a genuine common-driver construction. It does not assert an
arithmetic identity between H and a leg/hypotenuse, does not project
independent residue witnesses into a false joint witness, and does not
claim that this restricted synchronized family covers all odd H.

## 5. Rational series with the missing coordinate retained

The finite product of a controller and the modular compiler register
recognizes accepted digit words. Its ordinary accepted-word generating
series is rational. More strongly, attach the matrices (1) to its edges.
The sum of actual H over accepted words of length k is a coefficient of

    output * (I-z T)^(-1) * seed,

for a finite block matrix T, using a 9-dimensional block at each product
state. Output uses the fixed core/cut functional only at accepting states.
This series is rational over Q(z), even though the output functional is
signed. It counts address WORDS, including duplicate padding, rather than
isomorphism classes of tournaments or distinct H values. A canonical
padding restriction can instead be built into the controller.

The exact integer vectors and graph sidecars are needed for optimization;
the rational word series alone does not supply the minimum or coverage.

## 6. Reproduction and finite controls

Run `python3 -B engine.py`. No imported packages or external datasets.
`results.json` records the four explicit matrices, direct graph witnesses,
and all experiment statistics. The engine checks 4096 digit updates and
1440 newly implemented graph/count comparisons (90 core cuts, A,B<4).
These are implementation controls of this session's compiler, not an audit
of inherited coverage claims. Each displayed minimizing witness is also
counted directly by subset endpoint DP.

Unrestricted moduli 2,4,6,8 have respectively 4,64,324,1024 populated
labels and 4,64,361,1024 Pareto representatives. Their maximum retained
word lengths are 1,4,8,8. These are finite outcomes, not asymptotic laws.

A regular driver requires both addresses positive and an even number of
1-bits in each address. The exact least H values modulo4 across the five
cores are 37 in class1 and23 in class3. Witnesses are respectively
(A,B)=(3,5) and(3,3), both with singleton core and empty cut. The theorem,
not a scan through a fixed bit length, licenses the word “least.”

The separate Pythagorean experiment forbids the pair11 and initial pair00,
requires both addresses positive, and requires (m,n)==(1,0) mod5. It
minimizes H in each odd class modulo4. The exact least values are H=97 in class1 and H=79 in class3.
The witnesses have (A,B)=(4,10), parameters(31,20), triple(561,1240,1361),
and (A,B)=(30,1), parameters(21,10), triple(341,420,541), respectively.
The search retains4248 representatives in1249 populated labels. Full
words and graph rows are recorded in `results.json`.

## 7. Loss ledger

Preserved: one shared word, actual core and cut, both exact addresses,
exact H, reconstructible graph arcs, and in the ternary driver an actual
primitive triple. Compression discards dominated graph identities and
exact values of the Euclid register from the search STATE, but retains a
word that reconstructs one witness. It does not preserve future arbitrary
graph attachments, minimum vertex count, or minimum triple height. Those
would require extra objectives and, for nonmonotone ones, a new argument.

## 8. Generalization for genuine substitution languages

The antichain theorem only needs componentwise-monotone transition maps
N^d -> N^d and a compatible finite quotient, not linearity. Its proof
is unchanged for maps whose coordinates are polynomials with nonnegative
integer coefficients, using residues modulo M as the finite quotient.
A monotone output likewise has its exact fibre minimum at a retained
representative. Signed outputs are allowed only when the signed
coordinates are fixed within each label, as in the flank compiler.

This matters for the concurrently developed morphic Pythagorean lane.
If two letters have nonnegative transfer matrices X,Y, maintain the PAIR
(X,Y), not just the transfer matrix of the current expanded word. The
substitutions

    A: O -> OI, I -> I,
    B: O -> OI, I -> OIO,
    C: O -> O,  I -> OOI

act on this matrix pair by the respective positive-polynomial maps

    (X,Y) -> (XY,Y), (XY,XYX), (X,XXY).                 (3)

Here transfer is taken in the convention M(uv)=M(u)M(v). A column-vector
chronological convention reverses every product. For the starting word
OI, the evaluated matrix is XY. A nonnegative initial/final evaluation
of XY is monotone in all matrix entries. Thus the abstract theorem gives
a terminating exact minimum for every finite-register constraint in THIS
substitution family as soon as that lane supplies such a positive
transfer and endpoint evaluation.

Composition order is essential: (3) composes the current word morphism
on the RIGHT by the new substitution. The word receipt must reconstruct
that order. For a chronological left-composition controller, reverse the
operation word or track the appropriate right-composed incidence matrix
modulo its modulus. Do not silently use the same update for both orders.

This extension does not assert that the infinite expanded-word language
is regular. It works directly on the substitution semigroup through a
finite-dimensional positive-polynomial carrier. Section8 is an abstract
proved interface; this packet's runnable experiments use the explicit
nine-coordinate flank machine of Sections2–6, not an unprovided transfer
from the other lane.

## 9. Concrete morphic inverse: the positive carrier is now supplied

The companion `../spinor/PROOF.md` §§2–3 supplies exactly the required
nonnegative transfer. For a transitive-base/apex cut word, use coordinates
(H,P,e,1), initial vector t=(1,0,0,1)^t, and letter matrices

    U = [[1,0,0,0], [0,2,0,1], [0,0,0,1], [0,0,0,1]],
    V = [[1,1,1,0], [0,2,0,0], [0,0,0,0], [0,0,0,1]].

The graph proof there gives the actual Hamiltonian count after each letter.
In its column-vector convention M(uv)=M(v)M(u). Maintain 32 nonnegative
entries of (X,Y), initially(U,V), with transitions

    A:(X,Y)->(YX,Y), B:(X,Y)->(YX,XYX), C:(X,Y)->(X,YXX).

Read the CHRONOLOGICAL substitution address backwards. The exact H output
is `(YX)[0,0]+(YX)[0,3]`, a polynomial with nonnegative coefficients in
all32 entries. Thus every hypothesis of Section8 is now met by an actual
graph family. No conversion of a signed matrix or conjectural positivity
claim is needed. Keep reversed words during search and reverse them for
the graph/substitution/spinor receipt.

**Exact morphic inverse theorem.** For every target N>=1 there is an
algorithm that terminates and decides whether a Berggren-lift tournament
has H=N, producing an actual substitution address and graph when it does.
The same holds after any fixed regular restriction on addresses or finite
synchronized congruence-register restriction.

**Proof.** Run the Pareto algorithm on the32 coordinates with full-state
residues modulo N+1 and the chosen finite controller. Accept labels whose
output is congruent to N. The graph output is a positive integer. Its
possible values in this class are N,2N+1,3N+2,... . If no accepted endpoint
exists, no exact witness exists. Otherwise compute the least H among the
finitely many retained accepted endpoints. By monotonicity and the
all-height domination theorem it is the true minimum. The exact target
is realizable if and only if this minimum equals N. QED.

For a restriction stated on chronological addresses, use a finite
controller for its REVERSE language. A property such as containing B is
reversal-invariant. Other order-sensitive properties must be reversed
explicitly. Euclid residues can be carried by the right-composed
incidence matrix modulo a fixed modulus, evaluated on the initial count
vector; using a forward spinor update on a reversed address is incorrect.

The executable `morphic_optimizer.py` implements these matrices and the
reachable-state-only algorithm. It enumerates no q^32 universe and imposes
no depth cutoff. Dickson termination is an existence/termination theorem,
not a useful worst-case runtime bound. The packet still makes no coverage
claim, even for its unrestricted morphic family.

### Concrete outcomes of the all-height morphic optimizer

For moduli2,3,4,8,16 the respective populated residue-label counts are
6,147,19,81,377, with8,244,20,105,633 Pareto representatives. These are
reachable labels only. Modulo16 the exact least H in odd classes1 through15
are respectively

    17, 3, 5, 408071, 9, 107, 3181, 1615.

Their chronological minimizing addresses are

    AAA, empty, A, BCB, B, BA, AAAB, BB.

Thus13 and15 are absent from the morphic image at ALL depths, despite
being admissible targets for larger tournament grammars. This is a strict
boundary on the actual Berggren graph lift, not an obstruction to Chan–Pak.
Modulo8 the least H in class7 is1615. Under the regular restriction that
B occurs, the exact minima modulo4 are9 in class1 and107 in class3.

The exact-inverse examples use a power-of-two modulus strictly larger
than N. The proof works with ANY M>N, not only N+1. Targets3,5,9 are
accepted;1,2,7,13,15 are rejected in this restricted morphic family.
The saved receipt includes the actual chosen modulus, fibre minimum,
word, triple, and full graph rows. Newly displayed witnesses of at most
13 vertices receive independent subset-DP controls.

A preliminary run including modulus5 was voluntarily stopped after it
was slower than useful for these small controls; no modulus5 minimum or
failure is claimed. This illustrates why the theorem deliberately has no
practical complexity claim. Both normal and optimized final scripts
produce identical saved result files.

### FINITE-EXACT COMPLETE SEMIGROUP CERTIFICATE

The all-depth modulo16 minima and resulting missing-target statements are
**FINITE-EXACT COMPLETE SEMIGROUP CERTIFICATES**, not bounded-depth census
claims. The saved `morphic_mod16_certificate.json` contains all633 final
32-coordinate representatives and their reversed witness addresses. Its
construction exhausted the queue without a depth, size, or iteration cap:
377 reachable residue labels,633 retained representatives,763 additions.
The construction-additions number is a run statistic, not an invariant
needed by the independent closure proof.

Run `python3 -B verify_morphic_certificate.py` to check the certificate
WITHOUT rerunning the optimizer. It checks all of the following:

1. Every one of633 stored words reconstructs its exact saved32-vector.
2. The initial matrix pair is dominated by a representative with the same
   full-state residue label.
3. All1899 representative/branch transitions are dominated by a saved
   representative with the same residue label.
4. Every retained bucket is an antichain, and output minimization gives
   the stated eight exact minima.
5. The selected chronological words reconstruct the saved graph counts,
   triples, and row masks through their literal cut words.

Items1–3, positivity, and induction on operation-word length ALONE imply
that every semigroup endpoint is dominated in its label. Together with
the nonnegative polynomial output, they certify the all-depth minima.
No trust in a maximum observed depth is required. The verifier uses
explicit exceptions rather than removable Python assertions, so its
checks remain active under `python3 -O`.

In particular every positive target H==7 (mod16) below408071 is absent
from this morphic family, including7,23,39. These are NOT global tournament
exceptions; the global Chan–Pak omissions remain7 and21. The all-depth
omissions13 and15 follow from their respective fibre minima3181 and1615.
The companion verifier receipt includes chronological addresses, triples,
and all graph rows for the eight minimizing witnesses. The largest such
graph is21 vertices; no direct exponential-DP recount was performed on it.

### Incoming collision: optimization is not source-preserving repair

The concurrent `central_modular_chanpak_20261110/compiler/INCOMING_HEIGHT.md`
packet, received at main206e445 (session source copy
`session256/sources/INCOMING_compiler_INCOMING_HEIGHT.md`), supplies an actual
same-chart collision: C3 core, cut1, leftA=1, rightB=37 and39 both yield
H611, with g-values15 and5. Its carry difference is(-10,2,0,0), and either
common appended bit destroys equality. Those incoming counts are used as
supplied, not independently re-audited here.

This strengthens the reason to retain a Pareto antichain: equal-H states
need not be comparable and must not be merged on scalar equality. Our
optimizer retains incomparable endpoints unless another feasible endpoint
in the SAME finite label dominates them. It minimizes H by replacing the
address witness when legal; it does not preserve a specified source graph
or solve exact equality-repair continuations from that source. The latter
problem still needs its signed difference or another lawful carrier.
