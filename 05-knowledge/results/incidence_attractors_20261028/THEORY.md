# Maps, exact finite carriers, and the proposed common theory

## A1. The six-label cycle and two prime triples share a graphic chart

Let sigma(n) be the sum of all positive divisors. For a positive integer x,
let Q(x)=(v_p(x) mod2)_p. For pairwise coprime atoms b_i,

    Q(sigma(product b_i)) = sum Q(sigma(b_i)) over F2.             (1)

This standard parity mechanism is already used in Session82/83 and in the
classical power-values-of-divisor-sums literature. The contribution here is
the explicit chart, its forbidden selections, and the Fano extension.

For the bank B={4,2,5,3,7,13,11}, take parity coordinates(2,3,7):

| atom | sigma | squarefree kernel of sigma | graph edge |
|---:|---:|---:|---|
|2|3|3|0--3|
|3|4|1|loop at0|
|4|7|7|0--7|
|5|6|6|2--3|
|7|8|2|0--2|
|11|12|3|second0--3 edge|
|13|14|14|2--7|

Here0 is a formal root vertex, not a prime or an arithmetic atom. Deleting
its row from the graph incidence matrix gives exactly the three parity
coordinates. A selected subset has zero syndrome iff its subgraph is even
at every vertex, with a loop contributing degree two. The reduced matrix
has rank3 because the columns2,7,4 are independent, hence its kernel has
2^(7-3)=16 elements.

The minimal nonempty zero subsets are

    {3}, {2,11}, {2,5,7}, {5,7,11}, {4,7,13},
    {2,4,5,13}, {4,5,11,13}.                                  (2)

The graph has two triangles sharing edge7, the parallel pair2/11, and the
loop3. This is a genuine incidence-preserving map, not an analogy based on
the cardinality of a triple. In particular `{2,3,11}` is the union of the
loop and parallel-pair circuits; `{2,5,7}` is one circuit.

The only coprimality conflict between distinct bank atoms is{2,4}. Thus96
of the128 subsets are legal. Four kernel elements contain both2 and4, leaving
the following12 legal products, including the empty product:

    1,3,22,66,70,210,364,385,1092,1155,2860,8580.                 (3)

The two legal circuits A={2,5,7}, C={4,7,13} have an illegal symmetric
difference{2,4,5,13}. Its factorwise total is
3*7*6*14=1764=42², but its actual product is520=2³*5*13, with sigma1260.
Equation(1) never licensed that noncoprime multiplication. This is the exact
failure of closure, and no additional sign coordinate repairs a wrong domain.

For arbitrary finite atom banks the same construction gives a binary matrix
restricted by pairwise-coprimality constraints. It need not be graphic, as A3
proves. Nor are all graphs/matroids claimed realizable by arbitrary choices
of prime sigma values: primality of any proposed realizing labels is a guard.

## A2. Exact contextual quotients: 86,87,96

**Universe and language.** States are the96 products of legal subsets of B.
For each named b in B, request multiplication by b. It succeeds iff gcd(n,b)=1;
on success update to nb, on rejection keep n unchanged. Observe every success
or rejection and whether sigma(n) is square. The initial square observation
is included. These rules forbid repeats and conflicting atoms automatically.

**FINITE-EXACT theorem.** The coarsest recursively exact deterministic state
has86 possible values. Its only nonsingleton classes are the ten pairs

    (154,308), (462,924), (770,1540), (2310,4620), (910,1820),
    (2730,5460), (2002,4004), (6006,12012),
    (10010,20020), (30030,60060).                               (4)

Each pair changes the chosen2 atom to4 and retains its other atoms. Initial
guard answers recover prime support: asking any bank prime detects its
presence, with2/4 having the same support. Thus these are the only possible
nontrivial identifications. For a fixed support containing2, let q be the
current sigma syndrome and R the set of unused prime-atom syndromes. A legal
future can reach square sigma iff q is in span(R). Distinct current syndromes
are distinguishable whenever either lies in that span: choose a subset
representing that syndrome. If both lie outside, neither can ever reach zero,
and all guards agree. These are exactly the ten pairs in(4).

The producer independently refines literal transition/output signatures to
a stable86-class partition. The audit evaluates all96*128 subset contexts
using additive divisor sums, obtaining the same partition, and checks update
closure explicitly. Any successful word reduces to an allowed subset, while
rejected letters leave the state unchanged and are already observed. Hence
this finite table covers arbitrarily long words in this particular language.
Every pair in different rows has an explicit distinguishing context; every
merged pair remains merged under the updates. This proves minimality in the
stated finite model, rather than merely showing one sufficient encoding.

Adding the **Boolean** observer F=U+S yields87 classes. Only(154,308) splits:
154 is squarefree with(F,U,S)=(6,3,6);308=2²*7*11 has(10,3,7).
The other merged pairs already have at least four prime factors, and future
legal multiplication cannot reduce that support, so balance stays false.
Adding tau instead splits every pair and yields96 classes. The all-integer
numeric-defect observer in Session84 is a different, stronger experiment.

These numbers are not vector-space ranks. They also do not assert a minimal
carrier for unnamed tournament complement, arbitrary arithmetic multiplication,
or the full divisor-height graph.

## A3. A Fano/Hamming arithmetic extension, with its non-graphic obstruction

Use coordinates(2,3,5), numbered by the nonzero bit vectors1,...,7. The
following distinct primes realize every nonzero vector:

| vector | chosen prime | sigma | squarefree kernel |
|---:|---:|---:|---:|
|1|7|8|2|
|2|2|3|3|
|3|5|6|6|
|4|19|20|5|
|5|89|90|10|
|6|59|60|15|
|7|29|30|30|

Every subset is now legal. Its product has square sigma iff its column sum
is zero. The parity-check matrix has every nonzero F2³ column once, so has
rank3, dimension-four kernel, and no word of weight1 or2. Each unordered
pair x!=y completes uniquely to the triple{x,y,x+y}; there are21/3=7 such
triples. The sum of all seven columns is zero, so their seven complements
give the weight-four words. Together with the empty and full words, these
exhaust16 kernel elements. This proves the weight enumerator

    1 + 7z³ + 7z⁴ + z⁷,

and identifies the classical Hamming[7,4,3] code. Its seven weight-three
supports are the Fano lines. The arithmetic realization and every square
test are independently checked by direct divisor enumeration of all128
subset products. Each prime in the table is also the least realizing prime
for that specified square class, verified through its displayed value; no
minimum over all coordinate systems is asserted.

**Why the graphic model cannot be universal.** This rank-three matroid is
simple and has seven elements. A simple graphic matroid of rank3 could have
at most six elements: remove graph isolates; its components have total
sum(v_i-1)=3 and at most sum binom(v_i,2)<=binom(4,2)=6 edges. Loops or
parallel edges would contradict matroid simplicity. Thus this Fano matroid
is not a graph cycle matroid. In particular it is not the cycle matroid of
Petersen, K5, K3,3 or any other graph. Binary incidence is a broader common
language than ordinary graph cycles.

### A literal marginal obstruction in the arithmetic code

Let C_t={x in F2^7: Hx=t}, for the displayed parity-check matrix H. Compare
the uniform laws on C_0 and C_1, with1 denoting the first nonzero coordinate
vector, not an integer source. Each law has16 squarefree arithmetic sources.
For any chosen r<=3 coordinate bits, the remaining at least four distinct
nonzero columns span F2³ (a two-dimensional subspace has only three nonzero
vectors). Therefore every assignment of the selected bits has exactly2^(4-r)
completions, independently of t. All their marginals through order3 agree.

Each coset also has eight even-weight and eight odd-weight words: translating
by an odd-weight codeword gives a parity-reversing bijection within the coset.
Thus the marginal joint law of(mu(n),lambda(n)) agrees as well, with the two
signs equal because every source is squarefree. Nevertheless square sigma has
probability1 on C_0 and0 on C_1. The separate sign-law assertion does not claim
equality of its joint distribution with every selected bit tuple.

This is a precise finite counterpart of the prior Vitali/guard correlation
mechanism: low-order observations and an additional sign law can miss a joint
predicate. It is not an all-bank lower bound or an identification of these
arithmetic sources with Vitali atoms. All64 coordinate subsets of size<=3 are
checked, independently from both the parity matrix and direct divisor sums.

### Frobenius and an explicitly chosen tournament

Identify the parity space with F8=F2[alpha]/(alpha³+alpha+1), using the
polynomial basis(1,alpha,alpha²). The powers of alpha have bit values

    1,2,4,3,6,7,5,

so the corresponding prime order is `(7,2,19,5,59,29,89)`. Frobenius x->x²
fixes label7 and cycles

    (2,19,59), (5,89,29).                                    (5)

The second orbit corresponds to beta=alpha+1 and its conjugates, of bit
values(3,5,7). They are a normal basis: their eight binary sums are distinct.
Consequently Frobenius rotates coefficient decorations
`(d0,d1,d2)->(d2,d0,d1)`. The first orbit(2,4,6) sums to zero and is NOT a
normal basis. This retains the prior warning about trace-zero orbits.

Transporting a subset by(5) preserves zero syndrome and subset size, hence
square sigma and the source signs mu=lambda=(-1)^size. It does not preserve
arithmetic weights: `{2,5,7}` maps to`{19,89,7}`, changing source70 to11837
and sigma144 to14400. This gives an explicit commuting square for the Boolean
predicate, not a conjugacy of multiplication or of aliquot dynamics.

Choose the Singer order above and orient i->j iff j-i mod7 is in{1,2,4}.
There are no ties between distinct points. This is a chosen Paley7 tournament
on the prime labels; the opposite residue set gives its converse. Frobenius
doubles exponents and preserves the orientation. Exact enumeration gives189
Hamilton paths and21 automorphisms. The Fano plane itself has168 automorphisms:
every line-preserving permutation extends to a linear map, and GL3(2) has
7*6*4 elements. Thus incidence alone does not select the tournament frame.
Neither the additive field marking nor its orientation is claimed canonical
for the natural numbers. No4p threshold or81/27 classification follows solely
from this seven-point construction.

## A4. A common legal continuation breaks the 4/13 forward fold

Session43's N=2^5*3^7*5^3*7^4*11*13² gives h_N(4)=h_N(13)=2, after which
all their h_N iterates agree. The divisor-height realization theorem and
5759-base census are inherited, and the census is replayed here.

The common coprime multiplier205=5*41 has sigma252. Hence

    sigma(4*205)=7*252=1764=42²,
    sigma(13*205)=14*252=3528, not square.                     (6)

All m<205 coprime to52 are checked; none separates the two square outputs.
Thus205 is the least positive common legal multiplier for this particular
test. Equality of the divisor-height futures does not survive adding this
arithmetic attachment experiment. The fold was lawful in its original
language; enlarging the observation language changes the appropriate quotient.

## A5. F=U+S and a literal Berggren square-sum insertion

Use the recovered count definitions on1<d<n:
F counts all such divisors, U the primes, S the squarefree divisors. U is
contained in S. The inherited classification is n=p,p³,p²qr with distinct
prime letters, for n>=2; it is not an identity of disjoint subsets.

The circuit{4,7,13} yields364=2²*7*13, with(F,U,S)=(10,3,7), sigma784
and s364=420. It is the least p²qr square-sigma source, checked over every
smaller integer. Concurrent Session84 independently supplies this same minimum
and the general prime-guarded parametrization; credit is shared with that
incoming work. The June HYP-2220 source91=7*13 has sigma112 and s91=21;
coprime multiplication by4 supplies the missing square-class7, making364.
Only these elementary identities are inherited from that OPEN synthesis.

For distinct odd primes r=14k²-1, r!=7, the same circuit yields

    n=28r, s(n)=28(14k²+1), sigma(n)=(28k)².                  (7)

Here k>=1, and the displayed primality guard is essential (k2 gives55,
composite). The endpoints differ by56; no infinite prime-supply claim follows.

The October19 absorber theorem is a stronger geometric bridge than matching
numbers. The classical Berggren matrices preserve Q(a,b,c)=a²+b²-c².
The same branch words enumerate primitive Q=0 triples from(3,4,5) and odd
positive Q=1 triples from(1,1,1). This is a branch-preserving address bijection,
not a linear identification of the shells. At the middle child(5,5,7), any
square-sum edge x+y=t² has the local replacement

    u=25t²-x, v=25t²-y;
    u--x--y--v, with u+v=49t².                               (8)

At(x,y)=(364,420), (u,v)=(19236,19180). All four vertices are distinct,
and the edge-square roots are140,28,140, with host root196. This gives the
promised literal map from the arithmetic circuit to a graph insertion gadget.
An actual Hamiltonian host still must contain the chosen host edge and satisfy
all other coverage constraints; (8) alone makes no global existence assertion.

## U. Proposed unification and exact boundaries

The organizing object should be a **typed algebra of guarded attachments**,
whose observable states are residual responses to all legal future contexts.
The ingredients are joint incidence (who shares a port/resource), transport
(how a frame or parity changes), and realization (which formal choices are
actual integers, paths, labels or geometric pieces). We propose classifying
these response algebras and their quotient maps, rather than identifying all
the underlying mathematical objects.

The coarsest exact quotient is the standard context equivalence:
`x~y` iff every allowed one-hole context has the same output and legal status
on x and y. Closure follows by putting another operation around the context;
any exact recursive carrier must refine it. It can have infinitely many
classes and gives no automatic finite algorithm. Session84 independently
states and proves this standard theorem, crediting the older THM-840.
[Courcelle–Weil](https://arxiv.org/abs/cs/0609109) supplies the established
graph-algebra framework: recognizability depends on the operation signature.
[Lovasz](https://arxiv.org/abs/math/0408232) studies connection matrices and
graph algebras. Their homomorphism representation theorems are not imported
without their additional hypotheses. Distinct response rows and their linear
rank remain different invariants, as the prior Mobius packet already proves.

Three concrete directions follow.

* **Guarded binary realizations.** Classify which sigma parity charts are
  graphic and when non-graphic minors appear; record the coprimality constraint
  alongside the binary code. The graphic chart/Fano extension is a completed
  small example, not an all-bank characterization.
* **Boundary response algebras.** Combine the finite arithmetic observer with
  the graph pairing and tournament joint-port calculi, but compare only the
  operations actually shared. Theorem A of Session70 needs exactly(2r-1)!!
  matching states for unrestricted named path connectors; its graceful overwrite
  model needs k! permutations and can discard2^k signs. The conflicting answers
  are explained by different contexts, not a failure of either theorem.
* **Address, character and realization as separate layers.** Berggren words
  have unique integer addresses; a branch character can nevertheless update
  in one bit. Fixed endpoint residues lose that bit, while a normalized frame
  restores it. Our Frobenius arithmetic marking is another exact finite
  transport with a weight-loss boundary. Neither proves growing memory for
  every target, and neither turns a prime-filtered subset into a subtree.

For Petersen->K5 retain one-/two-run cluster pairings. Their exact local
degree rule produces six2-factors, each two5-cycles; connectivity is an extra
global test. For Petersen-minus-one->K3,3 retain the three suppressed-edge
marks: two of the six quotient Hamilton cycles lift to spanning9-cycles.
For graceful trees retain component incidence and colors; their overwrite
sign quotient does not solve existence for every prescribed abstract tree.
For the square-sum problem retain actual endpoint labels and future arrivals.
These inherited controls are replayed or explicitly cited in the receipts.

The12-vertex81/27 pair is still the hostile control for treating arithmetic
summaries as a universal tournament carrier: equal H=462915 gives equal
mu(H),lambda(H), but their counted complement responses are respectively
(-2916,-3321) and(297,-1323). Those figures are inherited from the pinned
duality packet, not reenumerated here. The current results do not classify
strong atoms or the full unnamed complement/re-marking system.
