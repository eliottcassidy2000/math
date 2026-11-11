# Session266: classify the frequency, then retain the missing state

**11 November2026. Final synthesis. Chan–Pak and Collatz remain OPEN.**

The session answers both readings of “how many times has an odd Hamiltonian
count appeared”: frequency at a particular order, and total frequency across
all orders up to a bound. It also finds an exact reason the requested missing
Collatz arcs display7 and21, and a common algebra behind tournament padding,
radix digits and the money-plant puzzle. Old proved inputs are inherited rather
than re-audited; new maps and constructors have bounded exact controls.

## 1. What is now completely specified, and in which scope

H(T) counts **directed Hamiltonian paths**, not cycles. Four different channels
are retained: ordinary isomorphism classes A_n(h), self-converse classes,
classes after identifying converse, and labeled tournaments. “Up to N” means
sum over sizes1 through N, with no auxiliary empty tournament included.

The deliverables have three separate completion scopes:

- **Every H through order8:** a complete inherited census, recast as an exact
  cumulative table with all ordinary, self-converse, merged and labeled
  channels. It covers7412 ordinary classes across sizes1..8, including explicit
  zero rows for7 and21. [Table and finite scope](fibers/FINITE_SCOPE.md).
- **Every order for each odd H<=73:** an executable exact compiler using the
  complete frozen strong-core inventory. The bound N is a query parameter,
  not an enumeration cutoff. [Compiler and proof](fibers/PROOF.md).
- **Any specified H:** an effective finite reduction. Enumerate the strong
  classes with H dividing that target, using the inherited bound
  |C|<=1+floor(3 log_5 H(C)), then use the same recurrence. This does not provide
  a practical closed inventory for every possible H, or prove that7 and21 are
  the only missing odd counts.

The older Session24 and opus S50 already established the strong-core method,
eventual polynomials and complete inventories below625. Those are credited
inputs, not discoveries of this session. The smaller executable input here
keeps automorphism and converse data together and is reusable at arbitrary N.

The actual law is simple once the right object is retained. Delete singleton
strong components to leave an ordered core word W=(C_1,...,C_k), of total size s.
All padding choices have the same H=product H(C_i), and the word contributes

    binom(n-s+k,k)                 at order n,
    binom(N-s+k+1,k+1)             across orders <=N.

A term is zero before its actual minimum size. Add these over the finite core
catalog. For labeled counts, weight the word by1/product|Aut(C_i)| and multiply
the order-n coefficient by n!, rather than multiplying ordinary classes by n!.
Converse requires its separate center-and-reflection formula.

Example: H63 has exactly two strong cores, both at order8 and converse-paired.
Thus at order N>=8 there are2(N-7) ordinary classes, but across all orders up
to N there are(N-7)(N-6). At N=30 these are46 and552; converse merging gives
23 and276. This distinction is the principal counting ambiguity in the request.

Two new structural consequences help turn tables into explanations. Ordinary
fixed-H frequencies are discretely convex even before polynomial onset.
Their second difference is exactly the birth count of strong classes whenever
H is irreducible in the realizable-count monoid. More generally, triangular
inversion over divisors recovers every strong inventory from the growth table.
A single scalar histogram does not recover its automorphism or converse sidecars.

## 2. An infinite family now has explicit degree and leading constant

For every h=7^a3^b5^c except7,21, the maximal number of nontrivial strong
components, and hence the ordinary multiplicity polynomial's degree, is

    r(h)=b+c+floor(a/2),

except for h=7*3^b, b>=2, where r=b-1. This is unconditional: all required
scalar atoms3,5,35,49,63,343 have inherited strong realizations. No assumption
about all other odd integers is needed.

The strong inventories at these six atom counts have total masses
1,1,2,8,2,103. The last comes from the old complete polynomial
P_343=12z^8+67z^9+16z^10+8z^11. From these inputs we give a closed formula for
the leading coefficient at **every**3,5,7-smooth h, not merely the degree.
For a=2m, it is8^m/(m!b!c!). For a=2m+1 outside the exceptional ray, it is

    [2c*8^m/m! + 1_(m>=1)*103*8^(m-1)/(m-1)!] / (b!c!).

On the exceptional ray it is2/(b-2)!. The full proof, equality patterns and
constructor are in [the smooth classification](smooth/PROOF.md).

This reveals a useful nonunique-factorization phenomenon:

    63² = 49*3⁴.

The left gives two irreducible strong cores; the right gives five. Accordingly
A_n(3969) has leading term n^5/3, even though the obvious two-core construction
only supplies quadratic growth. Similarly A_n(1225) has leading term4n³,
and A_n(16807) has leading term824n². Equal scalar counts can hide sharply
different amounts of future padding freedom. Among factorizations into the
six scalar atoms, the exact maximum ratio of longest to shortest length is5/2.

This does not supply all lower coefficients for arbitrary smooth h. Strong
cores with reducible scalar H contribute lower-degree terms and still need
inventory data. The leading classification is complete in its stated scope.

## 3. The two deleted Collatz arcs really do encode21 and7

Use C_+(x)=x/2 for even x and3x+1 for odd x. On negative integers, negation
conjugates the map to C_-(x)=x/2 for even positive magnitudes and3x-1 for odd ones.
The previous multiplication-table atlas assigns the sign-s source x the code

    h=4x+s.

The user's positive cut5→16 therefore has code21. The negative cut−2→−1,
viewed as the minus magnitude edge2→1, has code7. This is an **exact map**,
not a matching numeral found after the fact. Odd-source and oriented even-source
decoders are both inherited from the table atlas.

Every odd h>=3 has one such signed source code. We can therefore attach the
frequency A_n(h) to that signed edge, retaining h=1 as the transitive boundary
symbol. The entire tournament classification becomes a weighted signed-edge
atlas. Ordinal composition has the exact code product

    (x,s)⊙(y,t)=(4xy+t*x+s*y, s*t).

It preserves multiplication of H while forgetting SCC order, size and internal
structure. Those are precisely the data needed to recover actual multiplicities.
The assertion that the eventual support of the atlas consists of all edges
except the two prescribed ones is equivalent to the full odd-spectrum
statement, so it is not asserted as solved.

The transported raw Collatz step on codes depends on h modulo8:

    1: (h+1)/2;  3: 3h-2;  5: 3h+2;  7: (h-1)/2.

In particular it sends **realizable H3 to impossible H7**. This is a tiny
explicit obstruction to implementing that next-step rule by a universal
operation on tournaments. A numerical code map need not be a realizable-count
map, even when its multiplication has a genuine tournament interpretation.

The cuts also have different topology. Removing5→16 cuts a feeder but leaves
1→4→2→1 intact. Removing−2→−1 breaks the negative two-cycle itself. The known
negative cycles through−5 and−17 survive. No claim is made that the three known
negative cycles or one known positive cycle exhaust their respective graphs.

For a shortcut parity word with L steps and r odd branches, its unique periodic
candidate is sB/(2^L-3^r), with B>0. The two signs share modular word data,
but positivity depends on the **ordinary sign** of the denominator. That
missing height/sign coordinate cannot be recovered from a finite residue table.
[All definitions, maps and proofs](collatz/PROOF.md).

## 4. Four-vertex tournaments give a completion measure, not a cycle bijection

The four ordinary tournament types have H values1,3,3,5. The two H3 types
are source-plus-triangle and triangle-plus-sink, exchanged by converse. Hence
there are four ordinary classes but three converse classes. Separately there
is one strong type and three nonstrong types. None of those partitions proves
a correspondence with Collatz's observed cycle counts.

A genuine bridge comes from retaining an observed orbit segment as a path.
Completing its missing chords yields each tournament class exactly H/|Aut|
times. At order4 the eight completions split1,1,1,5 across the four types.
So a Collatz path does not select one tournament until its three absent chord
orientations are specified.

Delete a path edge instead, leaving ordered path-block lengths lambda. The
completion fiber is its ordered path-cover response C_lambda/|Aut|. In order4,
a gap after the first vertex gives weights4,2,2,8; a middle gap gives6,2,2,6.
The location of the missing arc changes the measure.

For every order, every tournament class has the same parity of this fiber:
it is odd precisely when adding the block lengths in binary has no carries.
This follows from the universal oddness of tournament Hamiltonian counts and
the binary multinomial identity. At a power-of-two order, every proper puncture
profile has even fibers. This is a concrete meeting point of missing arcs,
path-cover algebra and radix carries, with the required chord sidecar explicit.

## 5. Tao's plants reveal the same discrete-integration algebra

Starting before the Day1 harvest, each mother earns one coin and creates two
regular plants that first earn on the following day. Regular plants persist
and earn one coin daily. Thus

    coins at end Day2 =4a+2b+c,
    coins at end Day3 =9a+3b+c,
    coins at end DayX =aX²+bX+c.

The Day10 total427 uniquely gives4 mothers,2 regular plants and7 initial coins.
For any base B>=2, with a,b,c<B, Day B reads the three base-B digits directly.
Day B is uniquely both injective and gap-free on the B³ allowed initial states.
Keeping Day10 while changing the coefficient bound does not change the
numeric evaluation base to B.

The daily population matrix is I+N with N³=0. More generally a persistent
hierarchy realizes the forward-difference coordinates of any polynomial with
nonnegative integer coefficients. Tournament singleton-gap counts obey this
same unipotent finite-difference algebra. Converse reflection adds a possible
period-two sector. The link is an explicit recurrence, not plants relabeled
as vertices.

Radix evaluation p(x)↦p(B) has kernel(x-B): carries identify states that agree
at one measurement but have different future trajectories. The six-atom
relations above do the same for core factorizations: equal H hides different
polynomial growth. These are precise examples of why the recursive state,
not merely its present numerical output, must be retained.
[Puzzle sources, arbitrary bases and full proof](radix/PROOF.md).

## 6. An electrical sidecar survives arbitrarily long padding

Normalize the weighted path network by H and bracket a strong core C with a
source and sink. Its effective resistance rho(C) lies between2 and|C|+1.
It equals2 exactly when every vertex starts as many Hamiltonian paths as it
ends. Define delta(C)=|C|+1-rho(C).

When every pair of nontrivial cores is separated by at least one singleton,
ordinary series addition gives the exact identity

    H R(first,last)=n-1-sum delta(C_i).

The gap lengths disappear from the defect. A triangle contributes2, the strong
four-vertex core contributes650/219, and the regular five-vertex core contributes4.
Thus two H15 families have different persistent network defects:4 for the
single regular core and1088/219 for the triangle/four-vertex pair. The latter
wins in frequency as n grows because it has two nontrivial components.

The resistance-resolved generating series gives a finite limiting defect law
for each fixed H. This is a new refinement of the multiplicity classification,
not a Collatz step or an inverse reconstruction of all tournament arcs.
[Proof and literal constructions](geometry/PROOF.md).

## 7. What remains open

The missing task is not another check of the old small census. It is to grow
the strong inventory efficiently, find a structure-preserving consumer for
the signed-edge atlas, and control the sidecars that scalar H forgets. The
full odd-count spectrum, general Collatz convergence and exhaustive cycle
classification remain open. [Next questions](OPEN_QUESTIONS.md),
[sources and inheritance](SOURCES.md), [reproduction](REPRODUCE.md).
