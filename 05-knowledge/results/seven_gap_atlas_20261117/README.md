# A gap atlas for 7p: local placement, exact families, and the missing multiplier

**Proved.** Every prime target 7p in the first five nondegenerate dyadic
gap bands has a literal witness and is the anchor of an infinite
physically realized ray in the gaps. The ray remains divisible by 7 and
coprime to 3 at a fixed period; the quotient need not remain prime.
At all sufficiently large scales, every interior interval in a gap
contains an actual strong one-apex count of any prescribed admissible
residue modulo any fixed modulus. In the opposite direction, no 7p
witness can be built by replacing a vertex of a quotient by a homogeneous
p-path module, even if paths can revisit the module. Together these
separate numerical multiplication by 7 from the graph operation that
creates the factor 7.

## 1. What is being missed?

Let I_s = [18·2^s, 28·2^s]. The inherited [dyadic-wrapper
obstruction](../two_cover_feedback_20261117/README.md) is a missing-scale
obstruction of that architecture, not nonexistence in the tournament
spectrum. The prime number theorem yields on the order of 2^s/s targets
7p in I_s, for growing s. The atlas classifies the first five bands
s=2,…,6; no finite table certifies all subsequent prime targets.
The earlier two bands just give 21 (impossible) and 49 (see below).

The [literal atlas](witnesses.json) covers **all 35** prime targets
in the five bands. Each has an explicit retained-root word; the
[replay](replay.py) computes its Hamiltonian total in three ways:
the root start-menu transfer, direct subset DP on the constructed
tournament, and the [odd-cycle packing polynomial of Grinberg and
Stanley](https://arxiv.org/abs/2307.05569). The polynomial is obtained
by counting odd cycles independently, then summing disjoint packs;
each packed cycle contributes 2q, and each uncovered vertex q.
The [receipt](receipt.out) lists its coefficients, the strong-core
order, and its complete number of proper nonsingleton homogeneous
modules for the displayed retained-root witness. Selected additional
witnesses use smaller literal cores, a single apex, or cyclic gluing;
no minimal-order claim is needed for them.

A one-apex tournament is formed by giving a transitive spine one new
vertex. Its word records outcontacts 1, incontacts 0 in spine order.
After deleting leading 0's and trailing 1's (ordinal singletons), a
nontransitive word begins 1 and ends 0 and is strong. The count is

    H(w) = 1 + sum over i<j with w_i=1,w_j=0 of (2^(j-i-1) + [j=i+1]).

An odd cycle uses the apex and such an inversion, with an even subset
of intermediate vertices in increasing order; two cycles cannot coexist.
This is also the q=1 case of the inherited one-apex slot polynomial.
The first-to-last inversion forces H ≥ 1+2^(n−2) at spine length n.
Therefore the one-apex test in the replay is **complete at every order**
for each *fixed target*: it only needs active lengths through
floor(log2(H−1))+2. Exactly 25 of the 35 targets have one-apex words
(62 active words, not 62 isomorphism classes). The ten exceptions to this
lens are 77,161,301,707,763,791,1211,1267,1631,1687. These ten all have
witnesses by the retained-root construction.

The number 217=7·31 shows the distinction sharply. It has a strong
eight-vertex witness, a nine-vertex retained-root witness with no proper
nonsingleton homogeneous modules, and a ten-vertex one-apex witness
100110010. In the last lens its pack polynomial is

    q^10 + 24q^8 + 94q^6 + 84q^4 + 14q^2.

In the nine-vertex lens it is

    q^9 + 22q^7 + 100q^5 + 86q^3 + 8q.

Both specialize to 217, encoding different cycle structures.
The [last session](../two_cover_feedback_20261117/README.md) proved
that NO input to its two-run wrapper followed by any number of
near-sinks hits 217, even with a non-dyadic two-cover input. This
is a counterexample to an operation, not a sporadic missing count.

## 2. Every atlas anchor belongs to a physical infinite family

All three formulas below are exact, not asymptotic.

**Retained root.** On a normalized root with weighted even and odd
start counts E,O and deletion count g, normalization is O=E+g.
Appending a near-sink (loses to everybody but the root) sends E to 2E
and H to H+2E, preserving normalization. Thus k more such vertices
give alpha·2^k+delta for alpha=2E and delta=H−alpha.

**Apex tail.** For a head v beginning and ending 1, set A equal to its
binary integer. Appending t ≥ 1 zero contacts gives

    H(v 0^t) = A·2^t + H(v)+1−A.

The additional geometric sum is A(2^t−1); the last head vertex
contributes the extra adjacency bonus 1. Every nontransitive
one-apex word is in one such ray after stripping its final zero run.
To classify an individual ray in the gaps: normalize A by powers of
two to [16,32). An interior coefficient in (18,28) stays in a gap
eventually; an exterior coefficient does not. At the boundary, use
the sign of the constant. Divisibility by 7 is decided by the three
values of 2^t mod 7, and nondivisibility by 3 by its parity.

**Cyclic gluing.** Replace the three vertices of a directed triangle
by A,B,1. Put (u,R)=(P_A(1),P_A(2)) and (v,S)=(P_B(1),P_B(2)).
The singleton caps both other blocks at two visits. There are three
one-run quotient words, one additional word for each possible double
visit alone, and one with both double visits. Therefore

    H = 3uv + (R−2u)v + u(S−2v) + (R−2u)(S−2v)
      = (S−v)R + (3v−S)u.

Ordinal singleton padding in A replaces R by 2^k R and keeps u fixed.
This is a separate two-cover reader, not a special address for the
dyadic wrapper. The inherited general identity is
[THM-3121](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).

**Assignments.** The atlas stores the family as well as the witness.
For 30 targets a retained-root realization has alpha in the *same*
band as the anchor. For 217, 721, 889 use an apex tail (respectively
306·2^k−89, 1270·2^k−549, 1322·2^k−433). For 1687 use a rooted tail
2862·2^k−1175 that enters the next scale after growing. For 1757 use the
cyclic ray 1788·2^k−31. In each case N is the literal k=0 value.
The replay verifies that alpha, and H(k0)/2^k0, lie in the indicated
fixed band, with k0=0 or 6. The normalized count stays between these
two values for all k≥k0, so every subsequent multiple of six is in
some gap. Because 2^6≡1 mod 21, it stays divisible by 7 and coprime
to 3. For the exceptional anchor 49 the rooted profile (H,g,E,O)
=(49,9,24,33) gives 48·2^k+1, with the same conclusion.

**One family, two genuinely different geometries.** The count
768·2^k−5 comes from both the near-sink extension of the seven-vertex
rooted seed (6,12,96,21,7,27,59), root 2, by k+3 vertices, and from
C3[A→TT_k, TT_2, 1] with A=(126,124,88,112,96,68,0), profile (5,256).
These two witnesses are nonisomorphic at every k. In the first,
the six old nonroot degrees are {2,2,3,3,4,5}+(k+3); new degrees
are 1,…,k+3, and the root has degree 2. Exactly one vertex has
outdegree k+7. In the second, the same degree k+7 occurs at the
old degree-five vertex of A and at the singleton quotient vertex.
At k=0 and 3 the counts are 7·109 and 7·877, neither of which has
a one-apex witness at any order by the complete word bound. Thus
even two distinct mechanisms agreeing all along one ray cannot
be collapsed into the simplest apex grammar.

Other explicit cyclic families are 812·2^k+49 (k odd; e.g. primes
239,3719,59399 after dividing by 7) and 1788·2^k−31 (k=0 mod 3;
e.g. 251,2039). The apex head 10011 with t zeros gives 19·2^t−13;
t=5 mod 6 lies in I_t and gives prime quotients 5557 at t=11 and
93262146997 at t=35. Primality of these examples is checked by
trial division; nothing here implies infinitely many prime quotients
on any one ray.

## 3. The entire gap has a local decoder, even at a fixed scale

**Theorem.** Fix any open interval J inside (18,28) and any integer
modulus M≥1. For every sufficiently large s, and every class r mod M
compatible with oddness, there exists a strong one-apex tournament
with exactly s+6 vertices, H/2^s in J, and H congruent to r mod M.
For fixed J, the construction works with a threshold
s ≥ C(J)+O(M²), uniformly in r.

**Scale lens.** At spine length n, set L_top=sum w_i 2^(-i) and
R_top=sum (1−w_j)2^(j−(n−1)), with zero-based indices. The weighted
inversion formula gives

    H(w) = 2^(n−2) L_top R_top + O(n),

uniformly: unallowed terms j≤i have total at most n/2 and the
adjacency contribution plus 1 is at most n. For n=s+5, the
normalized value is 8 L_top R_top + O(n/2^s). With first bit 1
and last bit 0, each factor can be pinned arbitrarily closely
to any prescribed number in (1,2), by a fixed head and tail.
Choose their product strictly inside J/8. Every middle word is
now allowed without leaving J, once s is sufficiently large.

**Residue lens.** For exact concatenation U,V use

    H(UV)=H(U)+H(V)−1+L(U)R(V)+[last(U)=1,first(V)=0],
    L(U)=sum U_i 2^(|U|−1−i),    R(V)=sum (1−V_j)2^j.

Write M=2^a m with m odd. Let d≥max(1,a) be a multiple of
the order of 2 mod m (take order 1 if m=1). In a block of
length 2a+d+3, allow either ones at a,a+1 or at a,a+1+d,
all other positions zero. The boundary bits are the same.
Both L and R are unchanged mod M on swapping these blocks:
the shifted weights are equal mod m, and both positions
are at least a places from both edges for the 2-adic read.
The block's H change is 1 mod m and 2 mod 2^a. To see this,
sum the outgoing contributions as if every later site were
zero, then subtract the single one-one interaction, and add
the number of adjacent 10's. The individual geometric sums
are unchanged mod M; the interaction changes from 1 to 2^d,
and splitting the run adds one adjacency. Thus the total
change is 2−2^d mod M. The concatenation identity makes it
independent of the rest of the spine. Its gcd with M is
gcd(2,M), so M/gcd(2,M)−1 independent slots suffice to select
every admissible r. They consume O(M²) bits since d≤a+m.
Ordinary zero filler sets the precise length s+5. With the
anchored initial 1 and final 0, the resulting tournament
is strong.

In particular, take r=7 modulo 7M. This supplies physically
realized H=7q in every interior slice at every sufficiently
large scale with q avoiding any prescribed finite set of primes.
It yields infinitely many distinct prime divisors across quotients,
not infinitely many prime *quotients*. The latter is a genuine
global selection problem: a modulus fixed while s grows does not
control the residual within that residue class.

## 4. A seed can re-enter a gap; its prime factor cannot be copied locally

**Induced re-entry lemma.** Given any T with H(T) divisible by 7 and
not by 3, and any open J inside (18,28), T embeds as an induced
subtournament in an infinite strong family of 7-channel counts
whose normalized position is in J. Ordinally append a literal
(17,448) block with rows (126,60,120,112,96,64,2) and j copies
of the (11,100) block with rows (24,1,3,6,14), obtaining D_j.
Its profile (u,R) has 7 dividing both entries and 3 not dividing u.
The cyclic wrapper C3[D_j→TT_k,1,1] reads u+2^k R. Choose j so
that R normalized into [16,32) lies in J: R is a fixed number
times 100^j, and log(100)/log(2) is irrational. Then take k large
enough and of the parity making |D_j|+k even. The odd-cycle
polynomial has all exponents congruent to the order mod 2,
so P(2)≡(−1)^order H mod 3. The output is therefore 2u mod 3 and
stays in the desired region. This does not preserve a prime
quotient; common factors of the smaller profile can be inherited.

**No local 7-amplification.** For *every* tournament B and fixed
quotient G, replacing just one vertex of G by B with all arcs
to outsiders uniform cannot give H(T)=7H(B). The quotient is
allowed to depend on B, and there is no bound on its order.

Proof: put u=H(B), n=|B|. The paths visiting B once give H(G)u.
Each additional feasible quotient word with r visits contributes
c_r(B), the number of ordered covers by r nonempty paths. By
cutting one full path at r−1 places, c_r(B)≥binom(n−1,r−1)u.
Also c_n(B)=n!≥2^floor(n/2)u: partition all permutations into
orbits under swaps of disjoint consecutive pairs, only one
choice of orientation in each pair can occur in a Hamiltonian
path. If n≤3, u=1 or 3 and 7u is forbidden. Otherwise H(G)>7
is too large and H(G)=7 is forbidden, leaving 1,3,5. Count
1 means G transitive, with no repeat. Count 5 leaves a budget
of 2u; any repeat exceeds it by the two cut bounds. Count 3
means a unique directed triangle, by the odd-cycle packing
formula, so its only nonsingleton SCC is that triangle. For
completeness, if a tournament has a unique triangle, every
outside vertex is uniform to it (mixed contacts would create
another triangle). Contract it: a triangle through the
contracted vertex would lift, and the remaining quotient
is therefore transitive. If the replaced vertex is in
the triangle, H(T)=u+P_B(2), as in the cyclic formula; else
there are no revisits. The required P_B(2)=6u is impossible:
complementary vertex subsets pair up the terms of P_B(2),
and all terms are odd, so 4 divides P_B(2) for n≥2.

**Consequences for any 7p witness.** Its nontrivial SCC is unique:
SCC counts multiply, every nonsingleton SCC is strong and has
at least three paths, and the only nontrivial factorization of
7p uses the impossible factor 7. All other components are
ordinal singletons. Moreover *no proper module of the strong
core* has a path count divisible by p. Indeed the proper
quotient stays strong, has at least three vertices and has
H(G)≥3, so u≤7p/3 and an odd multiple of p there must be p
itself. The local 7-amplification theorem excludes this.

The synthesis is that the gap is a failure of one scaling law,
while the factor 7 in a minimal witness is a property of the
interaction between blocks, not a count carried inside one.
Scale and residue can already be separated by physical edits.
What remains is to retain that freedom while paying the *exact*
integer residual. An exact target interval with actual spacing,
not another residue net or another sparse exponential ray, would
make the next genuinely different proof step.

Replay: run `python replay.py`. Normal and optimized Python
have identical receipts. The finite witness scope is explicitly
the five bands plus the separately displayed examples. Chan–Pak
remains open.
