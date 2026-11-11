# Summand births, two-seed interiors, and a delayed parity frontier

Session259, 2026-11-11. Status: **PROVED in this research note**, with a separately bounded **FINITE-EXACT** companion. No numbered canon promotion is claimed. Novelty means a new development in this session, not an external priority claim.

## 0. Inheritance and concept board

Inherited definitions and mechanisms:

- `01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`: unordered paired-parent cospans, swap-fixed doubles versus squares, distinct closure of `{2,3}`, and its synchronous dyadic frontier.
- `01-canon/theorems/THM-2433-operation-fibre-deletion-incidence-and-startup-scar.md`: universal deletion convolution and the constant eventual additive scar. The deletion formula below is an application, not a new theorem.
- `01-canon/theorems/THM-362-natural-operation-graph-shadows.md`: the existential shadow loses co-parents. Its transitive tournament is NOT the object studied here.
- `07-reflections/the-summand-alternation-is-an-involution-parity-law-kps-S131.md`: historical suggestion to promote the swap-fixed locus, not only count its parity.

Hostile controls: deleting equal parents removes `4=2+2`; the `{2,3}` final hole is6, not4; synchronous node14 is late; one Fibonacci spine is not closure. The least-used sidecar here is the gap between the largest even and largest odd available values. It reveals a very long moving-boundary transient invisible in final closure.

Live board: paired-parent matching / numerical-semigroup interior / parity maxima / finite deletion moments / expression-versus-circuit cost. Anchor: vertex births. Niche: exact closure for arbitrary two seeds. Wildcard: exponentially large frontiers with only linear-time parity relaxation. No implication for LRC(14), prime gaps, or the open twin-center ancestry assertion is made.

## 1. Birth is a matching, not an ordinary adjacency increment

At numeric horizon N retain all positive vertices at most N and every cospan `{a,b}->a+b<=N`, with `a<b`. When n arrives, its complete new target fibre is

    {{1,n−1},{2,n−2},...,{floor((n−1)/2),ceil((n+1)/2)}}.

Each old vertex occurs in at most one pair in this fibre. For odd n it is a perfect matching on `[1,n−1]`; for even n the half-vertex n/2 is unmatched. Restoring equal parents attaches the diagonal singleton `{n/2,n/2}`. This is a statement about paired inputs; it is not the simple shadow.

The inherited count is r(n)=floor((n−1)/2), and total cospans through N are floor((N−1)^2/4). A target n first exists exactly when n is admitted, since both positive parents are strictly smaller. Parent cutoff and target cutoff are different filtrations: admitting parent N also licenses sums up to2N, but these are not yet vertices in the target-cutoff graph.

### A matching's exact failure budget

**PROVED.** For a fixed n, deleting old vertices destroys all strict witnesses iff it hits every matching edge. The minimum number of deletions is r(n), and there are exactly 2^r(n) minimum deletion sets, one endpoint chosen independently from each pair. In the weak fibre at an even n, the half-vertex is additionally mandatory, so minimum size is r(n)+1. Counts refer to minimum subsets of the relevant parent vertices, not arbitrary supersets.

Under independent retention of each old vertex with probability p, the probability of no strict witness is `(1−p^2)^r(n)`. The weak even version is `(1−p)(1−p^2)^r(n)`. Independence is within one target because its parent pairs are disjoint; events at different targets are NOT independent. A union bound nevertheless shows that iid positive-density retention almost surely has only finitely many targets with no strict witness, because the displayed failure probabilities are summable. This does not model primes or establish anything about their dependence.

## 2. Exact closure from any two distinct seeds

Let 0<a<b, let

    H(a,b)={ia+jb:i,j>=0},
    C(a,b)=Cl^<({a,b}).

Here H is an additive monoid including0. Then

    C(a,b)={a,b} union (a+b+H(a,b)).                         (S1)

**Proof.** Every nonseed construction must contain both seed types: the first nontrivial permissible operation is a+b. More formally, the right side contains the seeds and is closed under sums of distinct values. The only two-seed combination is a+b; any combination involving a translated-monoid element remains in that translate. Conversely, a+b is available, and adding a or b repeatedly is always legal because every running total is larger than both seeds. This constructs every ia+jb with i,j>=1. Thus equality holds, including numerical coincidences between different coefficient vectors.

This is stronger than a cofinite-existence argument: it identifies exactly which semigroup ideal the strict-parent rule selects. Distinctness does not mean each seed is used only once. It forbids equal VALUES at a binary gate, not repeated leaf labels across a derivation.

### Coprime case and the final hole

Assume gcd(a,b)=1 and 2<=a<b. Then

    largest positive hole F=ab,
    conductor c=ab+1,
    number of positive holes h=(ab+a+b−5)/2.                 (S2)

For completeness, the ordinary two-generator monoid has final hole ab−a−b. For any integer n select the unique j in `[0,a−1]` with jb congruent to n modulo a. For n>ab−a−b, n−jb is a multiple of a greater than −a, so it is nonnegative. The value ab−a−b cannot be represented: reducing a supposed representation modulo a forces its b-coefficient to be at least a−1, already too large. Translation by a+b gives final hole ab in (S1), and the two seed exceptions do not fill it.

The ordinary monoid's number of positive holes is

    sum_(j=0)^(a−1) floor(jb/a)=(a−1)(b−1)/2.

Indeed its least representative in residue jb modulo a is jb, and the nonzero residues permute1,...,a−1. Translating introduces a+b−1 additional positive holes below the translate. Restoring BOTH seeds removes two, giving (S2). Thus the correct numerator is −5, not −3.

Examples:

- `{2,3}`: holes `{1,4,6}`, h=3, F=6.
- `{3,5}`: holes `{1,2,4,6,7,9,10,12,15}`, h=9, F=15.
- `{2,b}`, odd b>=3: holes are1 and the even numbers4,6,...,2b. Hence h=b, F=2b. Every odd n>=b is present, and every even n>=2b+2 is present.

### General gcd and the seed1 boundary

For d=gcd(a,b), A=a/d and B=b/d,

    C(a,b)=d C(A,B).                                        (S3)

This follows because multiplication by d preserves sums and distinctness in both directions. If A>=2, the last missing MULTIPLE OF d is dAB=ab/d, and the number of missing positive multiples of d is `(AB+A+B−5)/2`. All nonmultiples of d remain missing forever. Calling this a cofinite subset of all positive integers when d>1 would be false.

If A=1, then

    C(1,B)={1} union [B,infinity).

For B>=3 the holes are2,...,B−1, conductor B and genus B−2. For B=2 there are no positive holes and the positive-integer conductor is1. This exceptional seed restoration is why formula (S2) must not be extended to a=1.

## 3. New connection to sum/product latency

The same coprime pair a,b>=2 has first mixed additive child a+b and last missing strict-closure value ab. Therefore the oriented sum/product transition

    a+b -> ab

has a precise closure interpretation: it links the entrance to the translated numerical semigroup to its final omitted integer. Its displacement

    ab−a−b

is the ordinary monoid's final gap, whereas the strict closure's final gap is ab. No assertion that multiplication causes additive holes is required.

Source: labelled two-seed cospan and its generated closure. Target: sum/product transition plus numerical-semigroup gap. Map: keep the same unordered parent pair. Preserved: both parent labels and their distinctness. Lost under just the displacement: the pair and the conductor separately. Sidecar: a+b or the unordered pair. Small test: `{2,3}` gives `2+3=5 -> 2*3=6`, latency1 and final hole6. For `{3,5}`, `8->15`, latency7 and final hole15.

## 4. Exact representation economics of the final closure

Let a,b>=2 be coprime, H be the finite hole set of C(a,b), F=max H=ab, h=|H|, and q=sum_(u in H)u. Let r_C(n) count distinct unordered parent pairs IN C landing at n. Targets need not themselves be in C for this definition, though closure forces every target with a witness to belong to C.

The inherited deletion calculus, specialized to this explicitly classified hole set, gives

    r(n)−r_C(n)
      = |H intersect [1,n−1]|
        − #{u<=v : u,v in H, u+v=n}.                       (S4)

The weak hole-pair correction includes its diagonal. This is necessary because a hole at n/2 was counted by the first term even though no strict diagonal edge existed. For n>2F this simplifies to

    r_C(n)=floor((n−1)/2)−h.                                (S5)

For every N>=2F, summing (S4) yields

    sum_(n<=N) r_C(n)
      = floor((N−1)^2/4) −hN +q+h(h+1)/2.                  (S6)

Proof: each deleted parent u contributes N−u to the first summed term; every weak unordered pair of holes contributes once to the second, and all such sums are at most2F. The total correction count is h(h+1)/2.

So a two-seed startup cost survives as an exact linear deficit in total witness storage, not merely as a finite list of absent vertices. For `{2,3}`, h=3,q=11, giving total `floor((N−1)^2/4)−3N+17` for N>=12. For `{2,b}` with b odd, q=b(b+1)−1, giving the analogous explicit formula with h=b.

The minimal existence certificate is radically cheaper than this full witness store. For any finite target cutoff, select ONE legal parent pair for each nonseed vertex of C. Values strictly increase along every dependency, so these selected pairs form a valid acyclic construction DAG with exactly `|C intersect [1,N]|−2` gates when N>=b. This is minimal when the output set must include every nonseed retained vertex and each addition creates one value. In contrast all legal cospans cost quadratically many records. This does not compute a shortest single-target addition circuit, nor does it preserve robustness or count all representations. One-parent-pair selection preserves constructibility and loses alternative witnesses and fault tolerance.

## 5. Synchronous growth has an arbitrarily long parity transient

Define `S_0={2,b}`, with odd b>=3, and

    S_(t+1)=S_t union {x+y:x<y in S_t},
    M_t=max S_t,
    T=(b+5)/2.

The eventual value set has last hole only2b. Yet its moving upper frontier has a parity defect lasting order b generations and reaching values exponential in b.

**PROVED.** For every 3<=t<=T,

    M_t=(3b+4) 2^(t−3)+2,                                  (S7)
    delta_t=b−2(t−3).                                      (S8)

The three largest members of the parity class containing M_t are

    M_t, M_t−2, M_t−4,

and the largest member of the opposite parity class is `M_t−delta_t`.
At t=T, delta_T=1, so S_T contains three consecutive top values. For all t>=T,

    M_t=(M_T−1)2^(t−T)+1.                                  (S9)

### Proof

Directly from the two seeds,

    S_2={2,b,b+2,b+4,2b+2}.

For b>=3 these are the five distinct values. At t=3 the top parity class is odd and its three largest elements are `3b+6,3b+4,3b+2`. The largest even value is `2b+6`, giving delta_3=b. One can check these maxima simply by the two parity classes of S_2: its even values are2 and2b+2, and its odd values are b,b+2,b+4.

Suppose a stage has global maximum M, its three top same-parity values M,M−2,M−4, and opposite maximum M−delta with delta>=3. Sums of two same-parity parents are even. The greatest such sum is2M−2, and2M−4 and2M−6 are also present. Sums
within the inferior parity class are at most2M−2delta−2<=2M−8, so they
cannot displace those three values. Sums of opposite-parity parents are odd; their greatest is2M−delta. Thus the new maximum is2M−2 and the opposite-parity deficit is delta−2. Retained old values cannot exceed the respective new maxima. This proves (S7)-(S8) by induction until delta=1. Notice that after the first step the global maximum is even; the argument uses parity classes, not a claim that the old maximum always had even parity.

When delta=1, the values M,M−1,M−2 are present. Their pair sums give2M−1,2M−2,2M−3, so three consecutive top values persist and the maximum follows (S9). No pair of distinct values at most M can exceed2M−1.

### Consequence and boundary

For b=101, static closure has last hole202, but the consecutive-top regime starts at generation T=53, when

    M_T=307*2^50+2.

Before that, the moving frontier keeps a strict even/odd imbalance; for3<=t<T the gap from the maximum to the second largest value is2. The large even and odd values approach one another by TWO per generation in deficit coordinates, while the whole frontier nearly doubles. This is an exact separation of absolute size growth from parity mixing.

This theorem concerns the top values. It does NOT assert that the full interval from the conductor to M_t has filled at generation T. The inherited `{2,3}` interval theorem is stronger in that special case. Nor is eventual top-step1 for EVERY coprime pair proved here.

## 6. A useful general frontier coordinate

For any finite synchronous state S with at least two values, put M=max S,

    D={M−x:x in S},
    g=min(D minus {0}).

Then the next maximum is `M'=2M−g`, and its exact deficit set is

    D'={d+e−g:d<e in D} union {M−g+d:d in D}.                (S10)

The second term retains old vertices. This is a lossless change of coordinates, not a probabilistic approximation. It turns exponential frontier motion into a subtractive update on small gaps, and makes the previous parity theorem natural: an opposite-parity deficit drops by2 while the frontier's local step remains2. Arbitrarily trimming D without a dependency buffer is unsafe because a deficit just beyond the cutoff can re-enter after subtraction by g.

Source/target/map: synchronous value set -> translated deficit set, x maps to M−x. Preserved: exact labelled state if M is retained, distinct pairs, and next-state evolution. Lost without M: absolute labels and the contribution of retained older vertices. Cheapest hostile: `{2,101}` versus `{2,3}`; both have normalized gcd1 and eventual cofinite closure but dramatically different onset times.

## 7. What is new, inherited, finite, and open

**New proved development here:** (S1)-(S3), their explicit final-hole/genus interpretation; the sum-to-product closure connection; (S6) applied to the classified general startup, with its certificate-storage interpretation; the exact parity-frontier family (S7)-(S9); and the lossless frontier coordinate (S10).

**Inherited or standard mechanisms, not novelty claims:** operation fibres and diagonal parity, matching structure implicit in the fibre definition, the deletion formula (S4), ordinary two-generator semigroup arithmetic, and independence of disjoint Bernoulli pairs. Proofs were included to make the extension self-contained rather than to audit old canon.

**FINITE-EXACT:** `experiment.py` checks all two-seed pairs1<=a<b<=16 against a numeric-birth closure through500, with expected scaled holes and genuses; checks (S4)-(S6) for the coprime a>=2 cases; constructs synchronous states only through t=8 for b in{3,5,7,9,11,15,21}; and checks the maxima and parity laws there. The b101,T53 value follows from the proof, not a giant enumerated state.

**OPEN:** exact synchronous stabilization time for arbitrary coprime a,b; whether three consecutive frontier values always eventually appear; a formula for when the ENTIRE interval from the static conductor to the moving frontier has filled; classification with three or more seeds; and useful minimal witness selections that simultaneously preserve robustness, generation depth, and small storage. These are genuinely different objectives; the cheap one-pair-per-vertex DAG solves only existence.

Reproduce from repository root:

    python3 05-knowledge/results/session259_arithmetic_graph_growth_20261111/summand/experiment.py

The scope is modest by design: the session's priority is new connections, not a large re-verification census.
