# Two semicircle guards turn a runner orbit into literal tournament flips

**Status: PROVED elementary geometric and arithmetic statements; FINITE-EXACT controls.**

This note gives a direct geometric bridge, rather than compiling an arbitrary finite-state machine into a tournament. Two fixed semicircle guards detect a runner's distance from the origin. At threshold 1/8, the four odd eighths become the four transitions of a Gray cycle. Integer speeds supply a regular-cover winding vector, and a gcd/2-adic rule describes every simultaneous transition. A boundary-aware sweep preserves isolated equality witnesses.

This is not a new proof of the general Lonely Runner Conjecture, a new spectrum theorem, or a realization in the repository's blue/black complement-line metagraph. The resonance condition below is an elementary consequence of the standard opposite-slope local-maximum argument; no literature-priority claim is made for it or for the geometric composition.

## 1. Sources and the connection being made

- [S168 navigation](https://github.com/eliottcassidy2000/math/blob/e92fc100d5e174b839e432deb31ccd1c86202f17/00-navigation/S168.md) and its [cellular-clock construction](https://github.com/eliottcassidy2000/math/blob/e92fc100d5e174b839e432deb31ccd1c86202f17/05-knowledge/results/session168_computational_structures_20261103/clocks/NOTE.md) distinguish scheduled physical arc flips from autonomous unmarked tournament dynamics and distinguish register period from profile observability.
- The earlier [source-aligned 7 by 13 tensor](https://github.com/eliottcassidy2000/math/blob/main/07-reflections/lrc-r5-source-aligned-relation-residue-7x13-spectrum-codex-20260816.md) is a different object: its seven-cell and thirteen-sheet coordinates are retained before contraction. Its physical-current/Fubini obligations are not resolved here. The numerical example with speeds 7 and 13 below is a fresh elementary clock example, not an identification of that tensor with this runner cover.
- The earlier [Lonely Runner spectrum map](https://github.com/eliottcassidy2000/math/blob/main/07-reflections/G-is-the-first-gap-lonely-runner-spectrum-conjecture-prior-work-map-opus-S116.md) distinguishes finite arithmetic tests from the open higher-dimensional spectrum assertions.
- Standard external context: Noah Kravitz, [Barely lonely runners and very lonely runners](https://arxiv.org/abs/1912.06034), defines maximum loneliness and develops the spectrum viewpoint. Opposite-slope maxima at times m/(v_i+v_j) are standard; see also Fan and Sun, [Amending the Lonely Runner Spectrum Conjecture](https://arxiv.org/abs/2306.10417). We use only the elementary definitions and supply the arguments needed here, not either spectrum conjecture.

**Connection contract.** The source is the actual phase vector (v_1 t,...,v_n t) on the n-torus. The target is a marked tournament on n+2 vertices, with two retained guards and n retained runner labels. The map preserves the exact threshold predicate and every scheduled single-guard-arc crossing. It discards the positions inside the four sectors. Retaining the runner-coordinate grouping and oriented event order recovers winding. Forgetting the guard marks can lose even the threshold predicate, including after passing to an unmarked tournament isomorphism class.

## 2. The two-guard detector

Fix 0<delta<1/4 and put two guards at positions

    a=delta,     b=1/2-delta

on R/Z. For a runner at phase x, orient the edge from a guard g to the runner exactly when

    (x-g) mod 1 lies in [0,1/2].                         (1)

At the two semicircle endpoints the guard wins the tie. This is an explicit pairwise convention, not a claim that the usual geometric orientation is unambiguous at coincident or antipodal points. Choose the guard edge A→B. Orient all runner-runner edges by any fixed tournament R, independent of time. These prescriptions define an ordinary tournament T(x_1,...,x_n) for every phase vector.

**Theorem 1 (exact domination detector).**

    ||x_i|| >= delta for every i

if and only if {A,B} out-dominates every other vertex of T: for every runner r_i, at least one of A→r_i or B→r_i holds.

**Proof.** The two closed outgoing semicircles are

    [delta, 1/2+delta] and [1/2-delta, 1-delta].

They overlap because delta<1/4 and their union is precisely [delta,1-delta]. This is the safe set for distance from zero. Runner-runner edges do not enter the test. ∎

Two guards are minimal within this architecture: one semicircle bit can accept a semicircle, its complement, the entire circle, or the empty set, none of which equals an arc of length 1-2delta strictly between 1/2 and 1. The minimality assertion concerns fixed semicircle-bit detectors, not arbitrary predicates or encodings.

The same construction works for arbitrary real phases and nonzero real speeds. The remaining arithmetic statements assume positive integer speeds and common initial phase zero.

## 3. At one eighth, every runner is a four-state Gray switch

Set delta=1/8. Write a runner's guard incidence as

    ( [A→r_i], [B→r_i] ).

On the open sectors the states are

    00 : (7/8,1) union (0,1/8),
    10 : (1/8,3/8),
    11 : (3/8,5/8),
    01 : (5/8,7/8).

At x=0 the state is also 00. Increasing phase traverses

    00 → 10 → 11 → 01 → 00.                            (2)

Each arrow reverses exactly one actual guard-runner edge. No runner-runner or guard-guard edge changes.

For n independent phase coordinates the entire marked state graph is

    C_4 □ ... □ C_4 = C_4^n ≅ Q_(2n),                  (3)

with 4^n marked tournaments. It is an induced subgraph of the labeled single-arc-flip graph: two members differ by exactly one arc precisely when their incidence vectors differ by one bit. No assertion of injectivity on unmarked classes is made.

The safe states are {10,11,01}^n. Their induced graph is P_3^n, with 3^n vertices. In particular the same four-state circle carries both the threshold predicate and the elementary physical flip mechanism. For uniformly distributed independent phases, the safe set has Haar measure (3/4)^n. The common-start runner trajectory need not be equidistributed in the full torus, so that measure is not an existence proof for it.

## 4. Closed boundaries matter: the exact event state is a union

For positive integer speed v, all guard-edge crossings in 0<=t<1 occur at

    E_v={ (2j+1)/(8v) : 0<=j<4v }.                     (4)

Let T^- and T^+ denote the marked incidence states immediately before and after an event shared by any number of runners. The exact state T^0 at the event satisfies coordinatewise

    bits(T^0)=bits(T^-) OR bits(T^+).                  (5)

**Proof.** Each affected runner crosses one of the four semicircle endpoints. The guard wins an endpoint tie, so its one changed incidence bit is 1 at the event. All unchanged bits retain their value. ∎

Equation (5) gives a literal single-arc schedule that visits the correct equality state: first perform all 0→1 guard flips; record this intermediate tournament; then perform all 1→0 guard flips. Within each group any order works. The intermediate graphs during a simultaneous batch need not arise at an actual nearby runner time, but the recorded midpoint does arise at that time. The total number of flips over one unit-time lap is exactly 4 sum_i v_i.

An isolated threshold witness can occur when one runner enters the safe region at 1/8 while another exits at 7/8. Sampling only open intervals or resolving simultaneous events in an arbitrary order can erase such a witness. The increasing-before-decreasing convention prevents that loss.

Thus a complete exact decision certificate consists of:

1. every event in the finite union of the E_v;
2. the exact closed event states from (5);
3. one rational midpoint of every complementary time interval.

This verifies existence for a given speed tuple. It is not a finite uniform proof for all tuples.

## 5. The regular cover retains the integer speed vector

Number the Gray states in (2) by 0,1,2,3. The grid graph on Z^n, with edges z↔z±e_i, maps to (3) by reducing each coordinate modulo four and then applying the Gray code. It is a regular graph cover with deck group 4Z^n.

Lift the chronological forward event word from z=0, increasing coordinate i by one whenever runner i crosses an odd eighth. Its terminal point at time one is

    (4v_1,...,4v_n).                                  (6)

Therefore the projected closed flip walk retains speed vector v as its deck translation, once its ordered runner coordinates and directions are retained. Different integer speeds can have the same starting and ending marked tournament while carrying different translations.

For n>=2 this grid is not the universal graph cover: it contains square cycles. If every product cube using edges from distinct runner coordinates is filled, including higher-dimensional cubes, while the within-runner C_4 cycles are left unfilled, the base is the cubical torus (S^1)^n; the corresponding filled grid is its usual contractible universal cover. The winding in (6), divided by four, is the standard torus winding vector.

**Important loss boundary.** The graph Q_(2n) alone does not dictate this torus. Taking *all* faces of the physical 2n-bit cube, including the within-runner Gray square and the corresponding higher-dimensional faces, gives the contractible full cube. The torus construction needs the pairing of the two guard bits belonging to each runner, and the instruction to retain each such four-cycle as a circle. The covering and winding are real mathematical structures on this marked model, not an intrinsic topology of the unmarked tournament-class graph.

### 5.1 Why the odd eighths give exactly four clock laps

The *labeled event clock* retains each crossing time and the set of runners crossing, but forgets which guard edge is flipping and its direction. Its period means translation invariance of this labeled point process on the real time line. The marked tournament trace retains every guard incidence, including the endpoint convention.

**Theorem 4 (symmetry jump at one eighth).** Let g=gcd(v_1,...,v_n). At delta=1/8 the labeled event clock has least positive period 1/(4g), while the marked tournament trace has least positive period 1/g. Thus precisely four event-clock laps are required to return the physical trace. For every other delta in (0,1/4), these periods are instead 1/(2g) and 1/g: precisely two laps.

**Proof.** One runner's four phase events are

    {delta, 1/2-delta, 1/2+delta, 1-delta}.

Their circular gaps alternate between 1/2-2delta and 2delta. If these lengths differ, their rotational symmetry group has order two, with least positive period 1/2. If they coincide, delta=1/8 and the four equally spaced events have symmetry group of order four. Scaling by speed v gives least event period 1/(2v) or 1/(4v), respectively. The event labels require a period common to every runner, whose least positive value is 1/(2g) or 1/(4g), by the integer gcd identity.

Each individual guard-incidence bit is a nonconstant semicircle indicator with least phase period one. A period of the complete trace must therefore satisfy v_i P integral for every i, and the least positive solution is P=1/g. The fixed runner-runner and guard-guard arcs add no restrictions. ∎

At eighths, write u=v/g, a primitive integer vector. Over one event-clock lap the open-sector Gray coordinates advance by u modulo four. The corresponding monodromy is

    z ↦ z+u on (Z/4Z)^n,

of exact order four because at least one coordinate of u is odd. At boundary times the OR rule in Section 4 must still be applied: the closed boundary state is not generally equivariant under a quarter-turn Gray relabeling. The regular chamber states and the boundary reconstruction together give the exact physical lift. Away from event times at a generic threshold, one event-clock lap (half a trace period) complements both guard bits in precisely the odd coordinates of u, giving order two.

This is a different four-lap mechanism from S168's signed-unit transport theorem. Here the multiplier comes from forgetting the edge identity and direction in an equally spaced geometric event clock; it is not inferred from Cayley–Dickson multiplication. The common mechanism is explicit retained sheet data over a simpler clock, not an identification of the two transport groups.

## 6. Every synchronization is controlled by gcd and 2-adic valuation

Write nu_2(v) for the exponent of two in v.

**Theorem 2 (exact odd-eighth synchronization).** For positive integers v,w,

    E_v ∩ E_w = E_gcd(v,w)  if nu_2(v)=nu_2(w),
                 empty     otherwise.                 (7)

The same statement holds for a nonempty finite collection: an intersection is nonempty exactly when all 2-adic valuations agree, and is then the event set for the gcd of the speeds.

**Proof.** Put v=ga and w=gb with gcd(a,b)=1. An equality of event times gives an odd integer divided by a equal to an odd integer divided by b. Coprimality forces a and b both odd, exactly the stated valuation condition. In that case the common times are precisely odd multiples of 1/(8g). The same reduction applies to any number of speeds. ∎

There is also an exact event-count formula that avoids inclusion-exclusion. In each valuation group write the speeds as 2^a m_i with m_i odd, and let D_a be the union of their sets of positive divisors. Then

    |union_i E_(v_i)| = sum_a 4·2^a sum_(d in D_a) phi(d).   (8)

Indeed, reduced denominators in that group are exactly 2^(a+3)d with d in D_a. Each contributes phi(2^(a+3)d)=4·2^a phi(d) distinct fractions. Different denominator groups are disjoint.

For speeds 7 and 13, the only common events are

    1/8, 3/8, 5/8, 7/8.

There are 4(7+13-1)=76 distinct events and 80 individual physical flips. The 91=7·13 common refinement and the four simultaneous events are explicit clock arithmetic. They are not the previously studied source-time 7 by 13 Fourier tensor or its missing temporal current.

## 7. A resonance obstruction to exact threshold equality

Let L(t)=min_i ||v_i t|| and ML(v)=max_t L(t).

**Theorem 3 (one-eighth strictness test).** If no pair of distinct speeds satisfies

    8 divides (v_i+v_j)/gcd(v_i,v_j),                   (9)

then every nonempty feasible set {t:L(t)>=1/8} contains a time with L(t)>1/8. In particular ML(v) cannot equal 1/8.

More locally, an isolated feasible time requires a pair satisfying (9), with one runner entering at phase 1/8 and the other leaving at phase 7/8.

**Proof.** At a feasible time, strict inequalities remain strict in a sufficiently small neighborhood. Equality runners at phase 1/8 enter the safe arc as t increases; equality runners at phase 7/8 leave it. If all equality runners are of one type, move slightly in the favorable direction to make every inequality strict. Thus a feasible time that cannot be locally improved has at least one entering speed v and one exiting speed w.

Write v=ga,w=gb with gcd(a,b)=1. Their simultaneous endpoint conditions imply

    vt = integer + 1/8,    wt = integer + 7/8.

Equation (7) gives t=k/(8g) with k odd. Therefore ak=1 mod 8 and bk=-1 mod 8. Since k is a unit mod 8, a+b=0 mod 8, which is (9). ∎

Condition (9) is necessary, not sufficient for tightness, and its absence does not by itself prove the feasible set is nonempty. For example speeds 1 and 7 satisfy it but are simultaneously at distance 1/2 at t=1/2.

The same opposite-slope argument gives the familiar general denominator restriction: if ML(v)=p/q in lowest terms and p/q<1/2, then for some pair i≠j,

    q divides (v_i+v_j)/gcd(v_i,v_j).                  (10)

To see it directly, at a maximizing time one active distance has positive slope v_i and another negative slope -v_j. Their phases are p/q and 1-p/q. With a=v_i/g and b=v_j/g, adding gives (v_i+v_j)t integral and hence the reduced denominator q divides a+b, since gcd(a,a+b)=1. This is a restatement of the standard local-maximum arithmetic, included to locate (9) within prior work.

## 8. Two controls: a tight orbit and complete unmarked blindness

**Tight control.** For speeds 1,2,...,7, the feasible times are exactly

    {1/8,3/8,5/8,7/8} modulo one,

and all four are isolated. Here is a short proof. If every ||jt||, 1<=j<=7, is at least 1/8, the eight points 0,t,...,7t have all pairwise circular distances at least 1/8. Their eight circular gaps sum to one, so every gap is exactly 1/8. Thus t must be an odd multiple of 1/8. Conversely every such t permutes the eight equally spaced points and works. This familiar tight family exhibits precisely the equality witnesses an interval-only sweep would miss.

**Unmarked-blindness control.** Take the internal runner tournament R to be transitive. The all-00 state is the transitive ordinal sum R→A→B. The all-11 state is A→B→R. They are isomorphic unmarked tournaments, yet {A,B} dominates no runner in the first and every runner in the second. Therefore not even the complete unmarked isomorphism class, and hence no Hamilton count or full path-cover profile, determines the threshold predicate without retaining the guards.

These states actually occur on one common-start trajectory with distinct odd speeds: at t=0 every runner has state 00, and at t=1/2 every runner has state 11. Both tournaments have H=1 and the same entire transitive path-cover profile. This is stronger than finding equal scalar counts in unrelated examples. The retained marked geometry, not the unmarked tournament, is where the loneliness information lives.

If strong tournaments are desired, replace either tournament T by the directed triangle substitution T→b→c→T with two new singleton vertices. It is strong because all three nonempty blocks lie on a directed three-cycle, and the two resulting unmarked strong tournaments remain isomorphic. The threshold test still applies to the retained original runner block and guards; it is not domination of the two new wrapper vertices.

## 9. Reproduction and scope

Run

    python3 check.py > results.json
    python3 -O check.py > results.optimized.json
    cmp results.json results.optimized.json
    python3 independent_intervals.py > interval_results.json
    python3 -O independent_intervals.py > interval_results.optimized.json
    cmp interval_results.json interval_results.optimized.json

The standard-library verifier uses exact rational arithmetic and explicit failures rather than removable assertions. It checks the detector on rational grids, (7) for every speed pair through 64, all tuples of at most four distinct speeds from 1 through 12, the total-event formula, chronological Gray lifts, literal one-edge changes, the closed boundary OR rule and its serialized schedule, deck translations, isolated-witness resonance, and the strictness implication. It also checks all 4^n marked cube vertices for n<=4 and literal Hamilton counts of the profile-blind controls through n=7. The general proofs above do not rely on these finite ranges. A separate checker imports no code from the first: it intersects the exact closed safe intervals for all 793 tuples, checks isolated witnesses directly from those components, and exhausts the event and trace-period symmetries for 28 rational-threshold/speed cases. Both scripts agree under normal and optimized Python execution.

The output separately retains the tight seven-speed family, speeds 7 and 13, the divisor family 1,7,13,91, dyadically separated speeds 1,2,4,8, and the earlier spectrum-map example 1,3,4,5,7,13,18. The script checks threshold 1/8, not that last tuple's exact maximum loneliness 3/23.

All tournament vertices and guards remain marked for the physical lift. Runner-runner arcs are deliberately fixed, so the entire tournament is not asserted to be the geometric circular tournament of all runner pairs. Simultaneous flips use an external schedule, and the torus filling uses coordinate sidecars. No autonomous tournament dynamics, blue/black metagraph identification, new universal Lonely Runner bound, physical source-time tensor transport, or complete Hamilton-path spectrum conclusion follows.
