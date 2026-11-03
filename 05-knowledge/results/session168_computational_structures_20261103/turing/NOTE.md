# Profile-invisible tournament computers: faithful registers, three-flip tape steps, and the size boundary

**Status.** Exact representation and simulation theorems, with literal finite tournament checks. The new contribution is a uniform, constant-arc-change compiler using the earlier equal-profile switch. Universality of Turing machines and FRACTRAN is classical, not a discovery here. The program, priority order, and microstep controller are supplied explicitly. Nothing below makes the unmarked tournament, its path-cover profile, or the repository's blue/black complement-line metagraph an autonomous computer. There is no claim about the open Hamilton-path spectrum conjecture.

S160 encoded a whole integer by replacing a potentially large transitive block, and separately showed why a finite vector of enabled FRACTRAN guards loses dynamical information. Here the complete registers, or the complete tape and head state, reside in the tournament itself. A fixed-capacity arithmetic instruction changes a bounded number of **actual arcs**, not an entire arithmetic-sized block. Every one of those changes can be invisible to the complete path-cover profile.

## 1. The five-vertex bit and its exact interface

Encode a tournament on vertices 0,...,4 by the lexicographic pairs 01,02,03,04,12,13,14,23,24,34, with a 1 bit meaning that the smaller vertex points to the larger. Put

    B_0 = mask 8,       B_1 = mask 10.

These are the order-five equal-profile edge already recorded in S167's small-order census. They differ only on arc {0,2}. Both contain the directed cycle

    0 → 4 → 3 → 2 → 1 → 0,

so both are strong. They are not isomorphic: each has a unique vertex of outdegree two; the outdegrees of that vertex's two outneighbors are {1,1} in B_0 and {1,3} in B_1.

For any tournament T let p_j(T) count partitions of its vertices into j nonempty directed paths, with order inside each path but no order on the collection. Set

    P_T(z) = Σ_j p_j(T) z^j,
    F_T(r) = Σ_j p_j(T) (r)_j.

Then F_T(r) counts r labeled path slots, allowing empty slots. Direct exact enumeration gives

    P_B(z) = 9z + 33z² + 31z³ + 10z⁴ + z⁵,
    F_B(r) = f(r) = r⁵ + 6r³ + 2r.                    (1)

The two independent algorithms in `check.py` enumerate all subset Hamilton paths and all permutation/cut descriptions, respectively. Thus (1) is a finite certificate, not an experimental conjecture.

For a bit word x of length K, use the ordinal sum

    E(x) = B_(x_1) ⊕ ... ⊕ B_(x_K),                  (2)

with every cross arc pointing from the earlier block to the later block.

**Lemma 1 (literal invisible cube).** The map x ↦ [E(x)] embeds the K-dimensional binary cube inducedly into the simple single-arc-flip graph of tournament isomorphism classes. Its outputs have 5K vertices, and

    F_(E(x))(r) = f(r)^K,         H(E(x)) = 9^K.       (3)

Changing bit j reverses exactly the literal arc {5j,5j+2}, with blocks indexed from zero. For fixed labels, the number of differing arcs between E(x) and E(y) is precisely their Hamming distance.

**Proof.** The strong components are exactly the five-vertex blocks, and their condensation has a unique linear order. An isomorphism therefore preserves every block position and its bit value. Within a block, a reversal can only change that coordinate if the result remains in the family. A reversal between blocks merges the entire interval between its endpoints into one strong component, so leaves the family. This proves inducedness. A path through an ordinal sum consists, in each labeled slot, of a path in the first block followed by one in the second, allowing either to be empty. Hence F_(T⊕U)=F_T F_U, proving the profile identity. For Hamilton paths H multiplies under ordinal sums, giving 9^K. The literal arc statement follows from masks 8 and 10. ∎

This reuses S166's SCC Cartesian-product argument and S167's smaller binary seed; it does not claim a new discovery of either input. The distance statement is exact within this encoded cube and for these fixed labels. It does not exclude a shorter path through other unlabelled tournament classes outside the cube.

## 2. FRACTRAN has bounded arc cost when its registers are retained

Fix an ordered list of positive reduced fractions a_i/b_i, i=1,...,m. Let p_1,...,p_s contain every prime appearing in those fractions and in the inputs under consideration. We initially restrict to P-smooth inputs

    N = ∏_(j=1)^s p_j^(e_j).

An outside-prime factor is invariant under the program, but omitting it would lose injectivity on arbitrary positive integers. It must be retained separately if those inputs are allowed.

Write α_ij=ν_(p_j)(a_i), β_ij=ν_(p_j)(b_i), and δ_ij=α_ij−β_ij. Applicability is exactly e_j≥β_ij for every j. The first applicable rule sends e to e+δ_i.

Choose capacities L_j≥1 and encode each register by a unary prefix

    u_(L_j)(e_j) = 1^(e_j) 0^(L_j−e_j),     0≤e_j≤L_j.

Concatenate these words in the fixed prime order and apply E. Put K=Σ_j L_j. The capacities and prime order are fixed parameters of the coding, not information guessed from H.

**Theorem 2 (bounded-capacity FRACTRAN compiler).** Every macrostate is represented injectively by a 5K-vertex tournament. All these tournaments, and all the microstates below, have the same complete profile f(r)^K. If the first applicable rule i has e+δ_i within capacity, its transition has a literal single-arc path of length

    D_i = Σ_j |δ_ij|.                                (4)

The number D_i depends only on the instruction, not on N or the capacities. The path is shortest in the encoded cube. If no rule applies, the arithmetic program halts. If its selected rule exceeds capacity, the simulator reports **overflow**, not halt, and does not choose a later fraction.

**Proof.** The unary code is injective. A threshold e_j≥β_ij is either automatic when β_ij=0, impossible when β_ij>L_j, or is exactly the value of the β_ij-th bit in that register. Thus the first-applicable rule is determined by the tournament together with the fixed decoder and ordered program.

If δ_ij>0, turn the first δ_ij zeroes into ones in increasing position order. If δ_ij<0, turn the last −δ_ij ones into zeroes in decreasing position order. These changes preserve unary-prefix validity even at intermediate times. The guards ensure no register becomes negative. The capacity test ensures enough zeroes exist. Apply this procedure in prime order. Exactly |δ_ij| bits change in register j, and no changed bit is visited twice, proving (4) and cube minimality. Lemma 1 supplies the actual arc reversals and the full-profile invariance. ∎

### What the controller must retain

At a macrostate, scan the rules in their given order and latch the first applicable index. Then finish its register updates before testing another rule. For a scanner uniform in the capacities, retain the s+1 register boundaries as marked cuts in the SCC order, or use the intrinsic triangle separators of §4 at an extra cost of 3(s+1) vertices. These cuts describe storage layout, not register contents. A rule index, register index, bounded counters up to the largest fixed guard or update threshold, and finitely many scan phases suffice; a marked cursor identifies the current block in the ordered SCC list. This is finite control with positional marks, not an unbounded arithmetic register hidden in the controller. The register magnitudes are the tournament word.

Locating a unary boundary may require traversing many blocks. Therefore (4) is a bound on **arc changes**, not on total scanner time. SCC order, program labels, cursor movement, priority, and the latch are specified operations of this marked model; they do not follow from a blue/black line or from the mere existence of a Hamilton path.

The latch is essential. In the example (3/2,5/3), the step 12→18 can first decrement the exponent of 2, giving an intermediate register value 6, and then increment the exponent of 3. Reselecting the program at 6 would start a new 3/2 instruction before the original one was finished.

### A fully instantiated finite run

With primes (2,3,5), capacities (3,3,3), and ordered program (3/2,5/3), the run is

    12 → 18 → 27 → 45 → 75 → 125 → halt,

with exponent vectors

    (2,1,0), (1,2,0), (0,3,0), (0,2,1), (0,1,2), (0,0,3).

Each macrostep changes exactly two arcs. All six macrostate tournaments and all ten successive microstep tournaments have order 45 and H=9^9=387420489. Their full profile, not merely H, is identical. At 12 both fractions are applicable; the first wins. Swapping them changes the first successor to 20.

In the same capacity box, `check.py` exhaustively controls all 64 exponent vectors: 45 have an in-capacity next step, 15 overflow, and 4 halt. It compares integer divisibility and integer multiplication with the register implementation and verifies every literal arc change. In particular overflow cannot silently become a halting observation.

## 3. A direct Turing tape costs at most three arc flips per binary step

Let a deterministic Turing machine have a finite state set Q of size q, including its halt states, tape alphabet Γ, and moves L or R. Choose an injective b-bit symbol code, where b=ceil(log₂|Γ|). Use a finite tape interval of L cells. For each cell, store:

- its b symbol bits;
- q one-hot head/state bits, one for each state.

Exactly one head/state bit is one across the whole tape. The code of the whole configuration has

    K = (b+q)L

bits, ordered first by cell and then by field. Its tournament is E of that word. No head position is lost on passage to unmarked tournament classes: the SCC order recovers the cell order, while the bits recover the tape, head, and state. The machine, field grouping, and symbol code remain fixed interpretation parameters.

**Theorem 3 (tape compiler).** A transition whose destination cell remains in the interval is represented by at most b+2 literal arc reversals on a tournament of order 5(b+q)L. Every intermediate tournament has the same complete profile f(r)^((b+q)L). For a binary tape the bound is three; more precisely the cost is

    2 + [the written symbol differs from the old symbol].       (5)

For an arbitrary alphabet, replace the indicator in (5) by the Hamming distance between the old and written symbol codes.

**Proof.** Flip the changed symbol bits at the old head. Turn on the new head/state bit in the neighboring cell. Turn off the old head/state bit. Because the head moves left or right, the last two coordinates are distinct and change from 0 to 1 and 1 to 0, respectively. The resulting word is precisely the successor configuration. There are no repeated changed coordinates, so this is also a shortest path within the encoded cube. Lemma 1 proves the tournament assertions. ∎

A microstate after the new head is added has two heads; it is not a legal machine macroconfiguration. With a nonbinary alphabet, partially written symbol codes can likewise be invalid until all their bits are changed. The controller retains its phase and marks the old head cell until the old bit is cleared. Its transition-table entry is finite data. This gives an exact simulation at macrostep boundaries, rather than falsely asserting that every intermediate arc flip is itself a complete Turing transition. A move beyond the represented interval is again **overflow**, not halt.

### A six-step, eighty-vertex demonstration

Take states A,B,H and these four rules, with H halting:

    A,0 → write 1, move R, enter B
    A,1 → write 1, move L, enter B
    B,0 → write 1, move L, enter A
    B,1 → write 1, move R, enter H.

Use cells −2,−1,0,1, initially all zero, with the head at cell 0 in state A. The macrostate trace is

| Step | Tape on −2,−1,0,1 | Head | State |
|---:|:---:|---:|:---:|
| 0 | 0000 | 0 | A |
| 1 | 0010 | 1 | B |
| 2 | 0011 | 0 | A |
| 3 | 0011 | −1 | B |
| 4 | 0111 | −2 | A |
| 5 | 1111 | −1 | B |
| 6 | 1111 | 0 | H |

Here q=3, b=1, L=4, and K=16. The seven macrostates are distinct 80-vertex tournament classes. Six machine steps use exactly sixteen arc flips, and every macro- and microstate has

    F(r)=f(r)^16,       H=9^16=1853020188851841.

The verifier constructs the actual 80-vertex tournaments and checks each changed arc. It does **not** enumerate Hamilton paths on eighty vertices: their profile equality follows from the proved module product. As an independent finite baseline it directly computes the full profiles of all four ten-vertex two-bit tournaments.

The verifier also exhausts all 192 valid four-cell macroconfigurations of this machine: 96 steps remain within capacity, 32 overflow, and 64 are halting states. This is a demonstration of the compiler, not a claim that this small machine is universal or a new busy-beaver bound.

## 4. Where unbounded computation enters

For a fixed finite tournament order and a fixed finite collection of mark types, the state space is finite. The compiled deterministic machine must halt, overflow, or eventually repeat. For the tape model there are at most qL|Γ|^L valid macrostates; for the register model there are exactly ∏_j(L_j+1). This yields a finite repeated-state decision procedure for the capacity-bounded model.

There are two legitimate unbounded interpretations, neither being an infinite computation in one fixed finite tournament:

1. **Finite-run realization.** Every finite Turing computation has some finite visited interval, hence has a realization entirely within one fixed profile fiber. Every finite FRACTRAN run has finite maximum prime exponents, hence has such a realization as well. The needed capacity is not promised to be computable in advance from a potentially nonhalting input.
2. **Growing marked tournaments.** Add new storage modules when capacity is exhausted. This is an explicitly enlarged dynamics allowing vertex insertion, not just arc reversal.

Here is a precise growing register implementation. Separate the s registers using s+1 directed triangles C_3, taking the ordinal sum

    C_3 ⊕ E(u_(L_1)(e_1)) ⊕ C_3 ⊕ ...
        ⊕ E(u_(L_s)(e_s)) ⊕ C_3.                    (6)

The three-vertex SCC separators and five-vertex bit SCCs identify register boundaries intrinsically, even when capacities change. To add capacity, insert a B_0 just before a register's right separator, orienting all cross arcs consistently with the ordinal order. No old arc changes. The represented exponent and arithmetic integer stay unchanged.

Before executing a selected instruction, insert

    g_j = max(0,e_j+δ_ij−L_j)

blank modules in register j. Since e_j≤L_j, one has g_j≤max(δ_ij,0), so any one arithmetic instruction adds at most

    5 Σ_j max(δ_ij,0)

vertices and then makes exactly D_i arc reversals. The number of added vertices per instruction is program-bounded, although each new vertex receives arcs to all old vertices. Thus this is **not** a bounded-neighborhood graph rewrite. A finite scanner can discover a missing blank at a boundary and insert one; it need not store e_j as an external integer.

For (6), writing K=Σ_j L_j and F_(C_3)(r)=r³+2r, the exact profile is

    F(r) = (r³+2r)^(s+1) f(r)^K,
    H = 3^(s+1) 9^K.                               (7)

Every fixed-capacity computation is still profile-invisible. A growth event changes the profile by a factor f(r), and multiplies H by 9, per new blank module. It would be incorrect to assert one constant profile through unbounded growth, since even the polynomial degree records tournament order.

The analogous tape extension appends or prepends one blank cell of b+q bit modules when the head would leave the represented interval. It adds 5(b+q) vertices, and subsequent transitions again use the bounded flip protocol. The stored finite interval is allowed extra blank padding. Coordinates are relative to its left endpoint, so the encoding is injective on allocated configurations in that convention; retaining an origin mark is required to distinguish absolute translations. It is not injective on infinite tapes modulo adding blank margins.

**Computability consequence, with its scope explicit.** Combine the uniform tape compiler with the classical existence of a universal Turing machine, or the growing register compiler with classical FRACTRAN universality. The resulting finite marked instruction system can simulate arbitrary computations on growing finite tournaments. Its halting problem is undecidable. This is a corollary of the inherited universal machine and the proved faithful simulation, not evidence that an unprogrammed tournament operation, a finite-order metagraph, or H alone is universal. The primitive system includes finite program control, ordered-module observation, positional marks, and module insertion.

Conway's original paper gives the first-applicable convention, prime-exponent registers, and a universal FRACTRAN program: J. H. Conway, *FRACTRAN: A Simple Universal Programming Language for Arithmetic* (1987), [primary paper, §§1,6–8 and Theorem 3](https://www.cs.unc.edu/~stotts/COMP210-s23/madMath/Conway87.pdf), [publisher DOI](https://doi.org/10.1007/978-1-4612-4808-8_2). The new tournament compiler above does not change those classical arithmetic semantics.

## 5. Strong tournaments and the observation obstruction

If strong output tournaments are desired at fixed capacity K≥1, apply the inherited S166 wrapper W to the core E(x) of order n=5K. It uses two fixed regular guards of order 2n+5, with cyclic orientation guard A → core → guard B → guard A. The result has 5n+10=25K+10 vertices.

The disjoint outdegree bands identify the core, so the class encoding remains injective and induced. Reversing a core bit arc remains one literal arc reversal. Replacing a module by a same-profile module preserves the whole profile, so all wrapped macro- and microstates have a common profile as well. Their common Hamilton count is **not** asserted to equal 9^K. The tape demonstration therefore also has a strong 410-vertex version; the 45-vertex FRACTRAN demonstration has a strong 235-vertex version. The small wrapper statement is inherited, not newly proved by the finite tests in this packet.

This construction supplies a sharp failure of profile-only state recovery. In the 45-vertex FRACTRAN family, 12 has a successor and 125 halts, while their full profiles agree. In the 80-vertex tape family, the initial and final configurations have the same profile but different tape contents and halting status. No function of the complete profile alone can classify these macrostates correctly. Equality of a compositional counting interface is therefore compatible with different computational states and different next actions.

The transition relation is the explicitly compiled one. Undirected one-flip adjacency includes many choices and does not select a program instruction. The simple flip-class graph used here remains distinct from the blue/black complement-line metagraph throughout.

## 6. Reproduction and inherited sources

Run, using only Python's standard library:

    python3 check.py
    python3 -O check.py
    python3 independent.py
    python3 -O independent.py

`certificate.json` records the seed invariants and profiles, all four full ten-vertex profiles, 256 four-bit fixed-label distance comparisons, aggregate counts for the 64 bounded arithmetic states and 192 bounded tape states checked by the script, the two demonstration traces, actual flipped arcs, and graph hashes. Explicit checks remain active under Python optimization. The independent program recomputes the small profiles by mandatory-cut counting, executes arithmetic and a signed-coordinate tape directly, and replays the certificate's listed arcs from their actual adjacency matrices. It also checks an explicit eleven-to-sixteen-vertex register growth event: all old arcs are preserved, and the following increment flips exactly one arc.

The all-size assertions rest on the proofs above. The finite certificates do not claim that a universal program was exhaustively tested or that an undecidability theorem can be proved by finite sampling.

- [S166 Cartesian flip graphs](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session166_planar_metagraph_20261103/cartesian_flip/NOTE.md): ordered SCC encoding, equal-profile products, and strong wrapper.
- [S167 ternary switch](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session167_recurring_structures_20261103/ternary/TERNARY_CONSTRAINTS.md): the order-five equal-profile edge 8–10 in the complete small-order census, and the distinction between literal labeled flips and quotient relabellings. Our binary seed needs no relabelling between steps.
- [S160 priority scanner](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session160_ordered_programs_20261103/priority_tilings/NOTE.md): first-applicable scanning and the failure of enabled-guard patterns as full states.
- [S160 encoded integer instructions](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session160_ordered_programs_20261103/instructions/NOTE.md): arbitrary ordered arithmetic updates by whole-block replacement, distinguished here from the bounded arc-change construction.

No literature-priority claim is made for the resulting encoding theorems.
