# Counting Hamilton paths by forgetting inactive endpoints

**PROVED.** Order the vertices of a tournament. At the cut after the first \(i\) vertices, call a prefix vertex active if *any* later vertex beats it. Suppose there are at most \(w\) such vertices at every cut. The full Hamilton count is computed by a finite set of nonnegative integer registers depending only on \(w\), updated once per vertex. The method allows one active vertex to be beaten by many later vertices. For fixed \(w\) the number of arithmetic operations is linear in the tournament order given the sequence of cuts and their incidences.

The inherited cutwidth automaton keeps pending backward arcs and path fragments. Here multiple future arcs into the same earlier vertex share one active label. This can bound the register count even when arbitrarily many arcs cross a cut. [THM-4878, graded spine grammars and the cut-state automaton](../../../01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md)

On the prefix consider all choices of arcs forming a spanning union of disjoint directed paths, allowing singleton paths. Each component exposes its head (needing a predecessor) and tail (needing a successor). Keep its ordered pair of endpoint labels, with a separate label for each active vertex and a single generic symbol for all inactive endpoints. Retain the **unordered multiset** of these pairs. If two components already have inactive heads, discard the union: neither head can ever acquire a predecessor, whereas a Hamilton path has only one global start. Thus there are at most \(w+1\) components. There are at most
\[
          \sum_{j=0}^{w+1}\binom{(w+1)^2+j-1}{j}
\]
possible multiset types at any cut. The register of a type counts partial path systems giving that type.

On introducing vertex \(v\), choose no incoming arc or choose one arc from the tail of an existing path into \(v\). Independently choose no outgoing arc or one from \(v\) into an existing head. Every inactive tail can feed \(v\), and no inactive head can be fed by \(v\). For an active endpoint consult its relation to \(v\). If both arcs use the same old component the result would be a directed cycle, so discard it. Otherwise concatenate the chosen paths through \(v\), or make it a singleton if neither was chosen. Update labels to the next active set and prune excess inactive heads. Distinct choices of component occurrences contribute separately, even if their recorded endpoint pairs coincide.

Every resulting arc set has a unique history, obtained by deleting the last vertex and its incident chosen arcs. For a full path, its restriction to each prefix is one of the retained unions: each head except possibly the global start must actually receive a backward arc from the suffix. The final count of single-component systems is therefore exactly the number of Hamilton paths, with both endpoint choices intrinsic. This gives a recursive grammar for Hamilton path constructions on these ordered tournaments in which the only historical information required is which exposed endpoints future vertices can still distinguish.

### A unique star grammar at width one

**PROVED.** For a vertex \(i\) beaten by a later vertex, let \(j_i\) be the last such vertex. Its full activity interval in cut coordinates is \([i,j_i)\). At width one these nonempty intervals are disjoint and may share boundary vertices. Thus the ordered tournament has a unique sequence of clean single cuts and stars: on an interval of length \(l\ge1\), the last vertex must beat its first, and any subset of the \(l-1\) interior vertices may beat its first. All other arcs go forward. Conversely any sequence of those pieces gives width at most one. Gluing occurs at the common boundary vertex of successive pieces, so their lengths sum to the number of cuts, one less than the nonempty tournament order.

Give each backward arc weight \(u\) and each cut weight \(z\). The ordinary generating function over **vertex-ordered** nonempty tournaments of width at most one is
\[
    \frac{1}{1-z-uz/(1-(1+u)z)}.
\]
Indeed a clean cut contributes \(z\) and a star contributes \(u z^l(1+u)^{l-1}\). At \(u=1\) this is \((1-2z)/(1-4z+2z^2)\), giving the counts at orders one through four as 1, 2, 6, 20. Differentiating with respect to \(u\) at one counts ordered-vertex, unoriented flip pairs whose union has width at most one. Mark a backward arc in the union and reverse it to obtain the unique other completion. This distinguishes the one-arc operation before projecting to its Hamilton heights. The all-size, fixed-height ordinary-class series in the prior packet uses a different quotient by relabelling.

### Deciding exact heights at arbitrary orders

For fixed \(w\) and fixed prescribed counts \(h<h'\), the vertex-order generating function of marked physical pairs with those counts is rational if the **union** of the two active sets at every cut has size at most \(w\). The coefficients count vertex-ordered tournaments with the lower orientation first. Different vertex orders of the same ordinary class can contribute. Positivity of the rational series is an exact existence test in this class.

To see this, supply \(w\) slots for the active labels. At each new vertex the finite alphabet records which stored vertices it beats in each orientation, which labels retire, and whether the new vertex receives a slot after retirement. Every omitted old vertex beats the new vertex in both orientations. Permit exactly one discrepancy over the entire word, affecting one earlier vertex at the arrival of its mate. A vertex slated to receive a later backward arc is marked active from its birth and may retire only on a step where it is beaten by the new vertex in at least one orientation. Reject a nonempty live set at the end. Assign the least available slot to each new active vertex. These rules give a unique word for every pair on the fixed vertex order because a label retires at its last backward incidence and the union has no other active labels. Allowing a label to stay on an earlier such step accounts for all possible later incidences. Use the common active labels for both path counters, permitting extra retained labels in one orientation. Pruning with the smaller set of definitely inactive heads is still sound.

Cap each register at \(h'+1\), which preserves all final comparisons since the updates use nonnegative sums. The two capped register vectors, bounded slot information and the edit flag give a finite automaton for the exact output pair. Counting accepted words by length gives the rational series by the finite adjacency matrix. A witness word recovers all physical arcs. This construction applies to arbitrary bounded-width inputs, although its state space can be very large as the height grows. [Ordinary isomorphism-class and marked-edge series via strong components](../chanpak_arc_fiber_grammar_20261122/README.md) · [Exact exponent languages](EXPONENT-LANGUAGES.md)

### Application to a backchain and all its neighbors

In the order defining a backchain, every nonterminal cut is crossed backwards by just the one boundary arc of its gap. Changing a single forward arc adds at most one additional active prefix vertex at each cut. Removing an original backedge can only reduce the active set. Thus the whole single-flip neighborhood of each [133-vertex exact lower witness](BACKCHAIN-LOWER.md) can be counted at width two. The proof permits arithmetic modulo any integer for preliminary screening and exact integer arithmetic for endpoint counts.

**FINITE-EXACT at exponent 61.** Each of the two order-133 hosts has 8,778 possible single-arc reversals. In each host, 107 lower the count, 82 preserve it and 8,589 increase it. The smallest increase is 4,835,850,847,907,644,178,432,000, attained by a unique arc in each host: vertices 70 and 94 with terminal widths 15,61 and vertices 92 and 116 with widths 61,15. Vertex labels start at zero. The respective numbers of neighbors in residue classes 0, 1, 2 modulo three are (6307,1240,1231) and (7871,458,449). The many edges leaving divisibility by three still miss the prescribed upper integer.

Now permute all seven factor words of the lower construction and independently reverse the two nonsymmetric width words. The two copies of width 2 are interchangeable. All 10,080 distinct width words give exact lower counts on 133 vertices. Reversing the entire width word takes the tournament to its converse up to relabelling, preserving both heights and the correspondence of physical arcs. This leaves 5,040 representatives. **None of their 44,241,120 reversals has the upper target**, even modulo \(2^{64}\). The low-six-bit comparison alone has 325,857 matches. This excludes 88,482,240 labelled edits across the original word presentations without asserting that the presentations are ordinary isomorphism classes.

For fast simultaneous counting, compute the ordinary prefix register column and the ordinary suffix linear functional at every cut. For each source \(u\) retain its label while continuing the unmodified prefix. At each potential mate \(v\), evaluate the changed transition and contract with the precomputed suffix at the next cut. All later arcs agree with the baseline. The two full hosts use multiprecision integers, and the exhaustive family test computes the exact indicated residues in unsigned 64-bit arithmetic. Independent subset Hamilton recurrences validate the altered-state computation on 119 small instances including genuine +2 edges.

There is also a shorter exhaustion. After the first positive gap of the second factor, both prefix registers in the [coalescence formula](BACKCHAIN-LOWER.md) are divisible by the first factor height. This common divisor is odd and greater than one. Both completions of a suffix-supported edit therefore share it. Applying the same argument after reversal of the width word forces a +2 candidate to have one endpoint strictly before the first such cut and the other strictly after its counterpart from the right. Among the 5,040 representatives this leaves 3,826,560 arcs. Of these 24,207 pass modulo 64 and none passes modulo \(2^{64}\). The all-arc census independently gives the unrestricted conclusion for this family.

From this directory, on a machine with GMP installed:

~~~sh
c++ -std=c++17 -O3 backchain_flip_sweep.cpp -lgmpxx -lgmp -o /tmp/backchain-flip-sweep
/tmp/backchain-flip-sweep
/tmp/backchain-flip-sweep --factor-bank
/tmp/backchain-flip-sweep --factor-bank-all-arcs
~~~

[Exact neighborhood and all-factor sweeps](backchain_flip_sweep.cpp)

### Which supersingular primes detect the wrong neighbor?

**FINITE-EXACT.** Among all 44,241,120 representative reversals from these *exact lower* witnesses, every choice of at most three of the 15 supersingular primes admits a neighbor matching the desired upper simultaneously. Four primes suffice to exclude every edge. There are 17 excluding four-element subsets. The one with the smallest product is
\[
                       \{13,29,31,47\},\qquad 13\cdot29\cdot31\cdot47=549289 .
\]
Removing one prime from this set in the displayed order leaves 332, 405, 793 and 643 surviving representative reversals, respectively. These are physical arc counts in a specified family. The prime 2 is included in the original 15 and every edge matches automatically at 2. [Arithmetic origin of the 15-prime list and unrestricted modular images](LOCAL-WORDS.md)

The best possible agreement by cardinality is **nine of the 15 primes**. Use the notation \(A=(2), B=(1,0,0,1)\) and the rows of the [prefix word table](BACKCHAIN-LOWER.md). Let \(U\) be the reverse of the 119,743 word and let \(V\) be the displayed 379,127 word. Integers 15 and 61 in the factor list denote one positive gap apiece. Both witnesses have exact lower height \(E_{61}-1\), and both orientations are strong.

| Ordered factors | Arc | Prime set at which the upper matches | Product |
| --- | --- | --- | ---: |
| \(V,U,A,B,15,A,61\) | \((11,105)\) | 2, 5, 7, 11, 17, 23, 29, 41, 59 | 21,120,361,570 |
| \(U,A,A,V,15,B,61\) | \((16,100)\) | 2, 3, 5, 11, 17, 23, 29, 31, 71 | 8,235,855,870 |

Run `/tmp/backchain-flip-sweep --ogg`. It contracts the exact transitions modulo the product of the odd primes, builds a histogram by agreement mask, then sums the histogram over supersets to decide simultaneous compatibility. For the reported witnesses the replay gives the full integer height of the reversed host. On the smaller bank of bounded primitive Farey directions the minimum excluding set has size two. On unrestricted canonical Farey words every finite prime set admits a simultaneous match. These differences identify which absolute constraints give the local tests their force. [Two-prime obstruction on absolute directions](FAREY-HYPERBOLA.md#a-least-set-of-supersingular-primes-after-fixing-absolute-scale)

The doubling construction for the exact lower counts has a specific boundary at prime 17. Both displayed edits occur after the end of the first positive gap. Prepending the doubling factors multiplies both prefix registers at that cut by \(P_j\), hence multiplies both completions by \(P_j\) relative to the corresponding graph with the enlarged last gap. For \(j\ge1\), \(1+2^{60}\mid P_j\) and \(2^{60}=-1\bmod17\). The two readings are therefore both divisible by 17 at those exponents. This explains why these nine-prime witnesses at 61 lose their upper match at 17 upon the first doubling, even though the lower count stays exact.
