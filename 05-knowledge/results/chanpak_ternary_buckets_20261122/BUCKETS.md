# An intrinsic Collatz bucket for each tournament vertex

**PROVED.** Write \(A_\sigma(x)=(3x+\sigma)/2^{v_2(3x+\sigma)}\) for odd positive \(x\) and \(\sigma\in\{-1,1\}\). Give **every** positive integer \(u\) the bucket
\[
 b_\sigma(u)=v_2(3\,\operatorname{oddpart}(u)+\sigma).
\]
For the shortened Collatz step, an even source moves within its bucket. An odd source makes one literal directed step from bucket \(b_\sigma(x)\) to bucket \(b_\sigma(A_\sigma(x))\). The intervening halvings preserve the new bucket. Thus a transition between consecutive valuations can be read on one edge of the original functional graph.

A tournament on \(n\ge2\) vertices may be labelled by distinct positive bucket values. To encode it, select one odd-origin step per unordered vertex pair, with direction specified by the two buckets. The unrestricted quotient on the buckets has both directions available for every pair. The recovered tournament therefore includes the selected transitions as data. Which paths it admits, how many orbit pieces it needs, and whether its bucket walk can close on an actual positive integer are three separate predicates.

The inherited [THM-4532, tournament clocks and the Syracuse clock](../../../01-canon/theorems/THM-4532-tournament-clocks-two-sheet-clocks-all-odd-tournaments-and-the-syracuse-clock.md) reads valuations as congruences; [THM-4849, rational Collatz cycles and circulant tournaments](../../../01-canon/theorems/THM-4849-circulant-tournaments-are-rational-collatz-cycles-shift-decimation-reversal-two-base-points-and-the-locking-theorem-behind-35-and-175.md) retains a cyclic word for a circulant. Here the observable is the valuation itself and source integers are positive. The earlier [cycle inverse](../chanpak_odd_square_spines_20261122/CYCLE-INVERSE.md) also distinguishes a finite parity history from an integer cycle.

## Hamilton paths are equal-mass residue classes

**PROVED.** Fix distinct positive labels \(b_1,\ldots,b_n\) and \(A=\sum b_i\). The number of odd residues modulo \(2^{A+1}\) whose first \(n\) consecutive odd-step valuations form a Hamilton path in a given labelled tournament \(T\) is **exactly \(H(T)\)** for either sign. Under natural density on odd initial integers, the probability is \(H(T)/2^A\). Conditional on those valuations listing all the labels exactly once, it is \(H(T)/n!\).

For any word \(w=(w_1,\ldots,w_m)\) of positive valuations, put \(a_i=\sum_{j=1}^i w_j\) and
\[
 C=\sum_{i=1}^m 3^{m-i}2^{a_{i-1}},\qquad
 3^m x=2^{a_m}y-\sigma C .
\]
The unique odd seed residue is
\[
 x\equiv 3^{-m}(2^{a_m}-\sigma C)\pmod{2^{a_m+1}}.
\]
The congruence first forces \(3x+\sigma\equiv0\pmod{2^{w_1}}\); after division by \(2^{w_1}\) the next iterate satisfies the remaining seed congruence and is odd. At the last step the extra high bit forces the final quotient odd. Induction proves all valuations are exact. Conversely, iteration gives the same congruence. Different words of the same length cannot specify the same seed when their valuations differ. Every permutation of the labels has the same \(a_n=A\). Counting those permutations that are paths proves the density and the exact integer assertion, without any distribution hypothesis on infinite Collatz orbits. The two signs give negatives of one another's seed classes at that modulus. A physical arc flip changes this class count by exactly its change in \(H\). [A 3-to-5 flip with common-terminal path exchanges](BUCKET-INTERCHANGE.md)

## Lift any finite word into a specified basin

**PROVED.** For either sign, any positive word and any fixed odd positive \(z\) prime to 3, infinitely many positive initial integers exhibit the word and reach \(z\) one odd step after it. To invert the word, its terminal odd integer \(y\) must have
\[
 y\equiv c=\sigma C 2^{-a_m}\pmod{3^m}.
\]
This single congruence suffices. Reducing first modulo 3 enforces the last backward division \((2^{w_m}y-\sigma)/3\); factoring it out enforces each earlier division in turn. All the quotients are odd. Set \(y=(2^e z-\sigma)/3\), choosing
\[
               2^e z\equiv 3c+\sigma\pmod{3^{m+1}} .
\]
The right side is a unit. The order of 2 modulo \(3^k\) is \(2\cdot3^{k-1}\), since \(v_3(4^j-1)=1+v_3(j)\) by binomial expansion or lifting the exponent. Thus \(e\) can be chosen arbitrarily large in one congruence class. All backward iterates are then positive. Distinct odd positions in the prescribed window give distinct integers for sufficiently large \(e\), since the coefficients of \(2^e z\) at two positions differ by a nontrivial ratio \(2^a/3^k\) with \(k>0\).

For the minus sign this realizes every Hamilton-path seed class inside each of the three known basins, and even with specified final anchor 7. For example the valuation word \(3,4,5,3\) comes from \(x=8706255902574489941491\) and reaches 7 after a fifth odd step of valuation 63. It also occurs with respective final anchors 1 and 17. The contrast between 7 and 21 is exact at the terminal condition: 7 is a permissible anchor, while 21 cannot follow an odd step of either sign because every odd image is a unit modulo three. Under the minus odd map, 7 goes to 5, and \(21\to31\to23\to17\). No classification of all negative Collatz cycles is required. [Known-cycle source, Lagarias](https://doi.org/10.4064/aa-56-1-33-53)

## The minimum number of contiguous pieces retains every Hamilton path

**PROVED.** Let \(r(T)=\max(1,\sum_v \max(0,d^+(v)-d^-(v)))\), where the positive part is taken at each vertex before summation. The \(\binom n2\) selected odd-origin steps can cover each arc exactly once in \(r(T)\) contiguous trails in accelerated time, in any chosen basin above. This number is minimal. Equivalently, in one positive trajectory exactly \(r(T)-1\) unselected connecting *odd steps* suffice and are necessary between the observed arcs. The number of halving steps within buckets is additional.

Every vertex with positive imbalance requires at least that many trail starts. If there is a nonzero imbalance, add edges pairing its unit deficits to its unit surpluses. The completed balanced weakly connected multigraph has an Euler circuit. Cutting at the added edges attains the lower bound. In the balanced case use an Euler circuit directly. The word-lifting theorem realizes the trails and intervening connectors. In particular for odd \(n\), a single trail exists exactly when the tournament is regular, since every imbalance is even. For even \(n\) at least \(n/2\) trails are required, since all imbalances are odd.

Every particular Hamilton path \(P:s\leadsto t\) can occur as a contiguous block in a minimum trail cover. Replace all its edges by one marked edge \(s\to t\), retaining any parallel edge. The imbalance at every vertex is unchanged. For \(n\ge4\) the underlying \(K_n\) minus a Hamilton path stays connected. Any cut separating at least two vertices on each side has size at least \(2(n-2)>n-1\), and a single vertex loses at most two edges. For \(n=2,3\) the remaining nonisolated multigraph is connected by inspection. Expand the marked edge in its optimal cover. Consequently \(H(T)\) also counts the distinct spanning vertex words that can appear as contiguous blocks over **all minimum** overlays of the same tournament in a fixed basin.

## Exactly when a cyclic overlay has translation symmetry

**PROVED.** Let \(T\) be a circulant tournament on \(\mathbb Z/n\) with connection set \(D\), and let \(S=\sum_{a\in D}a\). For a divisor \(e\mid n\), a cyclic Euler tour fixed up to rotation by the translation subgroup of order \(e\) exists precisely when \(\gcd(S,e)=1\). In that case the fixed tours are in bijection with the Euler tours of the quotient directed multigraph, taken up to cyclic rotation with parallel arc orbits distinguished.

Indeed the subgroup consists of multiples of \(n/e\). Each quotient vertex has one outgoing arc per \(a\in D\). Lifting any quotient Euler tour gives total displacement \((n/e)S\) modulo \(n\), independent of the tour. Repeating the lifted word traverses every arc upstairs exactly once iff that displacement generates the subgroup. Conversely a fixed cyclic tour upstairs carries the free action on arcs to rotations of its cyclic order. A fundamental block of length \(n|D|/e\) visits each quotient arc orbit exactly once, and its displacement must generate the subgroup. Different choices of the lifting vertex give rotations of the same tour when it generates. In particular the number of full-translation-fixed tours is \((|D|-1)!\) when \(S\) is a unit and zero otherwise.

At seven buckets the three regular tournament classes have \(175,171,189\) Hamilton paths, respectively, by direct enumeration in one-based rows 274–276 of the [seven-vertex catalogue](../../../session212/cycle-levels/classes7.rows). All these Hamilton paths extend to 21-edge Euler observations by the preceding theorem. The cyclic difference set \(\{1,2,3\}\) has two full-translation-fixed tours; \(H=175=7\cdot25\). The Paley difference set \(\{1,2,4\}\) has none because its sum is 7. Nevertheless its affine group \(x\mapsto ax+b\) with \(a\in\{1,2,4\}\) acts sharply transitively on its 21 arcs. It acts freely on Hamilton paths, giving \(189=21\cdot9\). Thus 7 and 21 occur here as symmetry orbit sizes even though neither is an attainable total height. A single reversal from any regular seven-vertex tournament changes its optimal number of orbit pieces from one to two, since the affected endpoints acquire imbalances \(+2,-2\).

## A positive periodic orbit pays an absolute valuation cost

**PROVED.** If a positive \(3x-1\) odd-step cycle of period \(L\) contains all projected arcs of a tournament \(T\) on distinct valuation buckets \(b_j\), repetitions and other buckets allowed, then
\[
 L\log_2(3/2)>\sum_j(b_j-1)\max(d^+_j,d^-_j).
\]
Indeed a cyclic bucket itinerary visits each bucket equally often as source and destination, so its number of visits is at least the larger of the two required degrees. Summing the excess valuations above one gives the right side as a lower bound. On an actual positive cycle, cancellation of \(\prod x_i\) gives \(2^{\sum b_i}=\prod_i(3-1/x_i)<3^L\). Combining proves the claim. Over all assignments of distinct positive bucket labels the best lower bound uniform in the labelling is obtained by assigning \(1,2,\ldots,n\) in decreasing order of the degree maxima.

For seven buckets, all degree maxima are at least three, so even a cycle with many extra observations would need **at least 108 odd steps**. A regular 21-arc Euler itinerary using labels 1 through 7 uses 84 halvings. It cannot close on a positive integer under either sign: for the plus sign \(\prod_i(3+1/x_i)\le4^{21}<2^{84}\), and the minus bound is stronger. A projected Euler circuit in a transient simply returns to the initial bucket at a different integer.

**FINITE-EXACT application.** In the two [118-vertex near witnesses](../chanpak_odd_square_spines_20261122/NEAR-WITNESSES.md), computing the degrees of the four literal tournaments gives the following label-independent numbers. The period bounds round the displayed inequality using integer comparisons of powers of two and three.

| Midpoint offset at exponent 61 | Orientation | Selected arcs | Minimum pieces | Minimum valuation excess over labellings | Necessary positive minus odd period |
|---:|:---|---:|---:|---:|---:|
|−844,710|lower|6,903|3,426|536,362|≥916,917|
|−844,710|upper|6,903|3,426|536,360|≥916,914|
|+184,814,588|lower|6,903|3,422|535,826|≥916,001|
|+184,814,588|upper|6,903|3,422|535,824|≥915,998|

These bounds concern periodic representations of those specific physical hosts. Their transient representations and their Hamilton-path residue classes exist unconditionally; the exact E61 flip and the full Chan–Pak conjecture remain open. Run `python3 -B collatz_buckets.py` from this directory for the seed classes, basin lifts, cyclic tours, seven-bucket Hamilton blocks and 118-vertex invariants.
