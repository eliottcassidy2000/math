# Physical inverse branches on infinite families

The baseline is [the signed descent reduction](../six_tutte_collatz_20261116/COUNTS.md): on count d, the three inverse branches are 2d−1, (4d−3)/3 on d≡3 mod6, and (4d+1)/3 on d≡5 mod6. These results realize each on explicit unbounded families, not all targets or an inductively closed bank.

## Exact one-apex switching table

For a word w, orient the old vertices as a linear tournament and put one apex pointing to precisely the positions carrying 1. Inherited state (h,A,p), initially (1,0,0), updates by
`1:(h,2A+1,1); 0:(h+A+p,2A,0)`.
Here h is the directed Hamiltonian path count, A the word as a binary integer, p the last bit. Each word in the table starts in 1 and ends in 0, so the tournaments are strong. The two words in each row give exactly the desired input and output at every k≥1:

| Lift | Source word | Target word | Input | Output |
|---|---|---|---|---|
| 2d−1 | 10100 1^k 00 | 1 0^(k+1) 100000 | 63·2^k+17 | 63·2^(k+1)+33 |
| (4d−3)/3 | 10 1^k 000 | 1^(k+2) 000 | 21·2^k−3 | 7·2^(k+2)−5 |
| (4d+1)/3 | 1000 1^k 0110 | 1000 1^(k+2) 00 | 81·2^k+5 | 27·2^(k+2)+7 |

For example k=1 gives 143→285, 39→51, 167→223. The identities follow by substituting the repeated-digit powers in the state updates; both repeated-digit matrices after the first digit have dependence only on 1 and 2^k. Direct subset DP in check_inverse.py checks literal orientations at small k. The operation changes the small core, rewiring only apex contacts outside the long repeated segment. It is not an insertion into an arbitrary witness of the same count.

## Reverse bridge: a preserved spacer, and when no finite bridge bank can repair

From the [two-apex theorem](../chanpak_target_ribbon_audit_20261116/TWO_APEX.md), gluing two strong one-apex words u,v with a single reversed arc between their apices and k middle chain vertices has H=F+G·2^k for k≥1, F=h(u)h(v) and G=(A(u)+1)(B(v)+1). Here B(v)=Σ_i(1−v_i)2^i. A singleton module is allowed for this padded formula: no odd cycles within it; the proof counting optional root paths is unchanged. All F,G are odd. Endpoint bounds give F≤G≤4F−3, with the second inequality strict unless both modules are singletons. Indeed each side P lies between h and 2h−1, which gives G≤(2h−1)(2j−1)=4F−2h−2j+1.

The following even lifts preserve the spacer coefficient while switching modules: (F,G)=(171,405)→(341,405) via words (11010,1000)→(empty,110101100); (243,525)→(485,525) via (100010,1000)→(100,1110100). Increase k by one. The outputs are strong, and the long middle interval is untouched other than that one new vertex.

But neither signed lift can be closed on an infinite padded ray by changing to a bounded library of these same modules, even permitting any fixed finite list of ordinal multipliers. To see this, absorb the multiplier q into F,G (inequality G≤4F−3 still holds, equality only for the singletons and q=1). Infinitely many successful k must use one target ray. Equality with (4d+b)/3, b=−3 or 1, then forces k'−k constant for large k of that subfamily, because 2^(k'−k) tends to a positive real constant but lies in a discrete set. Odd coefficients force k'=k+2, G'=G/3, F'=(4F+b)/3. Since F'≤G', this requires 4F+b≤G. For b=1 impossible; for b=−3 equality could arise only from F=G=1, incompatible with 3|G. The one-apex switches above escape by changing intercept signs, not merely choosing a larger spacer.

An **unpadded** two-apex bridge has a different intercept and evades this obstruction: write U=11 0^s (A₀=3) or U=101 0^s (A₀=5), s≥1; then h(U)=A₀·2^s−1, P(U)=A₀·2^s+1. With V=1^m0 one has h(V)=Q(V)=2^m+1. So H=2A₀(2^m+1)2^s−1. Changing (m,s)=(3,s) to (1,s+2) gives **18A₀·2^s−1→24A₀·2^s−1**, exactly the minus lift. This keeps total order: two right-chain ones are removed and two left-chain zeros inserted. The operation demonstrates why the zero-spacer boundary is not interchangeable with the positive-spacer theory.

## Opaque one-port inputs cannot be the universal gadget

Fix input count d. Any fixed tournament template, substituting the witness T at a single module port, reads a nonnegative combination Σ_j a_j P_j(T), where P_j counts ordered covers by j nonempty directed paths. This follows by partitioning Hamiltonian paths by the maximal consecutive runs in T. If a bounded library of such templates must produce a fixed target from *every* witness to d, ordinal-pad T by arbitrarily many sinks without changing d. For j≥2 the corresponding P_j grows unboundedly: retain a fixed whole path of T and divide the padded transitive chain into j−1 nonempty chunks (for j=2, append any proper prefix of the padding to the path of T and take the rest as the other run). Thus any template with a positive coefficient beyond a₁ eventually overshoots, and the others give multiples of d. For targets not divisible by d, no such bounded library works, even if the template is selected after seeing the witness. The result leaves open marked interfaces, unbounded libraries, and stripping the padding before rewriting the actual core.

## A repeatable marked inverse: two of the escape seeds

Fix a root x in T. Let g=H(T−x), and let r,s respectively be the even and odd path-order sums of H(T−V(P)) over every directed path P starting at x (including the singleton in s, and with H(empty)=1). Append a vertex z beaten by all old vertices except x, which it beats. The new pointed state is
`(h,g,r,s) → (h+2r,g,r+s−g,r+s)`.
Indeed, a new odd directed cycle is exactly z followed by an even-order path from x, so the inherited [odd-cycle packing identity](https://arxiv.org/abs/2307.05569) contributes precisely 2r. A root path containing z must end there. Appending z is possible to every old root path except the singleton at x; it reverses its parity. Paths not containing z are counted with neutral sink padding in the complement. Removing x from the new tournament also leaves only sink padding.

Thus r' doubles on each subsequent append, and the tail has intercept h−2s+2g. The state conditions 2r=h−1 and s−r=g are invariant, giving **indefinite** 2h−1 on the very same rooted tournament interface. The following literal masks, bit j meaning arc i→j in row i, work:

| Start | Root | Masks | (g,r,s) | Repeatable counts |
|---|---:|---|---|---|
| 13 | 0 | (12,21,24,18,1) | (3,6,9) | 1+12·2^k |
| 23 | 0 | (24,45,57,48,34,1) | (5,11,16) | 1+22·2^k |
| 41 | 3 | (8,65,3,6,79,31,45) | (9,20,29) | 1+40·2^k |
| 59 | 0 | (48,77,121,113,98,66,1) | (11,29,40) | 1+58·2^k |
| 75 | 0 | (24,77,121,112,98,67,1) | (11,37,48) | 1+74·2^k |
| 95 | 0 | (12,77,120,112,99,67,1) | (11,47,58) | 1+94·2^k |

The branch seed 27 has one exact step 27→53 at root 4 of (62,124,120,112,32,64,17), whose profile is (17,13,22); it does not iterate at that root (the next output is 89, not 105). Contrastingly 13 and 41, with the displayed rooted witnesses, can keep their root forever. This is an exact distinction in continuation data, not a statement about all possible roots at count 27.

A useful residual identity is s=(h+g)/2+e, where e sums H(T−V(P)) over odd-order root paths of size greater than one whose final vertex is beaten by x. All other nonsingleton odd root paths close an odd directed cycle at x. Thus the tail intercept equals g−2e: r controls the immediate step but this distinct defect controls repeatability.

## Why finite dyadic pencils do not solve the recursive problem

A directed library whose states are a·2^k+b with positive integral a and whose edges match affine count identities on all sufficiently large k cannot provide arbitrarily many signed steps without a new generator. The normalized odd part of a is unchanged by an even inverse and divided by 3 by either signed inverse; translates of k only change powers of two. A finite bank therefore runs out of 3-adic valuation. An actual recursion needs either a physical replenishment of that coefficient, an unbounded family of cores, or a more expressive interface.

### All-order obstruction at 27 for this particular operation

In a pointed tournament satisfying 2r=h−1, delete the SCCs strictly before the root component. Their Hamilton product divides h and r because no rooted path can reach those components, so it divides h−2r=1. Those components therefore are singletons (a nontrivial strong tournament has a directed odd cycle and H≥3), and stripping them leaves the full profile unchanged. The next near-sink append makes this trimmed tournament strong. Thus a purported indefinitely repeatable H=27 core would append to a strong H=53 tournament. By Busch's strong order bound as used in the [multiplicity archive](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md), this has at most eight vertices; 53³<5⁸. Its precursor can therefore be checked exhaustively through order seven. The scan from the complete archived order-seven representatives and all their restrictions, checking every possible root, finds no (h,g,r,s) with h=27, 2r=h−1, s−r=g. This excludes repeatability on *any* H27 witness in this near-sink framework; it does not exclude a change of port or a different operation.

For a pointed strong core, the stable intercept condition also bounds the port budget: if d is the root outdegree then binom(d,2)≤e=(g−1)/2 because every two out-neighbours make one open odd rooted path of length three (orient their internal arc). Thus the ready port can touch only O(√g) vertices outward. For fixed g this family is thin: the nontrivial SCC word of T−x has finitely many possibilities by multiplicativity and the strong order bound. If it has t singleton SCCs and more than one SCC, strongness forces x to send to the first and receive from the last. Those two vertices and any even subset of intermediate singletons give odd cycles through x, so H≥1+2^(t−2) for t≥3. The choice of placements for singleton runs and the ≤O(√g) outgoing contacts gives only polylogarithmically many such pointed strong witnesses below X, for each fixed g. Ready interfaces are a genuine recursive gain, but cannot be made universal by fixing a deletion count.

### A harder wall: no retaining extension from 53 to 105, at any order

Both counts are attained (the explicit 27→53 above; 105=3·35 by ordinal sum). But **no H105 tournament contains any induced H53 tournament**. This shows a real inverse sometimes has to alter or discard old vertices or arcs, not merely add more vertices.

For the finite-to-infinite reduction, 53 prime forces a unique nontrivial strong component A of a proposed source, with H(A)=53. In a target it remains within a strong component C. Monotonicity by odd-cycle packs implies H(C)≥53; SCC factorization makes H(C) a divisor of 105, whose only divisor ≥53 is 105. Busch now gives |C|≤9. The archived inventory and the same order bound put a strong H53 core at order 7 or 8 (the lower orders are checked via restrictions of the complete order-seven atlas). Any putative containment therefore occurs within an augmentation of an order-seven or order-eight source. The archived order-eight H53 representatives (all 14, also including decomposable ones) have zero successful one-vertex augmentations. From order-seven H53, the two-vertex scan tests all assignments, pruning only intermediates with H>105 by monotonicity; 7,168 surviving second-stage masks yield zero. Smaller targets can be padded by singleton sources to order nine for these tests without changing H. This proves the stated all-order noncontainment, not merely an obstruction to a preset root. Replay from the archive:

```sh
g++ -O3 anchored53.cpp -o /tmp/anchored53
/tmp/anchored53 ../actual_hamilton_multiplicities_20261021/session24-multiplicities/n7_representatives.txt
g++ -O3 h27_chain.cpp -o /tmp/h27_chain
/tmp/h27_chain ../actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv
```

There is a cost-sharp repair of the first impossible retaining step. From the specific H53 extension of the displayed H27 witness, reverse the old arc 0→1 and append a new vertex whose outgoing mask is 181. The final rows are (188,509,248,496,32,192,401,16,181), and have H105. Thus **exactly one old arc reversal is needed and sufficient** if the source vertices are kept and new vertices may be added. The lower bound applies to every choice of H53 source, and this example attains it with just one extra vertex.

### One step further needs another choice of host

There is a first-hit H105→209 port on a different tournament, with rows (6,12,24,49,99,71,15,383,127) and root 1. Its state is (105,33,52,77). But no H105 host is persistent for the same near-sink operation. As above, trimming any leading singleton components and taking a first hit gives a strong H209 output, of order at most 10 since 209³<5^10. Hence the trimmed input is of order ≤9. The [order-nine replay](ready105.cpp) extends every order-eight atlas class with H≤105 (790 classes; no other parent can produce H105 by monotonicity). It checks every root on all 381 labeled H105 outputs, including repeats, and finds no persistent root. Leading singleton padding includes all smaller inputs without changing the profile. So a first landing again differs from a reusable port even after the first necessary rewiring.
