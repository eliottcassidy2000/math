# A module-inflation repair of an actual eleven-matrix orbit gap

**Proved.** For every integer k≥6, all five successive odd counts
`2^k+1, 2^k+3, 2^k+5, 2^k+7, 2^k+9` occur as Hamiltonian-path counts of strong tournaments.
In particular the largest pure power-plus-three omission 262147=2^18+3 in the Wave-20 eleven-block census has an explicit 19-vertex witness. Other omissions covered here include 2051=2^11+3 and 2053=2^11+5. This is not a proof that the unaugmented matrix orbit contains any of those omitted values. It is an exact repair by an external physical family.

## A two-visit local law

Given a tournament R with marked vertex v, replace v by a transitive module of n≥1 vertices, each having the same orientation to every exterior vertex as v. Each Hamiltonian path projects to a word on the old vertices, with each exterior vertex used once and one symbol v for every nonempty consecutive clone run. Neighboring symbols must be distinct and joined by a directed arc in R. Conversely, for a projected word with j appearances of v, the number of lifts is the number of allocations of n clone vertices to j nonempty ordered bins: every bin traverses its assigned transitive vertices in the unique increasing order. Thus the weight is j! S(n,j), with S a Stirling number of the second kind. This is the special substitution/slot tensor of the repo's marked module constructions.

If v has exactly one outneighbor w, every appearance of v except possibly the terminal one must be followed by w, which is available once. Thus at most two runs occur. In the two-run case the skeleton has shape
`(list P), v, w, (list Q), v`,
where P and Q partition the exterior vertices other than w, and Q is nonempty because w does not beat v. Since every exterior vertex except w beats v, the two-run coefficient is
`C2 = sum_{∅≠Q ⊆ B\{w}} H(B\({w}∪Q)) St(B[{w}∪Q],w)`,
where B=R−v and St is the number of Hamiltonian paths starting at w. The one-run coefficient is C1=H(R). The exact count is
`H(blowup) = C1 + C2(2^n−2)`.

## Five literal generators

Outgoing bitmask rows are indexed from vertex zero. In each graph mark vertex zero.

| Core rows | C1 | C2 | Output |
|---|---:|---:|---|
| (2,4,1) | 3 | 1 | 2^n+1 |
| (2,12,17,37,75,23,47) | 67 | 32 | 2^(n+5)+3 |
| (2,12,17,37,43,71,31) | 69 | 32 | 2^(n+5)+5 |
| (2,12,17,69,43,15,55) | 71 | 32 | 2^(n+5)+7 |
| (2,4,25,19,35,79,31) | 41 | 16 | 2^(n+4)+9 |

The two-row path-count tables are computed by exact subset DP in replay_power_three.py. These seven-vertex cores are strong, and replacing a vertex of a strong tournament by a nonempty tournament preserves strong connectivity. For a given k use n=k in the first row and n=k−5 in the middle three, and n=k−4 in the last. No maximum k is built in. From a prescribed target N, testing the five differences N−d, d∈{1,3,5,7,9}, for a power of two produces the exact module length without scanning tournament sizes.

The +3 and +5 cores differ solely by reversal of the directed triangle on exterior vertices 4,5,6. In this pair the reversal changes C1 by two and preserves the whole C2 histogram by size of Q; the table is an exact boundary check, not a general triangle-reversal invariance claim.

The bridge is not a new count-state congruence or density argument. The small-state orbit can fail at a target even though this anchored substitution reaches it. But adding this family does not by itself prove that the omitted targets have a finite cover by such rays or that the original orbit is cofinite.

## A second infinite repair, through the actual omission 38933

There is also a strong tournament with `H=76·2^n+21` for every n≥1. Use the +7 core above, with its vertices 0 and 5 replaced by transitive modules of sizes n and 2 respectively. In a projected word, let i and j denote their nonempty run counts, and let C_{i,j} be the number of valid words before allocating the clones. The exact array of nonzero counts, rows i=1,2 and columns j=1,2,3, is

    [71  51  10]
    [32  22   4].

The lifted weight factors as i! S(n,i) times j! S(2,j), because the two allocations concern disjoint labeled modules. The second factor is 1,2,0 for j=1,2,3. The effective one-visit count is 71+2·51=173, and the effective two-visit count is 32+2·22=76. Therefore the count is 173+76(2^n−2)=76·2^n+21. The value 38933 occurs at n=9, and the replay counts that actual 16-vertex witness directly by subset DP. The same script computes the displayed array with a word-level recurrence, independent of the literal inflated counts.

## Power-plus-twenty-one, through the actual omission 65557

The two-visit identity also gives `2^k+21` for every k≥8. Use the nine-vertex tournament with outgoing rows

    (130,388,409,275,419,463,415,264,1)

and mark vertex 8, whose only outneighbor is vertex 0. The exact two-visit profile is (C1,C2)=(277,128), giving H=277+128(2^n−2)=2^(n+7)+21. Thus the actual Wave-20 miss 65557=2^16+21 has a concrete 17-vertex witness. Two companion nine-vertex cores, both marked at 8 and with C2=128, are:

| Core rows | C1 | Exact ray |
|---|---:|---|
| (130,420,265,403,391,477,415,260,1) | 269 | 2^k+13 |
| (130,388,465,263,427,399,443,264,1) | 275 | 2^k+19 |

Consequently [2^k+19, 2^k+21] is a second exact two-odd window at all k≥8. It does not fill the intervening offsets between this window and the five-odd window.
