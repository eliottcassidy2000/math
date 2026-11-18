# A repeatable paired recursion: the same word gives p and 7p

The earlier [rotation bank](PAIRED_ROTATION_BANK.md) gives matched readers for a common input. This construction goes in the other direction: start from two seeds, then apply the **same arbitrary sequence of rooted attachments** to both. The exact factor seven survives at every stage. This is a graph operation on actual tournaments, not a recurrence on counts without witnesses.

## The rooted state and nonuniform join

For a rooted tournament (A,a), write (h,e,t) for its ordinary count, the number of spanning paths ending at a, and the number of ordered two-path covers whose distinguished path starts at a (the other path may be empty). Attach a rooted tournament (B,b) by an ordinal join A→B except for the single reversed arc b→a. Retain b as the new root.

For B write (g,f,s) for the same three reads; let S be its number of full paths starting at b, R the number of two-path covers whose distinguished path ends at b, and V=f_(B−b)(2), allowing two optional paths (V=1 when the deletion is empty). Then
```
(h',e',t') = (g h − S e + R t,   f h − [|B|=1] e,   s h + V t).
```
In the first read, the paths lost from the ordinal join end the A-part at a and start the B-part at b. A new path has an optional A-prefix, a nonempty B-prefix ending at b, a nonempty A-suffix starting at a, and an optional B-suffix. All transitions except the one reversed arc are forward between blocks, so these are exactly the new paths. Paths ending at the new root cannot use the reversed arc. A lost old ending path exists only if B is a singleton, giving the second formula. For the third, the marked path starts at b. Without a backward crossing all A vertices lie in the ordinary slot; with a crossing the first edge of the marked path must be b→a. The remaining B vertices fill two optional tails, giving the factor V. The old cross-arc into b is not usable in this rooted read.

Every attachment is a homogeneous linear map on the state, so any exact scaling relation between seeds survives arbitrary depth. Conversely three literal one-step probes recover the full state from ordinary Hamilton counts: a singleton gives h−e+t; sink-rooted TT2 and TT3 give h+2t and h+4t. Thus equality of scaled future outputs for all words in this grammar is **equivalent** to scaling these three state coordinates. In particular the third reading is not an optional heuristic.

## Literal scaled seeds and exact target search

The lower seed has outgoing masks (8,5,9,18,7,31,63,127), root 4 and state (11,1,27). The upper has masks (8,17,11,50,69,23,47,127), root 5 and state (77,7,189). The replay confirms both directly. The upper factor is exact; a seven-rotation on these upper witnesses is not asserted.

For prescribed p, search the state graph with a fixed finite attachment palette. Always e≤h, and h'≥h because R≥g≥S and t≥h by splitting a path immediately before the root; in fact h'≥t because gh−Se≥0 and R≥1. Consequently any state from which p is later reached has all three coordinates bounded by p. At the terminal state only h=p is necessary; t may exceed p. A finite search with parent pointers therefore decides whether p occurs in this particular grammar and reconstructs both actual tournaments, without enumerating their isomorphism classes.

With the palette of all rooted tournaments through order four and the seed above, the replay selects words for the prescribed primes 100,003 and 100,019:

- **100,003:** append rooted blocks ((12,9,2,4),1) twice, then ((6,4,0),2), then ((12,9,10,0),2). The upper has exactly 700,021 paths.
- **100,019:** append ((12,9,2,4),1), singleton, ((4,1,2),2), ((2,0),0), then ((12,9,2,4),0). The upper has exactly 700,133 paths.

One-step identities and independently computed subset counts verify both seeds and multilayer controls (for example 227/1589 and 1289/9023). For the larger target words the proof computes both counts by the exact transfers. Both ordinary targets are certified prime by trial division.

## Two ways to use the group analogy

With B transitive and i vertices before its root and j after, the move is
```
h'=h+2^i t−[i=0]e,
e'=[j=0]h−[i=j=0]e,
t'=2^j(h+2^i t).
```
In the slice t=h, a source-rooted transitive block of order m≥2 followed by a sink-rooted one of order n≥2 returns to that slice. With d=h−e the pair update is (d,e)↦(A(d+e),2d+e) for A=2^(m+n−1). Within that branch, h'= (A+2)h−e and d'=Ah give a literal backwards test on a proposed terminal profile. This gives a family of paired recursions with a finite inverse search, rather than just a bank of rays.

More generally a finite attachment palette gives a rational generating function for the sum of the ordinary counts over words of a given total attached order. If M_B is the transfer matrix in the first formula, the series is the first-coordinate read of (I−sum_B z^|B| M_B)^(-1) times the state. The paired seed multiplies the whole series by seven. This weighted series counts with word multiplicity and is not the set of distinct attained heights. The constructive step is the word-to-graph map, not the pole data alone.

**Scope:** deciding reachability in the stated grammar is not proving it reaches every admissible prime. The general task is to choose attachments from the prescribed height, or add complementary grammars to cover what this one misses.

Replay: `python3 rooted_pair.py` alongside the earlier packet scripts.

## Let the target also choose whether to keep the current root

There is a four-coordinate version. Add d=H(A−a), with the convention H(empty)=1. When migrating, the fourth reading is d'=h·H(B−b), since deleting b removes the reversed edge. If instead the root remains a, the exact state update is
```
h'=gh−Se+Rt,
e'=d f,
t'=t f_B(2)−d s,
d'=d g.
```
A path ending a must use b→a as its last edge, after all of A−a and a path of B ending b; this gives df. For the marked-start two-slot read a is already at the beginning of the distinguished path, so the backward arc cannot be used. The old join supplies t f_B(2); it loses precisely the cases in which the marked path starts with a→b, leaving all A−a in the ordinary slot and a marked start read s on B. Removing a restores the ordinal join of the remainder and B. Thus keeping and migrating are both homogeneous linear moves on a common state. Also d≤h because every Hamiltonian path of A−a extends by inserting a. The finite target search argument therefore extends to a palette of all specified port/policy pairs in four coordinates.

Literal seeds, both rooted at 0, are (4,9,18,21,35,15,63,127) and (4,9,66,53,39,135,59,95); their states are (29,11,31,13) and seven times that. One exact mixed word for the prescribed prime 30,013 is:
1. Attach (6,4,0) at root 2, and migrate.
2. Attach (2,4,1) at root 0, and keep.
3. Attach a singleton, and keep.
4. Attach (12,9,2,4) at root 1, and keep.

The lower states are (29,11,31,13) → (153,29,153,29) → (889,29,1749,87) → (2609,87,3411,87) → (30013,87,108543,435). The upper follows the same choices and ends at (210091,609,759801,3045). An independent subset DP checked the literal graphs at all ten stages, through order 19, not just the recurrence. For a compact replay run `g++ -O3 rooted_pair_literal.cpp -o rooted_pair_literal` then `python3 rooted_pair.py --mixed-fixtures | ./rooted_pair_literal`. The separate raw graphs are specified by these fixtures. No universal target coverage or genuine sevenfold automorphism is assumed for these seeds.

From the same seeds, a second exact word gives **30,011 and 210,077**: keep the root after attaching, in order, a singleton, TT3 rooted at its sink, two copies of TT2 each rooted at its sink, and C3 rooted at 0. The final lower state is (30011,13,71949,39). Run the same literal harness with `--fixtures30011`. These adjacent prime targets make different physical choices, while both lift on the same pair of seeds.
