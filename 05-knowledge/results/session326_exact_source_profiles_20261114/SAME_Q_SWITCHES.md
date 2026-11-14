# Two sources with one slot budget: exact target switches

The previous session's bank gives a literal strong output
`99·2^r Q+cH` for each odd digit c=1,…,197 with the phase r listed in its certificate. If two sources share the same slot budget Q and have distinct Hamiltonian counts a,b, write g=gcd(a,b). Whenever b/g and a/g are available digits, swapping the source and its digit preserves the absolute target, not just a residue. Match the phases by adding singleton SCCs to the respective sources. For k at least the larger phase, both cuts then produce

    N_k = 99 Q·2^k + lcm(a,b).

The source H and digit genuinely change. This depends on an actual equal-Q interface, not comparing two H-only witnesses.

## Literal profile pairs

| (a,b) | Q | Source for a | Source for b | digit phases |
|---|---:|---|---|---|
| (5,19) | 256 | (2,4,9,3), followed by three singleton SCCs | (2,20,9,3,45,15) | 19→0, 5→0 |
| (79,91) | 1184 | (6,12,24,97,11,23,55) | (8,5,17,102,43,71,23) | 91→0, 79→2 |

In the first row the order-four core at H5 has Q32 before the three pads. This is inherited from the dyadic bank; the matched-digit application is new here. The second row has two strong sources of the same order (seven), equal Q and different H before any padding. All nontrivial masks in the table are strong. The two primitive rails are 25344·2^k+95 for k≥0 and 117216·2^k+7189 for k≥2. Both demonstrate prime targets: 25,439 (k=0) and 7,509,013 (k=6) respectively, verified by trial division. There is no assertion of infinitely many primes on either rail.

## Multiplicative row, not just one witness pair

Take an ordinal sum of d modules, each chosen from the same fixed (a,b,Q) pair. The feasible source row has a common Q^d and H_i=a^(d−i)b^i for i=0,…,d. Adjacent sources can always switch digits b/g and a/g as before. The exact common outputs at all adjacent seams are

    99 Q^d·2^k + a^(d−i)b^(i+1)/g   (0≤i<d).

Thus the selectable H bank grows at a literally fixed Q for that row, and switching costs only a known neutral phase. These rows do not fill intervals. Fixed a,b,Q give only two unbounded size parameters and the discrete positions in each row; the construction leaves arbitrary-target entry open.

## Rigidity of the round clock comparison

The [small Q atlas](PRIMORIAL_SWITCH.md) also shows that the sole strong class with Q=100 or Q=200 is the Q=100 order-five core (2,20,9,3,13), which is self-converse. Since every nontrivial component contributes at least two factors of 2 to Q, two such components cannot multiply to 100 or 200. Hence every Q=100 tournament is that core, and every Q=200 tournament is the same core with one singleton SCC before or after. In particular the response H=11 really is forced, not merely a feature of one realization. The singleton shift doubles the clock without changing H, whereas switching in the table changes H at fixed clock.
