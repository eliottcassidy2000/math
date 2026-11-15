# One changed spine edge, and an absolute inverse for a two-ended tail

## 1. A literally minimal edit preserving the full right interface

Changing only one apex arc on a fixed spine cannot retain S, since every vertex's weighted starting sum is positive. The [two-contact switch](TWO_SIDED_RECONTACT.md) therefore has minimum size among apex-only changes. Allow one interior spine edge to move instead, with **the apex mask held fixed**. The following two five-vertex spines differ in exactly the edge 2→0 versus 0→2:
```
A = (16,1,3,7,14),    B = (20,1,2,7,14),
common apex outgoing mask = 22 (vertices 1,2,4).
```
Direct subset DP gives
```
spine A: (h,Q_2,N,S,d)=(9,84,25,41,7)
spine B: (h,Q_2,N,S,d)=(9,84,23,41,7).
```
Thus a **single physical arc reversal**, leaving every apex contact and the continuation data fixed, changes H by two. It is as small as a nontrivial graph edit can be. The unmarked end weights differ, so this is not equivalent under arbitrary left contexts. Both spines are strong and the initial apex has both kinds of contact. A common final painted 0 gives strong seeds at 71,73 and common (h,S,d)=(9,82,0). The exact common-right family is
```
71 + 82t + 9p(t) + {0,2},  t≥0,
```
where p is the [binary suffix correction](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md). This is an arbitrary common-tail +2, not just matching one height. Since gcd(82,9)=1, even either arm has full projected odd-compatible residues by the scalar transfer proof of the [fixed-cycle modular theorem](../chanpak_cross_problem_theorems_20261114/FIXED_CYCLE_CONGRUENCES.md). This application uses only the terminal-zero cut state; no fixed-cycle claim is needed for this unguarded left end. Its harmonic word capacity is still too small for all absolute targets; residue completeness is not the missing step.

## 2. Select a two-ended tail from a prescribed N

The same new packet has a different option, a two-sided recontact keeping *both* interfaces. For the explicit regular five-core switch there, the two tail addresses x,y≥0 give
```
N_i = 41+2i +15[p(x)+p(y)] +63x+42y+120xy+9(x mod 2)+6(y mod 2),
i∈{0,1}.
```
The mixed term is a real interaction: two independently sized tails do not contribute independently to the scalar count. Nevertheless there is an exact O(sqrt(N) log(2+N))-candidate inverse over both tails, not a binary-only search on a predetermined side.

First, N_i≥120xy, so at least one address is ≤floor(sqrt(N/120)); try each such small address on either side. Include x=0 and y=0; they are allowed. With one side fixed, subtract its contribution to get an equation
```
Y = B z + h p(z) + e(z mod 2),
```
where h=15 and B=42+120x or 63+120y, and e=6 or 9. Reject Y<0. For q=2^r≥2 write z=qu+s, 0≤s<q, and use the exact binary split
```
p(qu+s) = p(s)+(q−1−rev_r(s))u+p(u)+δ_s(u mod2),
δ_s=1−highest_bit(s).
```
Enumerate the q values of s. Fixing each also fixes z mod2. Put
```
A=qB+h(q−1−rev_r(s)),
Z=Y−B s−h p(s)−e(s mod 2).
```
For Z≥0, the remaining u must be selected from the interval
```
max(0,ceil((Z−h)/(A+h))) ≤ u ≤ floor(Z/A),
```
then tested by exact equality. These inequalities follow from 0≤p(u)+δ_s(u mod2)≤u+1. Choose q to be the first power of two ≥max(2,sqrt(hY)/B). The length is at most hY/A²+h/A=O(1), and the number of candidates on the fixed-side choice is O(1+sqrt(N)/(x+1)) or the symmetric version with y. Summing this harmonic bound gives the asserted total. Each candidate uses O(log(2+N)) simple bit/arithmetic steps. This returns a physical pair of tails and the appropriate arm or certifies that this particular two-tail family misses the target. The count of candidate tests is not a witness-density assertion.

Replay: `python replay_single_edge_target.py` checks literal one-edge profiles and common guard strength, exhaustively compares the inverse with both-tail forward generation up to a small numerical threshold, and performs selected larger inverse controls.
