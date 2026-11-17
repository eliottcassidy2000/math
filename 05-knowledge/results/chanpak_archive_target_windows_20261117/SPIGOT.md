# A primitive physical target spigot: radix 128 with one deterministic carry

This is an exact target inverse and nearest-window construction inside a new infinite family of strong tournaments. It includes prime witnesses, but the family's sparsity forbids a universal absolute selector from this bank alone.

## Random archive routes and what actually transfers

An actual random sample of ten non-Chan/non-tournament titles from math/07-reflections, using Random(309071) on the sorted population, included:
* triangular-tower-overlaps-as-repunit-carriers.md: use weighted support carriers instead of assuming a visually special interval is arithmetic.
* lrc14-random031-certificate-spigot-codex-20260629.md: emit a digit only after retaining the exact pending carry; a visible output does not replace hidden boundary debt.
* nc2-spine-audit-fleet-convergence-and-the-remaining-composition-deathstar-S97.md: a separate no-carry and valuation perspective from the Gaussian-moment project.

The transfer here is a literal radix-and-carry identity on physical marked blocks, not a claim that those subjects imply Chan–Pak.

## Two tagged physical digits

In a tournament R tag a set of vertices. Let (J,W,K,E) be its state: J=H(R), W counts ordered path pairs where the second path is nonempty and starts tagged (the first is optional), E counts full paths ending tagged, and K=J−U, where U counts covers whose first path ends and second starts tagged. Append a tagged block by ordinal sum. If its local data are h=H, p=P₂, w=weighted starts, z=weighted ends, e=ordinary ends, s=ordinary starts, u=both tag endpoints, the new state is
    (hJ, pW+Jw, hK−Ju−Es−Wz, Je).
Indeed original/original endpoint location contributes hU to new U, new/new contributes Ju, old end/new start contributes Es, and new end/old start contributes Wz; the contributions are disjoint. For W split on where the second path starts. This is an exact path allocation identity.

Use two literal outgoing-mask graphs:
    digit A: (0,41,11,1,15,29), tag mask 12;
    digit B: (28,17,18,22,0,31), tag mask 43.
They share h=5,p=128,z=41,e=0. Other values (w,s,u,K) are (34,1,8,−3) and (106,5,32,−27), respectively, by independent literal subset DP. In either case w−3u=10. Start with a tagged singleton then a fixed A, which gives (J,W,K,E)=(5,162,−45,0); every additional digit has
    J'=5J,  W'=128W+(34 or 106)J,  K'=5K−41W−uJ,  E'=0.
The carry obeys
    3K'+W'=5(3K+W)+10J.
After m digits this gives 3K+W=5^m(27+10m), independent of the choices. Crucially gcd(5,128)=1, unlike normalized integer radix models with a persistent common factor.

Now append a vertex q beating exactly the tags in the entire ordinal stack, and append a uniformly attached TT2 that beats only q externally, losing to all original vertices. This final block has a single external exit. Its full count is W P₂(TT2)+K H(TT2)=4W+K: paths with one final-block run have 2W+K exterior itineraries, and two runs have W itineraries with 2 ways to allocate the transitive two vertices. Both counts follow by separating whether the final run ends the full path or exits via q. The graph is strong because q beats the first singleton, which reaches every original block and the TT2, and any vertex can reach q through the TT2. Thus
    3N=11W+5^m(27+10m).
Had we used an outer singleton the coefficient of W here would be 5, introducing a common factor. The two-vertex read is essential to the prime-producing version.

## Exact target placement and neighboring picks in this family

For digit word epsilon=(e₀,...,e_(m−1)),
    W=162·128^m +170(128^m−5^m)/123
         +72 Σ_i e_i 5^(i+1)128^(m−1−i).
All weights are integers and strictly superincreasing from left to right: the ratio of each successor to its predecessor is 5/128, so every weight exceeds the sum of all later weights. In particular lexicographic order of words equals exact Hamilton count order.

The bounds 162·128^m ≤ W < 167·128^m, and 5^m(27+10m)≤37·128^m, imply
    594·128^m ≤ N < 625·128^m.
These layers are disjoint. Given the prescribed absolute N, its only possible length is floor(log_128(N/594)); reject if outside the layer bounds. Calculate W exactly from N and the carry. Working backwards, the last digit is uniquely determined, if possible, by
    W_m−5^m(34+72e_(m−1)) ≡ 0 (mod 128).
The residues differ, because 5 is a unit and 72≠0 mod 128. Divide the resulting difference exactly by 128 and repeat; check the terminal W=162. This is an exact decoder in O(log N) arithmetic steps, counting arithmetic on integers of O(log N) bits. A successful decode constructs the tournament on 6m+10 vertices (tagged singleton, fixed A, m digits, q, TT2).

Even when equality fails, an exact predecessor/successor in this bank can be selected in the same order of work. For the preceding value in a layer, subtract the known all-zero count in W units, rounding the rational budget down, and greedily take each descending weight if it does not exceed remaining budget. Superincreasing weights prove this is optimal; the complement method gives the successor. At a gap between layers use the extreme in the adjacent layer. This explicitly places the adjoining **physical** reachable window about an arbitrary target, not merely a residue address. Its gaps are real: each layer has exactly 2^m values, hence the bank has Θ(X^(1/7)) values up to X, far from covering all primes.

Replay replay_spigot.py verifies both literal profiles, strong raw full graph counts at m=0,1, all exact inverses and all between-word brackets through m=6, and trial-certified prime examples 76717 (word 0) and 160923394019 (word 0110). Invariance and decoder completeness do not extrapolate from this range.

## Why a bigger palette of the same type is still insufficient

More generally, assume a fixed finite palette of tagged modules has common h,p,z,e (ordinary count, two-slot count, weighted marked ends, full marked ends). Initialize with J>0,W>0 and E/J=e/h. Put rho=z/(p−h); the transfer gives W'=pW+wJ and normalized carry (K+rho W)/J accumulating one fixed drift per chosen tile. At a fixed return order k≥1, the full count is (2^k−rho)W+J times that normalized carry. Since p≥2h and z≤p−h, the first term is of order p^m after m steps and the carry term is O((m+1)h^m); in particular targets ≤X use length O(log_p X).

If z>0 every palette tile has a tag, so by splitting paths and excluding second-empty covers its digit h≤w≤p−h. There are at most p−2h+1 distinct w digits, strictly smaller than p. The sequence of w digits determines W. Multiple physical modes at the same digit can still vary the carry, but that carry is just a sum of fixed drifts, hence offers at most polynomially many choices in m (determined by the histogram of modes). Consequently this entire fixed-profile strategy gives at most O((log X)^d X^s) distinct counts for some fixed d and s<1. If z=0 there are no tags in the repeated palette, so the w alphabet has one digit and the same conclusion holds. This gives a general reason to change the profile or the grammar, not just add more modes at the same fixed margins.
