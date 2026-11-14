# Three physically different routes repair the known eleven-block head

**Finite-exact; not a new global coverage range.** The named eleven are
the strict Wave-20 blocks J followed by empty,O,OO,OZ,OOO,OOZ,OZO,OZZ,
OOOZ,OZOZ,OZZZ, with the initial run capped at four. Their inherited
[exact coverage through 2^26](../chanpak_waves_20261113/README.md)
has 408 holes, last 418429. The first coordinate alone determines
membership, but not the remaining choice of reader.

1. Recount the capped head by strict matrix blocks. In a physical state
   0<=d<=H<=W. Each allowed block has H'>=H and H'−d'>=W; also
   W'>=H+W. Thus after two blocks H''>=H+2W>=3H.
   A fixed scalar cap gives a terminating exhaustive word search,
   even if one step is neutral. The possible first-run states are
   (g(z),z+1,1−p(z)) for 0<=z<16. This reproduces the 408 holes.
2. Lift the initial cap, but allow at most two of the same reader
   blocks. For a transitive puncture with root contact word z,
   W=z+1: choose the minimum of a nonempty part contacted by
   the root, and then any subset of later old vertices; also count
   root alone. The alternating sum is 1−p(z), since only a
   contacted last old vertex contributes an uncancelled even-root
   term. The full count is g(z). A fixed reader row (a,b,c) makes
   the exact scalar expression a*g(z)+b*(z+1)+c*(1−p(z)).
   For a nonempty reader b>=1 and a+c>=0, so once z+1
   exceeds the cap its count cannot fall back into the head.
   With no reader strip trailing O's to get g(2v)>=v+1.
   These bounds make the search exhaustive on the absolute target.
   Direct source, up to one reader, and up to two readers fill
   respectively 207, 299 and 363 of the 408 holes. The remaining
   45 have maximum 139201.
3. Use the distinct [13-port-chain alphabet](../kps_s18_ternary_chains_digraph_weights_20261109/README.md),
   not another run of the same semigroup. Its inherited exact
   certificate supplies all odd counts from 961 through 200000.
   The below-961 matrix closure leaves only 7,21,67,163,173
   among these 45. Finally, literal one-apex extensions of
   nontransitive punctures supply the last three admissible counts;
   their rows and independent path enumeration are in head_hybrid.py.

This is a constructive model of how the strategies compensate:
the strict orbit handles its normal forms, the free source fixes
the initial-interface bottleneck, and the other port or a nontransitive
puncture fills the residual. It gives no uniform bound on the
tail of the residual beyond the recorded census.

Run g++ -O3 -std=c++17 orbit_head.cpp -o orbit_head then
./orbit_head > orbit_head_misses.txt, and run python3 head_hybrid.py.
Pass a flag to the latter to print the 363 literal word certificates.
