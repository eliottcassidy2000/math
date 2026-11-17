# A common-block twin-prime pair inside the sevenfold spectrum

Changing the two pointed exterior tournaments removes the mod-3 obstruction for the [slope-128 symmetric pair](PLUS_RAY.md#symmetric-gap-and-why-it-cannot-be-a-prime-twin-construction). The new pair has the same nonempty block B and the same contacts, and produces strong counts

    7(129P−h),  7(129P+h),       P=P₂(B), h=H(B).

The two ladders have prime quotients at common heights, but no infinitude of prime pairs is asserted.

## Exact profiles and literal prime pairs

Use the one-exit attachment: exterior D−q beats B, and B beats q. For the [two-run law](PLUS_RAY.md#an-arbitrary-module-plus-reader), the resulting count is `W P+(G−e+J−W)h`. The exterior masks and checked profiles are:

| sign | outgoing masks, q=8 | (G,J,e,W) |
|---|---|---|
| − | (286,276,8,402,420,463,31,327,68) | (815,199,118,903) |
| + | (466,388,393,403,294,463,158,16,192) | (817,221,128,903) |

The common coefficient 903 is 7·129. The fixed terms are −7 and +7. The replay checks both assembled graphs directly in the following cases:

| B | (P,h) | prime quotients |
|---|---|---|
| TT₃ | (8,1) | 1031, 1033 |
| TT₅ | (32,1) | 4127, 4129 |
| masks (24,1,3,6,14) | (100,11) | 12889, 12911 |

Since 129=3·43, modulo 3 the two quotients are −h and +h, so both avoid 3 exactly when 3∤h. Changing only the height of TT_m already gives a complete **physical gap-14 family** with center 129·2^m. Exact twin primes in this family need m odd (otherwise 5 divides an endpoint), and if both sevenfold counts are to have exactly one factor of 7 then m must not be 1 mod 3. Together these allow m≡3 or 5 mod 6; m=3 and 5 are the first examples. These are necessary phases, not sufficient prime tests. A prescribed twin pair x,x+2 requires h=1 and hence a transitive block, so it must have x+1=129·2^m, not just a matching residue. This is not the universal selector.

## A physical four-arc path, but not a monotone one

Relabel the plus exterior (list old vertices in new order) by (3,1,7,4,2,0,6,5,8). The exteriors now differ on exactly four arcs:

    {0,4}, {0,5}, {4,7}, {5,7}.

All four avoid the distinguished vertex and B. Both exteriors contain the common path 8→2→3→7→0→1→4→5→6. This path, any Hamiltonian path of B, and the two connecting arcs (B→8 and 6→B) make a spanning directed cycle. Therefore reversing the four listed arcs in any order stays within strong tournaments, for every B. Four is minimal **between these two marked exteriors**, by an exact check of all 8! root-preserving relabelings, not a lower bound on all witnesses for the counts.

The additive gap alone does not permit a strictly increasing path. With r=P/h≥2, one first exterior flip from the negative endpoint has a count in the closed band between the two endpoints only for the following neighbor choices:

| flip | allowed r |
|---|---|
| {1,6} | [74/11, 8] |
| {3,6} | [48/23, 55/23] |

These intervals follow by checking each of the 36 one-arc neighbors: their coefficients (W,K) satisfy `−7≤(W−903)r+K≤7` only in the displayed cases. At r=8 the first choice remains at the lower count, and above r=8 neither is in the band. Thus for r≥8 there is **no strictly increasing path of exterior single-arc edits with the attachment fixed** from the lower to the upper endpoint, even if longer detours are allowed. The four-edit physical path must overshoot or go backwards in count. In particular this distinction holds at both explicit twin-prime heights. This does not constrain changes inside B or changes of its contacts.

## Upgrade to the old factor-seven path

The [published five-arc multiplicative route](PLUS_RAY.md#the-same-block-gives-the-unscaled-count-five-arc-reversals-away) is optimal between its marked exterior endpoints even after relabeling: the marked vertex has outdegrees 1 and 3, so at least two incident edits are necessary, and the final deletion has three edge-disjoint directed triangles on vertex sets {1,6,7}, {3,4,7}, {4,5,6} while the initial deletion is transitive. At least three further edits internal to the deletion are necessary.

Among the 5! orders of the five fixed differing arcs, exactly 29 give a strictly increasing path for **every** B. Each increment is ΔW·P+ΔK·h. Such an increment is positive on every B exactly when ΔW≥0 and 2ΔW+ΔK>0: P/h≥2, the singleton attains 2, and transitive blocks let the ratio grow without bound. The best order can be taken as 1,4,3,2,5 relative to the previously displayed order; it raises the guaranteed smallest increment from 26h to 60h. The replay checks the 32 exterior reads and all orders. This maximin statement is about these five edits, not all possible factor-seven paths.

Replay: `node replay_twin_seam.js`.
