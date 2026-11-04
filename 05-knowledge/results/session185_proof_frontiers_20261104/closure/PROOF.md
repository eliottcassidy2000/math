# Actual strong-component escapes from both blocked +4 branches

## Scope and result

The S184 SCC-at-most-four theorem completely classified the ordered root-attachment return grammar starting at all-marked TT3, and excluded +4 from either first-generation sibling. This packet constructs both missing +4 returns with actual strong suffix components of orders seven and eight. Each new component additionally has an actual +2 and +4 successor, by replacement with another strong component. These are exact all-clock finite constructions, **not reusable inductive closure**.

We inherit the physical signature (h,p,F,G,J,e,s), ordered-component multiplication, and root attachment from S182/S184. A prefix TT_t beats suffix B, ports W beat root r, and r beats the prefix and nonports. For t>=1 and W meeting the terminal SCC the full tournament is strong, and

    H_t(B,W)=c+2^t F,  c=h−J.

A return from x to y has H_(t+q)(y)=3H_t(x)+d exactly when

    2^q F_y=3F_x,  c_y=3c_x+d.

Only clocks with t>=1 and t+q>=1 are legal. All adjacency masks below give outgoing neighbors of vertex i in row i, and port masks specify actual subsets of those names. No native operation, old-name preservation, or Hamilton-complement construction is asserted.

## Literal strong carriers

The full seven-coordinate signatures, adjacency masks, and port masks are in `witnesses.json`; `verify.py` independently recomputes each coordinate from the literal directed graph.

| Address | Outgoing rows | Port mask | h | F | J | c |
|---|---|---:|---:|---:|---:|---:|
| 24 | 2,20,9,3,45,79,31 | 42 | 35 | 252 | 70 | −35 |
| 44 | 2,4,33,23,71,155,47,95 | 52 | 89 | 504 | 118 | −29 |
| 242 | 2,24,3,85,101,15,167,63 | 138 | 137 | 756 | 240 | −103 |
| 244 | 2,24,67,149,165,15,59,103 | 50 | 281 | 1512 | 382 | −101 |
| 442 | 2,4,201,83,39,15,179,59 | 104 | 171 | 756 | 256 | −85 |
| 444 | 4,9,98,21,7,27,59 | 74 | 65 | 378 | 148 | −83 |

All six suffixes are strongly connected. The old sources are address2=(S1,S1,C1), (F,c)=(21,−13), and address4=(S1,S0,S1,S0,S1), (F,c)=(21,−11). The six new edges have clock shifts

    2→24: −2; 4→44: −3;
    24→242: 0; 24→244: −1;
    44→442: +1; 44→444: +2.

The two arithmetic equalities above verify the identities for every legal clock, rather than just the tested clocks. Both descendants of each new carrier are therefore simultaneously available. A common original root clock t>=4 makes every listed branch legal. At t=4 the literal counts are

    107→323→973→{2921,2923},
    107→325→979→{2939,2941}.

This does not provide the whole third generation: addresses222,224,422,424 are not supplied by this packet. Nor is a return asserted from any terminal address listed here. The result breaks S184's necessary SCC-size boundary, but does not prove unlimited returns or odd-spectrum completeness. In particular clock decreases consume a finite clock budget along these particular finite paths.

## Why the old normalization bound cannot simply be enlarged

There is an all-length failure of any *uniform* bound on normalized v2(F), already within the SCC-at-most-five grammar and even when c is coprime to3.

Let D be the strong five-vertex tournament with rows (24,1,3,6,14), ports {1,2,3}. Its exact signature is

    (11,100,64,49,28,9,6).

Let V be a transitive 20-vertex word whose port bits, from first vertex to last, are the binary digits of A=671088 read low to high. Its exact signature is

    (1,1048576,671088,60357,38626,1,0).

These numbers can be obtained by twenty applications of the singleton multiplication law, with no graph-order search. V's last bit is1. Therefore each B_k=D→V^k, k>=1, begins with a marked SCC and ends with a port; there is no initial unmarked prefix left to normalize away.

Because 25A=16(2^20−1), the ordered-sum F recurrence yields

    F(B_k)=64+100A(1+2^20+...+2^(20(k−1)))
          =64·2^(20k)=2^(6+20k).

Modulo3, the V signature is (1,1,0,0,1,1,0). Each append leaves h=11 and changes c by −11≡1, independently of the remaining coordinates. The seed has c=11−28=−17≡1. Hence

    c(B_k)≡1+k (mod3).

For every k not congruent2 modulo3, c(B_k) is coprime to3 while v2(F)=6+20k is arbitrarily large. Thus the SCC≤4 normalization strategy relying on a bounded dyadic valuation of a marked-first word cannot extend to all SCC≤5 words. This does not rule out a target-dependent bound using the exact constant c, or settle the existence of particular targets with oddpart(F)=63. The companion growing-carrier work studies such a bound.

## Verification, novelty, and limits

Run `python verify.py` and `python -O verify.py`; outputs are identical. The verifier uses explicit exceptions, checks all nine literal source/new states by independent subset DP, all six new suffixes for strongness, and every full root attachment at clocks1 and2 for strongness and its literal Hamilton-path count. It checks all eight edge identities, the common-clock example, and thirty instances of the all-length valuation lemma. The mathematical recurrences above, not finite tests, establish every-clock and every-k assertions.

Discovery enumerated unlabeled tournaments by adjoining a vertex and canonicalizing among equal-score permutations, yielding the usual 12,56,456,6880 types at orders5,6,7,8, then evaluated actual port subsets. The large discovery banks are not part of the proof or publication deliverable. No exhaustive-minimality claim is needed here; any such claim requires its separate census audit. Targeted current-main checks at ad09a24c8fd25eca9ea320d1b240c04990ada96b found no matching B7 row/signature and S184 still states the SCC≥5 boundary without such an escape. Prior S180 paired-family witnesses already existed: this is not claimed as the first pair of clock-uniform digits in all architectures.

The global nonrealizability of7 and21 is inherited, not reproved. Realizability of every other odd count remains open in this work. All new return targets are exact physical tournaments; none is a proof of prime or7p coverage, recursive state closure, or a universal completeness theorem.
