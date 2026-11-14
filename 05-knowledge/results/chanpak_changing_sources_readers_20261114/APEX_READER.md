# A single seam can change both source and binary reader

Let B be any nonempty tournament, h=H(B), Q=P_B(2) the ordered two-path-cover count allowing empty slots, and W=h+Q. This is the same physical source used by the 67/69 port pair: take C3[r,w,B] in forward cyclic order. Its rooted state is (W,Q/2,h+Q/2) in (H, even, odd) path-sum coordinates. The rooted parity sums follow by either stopping at r or taking r→w and a path through a subset of B. The two parity sums over B both give Q/2: expand the Hamilton path odd-cycle packing count separately on the two sides. Every odd-sized block can change sides, with opposite signs, and B is nonempty. Counting odd cycles through r then gives H=W.

For a reader R with distinguished entry, place the source before R and reverse the single cross arc from r to the entry. Write η=H(R), and e,o for the weighted sums of even- and odd-order paths ending at the entry, using Hamilton path counts of the complementary induced tournament as weights. The seam outputs exactly
[
  eta W + 2o(Q/2) + 2e(h+Q/2)
    =W(eta+e+o)-h(o-e).
]
Proof: a cross odd cycle must traverse the reversed cross arc once and uses a single forward cross arc. Its source run starts at r and its reader run ends at the entry; their lengths have opposite parity for an odd cycle. Thus neither run pair tries to use the reversed arc as its forward arc. The remaining packings factor between the disjoint source and reader subsets, and the cross cycle has weight 2. Packings without any cross cycle produce ηW. No statement about a second port or continuation is needed.

## Binary-address reader and the inverse

Let A(y) be the tournament on a transitive spine with an apex whose outgoing bits from first to last are the binary digits of y (leading zeros harmless). Use R=A(y) reversed, with its apex as entry. In A(y), the weighted paths from the apex have a transitive complement and the possible first spine vertex is any 1-bit. For any such vertex except the last, choosing later vertices divides evenly by parity; the last gives exactly one even path. Including the apex-alone path gives e=ceil(y/2), o=1+floor(y/2). The path count g(y) is governed by
[
 g(0)=1,qquad g(2u+1)=g(u),qquad g(2u)=g(u)+u+(umod 2).
]
Appending a terminal 1 adds a sink relative to the spine paths in apex odd cycles (it cannot return to the apex); appending a 0 allows the old even apex-rooted paths, each giving a new odd cycle of weight 2. Thus (H(A(y))=g(y)).

The actual count of the one-seam construction is
[
 N=W[g(y)+y+1]-h[y	ext{ even}].
]
For y odd, N is the product of W and the odd integer g(y)+y+1>1. The even y lane gives an exact shifted divisor criterion for a **given** physical B: W must divide N+h and the quotient must equal g(y)+y+1 for some even y. The test is not a proof that every N has such a physical source.

## Dyadic checksum phase: changing both physical sides

Choose a width j≥1 and an even suffix v with (0le v<2^{j-1}). Let β be its j-bit reverse complement, set (M=2^j+eta) and (J=g(v)+v-M+1). Put (y=(u2^k-1)2^j+v) for u,k≥1. The prefix (u2^k-1) is obtained by appending k ones to u−1. In the recurrence for g, each suffix zero contributes the current prefix plus a parity correction. Summed across the suffix the prefix coefficient is β and the parity correction from the **last high-prefix bit** is 1 because the high prefix is odd and the first suffix bit is 0. All other suffix contributions are g(v)−1. Therefore the graph realizes
[
 N=K2^k+O,qquad K=WMu,qquad O=W[g(u-1)+J]-h.
]
Fix a common K and choose the *source* W and *reader suffix* M from the divisibility relation WM|K, then choose u=K/(WM). The remaining comparison must be made in the absolute offset O, not just modulo a power of two. This was inspired by the earlier two equal-weight carriers in [dyadic critical-run and cyclic-checksum extraction](../../../01-canon/theorems/THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection.md): strip an inert run and hold the shared carrier fixed, but keep the exact offset rather than merely its parity. There is no probabilistic-independence assumption here.

There is a short exact inverse for a fixed (B,j,v) when the clock exceeds the cofactor, (ule M2^k). First require W|(N+h). For each k≥1 with (M2^kle (N+h)/W-J-1), set (A=M2^k), (L=(N+h)/W-J) and try (u=lfloor (L-1)/A
floor). The recurrence gives (1le g(u-1)le u); hence in this guarded lane dividing after subtracting one recovers u exactly. Retain candidates only if (1le ule A) and (L-Au=g(u-1)). Without the guard the quotient may absorb a carry. In that case enumerate only (L/(A+1)le ule(L-1)/A), by the same inequalities. This is an inverse decision certificate, not a local-residue argument.

### Exact physical pairs and triples at every height

For j=1,v=0, (M,J)=(3,−1). B=TT2 gives W=5, phase u=3, offset 9; B=TT1 gives W=3, phase u=5, offset 11. Their common slope is 45. Thus (45cdot 2^k+{9,11}) for every k≥1, beginning 99 and 101, has two physical witnesses with both source and reader changed.

A three-consecutive-odd window is
[
  1913625cdot 2^k+{67913,67915,67917},qquad kge 1.
]
The rows below are literal outgoing bit masks before optional ordinal padding by one terminal singleton (which keeps h and doubles Q). All graph orientations use the formulas above.

| Outgoing rows of source module B before pad | Pad | h | Q after pad | j,v | M | u | O |
|---|---:|---:|---:|---:|---:|---:|---:|
| transitive triple | 0 | 1 | 8 | 4,6 | 25 | 8505 | 67913 |
| 36,33,26,99,75,84,7 | 0 | 125 | 1576 | 3,0 | 15 | 75 | 67915 |
| 24,73,11,32,78,87,13 | 1 | 123 | 2712 | 3,0 | 15 | 45 | 67917 |

All these witnesses are strong: the cyclic source is strong, and the canonical nonzero even address has first bit 1 and last bit 0, giving a spanning directed cycle through the apex of the reader before reversal. The seam joins two strong tournaments with an arc in each direction. The windows are unbounded in height but sparse by themselves. For a fixed finite source/suffix bank and any fixed required absolute offset O, the odd part of u is bounded: absorb all factors 2 of u into k, then from the recurrence (g(u-1)ge (u+1)/2) for odd u≥1. Thus any attempt to keep a constant finite offset menu in this bank is a finite collection of dyadic rails. Universal prescription needs additional freedom, e.g. varying offsets that grow with phase or an independently changing source/suffix bank of unbounded breadth.

## Short obstructions and a sharp gap

1. **Same-W source switch loses the unit gap.** For B of order at least two, Q is divisible by 4. Indeed, a term from a packing of c odd cycles and z leftover singletons contributes (2^{2c+z}); if c=0 then z≥2, otherwise 2c≥2. Thus if two modules have the same W, their h differ by a multiple of 4, and their count difference under a fixed reader is a multiple of 4. The exceptional singleton has W=3, unique, since each larger module has h≥1 and Q≥4. The obstruction is sharp: outgoing rows (50,92,105,17,68,26,41) and (12,73,34,84,7,27,53) give (h,Q)=(139,1664),(143,1660), W=1803; every even-address reader sees gap 4 and every odd-address reader identifies the two counts.

2. **Most ordinal reader inflations cannot create prime targets.** If the entry is not in the first SCC, the nonempty ordinal prefix precedes the entry's SCC. For every fixed terminal segment of the reader path in its own SCC, sum over the possible prefix subsets that join it. The signed sum is (sum_A(-1)^{|A|}H(P[A])H(P-A)=0) by assigning each odd-cycle/singleton block to the path or its complement. Later SCCs lie entirely in the complement. Thus o=e and the output is divisible by W, with a multiplier ≥2 from η and the nonempty rooted path sum. If the entry is in the first SCC, all later SCCs multiply the seam count by their Hamilton path count; for a prime output they have total count 1 and can be removed. Hence reader search for prime outputs only needs an initial SCC containing the entry with a nonzero imbalance (a singleton is allowed). Even the seven-fold prime regime lets any discarded tail supply only an actual factor of the prescribed target.

The reproducible arithmetic/graph checks are in `replay_apex_reader.py`. The unrestricted prime/7p target problem remains open.
