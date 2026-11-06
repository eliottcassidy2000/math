# The free edge-sector input cannot always survive scalar extraction

Status: **PROVED with complete finite certificates and a cited all-order size cutoff; independently audited.** The author checker compares two independent counting algorithms; a separate auditor recomputed all 204,800 final graph counts by direct subset DP. The inference to all tournament orders uses the cited Busch lower bound. No eventual spectrum theorem is claimed.

## 1. Result

Use the literal E(u) from Session202 `positive-overlap/EDGE_SECTOR.md`: its named arc is r->c, its private names form an increasing chain, all private names beat r, and c beats the j-th private name exactly when the corresponding binary digit of u-1 is one.

There is **no tournament T containing an induced copy of E(u)** with the indicated Hamilton-path count, in each of the following cases:

| u | Order of E(u) | H(E(u)) | Forbidden count of T |
|---|---:|---:|---:|
| 9 | 6 | 17 | 55 = 6u+1 |
| 17 | 7 | 33 | 103 = 6u+1 |
| 19 | 7 | 33 | 113 = 6u-1 |
| 35 | 8 | 61 | 209 = 6u-1 |
| 35 | 8 | 61 | 211 = 6u+1 |

In particular the two requested scalar extractors cannot both be realized merely by adding vertices to every E(u), even when every new edge orientation and the number of added vertices depend on the complete target. This strictly exceeds a failure of finitely many fixed graft templates: all vertex-addition recipes preserving the original induced input fail already at u=9.

The conclusion does **not** exclude surgery deleting old vertices or changing old edges, rebuilding from the binary recipe, using a different marked-sector carrier, or construction of 55/103/113/209/211 without the specified induced parent. The requested unrestricted extraction remains OPEN.

## 2. All-order reduction

Two elementary facts are needed.

**Induced monotonicity.** If A is induced in B, then H(A)<=H(B). For a one-vertex extension A+x, insert x into an old Hamilton path: at the beginning if x beats its first vertex, at the end if its last vertex beats x, or at a transition from a vertex beating x to one beaten by x. At least one such position exists. A deterministic choice gives an injection because deletion of x recovers the original path. Iterate for more vertices.

**SCC product.** Strong components of a tournament are totally ordered, with all edges between consecutive components pointing forward. Every Hamilton path concatenates one path in each component. Thus H(T) is the product of their path counts.

The graph E(u) is strong for every u>=2. Hence, if it is induced in T, it lies wholly within one strong component C. Induced monotonicity gives H(C)>=H(E(u)), and SCC multiplicativity makes H(C) a divisor of the target N.

In each row of the table the largest proper divisor of N is less than H(E(u)):

    N=55: 11<17; N=103: 1<33; N=113: 1<33;
    N=209=11*19: 19<61; N=211: 1<61.

Therefore H(C)=N in every case, even if T has arbitrarily many additional components and arbitrarily much source/sink padding.

The cited bound of Arthur H. Busch states that a strong q-vertex tournament has

    H(C) >= 5^((q-1)/3).

Consequently 5^K>N^3 excludes every q>=K+1. The exact integer cutoffs are

    5^8  = 390625  > 55^3  = 166375;
    5^9  = 1953125 > 103^3 = 1092727;
    5^9  = 1953125 > 113^3 = 1442897;
    5^10 = 9765625 > 209^3 = 9129329;
    5^10 = 9765625 > 211^3 = 9393931.

In every row K is precisely two greater than the order of E(u). Thus the strong component C can contain at most two vertices outside the specified induced E(u). Arbitrarily long sequences of births cannot escape this finite cutoff.

## 3. Complete finite certificates

The outgoing rows of the four fixed parents, in the literal names r,c,p_0,..., are

    E(9)  = (2,4,57,51,35,3),
    E(17) = (2,4,121,115,99,67,3),
    E(19) = (2,36,121,115,99,65,3),
    E(35) = (2,68,249,243,227,195,129,3).

For an n-vertex parent, add vertex n with outgoing mask b on the old n names, then add vertex n+1 with outgoing mask c on the first n+1 names. All

    2^n * 2^(n+1) = 2^(2n+1)

possibilities are enumerated, without strongness or isomorphism filters. The exact universe sizes are 8192,32768,32768,131072, totaling 204800 graphs. None has its indicated target count. Graphs with zero or one new vertex are also ruled out: a universal sink can pad them to exactly two new vertices without changing the path count.

The verifier compares the entire ordered vector of counts by two independent algorithms:

1. A subset-DP endpoint and two-path-cover formula for the final insertion. If s_j and e_i count full paths starting at j and ending at i, and Q_ij counts ordered disjoint nonempty paths spanning the old graph with first path ending at i and second path starting at j, the exact extension count for outgoing mask c is

       sum_(j in c) s_j + sum_(i not in c) e_i
         + sum_(i not in c,j in c) Q_ij.

   This follows by deleting the final new vertex from a Hamilton path and distinguishing its first, last, and internal positions. Q_ij is computed by subset convolution of independently computed start and end path arrays.
2. Independently enumerate every permutation of the final n+2 names, discard only permutations violating an old fixed arc, and distribute each remaining path word to every completion of its forced new-edge bits. This algorithm has no subset-DP or two-path-cover step.

Direct full-subset-DP counts are also checked for every one-birth graph and for deterministic final graphs with combined mask divisible by 997. The latter controls are supplemental; the complete comparison of algorithms 1 and 2 covers all 204800 two-birth graphs.

The retained histograms have respectively 261,971,1121,4095 distinct H-values. The observed ranges are 17..653,33..2635,33..2783,61..11219. These are support summaries for the finite extension universes, not unrestricted spectrum intervals.

## 4. Reproduction and citation scope

From this directory:

    g++ -O3 -std=c++17 certify_induced_obstruction.cpp -o /tmp/certify_induced_obstruction
    /tmp/certify_induced_obstruction > CERTIFICATE.txt 2> WORD_COUNTS.txt

The checker regenerates hist_u9.tsv, hist_u17.tsv, hist_u19.tsv, hist_u35.tsv. It uses explicit runtime exceptions; compilation with -DNDEBUG does not remove checks. All counts fit in unsigned 64-bit integers: at most 10! Hamilton words appear in any single graph.

External dependency: Arthur H. Busch, Electronic Journal of Combinatorics 13 (2006), N3, https://doi.org/10.37236/1141. The lower bound is inherited from the audited Session191 `grammar/INDUCED_CORE.md` and `grammar/ADAPTIVE_ENDS.md`; this packet does not claim an independent proof of Busch's theorem or a priority claim for the finite-cutoff method.

The new contribution is that the Session202 freely controlled E(u) sector itself has all-order induced-containment obstructions at both signs of the proposed scalar extractor. The strongest surviving route must change the induced carrier somewhere, rather than merely letting a graft vary with the target.
