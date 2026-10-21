# Exact Hamilton and converse-orbit counts for the exception towers

Status: **PROVED by the bijections below; FINITE-EXACT controls supplied.**
The ordinary walk-content transfer is inherited, not new. The contribution of
this packet is its converse-fixed refinement and efficient application to the
three structural exception families. No literature priority is claimed.

## Inheritance and scope

Closest proved mechanism: THM-3121,
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`,
already proves the full ordinary substitution transfer, the C3 diagonal law,
and H(C3[C3])=3159. Its unordered cover coordinate pc(r) becomes our ordered
coordinate c_r=r!pc(r). The same source explicitly warns that scalar H loses
block-run counts. The corrected near miss is treating bare repeated quotient
words as freely permutable before decorating them with actual disjoint paths.
Our ordered convention avoids that division altogether.

Least-used relevant sidecar here is the reversal-fixed cover vector d, not
another scalar Hamilton statistic. Hostiles: transitive quotients (no repeated
block visits); even numbers of paths under an odd-order anti-involution;
noncanonical block decompositions when claiming full automorphism groups.

External context, checked against primary abstract pages 2026-10-21:
Grinberg–Stanley, https://arxiv.org/abs/2307.05569, relates the Rédei–Berge
symmetric function to Chow's path-cycle function and Hamilton-path theorems.
Irving–Omar, https://arxiv.org/abs/2412.10572, consolidates the matrix-algebra
perspective. Neither abstract was used as a proof of our particular fixed
transfer. No claim that ordinary path-cover identities originate here.

Concept board: maximal block runs; ordered path covers; reflection-fixed
covers; canonical wreath groups; ordinary path orbits; converse pairing.
The signed-subset exception classification supplies the FAMILY to study,
not these Hamilton counts. These formulas are valid independently of that
classification being exhaustive.

## 1. Ordinary transfer

For nonempty B let c_r(B) count ordered lists of r nonempty vertex-disjoint
directed paths covering every vertex. Paths are ordered vertex sequences.
Set c_0=0 and c_r=0 outside 1<=r<=|B|. Thus c_1=H and c_n=n!.
For a loopless tournament A, with adjacency A_ij=1 when i→j, put

    W_A(x)=1^T X (I-A X)^(-1)1,   X=diag(x_i).

This counts nonempty directed WALK words, including repeated letters but no
consecutive repetition. For nonempty factors B_i,

    c_k(A[B_i]) = sum_r [x^r] W_A(x)^k product_i c_(r_i)(B_i).       (1)

Proof: read the ordered global paths in list order, splitting each at each
block change. In each block, list its resulting runs in this encounter order.
The compressed global paths form an ordered k-tuple of quotient walks, and
the run list is an ordered local cover. Conversely decorate successive
occurrences of each quotient letter with successive local paths, then join.
All interblock arcs exist by substitution. Splitting and joining are inverse.
There is no missing factorial and no SCC assumption.

For C3 with arcs 0→1→2→0,

    W=(x+y+z+xy+yz+zx+3xyz)/(1-xyz),
    H(C3[B])=3 sum_(r>=1)(c_r^3+c_(r+1)c_r²+c_(r+1)²c_r).          (2)

These are precisely THM-3121 in ordered coordinates, not a rediscovery claim.

## 2. Fixed-cover transfer: the new retained coordinate

Let B have odd order and let tau_B be an anti-involution. It has exactly one
fixed vertex: two fixed vertices would have their arc both fixed and reversed,
and parity forces at least one. On ordered covers apply tau_B to all vertices,
reverse each path, and reverse the LIST of paths. Let d_r(B) be its fixed count.
Then d_even=0: the unique fixed vertex cannot belong to a pair of exchanged
disjoint paths. The fixed list with 2l+1 paths has a central invariant path.

Take C3[B] with tau(i,v)=(-i mod3,tau_B(v)). Put

    Wt(X,Y)=(X+2Y+2XY+Y²+3XY²)/(1-XY²),
    Vt(X,Y)=(1+Y+Y²)/(1-XY²).

Then the exact full fixed-vector recurrence is

    d_(2l+1)(C3[B])
      = sum_(s,r) [X^s Y^r] Wt(X,Y)^l Vt(X,Y)
                    d_(2s+1)(B)c_r(B).                           (3)

Here 0<=s<=(|B|-1)/2, 1<=r<=|B|. All even entries vanish.

Proof with occurrence ordering explicit: choose the first l global paths;
the last l are forced by reflection. The central path's quotient word has the
form u,0,alpha(reverse(u)), where u is an arbitrary possibly empty walk whose
last letter enters0. Because C3 is a deterministic directed cycle, the halfword
u has content generating series Vt: lengths modulo3 give 1,Y,Y², and a
complete cycle gives XY². The full central quotient word contains an odd number
of0s and equal numbers of1s and2s. Occurrences in the first l words contribute
one X per0 and one Y per1 or2; that series is W(X,Y,Y)=Wt.

Reading ALL global runs in order, the block0 list is a reversal-fixed local
(2s+1)-path cover, giving d_(2s+1). For the swapped block pair, choose the local
cover of block1 in encounter order, giving c_r; the block2 cover is uniquely
forced in reverse order after applying tau_B. In particular both local copies
are covered exactly once. Conversely these choices uniquely decorate and join
all words. This is a bijection, not a heuristic reflection factor. It also
explains why there is no extra power of2 or r!.

Checks: d(K1)=(1); d(C3)=(1,0,2); d_1(C3[C3])=33. The companion independently
enumerates every symmetric vertex permutation and every symmetric cut pattern
at orders3 and9, matching the entire fixed vector, not just its first entry.

## 3. General fixed outer template, including R5 and P7

Let alpha be an anti-involution of odd A, with fixed vertex0 and p swapped
pairs. Give letter0 weight X and each letter in pair j the SAME weight Y_j.
Write Z=diag(z_i) for this assignment, and define

    V_A=1+1^T Z(I-A Z)^(-1)A e_0.                                (4)

This counts the empty halfword and all halfwalks ending immediately before0.
The endpoint factor A e_0 is essential. Set W_A(z) by the same substitution.
The general version of (3) is

    d_(2l+1)(A[B])
     = sum_(s,r_1,...,r_p) [X^s product Y_j^r_j]
            W_A(z)^l V_A(z) d_(2s+1)(B) product_j c_(r_j)(B).     (5)

Proof: identical encounter-order decoration. Every anti-palindromic quotient
word has unique central letter0, since alpha has just one fixed letter.
Repeated0s in the first half are allowed and contribute X. Reflecting a valid
halfwalk creates valid arcs because alpha reverses arcs. Conversely a fixed
compressed word necessarily has precisely this form. This verifies the proposed
generalization, with variable conventions stated unambiguously.

Use R5=Z5{1,2}, P7=Z7{1,2,4}, and alpha(i)=-i. Their base (a,e) are (5,1) and
(21,3), respectively. These may also be checked by enumerating all permutations
of the at most seven base vertices; see independent companion control.

## 4. Orbit normalization and metagraph meaning

Write a(B)=|Aut(B)| and e(B)=|C_Aut(B)(tau_B)|. When the displayed copies of B
are canonical maximal strong modules, the full ordinary group is the wreath
product and

    a(A[B])=a(A)a(B)^|A|,
    e(A[B])=e(A)e(B)a(B)^p.                                     (6)

For the second formula choose a quotient automorphism commuting with alpha.
At its fixed block the local coordinate must commute with tau_B (e(B)
choices). In each swapped block pair one local automorphism is arbitrary
(a(B) choices), and the other is uniquely forced by the commuting equation.
No extra quotient choices occur. These hypotheses hold for the strong towers
and the prime outer templates C3,R5,P7 considered here, including singleton
bottom blocks. Formula (6) is not asserted for arbitrary chosen partitions.

The automorphism group of any tournament has odd order, since an involution
would exchange two vertices and reverse their arc. In the extended group of
a self-converse tournament all anti-involutions are conjugate: Sylow2 subgroups
have order2 and are conjugate; conjugation can be taken in the ordinary group
by multiplying an anti conjugator by its source involution.

Aut acts freely on Hamilton paths. A converse-fixed ordinary path orbit has
exactly one stabilizing anti-involution per path, because any two stabilizers
would yield an ordinary automorphism fixing the ordered vertex sequence.
There are a/e anti-involutions, each fixing d_1 paths; a fixed ordinary orbit
contains a paths. Double counting gives

    h=H/a,    b=d_1/e,    number of converse-paired ordinary orbits=(h-b)/2. (7)

By the inherited sigma-twin theorem, THM-4731,
`01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`,
Statement 3, the fundamental black multigraph at a self-converse node has
degree (h-b)/2, with loops counting twice. Here its tc is h and its blue
multiplicity bm is b. Equation (7) independently establishes the corresponding
converse-paired Hamilton-orbit count. This degree counts incidences with
multiplicity, not distinct neighbors, and does not reconstruct adjacency.

For T_k=C3^[k], |T_k|=3^k, k>=0 (T_0=K1),

    a(T_k)=3^((3^k-1)/2),
    e(T_k)=3^((3^k-1-2k)/4).

These solve a_k=3a_(k-1)^3 and e_k=e_(k-1)a_(k-1); both exponents are integers.

## 5. Exact results

| tournament | H | d1 | h | b | paired-black degree |
|---|---:|---:|---:|---:|---:|
| C3 | 3 | 1 | 1 | 1 | 0 |
| C3[C3] | 3159 | 33 | 39 | 11 | 14 |
| R5[C3] | 197030475 | 7857 | 162165 | 873 | 80646 |
| P7[C3] | 123156134945235 | 5757723 | 2681562805 | 71083 | 1340745861 |
| T3 (27 vertices) | 417382500592116023859 | 9782085285 | 261792936934433 | 40255495 | 130896448339469 |

The additional R5[C3^[2]] order45 control has H=18298233974264660906837453770738731135711615, d1=1836100215435297226725, h=1049576450383154097794041024261123, b=93283555120423575, and paired-black degree524788225191577002255242951918774. Its complete normalization is in r5-45.json; rerun check_r5_45.py normally and with -O.

Full exact values at81 and243, including H,d1,a,e,h,b and paired degree, are
in receipt.json. Complete c,d vectors through81 are included. At243 the
implementation computes the full d243 vector transiently but retains only d1
in the output, alongside H and orbit normalizations. It does not compute the
entire c243 vector.
The ordinary15 count is independently checked by vertex Held–Karp. Fixed H
at15,21,27 is independently checked by a pair-mask vertex DP: choose one vertex
from every tau-pair, order those chosen vertices into a directed path ending
before the unique fixed vertex, and reflect. This uses no quotient transfer.
The21 ordinary count is computed by general quotient-word transfer, not claimed
independently verified by a full vertex Held–Karp census.

## 6. Explicit complexity, receipts, and loss ledger

For C3[B], m=|B|, let P_k=W^k truncated coordinatewise at m. With P_0=1,

    (1-xyz)P_k=(x+y+z+xy+yz+zx+3xyz)P_(k-1).

Visit coefficient triples in increasing coordinate order. Each new coefficient
is seven previous-layer terms plus one smaller coefficient in its own layer.
Two arrays of (m+1)^3 integers suffice. Computing all k<=3m and contracting with
the input c vector uses O(m^4) integer arithmetic operations and O(m^3) integer
storage. The fixed transfer analog uses O(m²) states and O(m³) operations.
These are arithmetic, not unit-bit complexity claims. All coefficients and
counts have polynomial bit length O(m log m) at each fixed-template step, so
the algorithm is also polynomial bit complexity with ordinary integer
arithmetic. Recursive tower costs are geometrically dominated by the final
transfer, giving polynomial time in order, not arbitrary-tournament counting.
Only c81 is needed for H243 by (2); only c81,d81 for d1(243) by (3).

For a fixed q-vertex outer template, denominator det(I-A X) has constant
support depending on q and constant term1. The same denominator recurrence
for W^k gives O(m^(q+1)) arithmetic operations and O(m^q) storage, constants
q-dependent. Equivalently endpoint walk DP computes the k=1 slice directly.
The fixed transfer uses p+1=(q+1)/2 variables. Thus all three exception families
admit exact polynomial-in-order Hamilton/orbit counting; the degree depends on
the fixed outer template. We do not implement a giant seven-dimensional full
cover array: the21 Hamilton slice needs only 4^7 bounded content cells.

Reproduction (from this directory):

    python3 verify_pathcovers.py --output receipt.json
    python3 -O verify_pathcovers.py --output receipt-O.json
    python3 verify_seeds.py > seeds.json
    python3 -O verify_seeds.py > seeds-O.json
    python3 check_independent.py > independent.json
    python3 -O check_independent.py > independent-O.json

Normal and optimized JSON pairs agree byte for byte. Validity gates use
explicit exceptions, not Python assertions. The tower receipt includes its
script SHA256. Inputs are explicit adjacency recipes; no outside numerical
census is silently inherited.

Loss ledger:
- Tournament substitution→quotient walk words forgets actual local paths;
  ordered local cover vectors restore exactly the missing run decorations.
- Ordinary c→Hamilton scalar loses repeat-visit capacity; retain full c for
  recursive composition. No scalar H-only composition claim is made.
- c alone loses converse-fixed information; d and the chosen anti-involution
  restore reflection constraints. All anti-involutions are conjugate, so the
  resulting fixed counts do not depend on which representative was chosen.
- Labelled paths→ordinary orbits loses labels; freeness and a restore counts.
- h→converse pairs loses fixed-orbit data; b=d1/e restores it.
- Degree statistics forget neighbor identities, edge multiplicities, and full
  metagraph adjacency. Those are not recovered here.
- Exception-family support is not multiplicity among ALL tournaments sharing
  a Hamilton count. This packet supplies exact values and orbit multiplicities
  within the classified family, not the global Hamilton-count distribution.
- No claim about the universal7,21 Hamilton gaps follows from these large
  examples alone, and no LRC conclusion is drawn.
