# Independent audit of the three-size uniform-module certificate

**Status:** PASS for the stated candidate proof and its stated universe. The root agent supplied the candidate; this lane independently inspected the combinatorics and ran an exact finite census. This is an audit report, not a priority claim or a replacement for the root's final theorem file.

## Audited statement and basis

Let a fixed tournament context have m external vertices and one nonempty uniform tournament module M. Each external vertex either beats every vertex of M or loses to every vertex of M. The observable is full H, or a path count with prescribed external start and/or external end. No port inside M is pinned.

Let K_j(M) count ordered j-tuples of nonempty vertex-disjoint directed paths covering M. The path decomposition into maximal M-runs gives

`F(M)=sum_(j>=1) A_j K_j(M)`, `A_j` nonnegative integers.

Each A_j counts outside skeletons with j dummy module runs. Empty prefix and suffix outside pieces are allowed; internal outside pieces are nonempty. The order of the module paths is already counted by K_j, so it must not be multiplied by j! again. For a nonempty transitive module,

`K_j(TT_k)=j! S(k,j)` for integers k>=1.

The fixed-size cutoff is sharper than j<=m+1: every internal outside piece starts at a vertex beaten by M and ends at a vertex beating M. These two types are disjoint, so it uses at least two external vertices. Consequently

`A_j=0 whenever 2(j-1)>m`.

## Why the positive support has no gaps above one

Suppose A_j>0 with j>=3. There are at least two internal outside paths. Write O for vertices beaten by M and I for vertices beating M. Each selected internal path begins in O and ends in I.

Two vertex-disjoint directed paths in a tournament have a directed merge preserving the order within each path: while both have unconsumed vertices, choose the current head that beats the other head. A transition within the same path follows its given edge; a switch follows the edge guaranteed when the previous head was chosen. When one path is exhausted, append the other. The first merged vertex is one of the original first vertices and the last is one of the original last vertices. Hence the merged path is again O-to-I.

In the skeleton list, replace one selected internal piece by this merge and omit the other selected piece. Reinsert the alternating dummy module slots. There are now j-1 slots. The prefix and suffix pieces remain unchanged, so prescribed external endpoints are preserved. The two selected internal pieces need not have been adjacent: every interface uses only the O/I type, so moving their interleaved vertices does not damage an interface.

Therefore `A_j>0 => A_(j-1)>0` for every j>=3. This is a positivity implication, not an injective counting estimate. No inequality relating the magnitudes of adjacent coefficients is asserted.

## Additional cutoff and divisibility refinements

The root supplied two refinements after the first audit; both pass the same structural check. With s prescribed distinct external endpoints (s=0,1,2), those endpoints require s vertices in the untouched prefix/suffix pieces. Thus `A_j=0` whenever `2(j-1)+s>m`. The full-count bound is sharp: in `C3[M,TT_(j-1),TT_(j-1)]`, all outside vertices can be used as j-1 middle pairs, each from the second block to the third, and the top coefficient is `((j-1)!)^2` (a matching of the two blocks, followed by an ordering of the pairs).

The internal outside paths can be permuted freely: they all have the same O-to-I interface, and their nonempty disjoint vertex sets make the action free. Hence every A_j is divisible by `(j-1)!`. In particular `A_3` is even, so positive E is at least12 and divisible by12. This argument also respects prescribed external endpoints because the prefix/suffix are fixed. It supplies divisibility, not a comparison between the sizes of different A_j.

## The three-size test and lower discrepancy

Write F_k=F(TT_k). Then

`F_1=A_1`, `F_2=A_1+2A_2`, `F_3=A_1+6A_2+6A_3`.

Thus `E=F_3-3F_2+2F_1=6A_3>=0`. If E=0, the support implication forces every A_j with j>=3 to vanish, and for every nonempty tournament M,

`F(M)=(F_1-2A_2)H(M)+A_2 P_M(2)`, `A_2=(F_2-F_1)/2`,

because `P_M(2)=2K_1(M)+K_2(M)`. The coefficient of H in this changed basis need not be nonnegative; only the K_j basis has that assertion.

If E>0, positivity gives, for every integer k>=1,

`F_k-[F_1+A_2(2^k-2)] >= (E/6)(3^k-3*2^k+3)`.

This is an all-k estimate from the proof, not a regression through finitely many values. The qualification k>=1 matters: at k=0 the naive displayed exponential expression is not the ordered nonempty three-cover count.

## Independent finite census

Run `python -B 05-knowledge/results/mrb_seven_cutoffs_20261113/curves/module_audit.py` (or the same with -O). Both runs pass **26,069 explicit checks**.

Universe: every labelled tournament on m=0,...,4 outside vertices, all 2^m uniform attachments to the module; 1,099 contexts total. Observers are full H, fixed outside start0 when available, and fixed outside start0/end1 when available: 3,293 observers total. Counts use direct subset dynamic programming, independently of the skeleton construction.

For every observer the audit computes F_1,F_2,F_3 and checks A_2 and A_3 are nonnegative integers, E is divisible by12, and the endpoint-refined cutoff holds. For every E=0 observer it checks the claimed arbitrary-module formula on TT4, C3, the strong-four module `[6,12,8,1]`, and the five-core `[22,4,8,3,14]`. Every module is explicitly checked to be a tournament before counting. For every E>0 observer it checks the discrepancy bound at k=4.

| Outside size m | Contexts | Observers | E=0 | E>0 |
|---|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 | 0 |
| 1 | 2 | 4 | 4 | 0 |
| 2 | 8 | 24 | 24 | 0 |
| 3 | 64 | 192 | 192 | 0 |
| 4 | 1024 | 3072 | 2904 | 168 |

The largest observed E is24. The cutoff correctly forces E=0 through three external vertices. Finite success supports the proof audit; it is not used to infer the infinite statement.

**Fixture correction during audit.** The first draft typed the four-vertex control as `[6,4,8,1]`, omitting edge1-to3. The broader path decomposition still passed on that directed graph, so a count-identity check alone did not detect the type error. It was replaced by `[6,12,8,1]`; explicit tournament-type checks were added, and both normal and optimized runs were repeated. All reported counts use the corrected fixture. The candidate theorem and proof were unchanged.

## Boundary of approval

The argument is approved for a single uniform nonempty module, an ordinary tournament outside it, and full/start/end path-count observers with all prescribed endpoints outside the module. It does not by itself prove the same criterion for nonuniform attachments, pinned internal module endpoints, signed differences of positive observers, several simultaneously changing modules, Hamiltonian-cycle conventions, or root-transport predicates. Each such extension needs its own nonnegative skeleton basis and a merge respecting its additional data.
