# Exact colour-resolved size of the strong-core cut

Status: PROVED. The total tiling recurrence is an equivalent reformulation of inherited THM-349 (with THM-354's SCC interpretation). The reflection-fixed recurrence and the complement-line cut formulas below are direct structural corollaries. A bounded source search did not locate these exact colour-resolved cut formulas; no exhaustive priority claim is made.

## Definitions

Use the fixed descending Hamiltonian path and bit convention in FOUNDATIONS_AND_MIXED_LANDING.md. All tilings, not isomorphism classes, are counted here. Set

T_0=B_0=1,
T_n=2^binom(n-1,2),
B_n=2^floor((n-1)^2/4)  (n>=1).

Here T_n counts all fixed-path tilings, B_n counts grid-reflection-fixed (blue) tilings, s_n counts strong tilings, and b_n counts tilings that are both strong and blue. A singleton tournament is strong, so s_1=b_1=1; there is no strong order-two tournament, so s_2=b_2=0. Empty tilings are auxiliary factors, not strong components.

Converse preserves strong connectivity, so strong/nonstrong is well-defined on merged nodes. It is unrelated to the abbreviation SC for self-converse.

## 1. First-component and mirrored-outer-component recurrences

For every n>=1,

T_n = sum_(k=1)^n s_k T_(n-k),

B_n = b_n + sum_(k=1)^floor(n/2) s_k B_(n-2k).

For the first identity, a fixed Hamiltonian path traverses the canonically ordered strong components in order, without returning to a component. Each component is therefore a contiguous interval of path labels. Choose the first component's size k, its strong fixed-path tiling, and an arbitrary tiling on the remaining n-k vertices. All cross-component arcs are forced. These choices are unique and reversible. This is THM-349 at x=1, with c_L=s_(L+1), together with the singleton term and THM-354.

For the second identity, a blue tiling is either strong or has distinct first and last components. Grid reflection induces an anti-automorphism reversing the strong-component order. Thus the first and last components have equal size k; the first is any strong tiling, the last is its uniquely prescribed reflection/converse, and the remaining middle is any blue tiling. Again all cross-component arcs are forced. When n=2k the middle is empty, contributing B_0=1. This proves the identity without enumerating classes or choosing anti-isomorphisms.

As formal power series (not analytic convergence claims), write

T(z)=sum_(n>=0) T_n z^n,  B(z)=sum_(n>=0) B_n z^n,
S(z)=sum_(n>=1) s_n z^n,  b(z)=sum_(n>=1) b_n z^n.

Then

T(z)=1/(1-S(z)),
b(z)=B(z)(1-S(z^2))-1=B(z)/T(z^2)-1.

## 2. Core/periphery support forces exact cut masses

For n>=3 every complement-line partner of a nonstrong tiling is strong. Indeed its Hamiltonian path begins and ends in distinct extremal strong components, so the initial vertex beats the terminal vertex. The path-preserving antipode flips this nonpath endpoint arc, closing the unchanged path to a directed Hamiltonian cycle. This is the inherited THM-333 apex-closing mechanism, now read on the complement metagraph.

Therefore the nonstrong nodes form an independent periphery: no nonstrong loops and no nonstrong--nonstrong lines. Every nonstrong tiling belongs to a distinct crossing line. Since antipode preserves blue/black colour, we obtain the exact table, counting literal line instances including multiplicities:

| Line location | All colours | Blue | Black |
|---|---:|---:|---:|
| Nonstrong--nonstrong | 0 | 0 | 0 |
| Nonstrong--strong | T_n-s_n | B_n-b_n | T_n-s_n-B_n+b_n |
| Strong--strong | s_n-T_n/2 | b_n-B_n/2 | s_n-b_n-(T_n-B_n)/2 |

Every loop lies in the strong core. The table counts each loop once, not twice. Its cut masses are not merely expectations or bounds: they are forced in the actual incidence graph. In particular s_n>=T_n/2 and b_n>=B_n/2. Mirror pairing also makes both black columns even; this follows independently from the inherited pairwise black-line multiplicity law THM-796.

These are stronger than global colour totals: a count-preserving rewiring that inserts any periphery--periphery line is impossible even if it repairs all scalar degrees elsewhere.

## 3. Exact values and independent audit

| n | s_n | b_n | All crossing | Blue crossing | All core | Blue core |
|---|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 1 | 1 | 1 | 0 | 0 |
| 4 | 5 | 3 | 3 | 1 | 1 | 1 |
| 5 | 50 | 14 | 14 | 2 | 18 | 6 |
| 6 | 903 | 59 | 121 | 5 | 391 | 27 |
| 7 | 30773 | 495 | 1995 | 17 | 14389 | 239 |
| 8 | 2032504 | 4026 | 64648 | 70 | 983928 | 1978 |

verify_strong_core_counts.py has no repository imports. For every tiling n3..7 it constructs the directed adjacency matrix, tests reachability from and to vertex0, tests the literal coordinate reflection, and directly enumerates complement pairs. All six measured columns agree with the recurrences, with no within-periphery line. An independent sorted-score/Landau audit of the existing census also verified the displayed table through n8. The companion script prints recurrence-derived values through n12 and explicitly distinguishes those from its direct n3..7 enumeration.

The old total strong counts already appear in THM-330; THM-340 is a related weighted composition identity, not the reflection-fixed recurrence. THM-349 and THM-354 are the more precise inherited routes. The new content here is the explicit coupling of those counts to a forbidden incidence region and an exact colour-resolved cut.
