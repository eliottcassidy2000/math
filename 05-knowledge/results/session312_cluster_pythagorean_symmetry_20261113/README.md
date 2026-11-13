# Session 312: Arithmetic memory and recursive tournament symmetry

- [Berggren–Price symmetry clock](ARITHMETIC.md): mixed-schedule group theorem, with the full path-slot polynomial fixed.
- [Cluster address chart](CLUSTER_CHART.md): positive Plücker coordinates, exact tropical exchange, and the boundary for full orientation.

# Two tests of the memory principle

## When the Price reset is invisible

The [mixed-schedule theorem](ARITHMETIC.md) shows that each Price step spends exactly one bit of dyadic precision. This is observer-dependent. If the palette is injective on projective residues modulo an odd integer $q>1$, **every mixed schedule is rigid and retains the entire seed residue**. The same holds for a mixed modulus $2^kq$ with $q>1$ odd, except that the dyadic residue is still retained only at the reduced precision.

For the odd assertion, compare two sibling triples under a putative rotation of one or two positions. All matrices are invertible modulo any odd prime $p\mid q$. For a shift $s$, the first matched child determines the other seed projectively as $T_s^{-1}T_0u$ in the ordered branch matrices $T_0,T_1,T_2$. Substitute this into the two other determinant equations. Up to nonzero scalar and sign, the necessary simultaneous zeros are:

| Branches | shift one | shift two |
|---|---|---|
| Berggren | $m^2+n^2,\ n^2-m^2$ | $m^2-n^2,\ mn$ |
| Price | $(m-n)^2,\ m^2-6mn+n^2$ | $m^2+2mn-n^2,\ m^2+4mn-n^2$ |

No pair vanishes at a nonzero column over an odd field. In the Price shift-one case substitute $m=n$; in shift two subtract to get $2mn=0$. The Berggren cases follow by addition or from $mn=0$. Hence only the identity rotation works, and one matched child recovers the full odd residue by inversion. Induct upward. For mixed modulus, odd alignment forces identity at every level even after the dyadic part becomes constant. The common-word determinant law still gives the dyadic precision loss exactly.

For the purely dyadic experiment, a fixed schedule with $r$ Price levels has exactly $2^{\max(k-r,0)}$ classes over all primitive positive seeds, since every dyadic projective residue occurs. In the chart with $m$ odd take $m=1$ and a positive even lift for $n$. In the chart with $n$ odd choose a positive even lift $m$ and $n=1+2^a m$ at precision $a>0$. Then $n>m$ and $\gcd(m,n)=1$. At precision one the union over all schedules and seeds has exactly $d+1$ classes at height $d$. For the proof, with no Price there are two distinct rigid classes, separated by their leaf color counts. A schedule with last Price at level $t$ depends only on $t$: the residual Berggren subtrees are fixed after the reset and everything above it is repeated inflation. The $t=1$ class is the pure Berggren class from first-coordinate parity zero: its $(1,1,0)$ tuple is a cyclic rotation of $(0,1,1)$. The remaining $t=2,\dots,d$ give different nontrivial group orders. For $d=0$ use the two leaf colors, so the $d+1$ formula is for $d\ge1$.

## What the cluster chart forgets

The [cluster chart](CLUSTER_CHART.md) separates three levels of observation. The valuation alone recovers the rooted ternary balls, so its isometry group is the iterated wreath tower of $S_3$: every permutation of the three children at every ball is allowed. The signed leading unit specifies the directed triangle between children and removes precisely the reflections, leaving the iterated $C_3$ group. Keeping the actual leaf colors takes the color stabilizer in that group. The real positivity of the Plücker chart does not distinguish these levels by itself.

The arithmetic minors exhibit a different finite clock. On $\mathbb P^1(\mathbb F_3)$ take $A=(1,0),B=(0,1),C=(1,1),D=(1,-1)$. Price's three branches send the four states respectively to $(D,A,A),(B,D,D),(C,B,C),(A,C,B)$. Only $D$ has all sibling minors units. From root $D$, the other three populations stay equal, and the $D$ population satisfies $g_{j+1}=3^j-g_j$, $g_0=1$. Therefore
$$
g_j=\frac{3^j+3(-1)^j}{4}.
$$
No two consecutive states of a Price ray are good, and no Berggren parent is good in its raw spinor chart. This gives an exact arithmetic alternation in a different observer; it is not a tournament readiness criterion.

These results suggest a more precise test for a count-preserving Chan–Pak repair: prove **which observer is allowed to forget which coordinate**. A scalar equality or chart identity alone is insufficient. A port-aware inverse must both separate the legal branches and lift them to physical operations; no such all-target inverse is proved here.
