# A two-arc switch inside one of the five cores

2026-11-10. **PROVED scoped construction**, using the inherited fixed-core
compiler. The existence of tournament pairs at every even distance was
already known in the workspace. This note isolates a concrete operation
inside the five-core bank; it is not a new all-even-gap theorem and does
not select an arbitrary absolute count.

Take S5a with outgoing rows (2,4,9,19,7). Keep both transitive flanks,
their binary cut addresses A,B, all core edges, and all apex-flank edges
fixed. Compare the two apex cuts O={0} and O={1}. In the positive compiler

    H = h+9(g(A)-1+g(B)-1)+84AB+iA+oB+s alpha+e beta,

their five varying coordinates (i,o,s,e,h) are respectively

    O={0}: (52,11,8,3,17),
    O={1}: (54, 9,8,3,15).

Subtracting gives the all-parameter identity

    H_{0}(A,B)-H_{1}(A,B)=2(B-A+1).                    (1)

Only two arcs change: the apex's comparisons with core vertices0 and1.
No root is silently replaced, and every other edge survives. In particular:

- A=B gives adjacent odd counts, differing by2, for every A>=0.
- B=A+d-1 gives a pair at distance2d for every integer d>=1.
- A=B+1 gives a nontrivial count-preserving two-arc switch.

The smallest-address choice A=0, B=d-1 gives the explicit lower count

    N_d = 9d+9g(d-1)+3((d-1) mod2)-3,

with the upper count N_d+2d. Both graphs have
6+bitlength(d-1) vertices (bitlength(0)=0) and are strong: the mixed apex
cut makes the core plus apex strong, while the nonempty right flank's
last vertex beats the apex and every core vertex beats the flank.
The parameter d fixes the difference, not the lower endpoint. For example
d=1 produces15/17 and d=2 produces27/31. This is precisely the absolute-
target information lost if (1) is presented merely as a gap generator.

## Context-dependent form of the same switch

Allow a fixed additional apex out-set R contained in {2,3,4}, and switch
R union{0} to R union{1}. The coefficients of A,B and parity remain as
in (1); the constant changes with R:

| R mask | 0 | 4 | 8 | 12 | 16 | 20 | 24 | 28 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Constant c_R | 2 | 4 | -2 | 0 | 0 | 2 | -4 | -2 |

The exact difference is 2(B-A)+c_R. Thus even this tiny switch needs its
retained apex context; a scalar rule attached only to the H9 source would
lose the constant. The two-cut formula is an application of the existing
cut interaction calculus, evaluated here on a specific bank core.

## The other H9 core has a different switch geometry

For S5b=(2,4,25,3,11), the singleton cuts at vertices1 and2 have positive
coordinates (60,9,8,1,17) and (66,15,8,1,23), respectively. These vertices
again have the same Hamilton start/end data, but now

    H_{2}(A,B)-H_{1}(A,B)=6(A+B+1)>0.

There is no count-neutral line in the nonnegative quadrant. Thus equal
ordinary endpoint data still loses the residual-path information that
distinguishes cancellation in S5a from strict separation in S5b. Neither
phenomenon follows from their shared uncoloured H=9 or L=84.

## Actual reset probe and its stopping reason

Set A=t+1,B=t in S5a. The two-arc surgery is then count-neutral. For
t=0,...,7 its counts are 77,321,725,1305,2063,2961,4019,5289. The paired
sources have different one-apex output supports, witnessed already by
135 versus133 at H77 and515 versus513 at H321. Hence the operation changes
future response while preserving current count.

Two separate exact tests assess ternary return. First, every ONE-APEX
attachment is enumerated: both H77 sources admit3H+2=233 and neither
admits3H+4=235; for every tested t=1,...,7, neither member admits either
target by one apex. Second, independently of this attachment test, the
rooted two-bin count F is evaluated at every vertex. Both H77 sources
are1-ready at vertex0 (F-H=1), and no other1/2-ready root occurs in any
of the sixteen sources. Thus the standard two-birth reset is available
there, while t=1,...,7 have no ready root. One-apex branch attainment and
two-birth readiness are different predicates; neither is inferred from
the other. These are finite tests on the named sources, not all-order
impossibility statements. The changed response need not hit the required
target. `reset_checks.json` retains all root defects, scopes and witnesses.

## Limits and next test

The switch changes future response while preserving H on its zero line, but
the zero line depends on the selected context. It is not a universal
reset between all realisers of an arbitrary count. To turn adjacent pairs
into a Chan–Pak proof one still needs a construction whose LOWER count is
an arbitrary requested N, or a closed target descent that retains the
flank addresses and the apex context. Equation(1) alone provides neither.

`check_switch.py` recomputes every signature from literal core paths,
checks all eight context constants, and independently recounts actual
children for 0<=A,B<=7 with subset Hamilton-path dynamic programming.
The identity for arbitrary parameters follows by the displayed coefficient
subtraction, not by extrapolating those finite tests.
