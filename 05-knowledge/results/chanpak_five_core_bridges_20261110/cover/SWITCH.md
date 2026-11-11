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

## Limits and next test

The switch can alter readiness while preserving H on its zero line, but
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
