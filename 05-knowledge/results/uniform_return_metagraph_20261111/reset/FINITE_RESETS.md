# A shorter141 reset and the exact limits of its prefix lifts

**PROVED literal construction and initial-deletion optimum; FINITE-EXACT
cut screens; PROVED response formulas for the specified lifts.**

## 1. The best initial one-vertex deletion retains H73

For the old eight-vertex source T_a(2), delete prefix vertex1 and keep old
labels(0,2,3,4,5,6,7). The actual remaining graph is T_a(1), with H73.
Adjoin one vertex with outmask53 in that seven-vertex order. The final rows
are

    (62,172,56,144,98,72,143,53),

and the exact trace is141→73→425. The target is strong, has order8, and
retains all seven old vertices and every arc between them. Outmask117
gives a second strong target, with rows

    (62,172,56,144,98,72,15,117).

All eight one-vertex cards of the specified T_a(2) have counts in the
multiset{9,17,33,49,59,59,73,73}. Thus73 is the maximum count after an
INITIAL single old-vertex deletion. Every subsequent insertion is
H-monotone, so this operation has the largest possible minimum count
among programs that begin with one deletion from this exact source and
then only insert vertices while retaining the remaining old arcs.

The optimum does not cover preparatory edge changes, earlier insertions,
multiple deletions, or arbitrary other H141 sources. It is a precise local
source-replacement optimum, with a complete eight-card certificate.

## 2. Why this shorter reset does not replace the uniform one

Let the remaining prefix vertex grow to TT_m in either new target. Its
positive run coefficients are respectively

    (A1,A2,A3)=(425,404,116),   (425,400,116),

and all higher coefficients vanish. Hence the all-m count laws are

    116·3^m+56·2^m−35,
    116·3^m+52·2^m−27.

At m1 both give425, but the desired uniform branch after deleting one
prefix from T_a(k) is204·2^m+17 with m=k−1. The positive3^m contribution
is an exact growth obstruction. This is why the longer but uniform
[middle-tail reset](UNIFORM_RETURN.md) remains necessary for the current
infinite-family result.

## 3. Nonuniform prefix cuts give real isolated hits

After deleting one prefix vertex from T_a(k), exhaust every one-birth
cut for k2..8. The3H+2 branch occurs at k2,k4,k5; the3H+4 branch occurs
only at k4 in this universe. At k4 its two cuts are106 and482, giving

    549→277→1651.

On the remaining three prefix vertices both cuts have pattern010.
This is a positive reason to investigate variable-incidence prefixes.
The simplest attempts to extend it uniformly nevertheless fail exactly:
inflating any one prefix vertex misses the required dyadic coefficient;
inflating any two equally has a6^m or9^m mode instead of pure4^m.
Inflating all three equally contains the last-two inflation as an induced
subgraph, whose positive9^m mode exceeds the intended8^m target growth.
The receipt retains every exponential coefficient, not fitted samples.

These statements exclude those seven equal-inflation choices for these
two actual targets. They do not exclude arbitrary block lengths, different
cuts, parameter-dependent sources, or the+4 family itself.

## 4. A lower-defect H141 source is also blocked until replacement

The [minimal-defect collar family](../metagraph/MINIMAL_DEFECT.md) gives
another H141 source, with rows(240,5,9,3,238,206,142,14), a three-port
root of defect6, and automorphism order3. Its complete strong frontier
under427 has sizes27,4,0 and contains neither425 nor427. The low-ratio
ear theorem applies since427<10·141, and both targets have all proper
divisors<141. Hence this actual source excludes BOTH targets at every
possible extension order. Smaller retained root defect is not a sufficient
source-selection criterion.

Deleting any one of its four transitive prefix vertices leaves H69 and
admits six one-birth425 cuts; deleting an in-port leaves H47 and admits
one such cut. No one-birth427 cut occurs in any of its eight deletion
cards. These last statements concern the explicitly complete eight-card
one-birth universes, not arbitrary replacement programs.

Reproduce with `python -B finite_cards.py`. The verifier rebuilds the old
source, all its deletion cards, both targets and every persistent arc,
the complete prefix-cut universes, and the marked module responses.
