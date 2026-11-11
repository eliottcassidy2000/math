# The missing residue phase requires a different address scale

Status: **PROVED in this note**, using the inherited fixed-core constructor. This note adds an infinite simultaneous-rebase family in residues 1/0 modulo three, complementing the earlier 0/2 family. It does not claim all-address coverage.

## The second C3 unit-switch family

For \(C\ge0\), put
\[
N(C)=3g(C)+12C+10.
\]
The actual C3 constructions
\[
F_{C_3,3}(0,2C+1)=N(C),\qquad
F_{C_3,1}(0,4C+3)=N(C)+2
\]
have equal rooted coordinates
\[
E=d=3,\qquad U=15\,2^c-3\operatorname{rev}_c(C),
\]
where \(c\) is the bit length of \(C\), with the zero-word convention used in the preceding notes. These formulas follow by substituting the two cuts in the fixed-core compiler; reversing their final one and two 1-bits gives the shared \(U\). At \(C=0\) the pair is 13/15.

Unlike the [earlier C3 pair](ADDITIVE_CARRY.md), whose lower count \(3g(C)+36C+30\) is divisible by three, this pair occupies residues 1/0. The all-prefix repeated-9 theorem therefore did not address this residue phase.

## No direct dyadic suffix rebase

**Theorem 1.** There are no integers \(B,t,s\) with \(B,t\ge0\), \(0\le s<2^t\), such that
\[
L(B)=N(2^t B+s).
\]
Thus appending any fixed binary word to the old address cannot align this second family with the H13 twins at that address.

**Proof.** For \(t=0\), the source-minus-replacement difference is
\(10g(B)+20B+6\beta+8>0\). For \(t=1\), the differences for \(s=0,1\) are
\[
10g(B)+5B+3\beta+8,\quad
10g(B)+8B+6\beta-4,
\]
both positive. For \(t=2\), the four differences are
\[
\begin{array}{c|l}
s&L(B)-N(4B+s)\\\hline
0&10g(B)-25B+3\beta+8\\
1&10g(B)-19B+3\beta-4\\
2&10g(B)-22B+6\beta-22\\
3&10g(B)-16B+6\beta-28.
\end{array}
\]
Use \(g(B)\le B+1-\beta\). All four expressions are negative for \(B\ge1\), except that the first requires the direct check \(B=1\), where it is \(-4\).

Appending additional bits cannot decrease \(N\), because
\[
N(2C)-N(C)=15C+3(C\bmod2)\ge0,
\quad N(2C+1)-N(C)=12C+12>0.
\]
This excludes every \(t\ge2\) when \(B\ge1\). Finally, \(L(0)=31\); equality \(N(C)=31\) would require \(g(C)+4C=7\), hence \(C\le1\), and neither candidate works. ∎

This is a scale obstruction for this restricted address grammar. It says nothing against other address maps, additional flank variables, or arbitrary representations of the same targets.

## A ternary address supplies the missing phase

Define the common defect under the address map \(C=3B\):
\[
D(B)=L(B)-N(3B)
=13g(B)-3g(3B)-4B+6(B\bmod2)+8.
\]
Binary recursion gives the exact identity
\[
D(2B)=D(B)+4(B\bmod2).
\]

**Theorem 2 (ternary seed and infinite phase-1 ray).** For every integer \(j\ge0\), let
\[
B=1002\,2^j,\qquad C=3006\,2^j.
\]
Then both H13 targets have simultaneous C3-source witnesses:
\[
\begin{aligned}
L(B)&=F_{C_3,3}(0,6012\,2^j+1)=45090\,2^j-4415,\\
L(B)+2&=F_{C_3,1}(0,12024\,2^j+3)=45090\,2^j-4413.
\end{aligned}
\]
These counts occupy residues 1/0 modulo three.

**Proof.** The finite seed has \(g(1002)=661\) and \(g(3006)=1531\), so
\[
D(1002)=13(661)-3(1531)-4(1002)+8=0.
\]
The seed and all its doublings are even, so the doubling identity preserves zero. Also
\(g(1002\,2^j)=661+(2^j-1)1002\). Substituting gives the stated targets and the two literal witnesses. ∎

The first targets are 40675 and 40677. This is a genuinely different residue phase from the repeated-9 family. Its present scope is an infinite ray, not an arbitrary-prefix theorem.

## A small exact carry system for the ternary map

The move \(C=3B\) is not opaque arithmetic: it has three binary carry states. Put \(G_r(B)=g(3B+r)\), \(r=0,1,2\), and \(\beta=B\bmod2\). The exact transitions are
\[
\begin{array}{c|cc}
 &2B&2B+1\\\hline
G_0&G_0(B)+3B+\beta&G_1(B)\\
G_1&G_0(B)&G_2(B)+3B+2+\beta\\
G_2&G_1(B)+3B+2-\beta&G_2(B).
\end{array}
\]
Together with the transitions for \(g(B),B,\beta\), this is a closed exact linear state system for the ternary defect \(D\). It preserves the target difference while recording the carry omitted by a direct dyadic suffix map. The all-zero input continuation in Theorem 2 is one cycle of this larger system; finding further synchronizing cycles is a precise continuation task.

## Controls

`residue_phase_probe.py` checks the carry identities on their stated finite universe, both literal C3 witnesses at the seed, and the infinite-ray formulas for a finite set of exponents. The proofs above supply the arbitrary quantifiers. No broad spectrum census is used.
