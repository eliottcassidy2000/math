# Eight exact source-rewrite languages

Status: **PROVED in this note**, using the inherited fixed-core formula and the binary recursion stated in [PROOF.md](PROOF.md). This is a partial rewrite grammar, not a complete coverage theorem.

Write
\[
L(B)=13g(B)+32B+6\beta+18,\quad T(B)=L(B)+2,\quad
\beta=B\bmod2,
\]
and let \(F_{C,O}(A,D)\) denote the literal fixed-core constructor's HP count. Here \(S_4=(2,4,9,3)\) and \(S_{5b}=(2,4,25,3,11)\).

**Theorem (eight iff rewrites).** In each row, the displayed source count equals the replacement count exactly on the listed address set. Each row has arbitrary integer quantifiers, not just a finite observed pattern.

| Rule | Source | Replacement | Exact equality set |
|---|---|---|---|
| U0 | \(T(B)\) | \(F_{S_4,8}(0,4B+3)\) | \(B=2^j-2,\ j\ge2\) |
| U1 | \(T(B)\) | \(F_{S_4,11}(0,2B+1)\) | \(B=0\), or \(B=2^j,\ j\ge1\) |
| U2 | \(T(B)\) | \(F_{S_4,4}(1,B)\) | \(B=4\), or \(B=2^j+4,\ j\ge4\) |
| L0 | \(L(B)\) | \(F_{S_4,0}(1,B)\) | \(B=2^j-1,\ j\ge1\) |
| L1 | \(L(B)\) | \(F_{S_4,1}(0,4B+2)\) | \(B\in\{0,3\}\), or \(B=2^j+3,\ j\ge3\) |
| L2 | \(L(B)\) | \(F_{S_4,4}(0,4B+2)\) | \(B=2\), or \(B=2^j+2,\ j\ge3\) |
| L3 | \(L(B)\) | \(F_{S_4,2}(1,B)\) | \(B=6\), or \(B=2^j+6,\ j\ge4\) |
| L4 | \(L(B)\) | \(F_{S_{5b},11}(0,B)\) | \(B=9\,2^j+2,\ j\ge3\) |

Therefore **simultaneous five-core failure of both twins is impossible on the union of these eight languages**. This assertion is target-preserving: each applicable row constructs an actual witness for one of the two targets. It neither covers every address nor proves the other twin is covered.

## Energy and the proofs

The binary energy
\[
f(n)=n+1-(n\bmod2)-g(n)
\]
has the exact recursions
\[
f(2q)=f(q),\qquad f(2q+1)=f(q)+q+(q\bmod2).
\]
Equivalently, if \(x_i\) are the binary digits of \(n\),
\[
f(n)=\sum_{i<j}w(j-i)x_ix_j,
\qquad w(1)=2,\quad w(d)=2^{d-1}\ (d\ge2).
\]
Both descriptions follow from the recursion for \(g\). In particular:

1. \(f(n)\ge0\), and \(g(n)\le n+1\).
2. \(f(n)=0\) iff the word has at most one 1-bit.
3. \(g(n)=n+1\) iff \(n=0\) or \(n=2^j\) for an integer \(j\ge1\).
4. \(f(n)=4\) iff exactly two bits are present and their distance is three: \(n=9\,2^t\), \(t\ge0\). Three bits already contribute at least six; among two-bit distances only three has weight four.

Subtracting the replacement from its source gives:

| Rule | Exact source-minus-replacement defect |
|---|---|
| U0 | \(8g(B)-4B+6\beta-16\) |
| U1 | \(8g(B)-8B+6\beta-8=-8f(B)-2\beta\) |
| U2 | \(8g(B)-7B+5\beta-12\) |
| L0 | \(8g(B)+6\beta-14\) |
| L1 | \(8g(B)-2B+6\beta-8\) |
| L2 | \(8g(B)-6B+6\beta-12\) |
| L3 | \(8g(B)-5B+5\beta-10\) |
| L4 | \(4g(B)-3B-\beta+10\) |

U0 was proved in Theorem 1 of PROOF.md. U1 follows from nonnegative energy. L0 forces \(\beta=1\) and \(g(B)=1\), hence an all-one word.

For U2, reduction modulo eight leaves \(B=8q+4\) or \(8q+7\). The odd branch gives \(8(g(q)-7q-7)<0\). The even branch gives \(8(g(q)-q-1)\), whose zeros are \(q=0\) and \(q=2^t\), \(t\ge1\).

For L1, the even branch \(B=2q\) gives \(4(2g(q)+q+2(q\bmod2)-2)\), which vanishes only at \(q=0\). The odd branch forces \(B=4r+3\), then reduces to \(8(g(r)-r-1)\).

For L2, the odd branch \(B=2q+1\) gives \(4(2g(q)-3q-3)<0\). The even branch forces \(B=4r+2\), then reduces to \(8(g(r)-r-1)\).

For L3, reduction modulo eight leaves \(B=8q+6\) or \(8q+7\). The odd branch gives \(8(g(q)-5q-5)<0\). The even branch gives \(8(g(q)-q-1)\).

For L4, an odd input \(B=2q+1\) would require \(2g(q)-3q+3=0\), hence \(q=2r+1\) and \(g(r)=3r\), which has no solution: \(r=0\) fails, and \(g(r)\le r+1<3r\) for \(r\ge1\). The even branch forces \(B=4r+2\), then reduces to \(4(g(r)-r+3)\). Its zeros satisfy \(f(r)+(r\bmod2)=4\). Energy is even, so \(r\) is even and \(f(r)=4\); therefore \(r=9\,2^t\), \(t\ge1\), giving \(B=9\,2^{t+2}+2\). ∎

## A recursive defect compiler

Every row belongs to the same exact family
\[
E(n)=a g(n)+bn+c(n\bmod2)+d.
\]
Reading the address from its least significant end updates its coefficient state by
\[
\begin{aligned}
E(2q)&:\ (a,b,c,d)\longmapsto(a,a+2b,a,d),\\
E(2q+1)&:\ (a,b,c,d)\longmapsto(a,2b,0,d+b+c).
\end{aligned}
\]
At the empty remaining word, equality holds exactly when \(a+d=0\). Thus every candidate suffix rewrite has an exact carry system. The state is not asserted finite; a finite closed collection, compatible with the paired tournament operations, would be substantially stronger than the eight languages above.

This is the operation-bearing connection with square-sum paired radix constructions: a compatible pair needs an explicit transition for every next digit and enough marked data to keep its target. It also connects graceful gap insertion's need for a marked baseline with an exact arithmetic remainder. The marked remainder here is \(E\), not just the residue of the eventual target.

## What the energy debt does and does not buy

The U1 replacement is defined on every address and obeys
\[
F_{S_4,11}(0,2B+1)=T(B)+D(B),\qquad D(B)=8f(B)+2\beta.
\]
Its defect has simple digit transitions:
\[
D(2q)=D(q)-2(q\bmod2),\quad
D(2q+1)=D(q)+8q+6(q\bmod2)+2.
\]

The H13 paired core switch has pure count gap \(3A-(A\bmod2)+2\), with all other rooted coordinates equal. These direct gaps are the positive even integers not divisible by six. The odd-step debt increment above is such a direct gap iff \(q\not\equiv2\pmod3\). The companion [portable-gap theorem](../dynamics/PORTABLE_GAPS.md) supplies **actual rooted switches of every positive even width**: factor the width as \(3^a v\), use a direct switch of width \(v\), and prepend \(a\) ordinal C3 blocks to scale its full state by \(3^a\). What remains missing is baseline compatibility: a switch of the right width has not thereby been attached to this particular S4 target representation while preserving its other root data.

Nor can the left flank of this same S4 chart absorb the debt directly. At fixed right address \(2B+1\), introducing \(A\ge1\) adds
\[
5(g(A)-1)+(64B+37)A+(A\bmod2),
\]
which is positive and exceeds the debt bound \(D(B)\le8B\). A successful actuator must change the representation more substantially than introducing that single marked flank.

The new positive objects are therefore (i) eight exact infinite target-preserving rewrite languages, (ii) a general recursive defect representation, and (iii) a sharply isolated missing operation: absorb an even energy debt by a rooted switch with a compatible baseline. [ADDITIVE_CARRY.md](ADDITIVE_CARRY.md) follows the defect compiler further to an all-prefix simultaneous rebase, using repeated hexadecimal 9 as an additive actuator.

## Finite receipt

`grammar_probe.py` checks all eight iff statements for \(0\le B<2^{16}\), exercises selected literal tournament witnesses, and checks the defect transitions. This finite universe is a control for the proofs, not a claim of exhaustive spectrum coverage.
