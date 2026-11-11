# A repeated-digit actuator for simultaneous rebasing

Status: **PROVED in this note**, using the inherited exact fixed-core constructor. This strengthens the eight regular languages in [GRAMMAR.md](GRAMMAR.md) to a family containing every finite binary word as a prefix. It still does not cover every address or preserve the target of the original unpadded address.

Let \(C_3=(2,4,1)\), and retain \(L,T,g\) from the earlier notes. Define
\[
E(B)=10g(B)-4B+6(B\bmod2)-12.
\]

**Lemma 1 (one defect, two actual rewrites).** For every \(B\ge0\),
\[
\begin{aligned}
T(B)-F_{C_3,1}(0,8B+6)&=E(B),\\
L(B)-F_{C_3,0}(0,16B+12)&=E(B).
\end{aligned}
\]
In particular, every zero of \(E\) gives simultaneous HP-3-source witnesses for both H13 twins.

**Proof.** The two C3 signatures are
\[
(3,9,-3,6,-2,3,-1,2,0),\qquad
(3,9,-3,9,-3,0,0,0,0).
\]
Thus their counts are \(3g(8B+6)+24B+20\) and \(3g(16B+12)\). Binary recursion gives
\[
g(8B+6)=g(B)+4B+4,\quad
g(16B+12)=g(B)+12B+10.
\]
Substitution proves both equalities. ∎

**Lemma 1a (a second complete unit-switch family).** The two C3 replacements have counts
\[
M(B)=3g(B)+36B+30,\qquad M(B)+2
\]
for every \(B\ge0\), even when their defect relative to H13 is nonzero. If \(b\) is the bit length of \(B\), with \(b=0\) for \(B=0\), both added-root states have
\[
E=d=3,\qquad U=39\,2^b-3\operatorname{rev}_b(B).
\]
Thus this is a pure count switch, portable under any common inherited \(O/Z\) continuation. The pair at \(B=1\) is 69/71, now represented by larger tournaments with HP-3 deletion sources.

**Proof.** The count formulas follow from Lemma 1. For a nonempty right address \(D\), the inherited one-flank rooted-state formula is \(E=d\) and
\(U=(d+o_+)2^{\operatorname{bitlength}(D)}-d\operatorname{rev}(D)\). Both replacement addresses are nonempty. Their reverse words are respectively
\(3\,2^b+\operatorname{rev}_b(B)\), and their lengths are \(b+4\), \(b+3\). Substituting \(d+o_+=3\) and 6 gives the same \(U\). Equality of \(E,U,d\) makes the difference in \(H\) invariant under both \(O\) and \(Z\). ∎

**Lemma 2 (hexadecimal 9 is an additive defect actuator).** For \(q\ge0\) and \(m\ge1\), set
\[
B_m(q)=16^m q+9\frac{16^m-1}{15}.
\]
This appends \(m\) copies of the hexadecimal digit 9 to the binary prefix \(q\). Then
\[
E(B_m(q))=10g(q)-4q+4m-6.
\]

**Proof.** One four-bit expansion gives
\[
g(16q+9)=g(q)+6q+4,
\qquad E(16q+9)=10g(q)-4q-2.
\]
For \(S_d(q)=10g(q)-4q+d\), the same identity gives
\(S_d(16q+9)=S_{d+4}(q)\). After the first appended 9, each further copy therefore adds exactly four to the residual defect. ∎

**Theorem (arbitrary-prefix simultaneous rebase).** Given any integer prefix \(p\ge0\), let
\[
q=4p+3,\qquad
m=q-\frac{5g(q)-3}{2},\qquad B=B_m(q).
\]
Then \(m\) is an integer at least two, and both H13 targets have literal five-core witnesses:
\[
T(B)=F_{C_3,1}(0,8B+6),\qquad
L(B)=F_{C_3,0}(0,16B+12).
\]

**Proof.** Every \(g(q)\) is odd, so \(m\) is integral. Also \(g(q)=g(p)\le p+1\), hence
\[
m=4p+3-\frac{5g(p)-3}{2}
\ge\frac{3p+4}{2}\ge2.
\]
Lemma 2 gives \(E(B)=0\), and Lemma 1 supplies both tournaments. ∎

For \(p=0\), this gives \(q=3,m=2,B=921\). For \(p=1\), it gives the binary word corresponding to hexadecimal \(7999999\). Every input pattern survives unchanged as a prefix, followed by `11` and a computable number of `1001` blocks. The padding is generally long: \(m\) is of order \(p\), rather than of order the bit length of \(p\). This is an exact constructive family, not an efficient completeness algorithm.

## Why a five enters the operation

The defect compiler's coefficient pair is \((a,b)=(10,-4)\). For the branch that keeps \(-a\le b\le0\), its \(b\)-values cycle
\[
-4\xrightarrow{1}-8\xrightarrow{0}-6\xrightarrow{0}-2\xrightarrow{1}-4.
\]
After normalization \(r=-b/a\), the cycle is doubling modulo one on
\(2/5,4/5,3/5,1/5\). Its four-step period is the order of 2 modulo 5. The corresponding suffix is `1001`, and the full coefficient state returns with additive drift 4 in its constant term. This is an operational appearance of the number five: a binary carry orbit, with an exact target defect and a realizable source rewrite.

More generally, if a binary word returns the nonconstant coefficient triple \((a,b,c)\) of a defect state while changing its constant by \(\lambda\ne0\), repeated copies form an additive actuator \(d\mapsto d+m\lambda\). A target-preserving source rewrite results whenever the original defect is canceled by a nonnegative integral number of copies. The theorem above gives one fully realized instance and a prefix preparation that always satisfies the sign and integrality conditions.

The rooted response supplies a further exact cancellation. On \(B=B_m(q)>0\), the H13 response minus its matching C3 response is
\[
\Delta U=6\,2^{\operatorname{bitlength}(q)}-10\operatorname{rev}(q),
\]
independent of the padding length \(m\). Indeed, the reverse of the repeated-9 suffix contributes
\(\tfrac35(2^{\operatorname{bitlength}(B)}-2^{\operatorname{bitlength}(q)})\), which cancels the growing term. Count alignment therefore does not discard the next-extension obstruction. In the theorem's preparation \(q=4p+3\), this response difference is negative; arbitrary subsequent continuation is not an equality-preserving replacement between the H13 and C3 representations.

## The more complex source has the cheaper future

**Theorem (all-prefix response boundary).** For every \(q\ge1\) with
\(m(q)=q-(5g(q)-3)/2\ge1\),
\[
\Delta(q)=6\,2^{\operatorname{bitlength}(q)}-10\operatorname{rev}(q)\le-6.
\]
Consequently, at every source rebase produced by the repeated-9 actuator, the H13 witness weakly lowers the HP count of every future common \(O/Z/J\) word relative to the corresponding H3 witness. Decreasing the displayed parent deletion count from 13 to 3 does not improve that future cost.

**Proof.** Write \(q=16p+s\), \(p\ge1\), and \(r=\operatorname{rev}_4(s)\). Concatenation gives
\[
\begin{aligned}
\Delta(q)&=(96-10r)2^{\operatorname{bitlength}(p)}-10\operatorname{rev}(p),\\
2m(q)&=2m(p)+5(r-9)p+K_s-5(p\bmod2)(1-\lfloor s/8\rfloor),\\
K_s&=2s-5g(s)+5.
\end{aligned}
\]
If \(r\le8\), then \(s\) lies in \(\{0,8,4,12,2,10,6,14,1\}\), and its respective \(K_s\) values are
\(0,-24,-12,-26,-6,-20,-8,-12,2\). Since \(m(p)\le p-1\), it follows that \(2m(q)\le-3p<0\). Thus this case is impossible when \(m(q)\ge1\).

If \(r\ge10\), then \(\Delta(q)\le-4\,2^{\operatorname{bitlength}(p)}-10\operatorname{rev}(p)\le-18\). If \(r=9\), then \(s=9\), \(\Delta(q)=\Delta(p)\), and \(m(q)=m(p)-1\). Strip such trailing 9s until one of the other cases applies or a one-nibble base remains. The bases with \(1\le q<16\) and \(m(q)\ge1\) are \(3,7,11,13,15\), with \(\Delta\)-values \(-6,-22,-34,-14,-54\). This proves the bound.

For either aligned twin, use the companion [continuation cone](../dynamics/PROOF.md) coordinates
\(\mathbf q=(H,U,U-d,H+U-E)\). H13 minus H3 has the exact difference
\[
(0,\Delta(q),\Delta(q)-10,\Delta(q)-10),
\]
which is coordinatewise nonpositive. The \(O/Z/J\) transitions are positive in these coordinates, proving every future count comparison. ∎

The equality bound \(-6\) is attained at \(q=3\) and \(q=57\). The result therefore gives a genuine boundary for this actuator, not merely a failed small search for an H3-dominating prefix. The source rebase lowers the deletion count but raises the relevant future response.

## Scope relative to the paired-square and graceful pointers

The old paired square-sum construction suggests preserving two witnesses under a digit grammar. Here the two H13 defects are literally the same function, so the same actuator cancels both; this is stronger than a scalar coincidence for one twin. The graceful gap-insertion pointer suggests retaining a marked obstruction that an operation changes predictably. Here the marked obstruction is \(E\), and each 9 changes its residual value by exactly four.

What remains missing is a way to apply such a debt actuator while retaining an arbitrary **original target**. Padding changes \(B\), and hence changes both target counts. The theorem rules out simultaneous failure on a language containing every finite pattern as a prefix; it does not rule out simultaneous failure at an arbitrary specified address. That distinction is essential for the Chan–Pak sufficiency step.

## Reproduction

`carry_probe.py` checks the two defect identities and the arbitrary-prefix construction, with a receipt that records exact universes and output sizes. Its finite controls are supplementary to the proofs.
