# An actual flank-block transfer and its target-repair gate

Status: **PROVED in this note** relative to the inherited H25 constructors and common-root amalgamation formula. The result identifies a physical block operation, classifies its preserved ports, and proves a sharp qualitative limit on target restoration. It is not an all-target Chan–Pak argument.

## 1. What is transferred from the old repair constructions

[THM-4479, strategy-cube distance to provability, Lemma 1.1](../../../../01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md) succeeds because an accepted first-descent block never encounters a repaired residue in its interior. Its exact heredity statement is the interface condition preventing rewiring from invalidating a block already chosen.

[THM-4850, paired square-sum radix construction](../../../../01-canon/theorems/THM-4850-ninja-pairs-of-the-square-sum-graph-start-at-42-and-the-base-49-numeration.md) supplies a different concrete interface: paired paths retain positional parity and their endpoint ports under every radix digit. The transferable requirement is to find a block move whose actual boundary response agrees, before interpreting its scalar effect as a usable digit.

For the H25 source, this requirement gives an explicit Diophantine block grammar rather than an analogy.

## 2. The H25 block grammar

Use the two literal sources and cuts from [the paired-source reset](../../h13_paired_descent_20261111/modular_cycles/PROOF.md):
\[
C_-=(50,8,3,5,46,14),\ O_-=38;\qquad
C_+=(50,40,35,5,46,8),\ O_+=35.
\]
For positive left/right addresses \(A,B\), write \(\alpha=A\bmod2\), \(\beta=B\bmod2\). Their counts are
\[
\begin{aligned}
F_-(A,B)&=25(g(A)+g(B)-1)+304AB+107A+15\alpha+122B+16\beta+48,\\
F_+(A,B)-F_-(A,B)&=-3A-\alpha+2\beta.
\end{aligned}
\]
The common forward ports are
\[
E=d=25,\quad
U=(304A+147)2^b-25\operatorname{rev}_b(B),\quad
K=304\,2^{a+b},
\]
where \(a,b\) are the two bit lengths and K is the ordered two-path-cover count after deleting the root. These are the ports retained by the actual additive-cube construction in [ADDITIVITY.md](../../h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md).

**Theorem 1 (complete fixed-port address transfer).** Let \(t\ge1\), \(M=2^t\), and \(0\le c,x<M\). The move
\[
(A,MC+x)\longmapsto(MA+c,C),\qquad A,C\ge1,
\]
preserves E,U,d,K, and the total number of vertices, if and only if
\[
304c+25\operatorname{rev}_t(x)=147(M-1). \tag{1}
\]
Moreover, every equality of E,U,d,K between two positive-address realizations with different left-flank lengths has this form, after reversing its direction if necessary.

**Proof.** Reversing the old right word gives
\(\operatorname{rev}_{b+t}(MC+x)=2^b\operatorname{rev}_t(x)+\operatorname{rev}_b(C)\). Substituting into U yields exactly (1). The total flank length is unchanged, hence so is K.

Conversely, equal K fixes the sum of flank lengths. If the new left flank is t bits longer, its right flank is t bits shorter. Equality of U forces the difference between the two right reversals to be divisible by \(2^b\), since 25 is odd. Their remaining b low bits are therefore equal, so the old right address is \(MC+x\). The equality then gives \(A'=MA+c\), with
\[
\frac{122(M-1)}{304}\le c\le\frac{147(M-1)}{304}<M.
\]
At equal left-flank length, equality of U already forces both addresses equal: a change in A contributes at least \(304\,2^b\), while the entire possible reversal change contributes less than \(25\,2^b\). ∎

The first solution is
\[
t=7,\quad c=61,\quad x=80,
\qquad304(61)+25(5)=147(127).
\]
It moves the seven-bit block `1010000` from the right address into the left block `0111101`, retaining the exact external forward response.

## 3. A uniform scalar effect, and why the pair toggles cannot cancel it

Let \(\gamma=C\bmod2\). For any solution of (1), put
\[
\begin{aligned}
P&=132(M-1)-25\operatorname{rev}_t(c)-304x,\\
S&=25(g(c)-g(x))+107c-122x+15(c\bmod2)-16(x\bmod2).
\end{aligned}
\]
Then the exact minus-branch scalar change is
\[
F_-(MA+c,C)-F_-(A,MC+x)=PA+10\alpha+16\gamma+S. \tag{2}
\]
In particular the arbitrary suffix C disappears except for its parity. Equation (2) follows from the binary concatenation formula for g; the coefficient of C is precisely the left side of (1) minus its right side.

**Theorem 2 (strict ordering throughout a port fibre).** Every nontrivial move in Theorem 1 strictly decreases H, even if either old or new core branch is independently toggled. Consequently the two-count sets belonging to different cut positions in the same E,U,d,K fibre are disjoint and strictly ordered by left-flank length. There is no nontrivial target-neutral rewrite inside this fixed-port H25 chart.

**Proof.** Equation (1) implies \(c<M/2\). Reducing it modulo two shows that \(\operatorname{rev}_t(x)\) is odd, hence \(x\ge M/2\). Therefore
\[
P\le-20M-132.
\]
Since \(g(c)\le c+1-(c\bmod2)\) and \(g(x)\ge1\),
\[
S\le132c-122x\le\frac{215}{76}M.
\]
Each plus branch has a strictly smaller H than its minus branch when A is positive. Thus the most favorable possible comparison is new-minus versus old-plus, which adds at most \(3A+\alpha\) to (2). It is at most
\[
(P+3)A+11\alpha+16\gamma+S
\le-\frac{1305}{76}M-102<0.
\]
The converse statement in Theorem 1 makes this a classification of the full port fibre, not only of a sampled list of moves. ∎

For the seven-bit move, (2) is
\[
-9906A+10\alpha-4768+16\gamma.
\]
At \(A=C=1\), the four literal 16-vertex realizations have
\[
\begin{array}{c|cc}
&-&+\\\hline
(A,B)=(1,208)&93803&93799\\
(A',B')=(189,1)&79155&78589.
\end{array}
\]
All have \((E,U,d,K)=(25,115181,25,155648)\).

## 4. A universal entry gap for positive continuations

The previous theorem decreases the target count. It remains necessary to determine whether a different continuation can repay that deficit.

**Theorem 3 (first positive count gap).** Let A be a rooted tournament with \(U_A>d_A\). For any finite construction program using O, Z, J and arbitrary fixed right-root amalgams, the resulting count either equals \(H_A\) or is at least
\[
H_A+\min(2H_A,U_A-d_A). \tag{3}
\]

**Proof.** Rearrange the inherited amalgamation formula as
\[
\begin{aligned}
H(A\prec_r B)-d_BH_A
={}&(E_A-d_A)(S_B-d_B)+d_A(H_B-d_B)\\
&+(U_A-d_A)(V_B-d_B).
\end{aligned}
\]
Every term on the right is nonnegative. If \(d_B\ge3\), the count increase is at least \(2H_A\). If \(d_B=1\) and B's root is not a global source, then \(V_B-d_B\ge1\), giving an increase at least \(U_A-d_A\). If B's root is a global source and \(d_B=1\), the count is unchanged and \(U\) is multiplied by \(K_B\ge1\), while d is unchanged.

The only other count-preserving operation is O, which doubles U and leaves d fixed. Z and J increase H by \(U-E\ge U-d>0\). Thus before the first positive increase, H remains fixed and \(U-d\) never decreases. The first increase satisfies (3), and every subsequent operation is nondecreasing in H. ∎

For each new 16-vertex source above, the threshold in (3) is 115156. Both old target counts lie strictly between the new count and that threshold. Therefore **no positive O/Z/J/right-amalgamation continuation can restore either old target from either new source**. This statement includes every fixed right gadget, not merely bounded searches.

This gives the precise transfer boundary from the old sparse-repair and square-pair constructions. A shared-port block is a real operation and composes safely, but the target-restoration gap must also be paid at its actual entry point. In this chart the gap cannot be paid by its own core toggle or by any later positive continuation. A viable all-target gate must add a source-changing or debt-carrying interface outside this chart; the all-prefix padding families alone do not supply that interface.

## Controls

`port_probe.py` checks the first physical four-corner example by subset path counting, verifies the block equations over a declared finite alphabet, and groups finite port fibres as a control for their proved ordering. The all-length claims above follow from the formulas and inequalities, not the finite search.
