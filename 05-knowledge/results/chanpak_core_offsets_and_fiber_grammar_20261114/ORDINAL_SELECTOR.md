# Exact ordinal offset selectors: nonlocal returns and valuation barriers

For a marked core \(B\), use the coordinates \(h,f,P,Q,S,T,E\) of
[the bilateral formula](../chanpak_movable_exact_windows_20261114/DIAGONAL_WINDOWS.md). Thus \(h=H(B)\) and \(f\)
counts all ordered complementary path pairs, allowing an empty path.
The apex beats precisely the contact set \(U\). \(P\) counts pairs whose
second path starts in \(U\); \(Q\) counts pairs whose first path ends outside
\(U\); \(S\) counts whole-core paths ending in \(U\); \(T\) counts whole-core
paths starting outside \(U\); \(E\) counts pairs whose first path ends and
second starts in \(U\). Put \(D=f-h-Q\), \(b=h-T\), and \(K=h+P-S-E\).

## Composition and the neutral-bank criterion

For ordinal sum \(A\to B\) with the union contact set,
\[
\begin{split}
h&=h_Ah_B,\quad f=f_Af_B,\\
P&=P_A f_B+h_A P_B,\quad Q=Q_B f_A+h_B Q_A,\\
S&=h_A S_B,\quad T=T_A h_B,\\
E&=h_B E_A+h_A E_B+P_A D_B+S_A b_B .
\end{split}
\]
Each old path proceeds through the blocks in ordinal order, and its path
choices in distinct blocks are independent. For the last identity, distinguish
whether the first path meets \(B\) and whether the second meets \(A\).
A family sharing the six coordinates other than \(E\) is called a neutral
bank here. Changes chosen independently within \(k\) such blocks have exact
\[
H-H^0=-\sum_{i=1}^k \frac{J}{h_i}(E_i-E_i^0),\qquad J=\prod_{i=1}^k h_i.
\]
This holds for the physical tournament at every pair of wing parameters,
because all other coordinates survive composition unchanged. It also holds
with inserted fixed marked cores by including their \(h\)'s in \(J\).
In particular **arbitrary transitive spacers and arbitrary contacts on them**
do not change the offset difference.

With \(d_i=(E_i-E_i^0)/2\), an integer by parity of the marked apex count,
the exact gcd of attainable differences in odd-step units is
\[
\gcd_i \left((J/h_i)\gcd\{d_i\}\right)
\]
when no fillers are displayed (include a singleton zero digit set for fixed
fillers). Varying a block alone proves exactness of the gcd statement.
In particular \(G=J/\operatorname{lcm}(h_i)\) divides all differences. If all
moduli \(h_i\) are pairwise coprime, the normalized residue support modulo \(J\)
is a Cartesian product, under the CRT and unit rescaling, of the digit supports
modulo \(h_i\). Its fraction is exactly \(\prod |\{d_i\bmod h_i\}|/h_i\).
For arbitrary \(h_i\), after dividing by \(G\), the normalized \(p^j\)-image
depends only on the blocks of \(p\)-valuation more than
\(\max_i v_p(h_i)-j\). Thus adding coprime moduli does not remedy local digit
omissions. In particular repetitions of \(h=9,45\) banks have
\(G=9^{a+b-1}5^{\max(b-1,0)}\) for \(a+b>0\).

## A nonlocal exact return with changing cores

Use these literal witnesses. Row \(i\) encodes outgoing neighbors by
bits; the final field is the contact mask \(U\).

| \(h\) | Common \((f,P,Q,S,T)\) | \(E\) | Outgoing rows; \(U\) |
|---|---|---:|---|
| 65 | (976,206,570,28,52) | 40 | (2,20,9,35,77,23,47); 17 |
| 65 | same | 42 | (2,20,33,7,45,75,31); 33 |
| 109 | (2848,371,1157,109,101) | 142 | (130,172,169,209,135,217,151,0); 130 |
| 109 | same | 144 | (0,9,147,101,171,199,23,75); 9 |
| 161 | (3104,681,2262,35,126) | 90 | (2,12,17,37,139,151,63,79); 130 |
| 161 | same | 92 | (2,12,17,37,139,151,63,79); 36 |

Select local \(E\)'s \((42,142,90)\) and \((40,144,92)\). Then
\[
(E'-E)/2=-109\cdot161+65\cdot161+65\cdot109=1.
\]
Each change alone has size at least 14,170 in the physical readout.
The combined change is exactly two, for *all* independent wing heights
and all insertions of arbitrarily marked transitive spacers, without
changing any of the six interface coordinates within a given placement.
The Diophantine mechanism is
\((a-c)(b-c)=c^2-1\), equivalent to \(c(a+b)-ab=1\):
with \(c=65\), \(a-c=44\), \(b-c=96\), the factors multiply to \(4224\).
More generally any such factorization with three neutral unit-digit banks is
an exact selector. The result does not prescribe an arbitrary absolute
baseline or show that these pairs overlap across placements.

## Absolute 3-adic obstruction for mixed repeated blocks

For ordered blocks in arbitrary contacts, writing \(h_{<j}=\prod_{r<j}h_r\)
and similarly for \(h_{>i}\), the closed formula is
\[
E=\sum_i E_i\prod_{r\ne i}h_r
 +\sum_i S_i b_{i+1}\prod_{r\ne i,i+1}h_r
 +\sum_{j<i}P_jD_i h_{<j}h_{>i}\prod_{j<r<i} f_r .
\]
Classify the first path's end block \(i\) against the start block \(j\) of
the second path. When \(i=j\), the shared block has the \(E_i\) rule. When
\(i<j\), no middle block could be covered, so \(j=i+1\). When \(i>j\), all
middle blocks are freely split. Similarly \(P=\sum_i P_i h_{<i}\prod_{r>i} f_r\)
and \(Q=\sum_i Q_i h_{>i}\prod_{r<i}f_r\). The formulas for \(h,f,S,T\) follow
from two-block composition.

Consider a word of \(k\) marked banks of types:
- **9**: \(v_3(h)=2,\ v_3(f)=1\), no divisibility assumed on endpoints;
- **45**: \(v_3(h)=2,\ v_3(f)=3\), and \(P,Q,S,T\) all vanish modulo 3
  (hence \(D,b\) do also).

Arbitrary fixed marked filler blocks, unconstrained, may be interspersed.
Assign a marked block type 9 costs \((c_i,d_i)=(2,+1)\), type 45 costs
\((c_i,d_i)=(1,-1)\), and exterior sentinels \(0,k+1\) cost \(c=0\). Put
\[
\Lambda=\max\left(2,\max_{0\le i<j\le k+1}
\left(c_i+c_j+\sum_{i<r<j}d_r\right)\right).
\]
Then every full readout has \(v_3(H)\ge\max(0,2k-\Lambda)\). The base valuation
is 2 per marked block. A free interior split changes that requirement by
\(-d_i\); an endpoint requirement can lose at most \(c_i\); a local \(E_i\)
term can lose at most 2. The two-boundary, one-boundary and non-boundary
intervals control \(E,P/Q,f\) as displayed in their exact formulas.
If an actual interval ends in a filler and has a marked block strictly
inside, move that end to the innermost included extreme marked block:
it replaces a cost \(d_i\) by \(c_i\ge d_i\), a safe upper bound on loss.
An interval having zero or one marked block has loss at most 2. Other
filler factors have nonnegative valuations and cannot weaken the bound.

For a pure type-45 tower this gives \(v_3(H)\ge2k-2\), even with arbitrary
fillers and all wing parameters. For a pure type-9 tower, if \(k\ge3\) and
\(P\) of the first and \(D\) of the last marked block are units modulo 3,
then \(v_3(H)\ge k-2\), **with equality exactly when** the product of \(h\)
for fillers outside the extreme marked pair and \(f\) for fillers inside
it is a 3-unit. Indeed, total loss \(k+2\) requires both extreme marked
endpoints and all \(k-2\) middle marked blocks freely split. That is a
unique term, with sign minus in the literal readout. Every other term,
including wing terms, loses at most \(k+1\). Thus modular control survives
arbitrarily many changing one-way fillers, rather than a fixed-order
numerical census.

## Replay

Run "python ordinal_controls.py". It checks each literal core and the
physical readout independently by a 3-by-3 endpoint-status transfer.
For every block count pairs indexed by first-path end status and second-path
start status (absent, valid, invalid). Composition takes the last nonempty
first-path end and the first nonempty second-path start. This replay does
not use the derived rank formulas in its direct count, and verifies sample
compositions with and without fillers, the Egyptian shift, mixed lower
bounds and pure-9 exactness.
