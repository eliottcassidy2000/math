# Observer-preserving signed lifts and composite-table families

Status: **PROVED in this note** for the arbitrary-parameter statements below. The exact table dictionary and common-tail lift are inherited; the marked-observer period and simultaneous composite-table corollary are the new synthesis. No convergence or complete cycle classification is asserted.

## 1. The recovered table map

The table convention in [the October 3 atom note, §3](../../tiling_modular_atoms_20261003.md) uses nonzero row and column indices \(1,\ldots,m-1\). For odd \(m\), its central block is
\[
\begin{pmatrix}d&c\\c&d\end{pmatrix},\qquad
d=4^{-1}\pmod m,\quad c=-d\pmod m,
\]
with least positive representatives. Explicitly,
\[
m=4n+1:\ \begin{pmatrix}3n+1&n\\n&3n+1\end{pmatrix},\qquad
m=4n-1:\ \begin{pmatrix}n&3n-1\\3n-1&n\end{pmatrix}.
\]
Primality is not needed for these formulas.

For \(\sigma\in\{-1,1\}\), write \(m=4n+\sigma\). The [current table packet](../../modular_centers_collatz_pell_20261111/PROOFS.md) records the exact raw map and its shortcut version. We use the accelerated odd map
\[
U_\sigma(n)=\frac{3n+\sigma}{2^{v_2(3n+\sigma)}}
\quad(n>0\text{ odd}).
\]
The inherited common-target lift is
\[
R_\sigma(n)=4n+\sigma,\qquad
U_\sigma(R_\sigma(n))=U_\sigma(n).
\]
On the table modulus this becomes \(m\mapsto4m+\sigma\). The source, sign, and chosen terminal are retained throughout.

## 2. A single marked table factor has an exact return clock

**Lemma 1.** Let \(m=4n+\sigma\), and let an odd \(q\ge3\) divide \(m\). Along the common-target lift, the j-th table modulus is
\[
m_j=4^j m+\sigma\frac{4^j-1}{3}.
\]
Then
\[
q\mid m_j\quad\Longleftrightarrow\quad
4^j\equiv1\pmod{3q}.
\]
Thus the exact positive return period is \(\operatorname{ord}_{3q}(4)\).

**Proof.** Reduce the formula modulo \(q\), using \(q\mid m\). Divisibility is equivalent to \(q\mid(4^j-1)/3\), hence to \(3q\mid4^j-1\). ∎

This works when \(q=3\) or contains a power of three; replacing the modulus \(3q\) by \(q\) would lose the valuation at three. Taking \(q=m\) turns any seed table, including a prime seed, into an infinite sequence of composite tables with the same next odd target. For example the plus-sheet source 1 has table 5; period \(\operatorname{ord}_{15}(4)=2\) gives sources \(1,21,341,\ldots\), whose tables \(5,85,1365,\ldots\) are all divisible by 5.

## 3. Exact least period for an entire marked head

Fix an actual positive odd trajectory prefix on one sheet,
\[
n_0\xrightarrow{a_1}n_1\xrightarrow{a_2}\cdots
\xrightarrow{a_r}n_r=h,
\qquad 3n_{i-1}+\sigma=2^{a_i}n_i,
\]
where every \(a_i\ge1\). Set \(A_i=a_1+\cdots+a_i\), \(A_0=0\), \(A=A_r\), and
\[
S=\sum_{i=0}^{r-1}3^{r-1-i}2^{A_i},
\qquad 3^r n_0+\sigma S=2^A h.
\]
For \(r=0\), take \(A=S=0\). Choose an odd positive observer modulus \(Q_i\) independently at every marked time \(i=0,\ldots,r\).

Define
\[
M_i=\frac{3^{r-i+1}Q_i}{\gcd(Q_i,3h+\sigma)},\qquad
M=\operatorname{lcm}_{0\le i\le r}M_i,
\qquad d=\operatorname{ord}_{M}(4).
\]

**Theorem 2 (marked observer lift).** A positive lift increment \(j\) preserves the exact r-step valuation head and all congruences \(n_i(j)\equiv n_i\pmod{Q_i}\) if and only if \(d\mid j\). Here the new endpoint is
\[
h_j=4^j h+\sigma\frac{4^j-1}{3},
\]
and the preceding states are recovered backwards using the original exponents \(a_i\). Every such state is positive odd, and \(U_\sigma(h_j)=U_\sigma(h)\), so the entire tail after the splice is unchanged.

For \(j=kd\), the source family has the one-state affine recurrence
\[
n_0^{(k+1)}=R n_0^{(k)}+K,
\quad R=4^d,
\quad K=\sigma\frac{(R-1)(2^A+3S)}{3^{r+1}}.
\]
The increment at the splice valuation is exactly \(2kd\); all preceding valuations remain \(a_1,\ldots,a_r\).

**Proof.** Backwards affine reconstruction gives the exact difference
\[
n_i(j)-n_i
=\frac{2^{A-A_i}(4^j-1)(3h+\sigma)}{3^{r-i+1}}. \tag{1}
\]
Since \(3h+\sigma\) is a three-adic unit and every \(Q_i\) is odd, the condition that (1) is an integer divisible by \(Q_i\) is precisely
\[
M_i\mid4^j-1.
\]
The simultaneous conditions are equivalent to \(M\mid4^j-1\), hence to \(d\mid j\). This proves both sufficiency and exact minimality, including the integrality needed by the head.

Each numerator \(3h+\sigma\) is positive and even. The differences in (1) are therefore positive even integers for \(j>0\); every reconstructed state is positive odd. Its edge equation has an odd successor, so its valuation is exactly the original \(a_i\). At the splice,
\(3h_j+\sigma=4^j(3h+\sigma)\), proving the common tail and the extra \(2j\) halvings. Substituting the endpoint lift into \(3^r n_0+\sigma S=2^A h\) proves the affine recurrence. Integrality of \(K\) follows from \(3^{r+1}\mid R-1\), already forced by \(M_0\). ∎

For \(Q_i=1\) at all marks, \(M=3^{r+1}\) and \(d=3^r\), exactly the inherited minimal head-preserving lift from [the October route-memory note, §4](../../collatz_atom_route_memory_20261003.md). The new coordinate is the family of retained observer residues; the gcd terms remove genuinely redundant factor constraints.

## 4. Any marked finite route has an infinite composite-table copy family

**Corollary 3.** In Theorem 2 choose
\[
Q_i=m_i=4n_i+\sigma.
\]
For every \(k\ge1\), each new marked table modulus satisfies
\[
m_i^{(k)}\equiv0\pmod{m_i},\qquad m_i^{(k)}>m_i\ge3.
\]
Hence every marked table in the new family is composite, with its original table modulus as a proper divisor. The exact accelerated prefix and shared terminal tail are retained.

**Proof.** Multiplication by 4 preserves the observer congruence, so
\(4n_i^{(k)}+\sigma\equiv4n_i+\sigma=0\pmod{m_i}\). Strict increase follows from (1). ∎

This is an explicit infinite-family compression: store one actual finite head, its chosen tail, the observer moduli, and one affine recurrence. It requires no construction of the full multiplication tables. Computing the least multiplicative order is a finite arithmetic task; no free factorization or universal computational speedup is asserted.

**Example.** Start with the plus route \(3\to5\to1\), and splice at \(h=5\) after its first step. Then \(r=1,A=1,S=1\), and the observers are \(Q_0=13,Q_1=21\). The two moduli are
\[
M_0=117,\quad M_1=63,\quad
d=\operatorname{lcm}(6,3)=6.
\]
The complete source family is
\[
n_{k+1}=4096n_k+2275,\quad n_0=3.
\]
Its first new route is
\[
14563\xrightarrow{1}21845\xrightarrow{16}1.
\]
The two marked table moduli are
\[
58253=13\cdot4481,\qquad87381=21\cdot4161.
\]
For every k, the valuation word is \((1,4+12k)\). The plus or minus sign and the original shared tail are part of the certificate. Using the minus-sheet seed \(5\to7\to5\) preserves that cycle's basin instead of changing it into a route to 1.

The simultaneous composite property concerns the marked odd states. The extra even states inserted into the expanded raw orbit at the splice are new states and are not claimed to inherit the old table factors.

## 5. The shortcut factor-three automaton

For the shortcut map, the table update is
\[
\widehat F_\sigma(m)=3^b(m+\sigma)/2,
\qquad b=((m-\sigma)/4)\bmod2.
\]
An odd shortcut step therefore lands in a table divisible by three. On the plus sheet this destination is at least 9 and composite. On the minus sheet the sole prime boundary is the destination 3 from input 1; every other odd input gives a composite destination.

After an odd step followed by h even shortcut steps,
\[
m_h=\sigma+\frac{m_0-\sigma}{2^h},\qquad
m_h\equiv\sigma(1-2^{-h})\pmod3.
\]
Thus the table is divisible by three exactly when h is even. Every next odd shortcut step resets this two-state phase. This is a literal finite factor automaton, not an assertion that an orbit visits prime tables with any specified frequency.

## Reproduction and boundary

Run `python3 05-knowledge/results/odd_tables_signed_carries_20261111/history/probe.py` from the repository root. The receipt records the finite route, observer, and lift universes. The general statements follow from the displayed identities.

The construction preserves a known head and tail and creates infinitely many certified copies. It does not place an arbitrary given source in one of those families. Prime/composite table structure can therefore be transported and compressed without resolving the separate global basin-coverage question.
