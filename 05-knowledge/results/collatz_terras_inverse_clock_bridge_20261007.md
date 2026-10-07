# From Terras count returns to synchronous inverse-clock debt returns

**Status:** PROVED in the Haar model, with FINITE-EXACT transition and
indexing controls. The alternatives are aligned merging or recurrent
bounded debt. Almost-sure merging, a return-time asymptotic, and universal
positive-integer convergence are not concluded.

**Reproduction:** run the matching experiment
    python 04-computation/experiments/collatz_terras_inverse_clock_bridge_20261007.py
and repeat with the interpreter option -O. Matching deterministic output
is saved beside this note.

## 1. What the incoming theorem proves

The inherited mechanism is
[THM-4569 — the Terras clock](../../01-canon/theorems/THM-4569-the-terras-clock-recurrence-of-two-collatz-orbits-is-unconditional.md),
specifically its elementary transition table and skipped-random-walk
dichotomy. The numerical box bounds and additional orbit-equivalence
statements are not needed here.

Let
\[
T(z)=\begin{cases}z/2,&z\text{ even},\\(3z+1)/2,&z\text{ odd}.\end{cases}
\]
At equal T-time suppose \(x_s=3^{j_s}y_s+c_s\), with \(c_s\in\mathbb Z_2\).
If \(p_s=y_s\bmod2\) and \(\epsilon_s=c_s\bmod2\), then
\[
j_{s+1}-j_s=\epsilon_s(1-2p_s).
\]
The translation updates in the four cases are
\[
\begin{array}{c|c|c}
(p,\epsilon)&j'&c'\\ \hline
(0,0)&j&c/2\\
(0,1)&j+1&(3c+1)/2\\
(1,0)&j&(3c+1-3^j)/2\\
(1,1)&j-1&(c-3^{j-1})/2 .
\end{array}
\]

For Haar \(y\), the successive parity bits are independent fair coins.
The disagreement bit is determined by earlier bits. Consequently the
nonzero increments of \(j\) are a simple symmetric random walk, selected
at predictable times. If only finitely many disagreements occur, the
future parity vectors agree forever, forcing the two current states to
be equal. This last implication follows directly because a length-r
parity vector specifies one residue modulo \(2^r\).

Thus, almost surely, either the aligned states merge, or \(j\) visits
every integer infinitely often. A finite initial parity condition only
forces a finite prefix; the same conclusion holds afterward.

This is a recurrence statement in a **count clock**. It does not identify
\(j_s\) with an odd-step valuation sum difference. The missing transfer is
proved below.

## 2. Exact alignment for the odd affine pair

Let \(Y\) be odd Haar, let \(v\ge1\) be fixed, and put
\[
X=3\cdot2^vY+1.
\]
Write \(a_i,b_i\) for the ordinary odd-map valuations of \(Y,X\), and
\[
A_k=\sum_{i=1}^k a_i,\qquad B_k=\sum_{i=1}^k b_i,\qquad A_0=B_0=0.
\]
The odd-map states \(U^kY,U^kX\) occur at T-times \(A_k,B_k\).
ROOT is not made absorbing in these definitions.

Let \(N_Y(s)\) count the odd states of the Y-orbit at T-times \(0,\ldots,s-1\);
define \(N_X\) similarly. Align
\[
y_s=T^sY,\qquad x_s=T^{s+v}X.
\]
The first v T-parities of X are those of 1, because \(X\equiv1\pmod{2^v}\).
In particular,
\[
x_0=3^{1+\lceil v/2\rceil}Y+d_v,
\qquad d_v=\begin{cases}1,&v\text{ even},\\2,&v\text{ odd}.\end{cases}
\]
At all later aligned times the exact coefficient identity is
\[
j_s=1+N_X(s+v)-N_Y(s).
\]
The initial Y parity is fixed odd; after that bit the Y state is full Haar
on \(\mathbb Z_2\), so the fresh-coin proof applies from time 1 onward.

Fix any integer ordinal offset \(h\). Returns to the level \(j=1-h\)
are precisely ties
\[
N_Y(s)=N_X(s+v)+h.
\]
On a nonmerging path, there are almost surely infinitely many such ties.

## 3. A one-half trial transfers count ties to inverse-clock returns

Consider a tie time \(\tau\ge1\), chosen from the past parity history.
Write \(k=N_X(\tau+v)\), so the next odd events of the two raw orbits
have ordinal indices \(k\) and \(k+h\).

If \(c_\tau\) is even, the two current parities agree. With conditional
probability \(1/2\), both are odd now. Their next odd-event times are
then \(B_k=\tau+v\), \(A_{k+h}=\tau\).

If \(c_\tau\) is odd, exactly one current state is odd. Whichever state is
even becomes odd one T-step later with conditional probability \(1/2\).
Indeed, after revealing the current Y bit, the next Y bit is fresh; the
other orbit's next parity is that bit xor a now-known disagreement bit.
The two next odd-event times then differ by exactly one after accounting
for the alignment shift v.

In both cases there is a trial chosen using only the past, with conditional success
probability exactly \(1/2\), and success gives
\[
\boxed{|v+A_{k+h}-B_k|\le1.}
\]

Choose successive level-return trials at least two T-steps apart. These
are stopping times for the parity filtration, and each trial uses only
its next two fresh bits. Conditional Borel--Cantelli, or the elementary
bound \(2^{-r}\) for failing the next r available trials, shows that
infinitely many available trials yield infinitely many successes almost
surely. This argument does not condition the coin law on eventual
nonmerging; it first proves that infinite trials with finitely many
successes have probability zero.

Successful disjoint windows correspond to strictly increasing matched
odd-event indices: both streams have advanced past the previous matched
event before the next window. We obtain the exact alternative

\[
\Pr\left(
\text{aligned merge}\ \text{or}\
|v+A_{k+h}-B_k|\le1\ \text{for infinitely many }k
\right)=1.
\]

For negative h the finitely many inadmissible indices \(k+h<0\) are omitted.
This holds simultaneously for all integer h by countable intersection.
It is the required bounded-return statement for the inverse clock, rather
than a renaming of \(j\).

## 4. The Mersenne debt indexing and the remaining obligation

For the two raw odd streams use \(h=0\), selecting the level \(j=1\).
The returned debt is \(v+A_k-B_k\).
For a paired stream beginning at \(U(Y)\), use \(h=1\), selecting \(j=0\).
If \(a_1\) is the first raw valuation and \(\widetilde A_k\) is the
cumulative valuation of that paired stream, then
\[
A_{k+1}=a_1+\widetilde A_k,\qquad
v+A_{k+1}-B_k=v+a_1+\widetilde A_k-B_k.
\]
This is the offset used by the lag-one paired debt, with its initial
constant retained. The incoming Mersenne representation
\(p=3q+2,\ q\in5+16\mathbb Z_2\), has forced parity prefix 1000;
the relation state moves from \((j,c)=(1,2)\) to \((2,1)\).
The checker independently confirms that boundary.

The recurrence half is therefore genuinely improved: outside the merge
alternative, the stated synchronous debt returns to a fixed finite band
infinitely often in the Haar model. No alignment or independence conjecture
is needed for this conclusion.

What remains is control of the translation coordinate along these returns.
A return to a fixed j-level does not bound the real size of c, and a bounded
valuation-sum debt does not make the affine identity hold. The relevant
remaining target is archimedean box recurrence or another uniform way of
turning these returns into successful merges. The merge-rate asymptotic
and the transfer from Haar-almost-everywhere to specified integer sources
are separate obligations.

The phrase “merges happen almost surely at equal T-time” in the incoming
statement is read as a statement about the clock of a merge, up to null
exceptional coincidences; it is not read as claiming probability-one
merging, which that same theorem leaves as a box-recurrence question.

## 5. Checks and scope

The local half-probability calculation is checked for 549 rational relation
states: \(-4\le j\le4\), translation numerators \(-30\) through 30, and
denominator \(3^{\max(0,-j)}\). Every one has exactly two successful
residues among the four possible two-bit source classes.

Independent integer replay uses odd \(Y=1,\ldots,127\), \(v=1,\ldots,6\),
and 60 aligned T-times per pair. It verifies 23,040 exact relation/count
states and 14,269 successful inverse-clock events, including the correct
ordinal offsets. Five malformed-type or range controls are retained.
Normal and optimized runs agree on 130,957 explicit checks.
These finite checks validate identities and boundaries; the infinite
recurrence assertion follows from the proof, not from the census.
