# Target-first compiler for the ordered five/eleven carry bank

The [subset-rank construction](f1.md) produces exact normalized carries,
and a separate calculation makes all baseline residues available. An
absolute target still needs a baseline in the **actual** winning window.
Here is a complete witness-or-rejection procedure at an absolute target N
for this bank, not just a check that the residue has some witness.

Take any positive numbers a,b of matched five-path and eleven-path
punctures respectively, in that order. Choose two nonempty transitive
wings with arbitrary positive binary addresses x,y. Compile the
common-apex union, with the low options, into (h,f,P,Q,S,T,K).
S,T are nonnegative endpoint counts, the negatives of signed responses.
For two successive nonempty compiled groups C,D,
\[
\begin{split}
 h'&=h_Ch_D,\qquad f'=f_Cf_D,\\
 P'&=f_DP_C+h_CP_D,\qquad Q'=f_CQ_D+h_DQ_C,\\
 S'&=h_CS_D,\qquad T'=h_DT_C,\\
 K'&=h_DK_C+h_C(K_D-h_D)+P_CQ_D+S_CT_D.
\end{split}
\]
The start/end formulas split a path at the ordinal cut, and the
path-count formula is the two-block case of the
[common-apex expansion](../chanpak_movable_reader_returns_20261114/INDEPENDENT_RETURNS.md).
Both six-coordinate local options have the same exposed data, so
compiling one reference graph suffices for the entire switch bank.

For given N:
1. Enumerate only a,b with 1024^a·3200^b≤N, because the positive
   two-wing cross term fxy is already at least f. This bounds both
   indices, and does not cut off large wing imbalances.
2. Test each integer d in [0,11a+10b] by the exact bounded-coin inverse
   for 11i+10j=d. Reject a d with no pair of actual subset ranks. For
   the saturated window only, restrict to a≥9, b≥10 and the central
   interval [90,11a+10b−90].
3. Solve the baseline exactly at
   \(B=N-2\cdot5^{a-1}11^{b-1}d\). The
   [two-wing inverse](../chanpak_movable_reader_returns_20261114/bilateral.md)
   tries the smaller positive address u≤sqrt(B/f), then peels
   low bits of the larger address until its remaining physical readout
   is strictly increasing. Binary search the larger address, and accept
   only equality. This considers both choices of the smaller wing.
4. Output the ordinal tournament with the prescribed wings and with
   i,j actual high choices in the two blocks. Their positions can be
   picked arbitrarily, since the gain depends only on the two ranks.

Every successful decode supplies a physical graph. Conversely every
graph in this bank has its a,b,d and smaller wing address in the
search. We need produce one rank pair per d: its exact summed gain
is the same even if some d has several representations.

For example the two opposite choices with a=b=1 and x=9,y=13 give
422,529,707 (prime) and 422,529,709 = 7·60,361,387 (prime cofactor).
The physical graph is computed by the compiler; the counts use the
proved transfer formulas and individually verified center profiles
rather than a full size-28 subset DP. Trial division certifies the
primality assertions. These two exact targets illustrate that the
decoder is not limited to targets with an obvious factorization.

Run **python3 compile_absolute.py N** for a target. It materializes the
outgoing masks for the first witness or reports rejection for this bank.
The **--saturated** option restricts to the proved central carry window.
Running with no target checks exact inverses at several absolute
targets, including strongly unbalanced wings in both directions.
There is no claim that an arbitrary admissible N avoids rejection.
