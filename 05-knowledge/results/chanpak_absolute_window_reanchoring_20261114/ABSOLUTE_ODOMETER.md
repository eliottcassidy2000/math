# An absolute carry test for neutral ordinal banks

**Proved algebraic/physical transfer.** This gives an exact obstruction
for arbitrary assemblies of binary neutral switches, and a sufficient
redundancy condition for a whole absolute window, if the required local
physical banks can be made. No assertion of such full local banks, or of
universal Chan–Pak coverage, is made.

Inspiration: the randomly sampled [Collatz Fourier experiments](../collatz_fourier_experiments_20261004.md)
separate real position from local digit phase and track an odometer
across compatible phases. We borrow that distinction, not a Collatz
claim. The physical input is the [ordinal offset law](../chanpak_core_offsets_and_fiber_grammar_20261114/ORDINAL_SELECTOR.md)
for marked cores with fixed six-coordinate interfaces; no signed
nonphysical choices are introduced.

## Separate a residue choice from its absolute carry

For independent neutral blocks with old path counts \(h_i>1\), let
\(J=\prod_i h_i\), \(W_i=J/h_i\), and \(D_i\) be the integer
local half-count offsets. After any fixed choice of transitive wings
or marked transitive spacers, the whole counts are exactly
\[
H^0+2\sum_i W_i d_i,\qquad d_i\in D_i.
\]
If some prime divides two \(h_i\), then \(J/\operatorname{lcm}(h_i)>1\)
divides every half-count difference, ruling out even a unit neighbor.
Assume henceforth the \(h_i\) are pairwise coprime. The target half-offset
\(t=(N-H^0)/2\) fixes a residue
\(r_i=(W_i^{-1}t\bmod h_i)\) in each block, since all other terms
vanish modulo \(h_i\).

For each residue build its *actual integer-height* set
\[
Z_i(r)=\{(d-r)/h_i:\ d\in D_i,\ d\equiv r\pmod{h_i}\}.
\]
The exact absolute target is possible if and only if
\[
(t-\sum_i W_i r_i)/J \ \in\
Z_1(r_1)+\cdots+Z_m(r_m),
\]
where the sum denotes all sums of one element per set. The quotient on
the left is an integer by CRT. The proof is direct substitution of
\(d_i=r_i+h_i z_i\). A successful choice gives literal local witnesses;
reading off just the residues removes this last condition.

In particular \(L\) consecutive odd totals require
\(\lvert D_i\bmod h_i\rvert\ge \min(L,h_i)\) in every block.
Thus **no product of binary neutral banks** with \(h_i\ge3\) can
have three consecutive odd totals: noncoprime blocks have no unit
step, and coprime blocks have at most two phases at each block.
Transitive spacers may move \(H^0\) but cannot alter this proof.

## The canonical CRT wrap is an exact clock

Now suppose each bank includes precisely the canonical digits
\(0,\ldots,h_i-1\). This is a hypothetical complete-residue test case.
Put \(u_i=(W_i^{-1}\bmod h_i)\) and let
\(d_i(t)=(u_i t\bmod h_i)\) in that range. The one actual representative
at phase \(t\) is \(F(t)=\sum_i W_i d_i(t)\).
Incrementing the requested residue gives
\[
F(t+1)-F(t)=1+J(B-C(t)),\quad
B=(\sum_i W_i u_i-1)/J,\quad
C(t)=\#\{i:d_i(t)+u_i\ge h_i\}.
\]
Indeed each wrap subtracts \(h_iW_i=J\), and
\(\sum_i W_i u_i\equiv1\pmod J\). Exact adjacent totals occur
precisely when \(C(t)=B\). Across \(t\bmod J\) the digits are independent
uniform variables by CRT, hence the exact fraction of adjacent phases
whose actual values advance by one is
\[
\Pr\!\left[\sum_i \operatorname{Bernoulli}(u_i/h_i)=B\right].
\]
A run of \(L\) consecutive totals requires \(L-1\) consecutive clock
steps all with \(C(t)=B\), not just individual local phase availability.

For instance \(h=(3,5)\) gives \(J=15,W=(5,3),u=(2,2),B=1\).
Eight of the fifteen cyclic steps are true unit steps. The actual
outputs in sorted order are
\[
0,3,5,6,8,9,10,11,12,13,14,16,17,19,22.
\]
They occupy all fifteen residue classes, but no single interval of
width fifteen has one point in every class.

This failure has a uniform height bound for all such canonical banks
with two or more blocks. The fraction of representatives in any
half-open interval of width \(J\) is the probability that
\(\sum_i d_i/h_i\) belongs to the corresponding interval of width one.
Condition on all but two blocks, and extend the two discrete uniforms
to independent continuous uniforms on \([0,1]\) by adding one uniform
mesh increment to each. An occupied interval of width one maps inside
an interval of width \(1+\delta\), where
\(\delta=1/h_i+1/h_j\). Integrating the triangular density for the
continuous sum gives maximum mass
\(1-(1-\delta)^2/4\) in that longer interval. Thus when both selected
moduli are at least three the original mass is at most \(35/36\).
In particular the full modulus-sized window cannot be made by shifting
the canonical chart to a different absolute baseline alone.

## A sufficient physical redundancy condition

For each block suppose every residue \(r\) has actual representatives
\(r,r+h_i,\ldots,r+s_i h_i\). Put \(Q=\sum_i s_i\).
The canonical representative of a target residue \(r\bmod J\) is
\(F(r)=r+c(r)J\) with \(0\le c(r)\le m-1\). The hidden quotient
choices can give every additional multiple \(aJ\) for
\(0\le a\le Q\). Consequently the physical chart realizes **every**
half-offset in
\[
[(m-1)J,\ (Q+1)J).
\]
In particular one extra height per residue in just \(m-1\) blocks
suffices for one complete modulus-sized absolute window. For example
in the arithmetic \(h=(3,5)\) example, adding \(3,4,5\) as extra
digits at the first block fills all fifteen consecutive offsets
\(15,\ldots,29\). This is an arithmetic certificate, not a claim
of physical existence at those small \(h\).

For a literal contrasting example, the known \(h=(65,109,161)\)
Egyptian unit banks have weights \((17549,10465,7085)\) and two
digits \(\{0,1\}\) each. The eight offsets include the adjacent
pair \(17549,17550\) from the signed relation
\(-17549+10465+7085=1\); that is the only unit pair. More generally
if signs \(\epsilon_i\in\{\pm1\}\) satisfy \(\sum_i\epsilon_i W_i=1\)
and each \(D_i=\{0,\ldots,\ell_i-1\}\) has \(\ell_i<h_i\),
any unit increment must change each local digit by exactly
\(\epsilon_i\): its residue must, and no wrapped representative is
available. Thus all physically placed consecutive intervals in
that chart are chains in the fixed digit direction, stopping at
the first exhausted local range. This explains why the existing
binary Egyptian exchange is a two-point repair rather than the
start of a scalable consecutive window.

The remaining obstruction to arbitrary targets is substantive:
construct actual local digit alphabets with the residue-wise extra
heights, or a moving absolute baseline whose required quotient
always belongs to the attainable sum. The criterion does not
establish either condition. Run **python3 replay_odometer.py** for
the two ideal arithmetic examples, the general finite carry checks,
and independent subset-path profiles of all six physical cores in
the Egyptian contrast.
