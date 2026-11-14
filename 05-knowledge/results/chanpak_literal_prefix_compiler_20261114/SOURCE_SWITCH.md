# Changing the physical source by swapping two corridors

Let \(M_q=D_2C_q\) be the physical step of length \(q\) with
next width 1 (width-power 2); let \(w=(47,56)\). The physical
length-7 step with terminal width 2 reads \(7w\) since
\(C_7(1,4)^t=7w^t\).

Given any physical source \(z=(a,b)^t\), append \(M_q\) with
\(q=10\) or \(14\), then append \(M_1,M_2\) in chronological
order \((2,1)\) or \((1,2)\), and finish with that length-7
terminal. The two cofactors are
\[
 p_{21}=wM_1M_2M_qz,\qquad p_{12}=wM_2M_1M_qz .
\]
The commutator of the transfer products gives
\[
 M_1M_2-M_2M_1=\begin{pmatrix}0&1\\-2&0\end{pmatrix},\qquad
 w(M_1M_2-M_2M_1)D_2=-2(56,-47).
\]
Thus the exact rows \((56,-47)C_{10}=(-1,0)\) and
\((56,-47)C_{14}=(1,0)\) give respectively
\(p_{21}-p_{12}=2a,-2a\). This works for every physical
source and uses no rank-one corridor. A commutator of different,
legal forward contexts manufactures the correction, rather
than trying to return a difference of two fixed states by a
common continuation.

The four exact cofactor rows are
\[
\begin{array}{c|r|r}
q&(2,1)&(1,2)\\\hline
10&(3805977,4534839)&(3805975,4534839)\\
14&(43557897,51898713)&(43557899,51898713).
\end{array}
\]
The \((2,1)\) branch forces a factor 3. The swap removes that
factor exactly when \(3\nmid a\), leaving the coefficient of
\(b\) unchanged.

In particular the genuinely physical single-gap source
\(z=(1,2^r)^t\), \(r\ge1\), gives an unbounded baseline
with a correction of exactly two in the cofactor. The literal
whole graph is specified by this gap word:
\[
(r,0^q,1,0^i,1,0^j,1,0^7,2),\quad (i,j)\in\{(2,1),(1,2)\}.
\]
For a prescribed total \(N\), test \(7\mid N\), then
\((N/7-A)/B\) for a positive power of two at least 2 in each
row \((A,B)\). On the swapped branch \(q=10,r=3\) gives
\(H=7\cdot40084687=280592809\); \(q=14,r=8\) gives
\(H=7\cdot13329628427=93307398989\). Both cofactors are prime
by trial division, and their companions differ by the specified
exact two and have a factor 3. This supplies actual target
witnesses but not arbitrary cofactor coverage or infinitely
many primes along the rails.

## Two independent widths, with an explicit obstruction to primitivity

The words \((1,r,0,1,0,s)\) and \((r,0,1,1,0,s)\), for
\(x=2^r,y=2^s\) and \(r,s\ge1\), have ordinary counts
\[
 H_+=(6x+5)(6y+5)+2,\qquad H_-=(6x+5)(6y+5).
\]
Their source vectors before the last variable-width step
are \((3+6x,12+12x)\) and \((5+6x,10+12x)\).
Their difference is \((-2,2)\); applying \(C_1\) produces
\((2,0)\) and the last free width \(s\) cannot amplify it.
Factoring \(N-2\) or \(N\) provides the complete inverse on
each chart. The first branch always has a factor 3 and
the second is a genuine product of two counts at least 17.
Thus genuine two-width elasticity has been constructed,
but not simultaneously with a primitive prime-producing
rail. The preceding commutator gives primitive cofactors
with one free width when \(a=1\).

## A full equal-retained-coordinate unit is impossible for just a free first width

For a physical matrix word \(P\), including an empty one,
\(|\det P|\le P_{22}\) by \(|\det C_q|\le g_q\) and
positivity of products. Suppose two fixed postcontexts
\(P_i\) following a free first width, with fixed dyadic
shifts allowed, produce states with equal rank-13 retained
coordinates and constant transverse unit
\(d=\pm(112,-94)^t\) at every width. The states are
\(c_i+x\beta\), with \(c_i=P_i(1,0)^t\) and the shared
\(\beta=2^{h_i}P_i(0,1)^t\). But
\[
 |\det(c_i,\beta)|\le\beta_2,\qquad
 |\det(c_1-c_2,\beta)|=112\beta_2+94\beta_1>2\beta_2.
\]
This is impossible. At a free marker after changed fixed
prestates with first coordinates \(a_1,a_2\), the same
argument gives the necessary inequality
\(a_1+a_2\ge112+94\beta_1/\beta_2\).
This does not rule out larger changed presources, nor the
general relaxed-prefix realization problem. The scalar
commutator above instead keeps one cofactor coefficient,
not the complete state.

Run replay_source_switch.py for an independent interval-packing
check, actual gap graphs, literal small Hamilton DP, and prime
and inverse checks.
