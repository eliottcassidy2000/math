# Exact Leech octad baselines and a second prescribed prime repair

This follows the [main packet](README.md). It gives literal lattice
vectors and a separate physical tournament witness, not a map
from one to the other.

## Orthogonal decimal baselines in the actual Leech lattice

Use the archived [Golay-code coordinates and mod-8 sum guard](../session97_overnight_arithmetic_carriers_20261029/leech/GOLAY_LEECH_SUM_GUARD.md).
The disjoint octads below partition 24 coordinates (printed 1-based):

- \(A=\{1,3,6,9,10,11,12,13\}\)
- \(B=\{5,7,8,14,15,18,19,20\}\)
- \(D=\{2,4,16,17,21,22,23,24\}\)

The three vectors \(100\mathbf1_A/\sqrt8\), \(100\mathbf1_B/\sqrt8\),
\(100\mathbf1_D/\sqrt8\) are mutually orthogonal Leech vectors of squared
norm \(10\,000\). Their sum has squared norm \(30\,000\).
Indeed their numerators are \(0\bmod4\), and their sums are
\(0\bmod8\). The lattice membership here does not require their
supports to be octads, but the indicated partition genuinely is a
partition by codewords and supplies the signed corrections below.

A balanced signed octad means entries \(\pm2/\sqrt8\), four of each.
Take \(a,b,r,s\) on supports \(A,B,R,S\) with the indicated plus sets:

| Vector | Support if new | Plus set |
|---|---|---|
| \(a\) | \(A\) | \(1,6,10,11\) |
| \(b\) | \(B\) | \(5,7,8,14\) |
| \(r\) | \(2,3,4,5,7,8,10,15\) | \(2,4,5,10\) |
| \(s\) | \(1,2,3,4,6,7,9,14\) | \(1,2,4,6\) |

All four are lattice vectors by the Golay and sum guards, of squared
norm 4. The inner products \(\langle a,b\rangle=0\),
\(\langle a,r\rangle=1\), and \(\langle a,s\rangle=2\) are checked
from the signs. The corrections \(a-r,a+b,a+r,a+s\) therefore
have squared norms \(6,8,10,12\). The first two are orthogonal to
the \(A\)-packet (their restrictions to \(A\) sum to zero); the last
two are orthogonal to the full vector. Thus there are explicit
Leech squared norms \(10006,10008,30010,30012\), exactly one
below the four prime tournament counts in the main comparison.
The +1 alignment is an exact comparison of separately constructed
objects, not an implication about arbitrary primes.

## Two physical switches require both a difference and a midpoint

The \((\pm4e_i\pm4e_j)/\sqrt8\) minimum vectors give signed
two-coordinate directions. On a cube of 24 chosen physical switches
one can associate these directions with two toggled coordinates, but
the lattice gives no information about the function on that cube.

In fact on a two-switch slice, writing the path count as
\(F(s,t)=c+\alpha s+\beta t+\gamma st\) for \(s,t=\pm1\),
a change along the antipodal diagonal \(--\) to \(++\) has gain
\(2(\alpha+\beta)\) and midpoint \(c+\gamma\). A prescribed
\(p,p+2\) needs **both** \(\alpha+\beta=1\) and \(c+\gamma=p+1\).
The direction tells us neither equality.

For a second target chosen before search, the 11-vertex source
\[
(1232,389,1977,979,770,1819,1462,560,1665,67,666)
\]
has \(H=30089\), and reversing the edges on pairs \((0,6),(3,4)\)
gives 30091. The first and second single flips yield 23331 and
38397. Both endpoints are prime, and the four-corner parameters are
\((c,\alpha,\beta,\gamma)=(30477,-3766,3767,-387)\).
The first-order effects sum to \(+1550\) and their mixed correction
is \(-1548\). These are disjoint switches, unlike the two edges in
the main packet. The upper count 30091 has no one-apex realization
at any order by the same spine-length cutoff and complete replay;
30089 itself has a converse pair of one-apex realizers. This is a
second exact placement, not a completeness result.

Run python3 octads.py to check the lattice vectors, prime counts,
two flips, and one-apex census. The code constructs the Golay span
from the archived generator and checks every vector against the
three independent guards.
