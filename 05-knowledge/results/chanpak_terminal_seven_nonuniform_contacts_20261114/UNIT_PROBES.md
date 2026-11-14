# Exactly which continuations recover a transverse unit

Use the known physical pair
\(z_1=(291,758)^t,\ z_2=(179,852)^t\) with difference
\(2(56,-47)^t\), and the suffix convention of
[the seven-axis note](AXIS_FRONTIER.md). A shared suffix gives
ordinary counts differing by exactly two iff its backwards row
\(R=(u,v)^t\) satisfies \(|56u-47v|=1\). Treat the width \(y\)
below as \(2^r,\ r\ge1\). The following is the full list, even
allowing arbitrary suffix length and resets.

| Steps in backwards build order | Counts from \(z_1,z_2\) |
|---|---|
| \(q=10\), arbitrary \(y\) | \(423939+505125y,\ 423941+505125y\) |
| \(q=14\), arbitrary \(y\) | \(4851811+5780875y,\ 4851809+5780875y\) |
| \(q=8,y=2\) | \(423939,\ 423941\) |
| \(q=12,y=2\) | \(4851811,\ 4851809\) |
| \((q,y)=(2,8),(4,2)\) | \(1209689,\ 1209691\) |
| \((q,y)=(1,4),(2,8),(4,2)\) | \(11985689,\ 11985691\) |

Physical chronological suffix order is the reverse of the left column.

**Proof.** Nonempty backward probes have \(u\ge3\) and odds \(v/u\)
in \([1,4/3]\), with the sharper odds bounds in AXIS_FRONTIER.
Consider the final backwards step of the probe, let its
input precede the width multiplication be \((a,b)^t\),
and put \(T=yb\ge2a\). The signed transverse readout on this
output is the dot product of \((56,-47)C_q\) with \((a,T)^t\).
For \(q=0,\ldots,16\) these rows are
\[
 (9,9),(-38,18),(27,-20),(-2,7),(-13,5),(12,-8),
 (-3,4),(-4,1),(5,-3),(-2,2),(-1,0),(2,-1),
 (-1,1),(0,0),(1,0),(0,1),(1,1).
\]
Thereafter use positivity and the recurrence. Since \(T\) is
even, the rows at \(q=1,5,7,9,11\) cannot give absolute
value 1. Rows \(0,2,3,6,13,15,16\) and later cannot by
\(T\ge2a\) or zero. At \(q=8,10,12,14\) the claimed
solutions follow from \(a=1\) forcing the initial \((1,1)^t\).
For example \(|5a-3T|\ge a\) at \(q=8\), with equality
only at \(T=2a\).

At \(q=4\), \(|-13a+5T|=1\) forces \(y=2\), since
\(y\ge4\) gives \(-13a+5T\ge7a\). The preceding probe
must have been made by \(q=2\): for the empty probe the
expression is \(-3\), and for every other type the odds are
at most \(11/9\), giving \(10b-13a\le-7a/9<-1\).
Write the \(q=2\) probe as
\(a=3U+3YV,\ b=3U+4YV\). The desired expression becomes
\(-9U+YV\), forcing \(Y=8\); for \(Y\le4\) use \(V/U\le4/3\),
and for \(Y\ge16\) use \(V/U\ge1\).
Now \(|8V-9U|=1\) holds for \(U=V=1\). If \(q=0\) made
the preceding probe, \(U=V\ge3\) fails. If \(q\ge2\) made
it, use odds at least \(7/6\) and \(U\ge9\), giving
\(8V-9U\ge U/3>1\). Thus the remaining way is \(q=1\):
write \(U=A+2ZB,\ V=2A+2ZB\). Then the expression is
\(7A-2ZB\). For \(Z=2\), this is at least \(5A/3>1\);
for \(Z\ge4\) it is at most \(-A\), so \((A,B)=(1,1)\).
The equation now gives precisely \(Z=4\).
This exhausts all continuations.

## A general boundary lift and a no-return result

A formal coefficient at terminal weight zero becomes physical
with two fewer zero gaps:
\[
 C_q(1,0)^t=C_{q-2}(1,2)^t\quad(q\ge2).
\]
Likewise \(C_q(0,1)^t=C_{q-1}(1,1)^t\) for \(q\ge1\).
Both are exact by the first three rows and the recurrence.
These explain why the \(q=8,12\) endpoints above are the
unattained-weight-zero baselines of the two longer rails, but
appear as ordinary counts in shorter words.

No nonempty common continuation can return any nonzero
difference vector to itself or its negative. Let \(d_q=\det C_q\).
The signed row recurrence for \(d_q\) gives
\(|d_q|\le g_q\), strictly for \(q\ne1\): the initial values
are \(0,-2,3\) while the \(g\) values are \(1,2,4\);
use \(d_{q+3}=d_q-d_{q+1}-d_{q+2}\) and the ordinary
positive \(g\)-recurrence. For a product \(M\) of the forward
matrices \(D_yC_q\), the bottom-right entry is at least the
product of the individual bottom-right entries, hence at
least \(|\det M|\). The top-left entry is \(\ge1\), strictly
larger for at least two factors. Therefore the trace is
strictly larger than \(|\det M|+1\), except the lone \(q=1\)
factor (where the sign is \(\det M<0\) and
\(\operatorname{tr} M=1+|\det M|\)). In either case the
characteristic polynomial does not vanish at \(1\) or \(-1\).
A shared positive transfer therefore cannot preserve the signed
difference even after a full loop of common contexts. This is
why the constructive next step must vary actual paired sources,
or use distinct forward contexts with an additional state, not
just append a clever fixed common reload.
