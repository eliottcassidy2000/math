# Two sliding cores give arbitrarily wide consecutive odd windows

**Theorem.** There are genuine, arbitrarily wide intervals of consecutive
odd Hamiltonian path counts, all carried by **strong tournaments at the
same order** and the same two-sided painted-spine interface. In fact for
every sufficiently large order \(n\) there is such an interval of at least
\((145/12)n-O(1)\) consecutive odds. This is growth of the width, not a
claim of arbitrary placement.

This is different from continuing a fixed window at infinitely many
heights. One counts how many times a core can slide across its matched
tape on either side: that integer is a genuine new degree of freedom.
Each slide changes the count and returns the same observable interface.
The two mechanisms below have coprime weighted sizes, so their *bounded*
integer moves eventually overlap.

## Two literal pumps

Each pump uses a strong core with one fixed mark mask, and two equal-length
transitive words \(z,w\). Write \(Cz\) for core then word; the apex is
common to everything. The notation and exact slide theorem are in
[Matched slides](MATCHED_SLIDE.md). Put a one-bit guard \(p\) on the left
and \(q\) on the right, where \(p\) is the last bit of \(w\) and \(q\)
the first bit of \(z\). The physical spines are
\(p\,w^j C z^{k-j}\,q\), \(0\le j\le k\). Their full profile
\((h,Q,S,d,t,f)\) is fixed at each \(k\).

| \(h\) | Outgoing rows | Mask | \(z\) | \(w\) | Absolute count step |
|---:|---|---:|---|---|---:|
| 13 | (6,12,24,1,11) | 2 | 101001000010 | 010000001101 | 10 |
| 29 | (70,8,2,5,15,31,62) | 49 | 0001000011010010111100100010 | 1001010001000011000010110011 | 6 |

Each row is independently replayed using exact local path counts and
ordinal transfer. The first left guard is marked 1 and the second right
guard is marked 0. Thus when joined with one apex the output is strong:
the apex reaches the first vertex of the spine, and the last vertex
can return to it. All intermediate vertices lie between these ordinal
guards; this also holds with any additional transitive spine between
the pumps.

## Absolute decoder and exact interval

Let \(A_k\) and \(B_l\) be the **least** local apex counts in the first
and second guarded families at heights \(k,l\), respectively. They are
explicitly computed by forming the guard words, counting the constant
cores once, and using the exact ordinal transfer; equivalently run
the formulas in the replay, comparing the two ends of each slide.
Put
\[
N_0=29A_k+13(B_l-29)+t_B S_A+f_B d_A.
\]
The last terms depend only on the matched profiles and hence do not
change across either pump. The complete attainable set of this ordered,
two-pump bank is
\[
N=N_0+2(145i+39j), \qquad 0\le i\le k,\quad 0\le j\le l .
\]
The loss of a factor through scaling is not ignored: the individual
half-steps 5 and 3 become \(29\cdot5=145\) and \(13\cdot3=39\).
Their gcd is still one.

Consequently for \(k\ge38,\ l\ge144\), every half-offset between
\(C=144\cdot38=5472\) and \(145k+39l-C\) is realized, and the adjacent
points just outside that particular interval are absent. This is the
bounded two-generator case of [the exact CRT amplifier](CRT_AMPLIFIER.md),
with \(h_1=39,\ h_2=145\) in its **arithmetic** role: the local steps
happen to be different from the two physical spine counts. In particular
at \((k,l)=(38,144)\) one gets 183 consecutive odds. For each target \(N\)
subtract the physical \(N_0\) and put \(y=(N-N_0)/2\). Solve
\(i\equiv145^{-1}y\bmod39\), \(j\equiv39^{-1}y\bmod145\) and distribute
the remaining carries in the available integer ranges. This returns
literal slide positions; rejecting a carry means no witness in this bank,
not no tournament.

The order is \(17+12k+28l\). Fix \(l=144\), increase \(k\), and, if
needed, insert up to eleven transitive singletons between the pumps.
Their count is one, so all differences and the shared full profile
survive, while the baseline is recomputed by transfer. This gives the
stated every-order lower bound on the number of consecutive odds.

**Remaining placement problem.** The growing absolute interval is still
short compared with how far its endpoints move as the tapes grow. For
these two pumps alone (bounded filler), the strong-order bound implies
\(k,l=O(\log X)\) whenever \(N\le X\), so only \(O((\log X)^3)\) distinct
counts through \(X\) can come from this construction. It cannot prove
cofiniteness as it stands. But it supplies exactly the resource missing
from the fixed-width construction: a **portable** window with arbitrarily
large radius. Try using it as an inner repair, with both endpoint
observables intact, while solving placement with a different outer
reservoir. The quantitative aim is to get the residual absolute distance
inside the window (linear in tape length), not just into its residue class.

## Why the most obvious placement reservoir is still insufficient

For one fixed \((k,l)\) take no filler. The outgoing signed start is
zero because the last guard is 0. An arbitrary common transitive suffix
has exact translation \(St+377\,p(t)\) in the
[binary-address convention](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md),
where \(p(t)\ge0\). This turns placement at prescribed \(N\) into the
following exact test for the central interval:
\[
N-N_0-2(M-C)\ \le\ St+377p(t)\ \le\ N-N_0-2C,
\qquad M=145k+39l.
\]
The parity is automatic, since \(S\) is even when the signed start
vanishes, and \(p(t)\) is even. The integer inequalities and the
preceding CRT inverse give an actual tournament whenever they succeed.

There is a quantitative reason to try a *different* outer reservoir.
The initial marked vertex gives \(S\ge 2^{12k+28l}\): for each subset
of the \(12k+28l\) transitive tape vertices, there is one path from that
initial vertex using precisely those vertices in order, counted with a
positive complement weight. The number of suffix addresses that can
contribute through \(X\) is at most \(1+X/S\). Multiplying by the central
width, bounded by \(145k+39l+1\), and summing over all admissible pump
heights gives a coefficient of \(X\) less than
\[
\sum_{k\ge38,\ l\ge144}(145k+39l+1)\,2^{-12k-28l}
<2^{-4470}.
\]
The last inequality follows by two geometric sums and their first
moments. The \(+1\) errors sum to \(O((\log X)^3)\), since the
pre-suffix strong count is at most \(X\) and Busch bounds its order.
This is a very thin family even after trying every such common suffix.
The growing repair radius is real, but an outer reservoir must be able
to shift the baseline **without paying this same exponential cost**.
An absolute-target proof needs both an amplification and an economical
placement theorem; this construction cleanly separates them.

## The next inverse problem, with a changing left reservoir

Place a **different** painted reservoir \(L\) to the left of the entire
two-pump bank, still with one apex. Let \((g,n,s,d)\) be its actual
ordinary count, apex count, total start and signed start. Write
\(a=377\), and fix the chosen bank's common unmarked end response
\((t,f)\). If its central portable interval is \([B_-,B_+]\) by twos,
the actual intervals in this new layout are
\[
a n+t s+f d+g(B_--a),\
a n+t s+f d+g(B_--a)+2g,\ \ldots,\
a n+t s+f d+g(B_+-a).
\]
Thus the exact upstream admission condition on a proposed \(N\) is
\[
a n+t s+f d\equiv N\pmod g,\qquad
B_- \le a+\frac{N-a n-t s-f d}{g}\le B_+.
\]
Here all counts \(a,g,n,N\) are odd; \(ts+fd\) is even (write
\(s,d\) and \(t,f\) as the sums and differences of even and odd parts),
so after division the parity of the required local count is correct.
This test is sufficient, and necessary for admission through the
**central interval**. It is not necessary to neighboring attainable
offsets outside that interval.

This form pays for the repair by multiplying its radius by the
reservoir's ordinary count \(g\), while spelling out the extra residue
debt. Merely taking \(g=1\) by lengthening a transitive word cannot use
that leverage. For a prime or prime-times-seven least-missing target,
choose \(g\) and the actual selective marks in \(L\) together, or use a
second left repair palette to adjust this residue without moving the
right signature. No unconditional supply of these \((g,n,s,d)\) is
claimed. The useful concrete prediction is that future search should
measure the **integer interval defect after division by \(g\)**, not
the raw gap \(N-N_0\) or the distance to a congruence class.
