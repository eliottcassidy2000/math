# 3-connected Z counterexamples in every rank from six onward

2026-11-27. The [previous session](../z123_exact_rank_threshold_20261127/README.md)
proved that every matroid up to rank five has a log-concave Z-polynomial.
Its rank-six examples were 3-connected after one extra coloop and a
second truncation, but its extension to higher ranks used disconnected
padding. Here one geometric block suffices at every rank, while keeping
3-connectivity. The mechanism is an exact incidence-cycle formula.

## Main construction

Let \(A_q\) be the simple rank-four affine packing matroid from
[Bérczi–Coron–Eur–Li, Definition 4.10 and Corollary 4.11](https://arxiv.org/html/2610.17135#S4.SS1),
for odd prime powers \(q\). For \(d\ge1,\ s\ge d+2\), put
\[
M_{q,s,d}=\mathrm{Tr}^{d}(A_q\oplus U_{s,s}),\qquad r=s+4-d.
\]
Here \(U_{s,s}\) is free. Let \(z_i=[t^i]Z_{M_{q,s,d}}(t)\).
At fixed \(s,d\) the gap satisfies
\[
z_2^2-z_1z_3 \sim
-\tfrac12\binom{s}{d}\binom{s}{d+2}\binom{d+1}{2}\,q^{17}.
\]
In particular, choosing \(d=2,\ s=r-2\) gives a 3-connected
counterexample in each rank \(r\ge6\) for all sufficiently large odd
prime powers. The Fano summand is no longer needed. This is the exact
rank threshold for 3-connected counterexamples to log-concavity;
it is not a minimum-rank theorem for failure of real-rootedness.

**Connectivity.** The truncated matroid is simple. In a partition
into two sets each of size at least two, if neither spans the
truncation, the ranks are unchanged and sum to at least \(s+4\);
the connectivity value is then at least \(d\). If one spans, the
other has rank at least two. For a singleton partition, removing one
point drops the original rank by at most one, so the complement
still spans and the connectivity value is one. Thus \(d\ge2\)
ensures 3-connectivity. This does not assert ordinary
\((d+1)\)-connectivity for larger \(d\).

## The interface is the cycle rank of a point–line graph

For a general loopless rank-four block \(A\) and loopless exterior
matroid \(E\), truncate \(A\oplus E\) to rank \(r\ge6\).
Write \(t=r-6\), and assume \(\mathrm{rk}(E)\ge t+2\).
Use \(W_i\) and \(f_I\) for single flats and strict flags; use
\(f_{123}(A)\) for strict flags inside the first summand. Then the
coefficient of \(f_{123}(A)\) in the coordinate-rank expansion of
\(z_3\) is
\[
\sum_{\substack{X\in L(E)\\\mathrm{rk}(X)=t}}
\bigl(f_{12}(E/X)-W_1(E/X)-W_2(E/X)+1\bigr).
\]
Each summand is the cycle-space dimension of the bipartite
point–line incidence graph of the contraction (a point is a rank-one
flat). The graph is connected, since two distinct points determine
a line, and every line contains points. In rank two it is a star.
In rank at least three, three independent points and their distinct
joining lines give a six-cycle. So the displayed coefficient
vanishes without truncation, and is positive for any positive
truncation depth. This concerns the formal flag coefficient, not
linear independence of all flag statistics on actual matroids.

**Proof.** By palindromicity and the defining weighted sum in
[Proudfoot–Xu–Young](https://arxiv.org/html/1706.05575), read the
coefficients of degrees \(r-1,r-2,r-3\). The corank-three and
corank-four linear KL coefficient is the difference between
hyperplanes and points, and for corank five the quadratic KL
coefficient follows by comparing degrees two and three in its
palindromic Z. This gives
\[
\begin{aligned}
z_1&=W_{r-1},\\
z_2&=W_{r-2}+f_{r-3,r-1}-f_{r-3,r-2},\\
z_3&=W_{r-3}+f_{r-4,r-1}-f_{r-4,r-3}
 +f_{r-5,r-2}-f_{r-5,r-3}\\
&\quad+f_{r-5,r-3,r-1}-f_{r-5,r-3,r-2}
 -f_{r-5,r-4,r-1}+f_{r-5,r-4,r-3}.
\end{aligned}
\]
All proper flats and flags here survive truncation unchanged. Only
the triple flags can read \(f_{123}(A)\); subtract \((1,2,3)\) from
their rank strings to get
\((t,t+1,t+2),(t,t+1,t+1),(t,t,t+2),(t,t,t)\), respectively.
Drop repeated ranks, then group through the first exterior flat.
The four contributions are edges, minus points, minus lines, plus
one. This is the inclusion–exclusion analogue of keeping boundary
information intact in the repo's [polygon and closing-tile work](https://github.com/eliottcassidy2000/math/commit/71cb222faec6dbd6a367b71d1832459fbcfc28e2):
the local cancellation leaves the cycle obstruction, not a raw
count of contacts.

For a free exterior each quotient in this sum is Boolean of rank
\(d+2\). Its incidence graph subdivides \(K_{d+2}\), with cycle-space
dimension \(\binom{d+1}{2}\). There are \(\binom{s}{d+2}\) exterior
flats of the required rank. Now \(W_3(A_q)\sim q^8/2\) and
\(f_{123}(A_q)\sim q^9\); all other affine flags are \(O(q^8)\)
by Corollary 4.11. The coefficient of \(W_3(A_q)\) in \(z_1\)
is \(\binom{s}{d}\) and every other term in \(z_1\) has lower
\(q\)-degree. Thus \(z_1\sim\tfrac12\binom{s}{d}q^8\),
\(z_2=O_{s,d}(q^8)\), and the sole \(q^9\) term in \(z_3\) gives
the fixed-rank asymptotic.

## Growing rank: a cutoff for each truncation depth

Hold \(d\ge1\) fixed and let \(s,q\to\infty\) with \(s/\sqrt q\to u<\infty\).
Set \(A_d=d(d+3)/2,\ B_d=d(d+1)(d+4)(d+5)/12,\ H_d=\binom{d+1}{2}\).
In this regime
\[
\frac{z_2^2-z_1z_3}{q^{17}s^{2d+2}}
\ \longrightarrow\
\left(\frac{A_d^2}{4((d+2)!)^2}
-\frac{B_d}{4d!(d+4)!}\right)u^2
-\frac{H_d}{2d!(d+2)!}.
\]
The positive cutoff is
\[
u_d^2=\frac{12(d+1)^2(d+2)(d+3)}
 {2d^4+18d^3+56d^2+54d-10}.
\]
The denominator is positive for positive integer \(d\);
\(u_d^2\to6\) as \(d\to\infty\), after taking the fixed-depth limits.
In particular \(u_2^2=360/83\), different from the
[earlier \(4/3\)](../z123_rank_six_uniform_spike_20261127/README.md)
for padding after truncation. Below the cutoff this index fails
log-concavity, above it this index passes. The equality case and
the other indices are not classified.

**Leading coefficients.** For total rank string \(r+e_1,\ldots,r+e_k\)
split off a nondecreasing \(A\)-rank string \(a_i\in[0,4]\).
The exterior deficits \(h_i=d+a_i-4-e_i\) are nonnegative and
nonincreasing. Its flag count is
\(\binom{s}{h_1}\prod_{i=1}^{k-1}\binom{h_i}{h_{i+1}}\).
Thus in the \(\binom{s}{j}\) basis, the coefficient polynomials of
\(W_3(A_q)\) in \(z_1,z_2,z_3\) have top terms
\(\binom{s}{d}, A_d\binom{s}{d+2}, B_d\binom{s}{d+4}\).
For \(z_2\) the top term comes from the two pair flags whose first
coordinate rank is 3. At the next coordinate, \(A\)-rank 3 or 4 gives
Pascal sums and the difference is \(\binom{d+2}{2}-1=A_d\).

For \(z_3\), the top term takes first \(A\)-rank 3 and first deficit
\(N=d+4\). The two pair flags beginning at \(r-5\) give
\(\binom N3-N\). The four triple flags, with masks \(333,334,344\),
give the signed sum
\[
 (\binom N2+N)\binom{N-1}{2}
 -(N+1)\binom N2-(N+1)\binom N3+N(N+1).
\]
Together these equal \(B_d\). Including the pair flags is necessary.
The coefficient of the empty \(A\)-flag in \(z_i\) has \(s\)-degree at
most \(d+2i-1\); each other nonempty flag except \(123\) and \(3\)
begins at \(A\)-rank at most 2, and has \(s\)-degree at most
\(d+2i-3\) (absent terms contribute zero). These bounds follow from
the earliest total rank \(r-(2i-1)\). The empty flag itself equals
one. Consequently the remaining terms are negligible in the
displayed normalization, using \(s=O(\sqrt q)\), and the asymptotic
follows by multiplying the leading terms. No uniform assertion for
depth \(d\) growing with \(q\) is made.

Run `node replay.js` for exact arithmetic checks of the top
coefficient identities at depths 1 through 20 and numerical sign
controls. The replay is supplementary: the identities above prove
the fixed-depth results. It also searches the simpler \(r=6,d=2\)
family for its first odd prime-power negative tests, without claiming
that this family beats the earlier Fano specimens in ground size.
