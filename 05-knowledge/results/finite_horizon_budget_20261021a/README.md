# Uniform scale coverage at each fixed arrival horizon

Status: **PROVED** scoped counting/energy bound and computable scale tails;
**FINITE-EXACT** independent forward-prefix controls; uniform coverage at
every horizon for a fixed scale remains **OPEN**.

## Inheritance and target

Use the native primitive map, fixed 627-point grounded frontier V, charge
alphabet, stopped arrival convention and energy of
[session18_native_annulus](../session18_native_annulus/README.md). The closest
inverse mechanism is the one allowed strip residue modulo 3 per parent and
branch in [sibling_annulus_guards_20261020i](../sibling_annulus_guards_20261020i/README.md).
The hostile control is N=7, depth 8, where the native increment exceeds 153N;
this packet does not certify that pair. The corrected near miss is counting
eventual arrivals as timely arrivals. We examine only prefixes through the
requested cut; a path that has not reached V is never counted as an arrival.
The useful sidecar is the actual accumulated division clock, not just its
parity or its modular phase.

Let P be the primitive map, E(M,d) the native stopped energy at scale M and
arrival depth d, and 2^K <= N < 2^(K+1). Define

\[
B(K,d)=\sum_{r=0}^d 2^r {\lfloor(K+8+r)/6\rfloor+r\choose r}.
\]

Then, for every integer N>=1 and d>=0,

\[
E(64N,d)\le 3|V|B(K,d)^2=1881B(K,d)^2.                 \tag{1}
\]

In particular, if the right side is at most 153*2^K, then
E(64N,d)-E(N,d)<=153N throughout that dyadic band. The same holds at all
earlier cuts because B is increasing in d. No convergence conjecture or
negative interaction credit is used.

## Proof: clock cost makes the inverse tree thin at a fixed horizon

For an actual primitive edge s -> u of clock c, deleting B(x)=4x+1 strips
gives 2^c u <= 3s+1 <= 4s. For a first-arrival word of r primitive edges,
total clock T, source s<256N and terminal u>=1, multiplication gives

\[
2^T u\le 4^r s<2^{K+9+2r},\qquad T\le K+8+2r.        \tag{2}
\]

For a fixed primitive parent u and branch a in {1,2}, the inverse is

\[
s=\frac{2^a(3u+1)4^j-(2^a+3)}9.
\]

Integrality allows exactly one j residue modulo 3, so its clock is
c_0(u,a)+6k, k>=0, with c_0>=1. After fixing the r branch choices and the
r nonnegative k values, the backwards word and source are unique if valid.
Inequality (2) implies

\[
\sum k_i\le\left\lfloor\frac{K+8+r}{6}\right\rfloor.
\]

Weak compositions with sum bounded by m number (m+r choose r). There are
at most 2^r branch words. Invalid, negative, out-of-band, and prematurely
stopped words can all be included in this upper bound. Thus each terminal
has at most B(K,d) distinct first-arrived sources at cut d. The r=0 term
counts the terminal itself.

Every source charge has Q<=3: the lower alphabet has period three and sums
to zero over a period; at most one upper-band letter occurs. The finitely
many reduced lower/upper cases verify Q<=3. The sign from the actual clock
does not change Q. The triangle inequality for the norm sqrt(Q) bounds a
terminal fibre of size h by Q(sum charges)<=3h^2. Summing over V proves (1).

This proof preserves sources, first-arrival depth, primitive branches and
integer clocks. It deliberately discards cancellation, source band guards
and mod15/mod35 restrictions in an **upper** count; retaining those filters
could only reduce the count. The literal checks still record compatible
source registers. No quotient representative is substituted for a source.

## Six-band certificate for an infinite scale tail

Set h=floor((K+9)/6). Increasing K by six increases each binomial parameter
m_r by one. For 1<=r<=d its growth ratio is

\[
\frac{m_r+r+1}{m_r+1}\le1+\frac d{h+1}.
\]

If d<=7(h+1), then B(K+6,d)<=8B(K,d), hence the energy majorant grows by at
most 64 while the scale budget grows by exactly 64. The condition persists
at larger K. Therefore checking (1)<=153*2^k at k=K,...,K+5 proves coverage
for **every integer N>=2^K**, at all cuts through d.

The least K for this particular sufficient certificate is:

| d | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K | 10 | 16 | 22 | 28 | 34 | 41 | 53 | 65 | 79 | 104 | 129 |

These are thresholds for this coarse proof, not sharp thresholds for native
energy. In particular d<=1 is now uniform for the entire intended N>=4096
range. The quantifiers are **for each fixed d, all sufficiently large N**;
exchanging them would leave precisely the long-horizon obstruction open.

There is also a simple closed formula, without a lookup table:

\[
\boxed{N\ge2^{8d+12}\quad\Longrightarrow\quad
E(64N,h)-E(N,h)\le153N\quad(0\le h\le d).}             \tag{4}
\]

Indeed, with m=floor((K+8+d)/6), the binomial bound gives
B(K,d)<=(d+1)2^{m+2d}. The elementary inequality
(d+1)^2<=(9/4)2^d holds for all d>=0: check d=0,1,2, then its ratio
decreases. Consequently

\[
B(K,d)^2\le(9/4)2^{K/3+16d/3+8/3}
\le(9/4)2^{K-16/3}\quad(K\ge8d+12).
\]

The exact constant inequality `(1881*9/(153*4))^3<2^16` proves (4).
Thus the certified horizon can grow logarithmically with the scale. The
six-band certificates are sharper for the listed small depths. Neither
allows d to grow arbitrarily while N remains fixed.

## Residual accounting without recycling

Keep the existing exact decomposition Delta E=A+I, A+M=153N, and
I=R'+J', where J' includes each selected negative pair exactly once after
external credit is moved from R to its bank. Equation (1) gives

\[
R'\le1881B(K,d)^2-A-J'.                              \tag{3}
\]

On a certified band this is R'<=M-J', the exact compensated criterion for
the original budget. It does **not** prove the stronger R'<=M. No term in
J' is spent a second time, and changing the ownership partition cannot
alter Delta E. This direct bound can cover early cuts while source-owned
blocks handle other cuts; neither certificate licenses filling a gap by
assuming eventual arrival.

## Reproduction and failure boundary

```text
python -B 05-knowledge/results/finite_horizon_budget_20261021a/check.py
python -B -O 05-knowledge/results/finite_horizon_budget_20261021a/check.py
```

The immutable receipt pins the grounded frontier; checks all reduced charge
cases, inverse guards, exact six-band tails through d=20, and independent
literal forward prefixes at (N,D)=(7,8),(64,3),(4096,1). Prefix enumeration
stops at D even for an unarrived source. Checks remain active under -O.
Boolean/float aliases are hostile inputs. The receipt records 20,285 body
checks; verification adds one immutable receipt comparison.

The next precise gap is a source-specific bound beyond this scale/depth
wedge. The new linear credit supply has an unconditional sign, but its
timely arrival count is not bounded uniformly. Neither finite computation
nor the inverse-word count resolves all-horizon 153N or Collatz convergence.
