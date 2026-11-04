# A reusable nested grammar must acquire unbounded depth or new cores

2026-11-03. **PROVED counting obstruction.** This extends the incoming
finite-fixed-suffix-bank argument in
[`session180_170807_returns_20261104/ports/PROOF.md`](../session180_170807_returns_20261104/ports/PROOF.md)
(pinned incoming commit `838b01cec4`). It does not rule out unbounded-depth
recursion, changing the bottom-core family, or changing the attachment type.

Fix K finite bottom tournaments, each strong and of order at most M. Starting
from any of them, repeatedly apply

    T_i = Z(TT_(s_i); T_(i-1), w_i),   s_i>=1,

where **any actual vertex** w_i of the preceding tournament may be chosen.
All created vertices are new; the only quotient used here forgets vertex
names when counting possible numerical Hamiltonian counts. The attachment
is precisely the zero-ready operation: TT_(s_i) dominates T_(i-1), and a
new root is beaten only by w_i. No restriction to one predetermined mark,
one finite field phase, or one fixed sequence of marks is imposed.

For at most d layers, the number N_d(X) of distinct obtainable Hamiltonian
counts at most X obeys, with L=floor(log_2 X),

    N_d(X) <= K sum_(j=0)^min(d,L) binom(L,j) (M+L+d)^j.       (D1)

In particular, for fixed d this is O((log X)^(2d)).

Proof. Split a Hamiltonian path of T immediately after w. This injects
Hamiltonian paths into the ordered two-slot covers counted by E_w(T,2),
including the case of an empty second slot. Thus E_w(T,2)>=H(T). The exact
zero-ready identity gives

    H(T_i) = H(T_(i-1)) + 2^(s_i) E_(w_i)(T_(i-1),2)
           >= (1+2^(s_i)) H(T_(i-1))
           >= 2^(s_i) H(T_(i-1)).

Consequently sum_i s_i<=L whenever the final count is at most X. For exactly
j layers there are binom(L,j) ordered positive length tuples with that sum
bound (the standard partial-sum subset bijection). Every intermediate graph
has at most M+L+d vertices, so there are at most (M+L+d)^j possible sequences
of actual mark selections. Multiplying by K and summing proves (D1). Any
isomorphisms, duplicate counts, or extra finite phase labels can only reduce
the represented set of numerical counts. The argument bounds actual graph
construction histories before taking any quotient; it does not identify
separate joint-port states.

Now suppose this grammar represented all words of length k in both maps
H->3H+2 and H->3H+4 starting from H_0. Those words have 2^k distinct values,
all at most X_k=3^k(H_0+2)-2. Subtracting the all-2 word leaves distinct
base-three words with digits 0 and 2. Therefore

    2^k <= N_d(X_k) = O(k^(2d))

would hold for every k, a contradiction for fixed d.

Even if d=d_k varies with k, this argument gives the necessary lower bound

    d_k = Omega(k/log(k+2))                              (D2)

for a fixed finite bottom bank and unrestricted mark selection. To see this,
every useful layer has s_i>=1, so layers beyond L cannot contribute. Replace
d_k by min(d_k,L). The right side of (D1) is at most

    K(d_k+1) ((L+1)(M+L+d_k))^(d_k),

and L=O(k), d_k<=L. Taking logarithms against 2^k proves (D2).
This is a necessary complexity bound, not a construction or a sharp bound.

## Implication for the present search

The actual 170807 source has a five-layer nested zero-ready decomposition
above its five-vertex core. The new finite returns show that its suffix and
port can change while all original names survive. They do not establish
closure of a fixed-depth grammar. Together, the exact one-exit formula in
[`structure/README.md`](structure/README.md), its arbitrary-core obstruction,
and (D1) direct the search toward an unbounded family of changing multi-port
contexts. Both the port-incidence information and the depth must remain
available; a fixed menu of scalar clocks and finite phase labels cannot
supply fully branching closure by itself.

This is a restriction on one proposed proof architecture. It is not an
obstruction to the conjecture that 7 and 21 are the only missing positive
odd Hamiltonian path counts.
