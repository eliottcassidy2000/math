# Outlet block placement: exact pair, prime inverse and an all-height palette obstruction

The complement-star family gives rooted sources (W_e), indexed by integers
(e\ge1), with root (u) after a forward spine of (m=2+floor(log₂ e))
vertices and with nonempty marked positions (S) that (u) beats. The last
spine vertex beats (u), and position (i) is in (S) iff bit (m-2-i) of
(e) is set. The first spine position is marked. If (J(e)) sums
(\omega_{j-i}) over set-bit pairs (i<j), where (\omega_1=2) and
(\omega_d=2^{d-1}) for (d\ge2), the inherited chart gives
```math
H(W_e)=1+2e-J(e),\qquad E_u(W_e)=e.
```
This is the complement-star inverse in the [preceding framed session](https://github.com/eliottcassidy2000/math/blob/bd79677a2da8e271fbe420cdd032d934e764e2ad/05-knowledge/results/chanpak_framed_inverse_20261115/COMPLEMENT_STAR_INVERSE.md). Throughout this note (P(B)) means **two slots with empties**,
the symbol (Q(B)) in the inherited arbitrary single-outlet replacement
formula in `chanpak_two_path_session/ONE_OUTLET.md`, §1. In particular,
```math
h=H(B),\qquad P=P(B)=\sum_{I\subseteq B} H(B[I])H(B\setminus I).
```
The empty tournament has path count 1. Both formulas and obstructions below
are for physically specified tournaments, not for arbitrary integer tuples.

## 1. Independent outlet block and exact repeated return

Insert any nonempty (B) after the entire spine and before (u), with
(W_e-u\to B\to u). Then
```math
N_B(e)=h+2Pe-hJ(e).                                     \tag{1}
```
This specializes the inherited general replacement theorem for a vertex
with one out-neighbor. Here is a direct path proof. There are at most two
visits to (B), since any departure must go immediately to (u) and can
occur at most once. The one-visit contribution is (h[H(W_e)+2e]):
contract the sole visit to the outlet vertex. With two visits, the
intervening segment starts at (u) and contains at least one spine
vertex. If its first vertex is (i\in S), the suffix after (i) can be
arbitrarily allocated to that increasing segment or the increasing
prefix before the first (B)-visit. All positions earlier than (i)
belong to that prefix. This gives exactly (2^{m-1-i}) choices, summed
over (i) to give (2e). The choices internal to the two (B)-visits
contribute (P-2h). Summing proves (1).

The resulting tournament is strong: the spine opens to (B) which
opens to (u), and (u) returns to the first spine vertex. At the same
root the new readout is (E'_u=Pe). Indeed a second appended singleton
still feeds (u) and can be amalgamated as an ordinal final vertex of
(B), leaving (h) unchanged and doubling (P). The change in (1) is
exactly (2Pe), so the one-outlet identity gives the specified (E'_u).
In particular, the pair (H-2E,E) has diagonal transfer
```math
(H'-2E',E')=(h(H-2E), PE).
```
Repeated block returns are thus exact, but this particular retained root
cannot realize a nontrivial translation in those two coordinates. For
arbitrary nonempty (B_1,B_2), both (H) and (P) multiply on the
ordinal sum (B_1\to B_2): each of the two ordered paths concatenates
its choices on the two parts, with either choice allowed to be empty.

Pure transitive padding aliases the original scalar spectrum:
the previous (N_t(e)=1+2^{t+2}e-J(e)) is exactly (N_0(2^t e)),
since shifting all set bits preserves (J). More generally any
singleton SCC of (B) can be removed and its effect absorbed by
doubling (e). This identity concerns scalar output, not the graph
shape, path length, or labelled rooted remainder.

## 2. A strong four-vertex outlet supplies a new physical pair

Use the fixed labelled outlet (B=(8,1,3,6)) in outgoing bit masks,
with ((h,P)=(5,32)). Addresses (x+2^k y), (k\ge8), with
```math
(x,y)=(106,17773)\quad\text{or}\quad(103,18613)
```
give exactly
```math
141775221\,2^{k-7}-56761,\qquad
141775221\,2^{k-7}-56763.
```
Indeed the intervening zero bit separates the low seven bits from the
high word; the cross term of (J) is (2^{k-7}y\operatorname{rev}_7(x)).
The two reversed words are 43 and 115, and the two sums (J(x)+J(y))
are 12710 and 12672. Both slopes are
(y[64\cdot128-5\operatorname{rev}_7(x)]=141775221). The intercepts
differ by two. This is an adjacent pair at every height, not a long
consecutive window.

Both have address popcount 12 and length (k+15); the same unlabelled
simple graph represents both in the framed chart. Follow the increasing
spine, then a fixed Hamilton path through (B), then (u). The backward
tile graph consists of the star joining the last gap to the twelve
marked gaps, the fixed graph from the internal arcs of the (B) frame,
and isolated gaps. The two pieces are disjoint: star leaves precede
the block, and its center is the last gap, outside (B)'s internal
gap set. The complements also have positive undirected Hamilton
counts. There are at most three internal block gaps and more than
four isolated gaps; the star center and the at most three internal gaps form an exceptional
set of at most four vertices. Five universal vertices suffice to
separate them in a Hamilton path. The rest of the graph is a clique,
so its vertices can be appended after the last universal vertex.

## 3. Target-first pruning in the prime and 7p lanes

For prescribed (N) and a specified physical profile ((h,P)), every
possible address satisfies both
```math
\frac{N-h}{2P}\le e\le \frac{N-h}{2P-h},\qquad
  2P e\equiv N\pmod h.                              \tag{2}
```
The bounds use (0\le J(e)\le e), since the sum of all lower
weights below a fixed upper bit (j>0) is (2^j). Put (d=\gcd(h,P)).
If (d\nmid N), reject; otherwise the congruence gives one progression
of spacing (h/d). Since (P\ge2h), there are at most
(1+Nd/[2P(2P-h)]) candidate addresses. An actual witness is obtained
only by testing (1), not by accepting (2).

Every nontrivial strong SCC has order at least 3 and (P\ge4h):
splitting each Hamilton path in all (n+1) places injects into
the two-slot count. Thus for any (B) with a nontrivial SCC,
(N_B(e)\ge 2P\ge8h). A prime target forces (d=1).
For (N=7p) with odd prime (p\ne7), it forces (d\mid7):
if (p\mid d), then (N_B(e)\ge8p>7p). For (N=49), the same
approach gives (d\mid7) (if (d=49), it is too large).
This is finite inverse pruning, not a supply proof for a physical
profile. For example, once (P^2>Nd/3), the interval/progression
contains at most one address by (2P(2P-h)\ge3P^2).

## 4. Unbounded order with a fixed finite SCC palette

For any given ((h,P)), its presentations through (X), counting
addresses rather than distinct values, have asymptotic mass
```math
\lim_{X\to\infty}\frac{\#\{e\ge1:N_B(e)\le X\}}{X}
  =\frac1h\log\frac{2P}{2P-h}.                        \tag{3}
```
The binary row Riemann-sum method is inherited. For completeness,
split (e=x+q a), (q) dyadic and (q\asymp\sqrt X). The interaction
is (J(e)=J(x)+J(a)+qa\operatorname{rev}_{\log_2q}(x)/q+epsilon),
where (epsilon\in\{0,1\}) corrects the distance-one cross term.
Thus for each residue (x), the count in (a) is
(X/[q(2P-h\operatorname{rev}(x)/q)]+O_{h,P}(1)),
because (J(x),J(a)=O_{h,P}(\sqrt X)) for all addresses that can
matter. The bit reversals permute the uniform grid, establishing
(3) with error (O_{h,P}(\sqrt X)) before division by (X).

Allow every outlet block whose SCCs have order at most seven, with
no bound on the number of SCCs. Remove singleton SCCs as above.
Each remaining distinct strong profile (i) has coordinates
((h_i,P_i)), and with its multiplicity (\alpha_i\ge0) the total
profile is ((\prod_i h_i^{\alpha_i},\prod_i P_i^{\alpha_i})).
This can overcount profiles or outputs, which is harmless for an
upper bound. The all-zero profile is represented by a singleton
and has mass (\log(4/3)). For the nonzero profiles the bound
(1/(2P-h)\le 4/(7P)) dominates the normalized presentation counts.
Its sum converges over the finite bank's unbounded multiplicities
by the geometric-product identity. Thus (3) may be summed, and the
upper bound for all distinct outputs divided by (X) is
```math
\log(4/3)+\sum_{j\ge1}
\frac{G_j-1}{2^j j},\qquad
G_j=\prod_i \left(1-\frac{h_i^{j-1}}{P_i^j}\right)^{-1}.   \tag{4}
```
This follows by expanding the logarithm as a positive series and
summing the nonzero multiplicity vectors term by term. It also
holds with any other finite strong palette and is an obstruction
whenever the right side is below one half.

The exact profiles through order seven are recorded in
`profiles.txt`, with counts 1, 1, 5, 16, 131 by order 3–7.
A labelled orientation census is sufficient because profiles
ignore labels; the replay enumerates masks, filters by forward
and reverse reachability, and computes the full induced-subset
path DP before evaluating ((h,P)). No order-eight or larger
**strong core** is claimed included. The rational calculation
in `check.py` bounds (4) strictly by .493. In that calculation
(G_j-1\le 4^{-(j-1)}(G_1-1)) follows from (h_i/P_i\le1/4) by
expanding over nonzero multiplicity vectors; the tail from
(j=4) is at most ((G_1-1)/3584). The all-zero contribution
uses the positive atanh series at (1/7) with an explicit geometric
tail. Consequently even this unbounded changing-block architecture
misses more than 1.4 percent of odd counts in lower density.
The numerics concern a small palette of local profiles **only**;
the exclusion is for infinitely many outputs and unbounded order.
The baseline singleton already has all fixed-compatible odd
residues by the inherited prefix programming theorem, so this
failure cannot be detected by fixing and checking a modulus.

## Reproduce

From this directory:
```sh
c++ -std=c++17 -O3 profiles.cpp -o /tmp/chanpak_outlet_profiles
/tmp/chanpak_outlet_profiles | diff -u profiles.txt -
python check.py
```
The Python verifier independently builds and counts 116 literal small
tournaments, then checks the pair formula and the rational mass cap.
The strong-palette proof does not constrain arbitrary B, changing the
marked host, or recontacting a new root. Those remain possible routes
to exact placement; Chan–Pak is not proved here.
