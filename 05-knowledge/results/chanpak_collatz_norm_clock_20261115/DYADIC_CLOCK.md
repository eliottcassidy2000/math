# An all-order physical dyadic clock, with an absolute error budget

**Proved.** For all sufficiently large odd N, there is a k ≥ c sqrt(log₂ N)
and a three-count tournament window on each side of N, each at distance
O(N^(3/8)). In each window H₀,H₁,H₂ differ by exactly +2,+2. The same
index j on both sides satisfies Hⱼ ≡ N (mod 3·2ᵏ). For fixed k the radius
is Oₖ(N^(1/3)). No assertion of a zero residual is made.

## Physical counter and its nested switches

Use the literal six-vertex center bank of the [target-address session](
../chanpak_target_address_20261115/README.md): for any positive addresses
x,y all three counts come from actual strong tournaments, and

    H_j = 432xy + 129x + 258y + 15[x odd] + 30[y odd]
          + 45(p(x)+p(y)) + 135 + 2j.

Here p(2u)=p(u)+u+[u odd], p(2u+1)=p(u), p(0)=0, hence p is even.
Use pre-addresses U,V and append 1 to both at the end: x=2U+1, y=2V+1.
Before that final append, always append complementary bits. Define
S=U+V, r=U mod 2, W=p(U)+p(V), K=W−U+r. After complementary parity
has been established, the exact two-letter law is

    append t to U, 1−t to V:
    S' = 2S+1,  U' = 2U+t,
    K' = K + t(S+1−2(U+r)).

Even before complementary parity, t=0 preserves K, by p's recurrence.
Set M=2^k, k≥2. When S≡−1 (mod M), this remains true and complementary
appends preserve K mod 4. Start every block with U≡0 (mod M), using k
zero letters at the start. For h=2,...,k−1 use, at will,

    1 0^(h−1) 1 0^k      or      0^(h+1+k).

The first 1 contributes zero mod M. Just before the second 1,
U≡2^(h−1) and r=0, so this 1 adds −2^h to K mod M; both choices end
with U≡0 mod M. The clocks do not change the actual final sum: every
letter is paired with its complement, and the two options have equal
length. Thus the bits of any desired K correction in 4Z/MZ are independent.
The control length B=k+Σ(h+1+k)+2 is O(k²).

Crucially, the last both-1 append escapes the frozen mod-4 sector.
Use one terminal two-bit word 00 or 10 on U, with complement on V.
Either keeps K mod M and ends in a prescribed u=0 or 2 mod M (for k=2
the same statement holds). Then U even, both final addresses odd, and

    H_j ≡ 51 + 2j + 45K − 213u − 432(2u+1)²  (mod M).

This follows by x≡2u+1, y≡−2u−1 and W≡K+u. In particular
H_j≡K+3+3u+2j (mod 4). Pick j=2N mod 3, then choose u to match N mod 4,
then set K using the independent clock bits because 45 is odd. This
hits N mod 3M in the designated center, while retaining both +2 switches.

## Putting the clock at the crossing, not after it

Let L=√N, choose g a power of two comparable to L^(1/3), and choose
the least integer pre-sum S≥(ceil(√(N/108))−2)/2 with S≡−1 mod g
(and 2S+2≥ceil(√(N/108))). Then R=2S+2 satisfies
0≤R−√(N/108)=O(g). Select k with q=2^B≤N^(1/24) and q|g.
This permits k≥c sqrt(log₂ N), for absolute c>0 and large N. When two
prefixes are held fixed on a length-q grid, the final guard replaces
only the last B bits of complementary g-blocks. The upper two
prefixes at the beginning of the guard have sum ≡−1 mod M for
large N: (S+1)/q is divisible by M, since qM|g. Thus the literal clock
works in every grid cell. One representative per cell suffices.

For U=S/2+d with d≥0, the polynomial part after both-1 is

    108 R² + (387/2)R − 1728 d² − 258d.

The omitted flags and constants are nonnegative. A full first
cell above d=0 lies above N, since d≤2q and q²=o(L).
A full cell at d=O(√(Lg)) sufficiently far out lies below N by
0≤p(U)+p(V)≤S, even for j=2. In the crossing band, previous
[complementary-block calculation](../chanpak_qrh_residual_20261115/COMPLEMENT.md)
gives W=C(S,g)+O(d+g+S/g+1) for one common C(S,g). Indeed the two
low blocks are complementary; their reversed coefficients sum to g−1;
subtract that common linear term, and the remaining term is bounded by
the imbalance plus the two low and high checksum sizes.
Thus neighboring representatives, with address gap at most 2q, have
output gap O(q·N^(1/3)+N^(1/3)); the actual sum and the chosen
residue are shared. Bisect by sign on each end of the three-count
window (no monotonicity needed), obtaining an entire window on each
side. This gives O(N^(3/8)), and O_k(N^(1/3)) for fixed k. The normalized
residual is a multiple of one and of size O(N^(3/8)/(3·2^k)), but this is
not an exact physical method for setting that quotient to zero.

**Capacity warning.** In a fixed-sum cell, programming more low bits
costs absolute address resolution. Ignoring the q multiplier in the
quadratic crossing would incorrectly claim a uniform N^(1/3) radius
for growing k. The spectral conjecture remains open.

Replay: python dyadic_clock.py.
