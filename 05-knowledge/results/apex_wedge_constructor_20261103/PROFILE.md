# Full path-cover responses and weighted binary charges

**PROVED scoped.** This extends the new edge construction using the incoming
[positive profile operator](../session159_constructions_20261103/profile_operator/NOTE.md).
The source is still a literal Boolean two-apex tournament. The extra argument
r is the number of ordered path slots, not a relaxed edge weight.

## P1 The all-slot contracted edge kernel

Let F_r(T) count ordered r-tuples of disjoint directed paths covering T,
allowing empty paths; F_1=H. Put Delta_r=F_r(T+)-F_r(T-).
For r>=1 and m>=1 define K(r) by zero diagonal and, for h=j-i>=1,
\[
 K(r)_{ij}=r^{m-1-h}(r-1)^{h-1},\qquad
 K(r)_{ji}=r^{m-1-h}(r+1)^{h-1}.
\]
Here 0^0=1. Define
\[
 b_j(r)=\frac{r^j}{2}\big((r+1)^{m-1-j}+(r-1)^{m-1-j}\big),
\]
and the alternating matrix S_m(r) with upper entries
\[
 S_m(r)_{ij}=\frac{r^{m-1-h}}2
       \big((r+1)^{h-1}-(r-1)^{h-1}\big).
\]
All displayed half-sums and half-differences are integer polynomials. Then
\[
 \Delta_r=2r\{b(r)^T(u-w)+w^T S_m(r)u\}.                 \tag{P1}
\]
Set Delta_0=0; for m=0 all Delta_r vanish. At r=1 this is exactly G1.

**Proof.** A cover affected by the reversal has the block a,b in one of its
r slots. Its base vertices choose the gap before that block, the gap after
it, or one of r-1 other increasing paths. If both distinguished gaps are
nonempty, write i for the first gap's maximum and j for the second's minimum.
For i<j, vertices outside [i,j] have r choices and those strictly inside
have r-1 choices. For i>j the inside vertices have r+1 choices. This gives K(r).

If only the before gap is nonempty with maximum i, its coefficient is
L_i=r^i(r-1)^(m-1-i). If only the after gap is nonempty with minimum j,
its coefficient is R_j=(r-1)^j r^(m-1-j). If both are empty the contribution
is (r-1)^m. Multiply the sum by r for the named slot. Subtracting the
opposite block cancels the last constant and gives a linear term
(L+R+K^T 1)^T(u-w) and the skew bilinear term.
The geometric sums in each column give
\[
 L_j+R_j+\sum_iK(r)_{ij}
 =r^j\big((r+1)^{m-1-j}+(r-1)^{m-1-j}\big).
\]
This proves P1. In particular Delta_r is divisible by 2r at every integer
r>=1; this is a uniform statement about these edge reversals.

## P2 A weighted binary charge and its Lucas recurrence

Put
\[
 L_h(r)=((r+1)^h+(r-1)^h)/2,\qquad
 L_0=1,\quad L_1=r,\quad
 L_{h+2}=2rL_{h+1}-(r^2-1)L_h.
\]
For a positive integer N with binary digits n_0,...,n_(k-1), set
\[
 C_r(N)=\sum_{i=0}^{k-1}n_i r^i L_{k-1-i}(r).
\]
This is a specific Lucas-type recurrence, derived from the two roots r+1
and r-1. At r=1, C_1(N)=ceil(N/2).

For the balanced chart p01s,p10s of A2, let k=bitlength(P),t=bitlength(R).
Then the complete profile difference is
\[
 \Delta_r=2r\big(r^k C_r(R)-r^t C_r(P)\big).              \tag{P2}
\]
Thus the former binary charge cancellation has an exact weighted lift;
the two length clocks r^k,r^t are indispensable away from r=1.

**Proof.** Substituting u-w=e_k-e_(k+1) into P1 gives
\[
 \Delta_r/(2r)=
 \tfrac{r^k}{2}((r+1)^t-(r-1)^t)
 -\sum_i p_i r^{t+i}L_{k-1-i}
 -\sum_j s_j r^{k+t-1-j}L_j .
\]
The telescoping identity
sum_(j=0)^(t-1) r^(t-1-j)L_j=((r+1)^t-(r-1)^t)/2
cancels the full suffix sum. Since R is the reverse-complement of s,
what remains is P2.

For every j>=2, the offset-zero unit family P=2^j,R=2^j+2 satisfies
\[
 \Delta_r=2r^{2j+2}.                                    \tag{P3}
\]
Only the second-last bit changes, so C_r(R)-C_r(P)=r^j.
For every j>=3, the offset-three family P=2^j+3,R=2^j+5 satisfies
\[
 \Delta_r=2r^{2j}.                                      \tag{P4}
\]
The tail changes from 011 to101; its charge change is
r^(j-2)L_2-r^(j-1)L_1=r^(j-2).
These identities are all-size, all-slot formulas, not scalar interpolation.

## P3 Transport through actual triangle substitution

Let G(X) substitute X and two singletons into X->c->d->X. The incoming
profile theorem gives
\[
 F_r(G(X))=rF_{r+1}(X)+r^2F_r(X)-rF_{r-1}(X).
\]
Therefore the edge-response polynomial transforms by the same operator:
\[
 \Delta'_r=r\Delta_{r+1}+r^2\Delta_r-r\Delta_{r-1}.        \tag{P5}
\]
The old source vertices, including both apexes, remain literal vertices.
At r=1 the new Hamiltonian jump is Delta_2+Delta_1.

| balanced P,R | original counts | full Delta_r | jump after G |
|---|---|---|---:|
| 1,3 | 131,133 | 2r^3 | 18 |
| 4,6 | 1351,1353 | 2r^6 | 130 |
| 6,8 | 2625,2627 | 2r^5(2-r^2) | -126 |

Thus even the sign of an improving move can reverse under a specified
strong attachment. The full polynomial predicts this constructively.
For the two infinite families, P3/P4 give after-G jumps
2+2^(2j+3) and 2+2^(2j+1), respectively.

The first row has a useful contextual realization. Reversing one edge of
a directed triangle to make it transitive changes F_r by 2r. Adding two
transitively ordered singleton components multiplies that by r^2.
Hence the strong seven-vertex pair realizes the same profile difference
as that five-vertex, nonstrong padded pair. This equality persists in any
one-hole context built from ordinal sum with fixed tournaments and G:
ordinal sum multiplies F_r by the fixed profile, while G acts linearly.
It is equality of differences, not equality of the individual tournaments.

Common leading 00 padding, invisible to J at r=1, multiplies the full
Delta_r by r per padded vertex. Enlarging the observation language from
ordinary edge jumps to path-cover responses therefore restores the missing
length clock. This explains exactly why C3's normalized-word quotient needs
a length sidecar for these additional operations.

[profile.py](profile.py) compares P1 with ordered subset convolution on all
word pairs through length4 and slots0,...,4 (1,705 checks), checks the weighted
charge identity, and directly counts the six displayed small G constructions.
