# Slot spectra: arithmetic carry, positive reconstruction, and a physical third-slot switch

Write Q_j(T) for ordered j-slot path covers of T, allowing empty paths. Thus H=Q_1 and the previously used Q is Q_2. This note applies the [Rédei–Berge odd-cycle expansion](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20) and the repo's [reflection-positivity margin study](../pythagorean_form_conic_graphs_rp_margin_20261018.md) to the actual absolute-target profile problem. The spectral statements below concern finite tournament polynomials, not Yang–Mills mass.

## 1. Positive spectral data and what finitely many slot observations remember

A packing of k disjoint odd cycles on s vertices contributes 2^k to H. In the j-slot sum every whole cycle and each unused vertex independently chooses a slot, yielding 2^k j^(n−s+k). Thus
```
Q_j = sum_(d=0)^D c_d j^(n−2d),
D=floor((n−1)/2),   r=n−2D∈{1,2},
c_0=1 and c_d are nonnegative even integers for d>0.
```
Here d=(s−k)/2: cyclic constituents are contracted before the palette is chosen. For a strong tournament all c_d>0. In fact c_d≥2(n−2d) for d≥1: take a nested chain of strong induced subtournaments of all orders from 3 to n using [Moon's vertex-pancyclicity theorem](https://doi.org/10.4153/CMB-1966-038-7), and at each order ≥2d+1 take a (2d+1)-cycle through its new vertex. Each of these n−2d cycles is new.

Set R(j)=Q_j/j^n. For arbitrary distinct positive slot parameters u_i, the matrix
```
K(i,l)=R(u_i u_l)=sum_d c_d u_i^(−2d) u_l^(−2d)
```
is a Gram matrix. This is a literal positive reflection kernel on the multiplicative scale variable. For a strong tournament its rank is min(number of scales,D+1) by Vandermonde; without strength the rank is at most the number of nonzero coefficients. After subtracting the empty packing, R(e^t)−1 is a sum of decaying modes of exponents 2d, and its leading decay exponent is exactly 2 for every nontransitive tournament (triangles occur). The rank and decay are in this explicitly diagonal finite model. They are neither a mixing estimate for changing the tournament nor a continuum mass-gap inference.

Let G(z)=sum_(d=1)^D (c_d/2) z^(D−d). If two order-n tournaments agree at Q_1,…,Q_k, their G-polynomials differ by an integer multiple of the monic product ∏_(j=1)^k(z−j²) in Z[z]. If k≥D they agree at every slot count. Otherwise, at the next slot one has the exact carry spacing
```
ΔQ_(k+1) ∈ 2 (k+1)^(r−1) (2k+1)! · Z.
```
This follows by evaluating the product at (k+1)². In particular a third-slot shift preserving order and both H,Q_2 is a multiple of 240 in odd order and of 720 in even order. The equal-H,Q_2 profile is not generally reconstruction-complete. This is a minimum spacing when a nonzero shift exists, not a theorem that every multiple is physical.

## 2. An exact carrying filter for the bounded profile progression

For a proposed order n≥3, h and q=Q_2, first subtract the strong baseline β_d=n−2d if strength is required, or zero otherwise. Put
```
A=(h−1)/2 − sum β_d,
S=(q−2^n)/2^(r+1) − sum β_d·4^(D−d).
```
Then a necessary condition for a physical tournament is the following **exact** coefficient-relaxation test:
1. A,S are nonnegative integers.
2. A≤S, A≡S mod 3.
3. A is at least the digit sum obtained by writing S in base 4 while allowing the top digit (at 4^(D−1)) to overflow; explicitly, it is floor(S/4^(D−1)) plus the base-4 digit sum of the remainder.

These conditions are sufficient for a *nonnegative integer coefficient vector above the baseline*, not for a tournament. To see it, identify S with the total value of pieces of weights 1,4,…,4^(D−1) and A with their number. Merging four equal pieces whenever a larger type exists minimizes cardinality and gives exactly the bounded-height base-4 expansion. Repeatedly splitting one piece into four of the next smaller weight increases the cardinality by three until all S pieces have weight 1. Thus no admissible sizes are skipped in the coefficient relaxation. This gives a precise absolute rejection test in N=L(Q_2−h)+Mh: evaluate each candidate (h,Q_2) at a proposed order before spending effort on realizing the packing counts.

Example: proposed (n,h,q)=(5,11,112) passes Q_2≥(n+1)h, the elementary mod-3 guard Q_2≡(−1)^n h, the 4-divisibility guard, and 2^n≤Q_2≤2^n(h+3)/4. After the strong baseline, A=1 and S=7, with A≡S mod 3. But the bounded-height digits of 7 are 13 in base 4 and need four pieces, so this profile is not physical for a strong order-five core.

## 3. A literal switch invisible to every two-slot host

Two strong order-eight tournaments with outgoing rows
```
(92,141,32,228,206,19,38,101)
(52,229,208,231,10,148,49,81)
```
both have H=345 and Q_2=5640. Their Q_3 are 45225 and 44505: the difference is exactly 720, attaining the even-order minimal nonzero spacing. The polynomial difference is carried by a multiple of (z−1)(z−4). In particular no host whose exterior sees the core only via (H,Q_2) can distinguish these. Such a host includes the [one-outlet two-core formula](../chanpak_quadratic_anchor_atlas_20261115/README.md).

To see the switch in an absolute H readout, put either core V into C3[V,TT2,TT2]. Classifying the cyclic run words gives
```
H(C3[V,TT2,TT2]) = 4 Q_3(V) + 5 Q_2(V) − 11 H(V).
```
Indeed, if L_j is the j-slot count with every slot nonempty, the cyclic run coefficients of L1,L2,L3 on V are 11,17,4, respectively. Each TT2 side has L1=1, L2=2 and higher L_j=0. Substituting L2=Q2−2H and L3=Q3−3Q2+3H gives the formula. The two literal H counts differ by 2880 with every two-slot datum of the core unchanged. Common singleton padding of V preserves equality of the first two observables between the two cores and multiplies the three-slot difference by 3 per added singleton. It does not shrink the absolute shift.

The principle is related to the chromatic port work's [boundary relation rather than palette size](../six_colour_interface_profile_20261015f.md): on changing the admissible consumer (here, the number of visits), matching a lower readout no longer guarantees equivalent behavior. The gain is a true additional actuator; the factorial tariff quantifies why it cannot immediately be used as a +2 finish.

Replay: `python replay_slot_spectral.py` checks both physical profiles and both literal cyclic-wrap counts by subset-endpoint DP, and validates the bounded-digit criterion on finite coefficient cones.
