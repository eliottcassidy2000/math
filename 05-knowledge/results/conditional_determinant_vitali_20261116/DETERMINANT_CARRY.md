# Paired arm carries choose exponentially many nonzero determinants at one prescribed gap

Use the [exact domino contraction formula](../unit_shields_wheel_followup_20261116/apex/SHIELD.md).
A spine letter 2a+b gives the win bits of the two apices A,B.
Write L=(L_a,L_b) for the forward binary addresses, R=(R_a,R_b)
for the reversed-complement addresses, and D=L_a R_b−L_b R_a.
Let Δ=H(A→B)−H(B→A).

## Grow both arms, carry one integer

Take a core c=0123, with common-letter arms U,V of the same length k
(all their letters are 0 or 3). Write X for the ordinary binary address
of U (3 becomes 1) and Y for the **reversed-complement** address of V.
Put r_k(x)=2^k−1−rev_k(x), so V has ordinary address r_k(Y).
The core's marked rows satisfy

    L(c)=(3,5), R(c)=(3,5), Δ(c)=0.

Within the common-letter arms Δ=0; the cross determinant from U to V
vanishes. The two other cross contributions in the domino formula
are −2X and +2Y; the contact adjacency terms vanish because the
first and last letters of c are common. Therefore

    Δ(U c V)=2(Y−X).

For a specified half-gap g≥0 select Y=X+g, with
2^(k−1)≤X≤Y<2^k. Both apices beat the first spine vertex and
lose to the last; hence both mutual-arc orientations are strong.
No other arc is edited. Swapping orientation negates the signed gap.

Let A0=2^(k+4)X+r_k(Y), R0=2^(k+4)Y+r_k(X).
The assembled marked vectors are

    L=(A0,A0)+2^k(3,5),  R=(R0,R0)+2^k(3,5).

In particular their determinant is exactly

    D = 2^(k+1)[−16·2^k g+r_k(X+g)−r_k(X)].       (1)

The gap only depends on the ordinary carry g; the determinant still
sees its individual bit carries via the reversed difference.

## Exponential nonzero capacity at an exact chosen gap

Let j≥1, k=2j+1 and g=(4^j−1)/3=1+4+...+4^(j−1).
For each independently chosen bit e_i in {0,1}, set

    X=2^(k−1)+∑_(i=0)^(j−1) e_i 4^i.

When adding g to X, each even position 2i has an adjacent zero
sentinel at 2i+1. If e_i=0, reversal records a new 1 at weight
2^(k−1−2i); if e_i=1, there is one local carry into the sentinel,
changing the reversed amount by −2^(k−2−2i). There is no interaction
between these j carries, and they never reach the highest bit. Thus

    r_k(X+g)−r_k(X)=−4g+6E,
    E=∑_(i=0)^(j−1) e_i 4^(j−1−i).

By (1), all 2^j bit packets give distinct strictly negative
determinants at the **same literal gap 2g and same width 2k+4**.
This is an exact encoding with an exact decoder: solve (1) for E,
read base-four digits, and accept iff they are 0 or 1 within width j.
The determinant values are a sparse digital set, not an interval.
The two endpoint counts are not controlled by this decoder.

## Fixed +2 also admits growing determinants, but carries less freedom

Set g=1. If the final run of ones in X has length τ, the same carry
calculation now gives

    D=−30·2^(2k)−3·2^(2k−τ),   0≤τ≤k−2

in the strong range X≥2^(k−1), X+1<2^k. Every τ in that range
occurs. A fixed +2 alone therefore gives only k−1 determinant
options per width in this particular arm chart. The exponential
menu above came from many independently movable carries, not
just more room around a single carry.

One can also preserve +2 through a compatible symmetric outer
wrapper. Take s of width k with address 2^(k−1)−1, and set
P=s 0123 rev(s). Its gap is +2, and
D(P)=−2^(k+1)(15·2^k+3) by the same τ formula at τ=k−1.
For any width-m common-letter u starting with 3, consider u P rc(u),
where rc reverses and complements. This satisfies the strong
[mirrored exterior program](../unit_shields_wheel_followup_20261116/apex/PROGRAM.md).
Its full L and R vectors are each the inner vector scaled by 2^m
plus **the same** multiple of (1,1). The cross terms in D cancel
since L_b(P)−L_a(P)=R_b(P)−R_a(P). Hence the outer determinant is
4^m D(P), independent of the other m−1 exterior bits. The resulting
positive magnitudes are injective in (k,m): their odd part
15·2^k+3 specifies k and the 2-adic valuation specifies m.
This is a separate sparse fixed-gap atlas, not exact placement.

The carry identity is analogous to the [67/69 exact contrasts](ARCHIVE_67_69.md)
in one important way: a literal path comparison is retained when
common bulk terms cancel. Its distinguishing new feature is that the
free choices leave not just two marked rays but explicitly many
nonzero signed determinants. A counting or Vitali argument cannot
turn that menu into an absolute landing law without controlling H.

## Exact carry automaton and the Fibonacci maximum

The displayed independent sentinels give a simple subfamily. In fact
the **entire** determinant menu of this arm chart admits a two-state
automaton, for any prescribed 0≤g<2^(k−1).

Number bits from right to left, and let c_i∈{0,1} be the carry into
position i when adding g to X, with c_0=c_k=0. The strong range makes
the last top bit of both X and Y equal one, so c_(k−1)=0 too. From
y_i=x_i+g_i+c_i−2c_(i+1), reversal gives

    rev_k(Y)−rev_k(X)
       = rev_k(g) − 3∑_(i=1)^(k−1) 2^(k−1−i)c_i.          (2)

This follows by collecting the two appearances of each intermediate
carry: its coefficients are 2^(k−1−i) and −2·2^(k−i), summing to
−3·2^(k−1−i). Thus (1) encodes the *whole carry word* injectively.
To decode a candidate D, extract r_k(Y)−r_k(X) from (1); add
rev_k(g), divide by 3, read the carry bits, and check their range
and the two endpoint zeros. At each low bit:
- If g_i=c_i, c_(i+1)=g_i is forced and x_i can be either bit.
- If g_i≠c_i, choose x_i=c_(i+1), which realizes either carry.

These conditions are necessary and sufficient. Write (u,v) for the
numbers of carry paths currently at zero and one. Starting at (1,0),
scan the k−1 low bits of g. The updates are exactly

    bit 0: (u,v)→(u+v,v),     bit 1: (u,v)→(u,u+v).

The final u is the number of distinct determinants. For k=2j+1
and g=(4^j−1)/3 the input bits from low to high are j copies of 10,
so the exact menu has F_(2j+1) entries, where F_0=0,F_1=1.
The earlier 2^j values are only the sentinel subfamily of this menu.

With m scanned bits, the sorted pair is bounded coordinatewise
by (F_(m+1),F_m); this follows inductively because an update replaces
one coordinate by the sum. Alternation maximizes the larger
coordinate, and the final update may always be chosen to grow the
zero coordinate. Thus F_k is the exact maximum over g in this chart.
There is a literal golden-ratio growth law here, not just a numerical
resemblance. It controls determinant options at chosen gaps, not
the size or signs of gaps between successive Hamiltonian counts.

One can fix both of the two highest arm bits to 1 instead of just
the highest. Choose k=2j+2, the same alternating g, and require
X,Y in [3·2^(k−2),2^k). The low 2j bits and their automaton are
unchanged. Now the first two spine vertices beat neither apex
(both apices beat both vertices), and the last two beat both apices.
Use distinct outer endpoints for two directed odd cycles through
A and B, respectively. Assign each other spine vertex to either
cycle or neither and impose even extra vertices on each. The
positive odd-cycle-packing formula gives a lower bound 3^(n−4)
for both counts, with n=2k+4; three-run assignment gives upper
bound 2·3^n. Thus at a fixed gap and comparable absolute size
scale H=Θ(3^n), this chart produces Θ(H^(log_3(phi)/2))
distinct determinant options, phi=(1+sqrt(5))/2. The theorem
does not assert those options land distinct counts, or cover
all targets at that scale.
