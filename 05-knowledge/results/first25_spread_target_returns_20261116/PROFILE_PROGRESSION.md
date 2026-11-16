# The exact new obligation: a spectrum of two-slot profiles

Use the target equation N=Cv+I_j u in TARGETS.md.
There is a complete decomposition of the possible right profiles,
not merely an integer lattice: for each fixed odd u, factor it
into an ordered word of strong counts u_i>=3. Choose an actual
strong tournament B_i with that count. Then and only then

    v = 2^t product_i P_{B_i}(2),   t>=0.

For necessity compress the SCCs of B in their unique ordinal order;
both Hamilton paths and two-colour assignments factor independently
along this order. Every singleton contributes (1,2). The converse
uses the actual ordinal concatenation. In particular the empty
nontrivial word describes only v a power of two, u=1. Since
P_B(2)>=2^|B|, exact membership of a proposed (u,v) always has
a finite bounded-order test; this says nothing about universal
coverage in the target progression.

## Corrections to the inverse sieve

The first test is N≡I_j u mod C. But v is even for nonempty B
(by pairing a subset with its complement), so actually
N≡I_j u mod 2C. If |B|>=2, then v is divisible by four.
Indeed the ordinary odd-cycle pack identity for the path-colour
polynomial gives a term 2^(|B|−|vertices in packs|+2s) for a
pack of s cycles at argument 2. The empty pack has exponent
|B|, and every nonempty pack has exponent at least two.
For even |B|>=4, every pack exponent is at least three, so
the necessary modulus there is 8C. The same polynomial satisfies
P_B(-1)=(-1)^|B| H(B), hence v≡(-1)^|B| u mod3. The candidate
u therefore must solve the stronger linear progression

    N ≡ (I_j + C(-1)^|B|) u (mod 3C),

along with its dyadic constraints. These filters rule out
unrealizable integers v before asking for a physical block.

## A density barrier for tiny composite sources

No bounded list of u's can close Chan–Pak using just the three
intercepts and variable height. More quantitatively, let every
chosen u be at most U, and bound the target by X. The number of
distinct outputs in this family is at most

    exp(κ (log U)^2 + O(log(U+2) + log log(X+3))),
    κ = 9 log(2) / (2 log(5)^2).

Here is the counting argument. Use the Busch lower bound
n_i−1≤3 log(u_i)/log5 on each nontrivial strong SCC.
Their total order is O(log(U+2)), since each u_i>=3; their
total internal arc slots are at most
9(log U)^2/(2 log(5)^2)+O(log(U+2)), using
sum(log u_i)^2≤(sum log u_i)^2. Enumerating the ordered
size lists costs exp(O(log(U+2))); orientations cost at most
two to the number of slots. The singleton dyadic factor 2^t
and the row height 2^s merge in Cv, so each resulting profile
mantissa contributes only O(log(X+3)) distinct outputs. For
example N≥(756·2^(s+t)−33/2)v_0 and v_0≥2u on any
nonempty nontrivial word, which bounds s+t. The empty word
also contributes O(log X). Fixed three choices of j change only
the constant.

Thus U≤exp(c sqrt(log X)) with c<1/sqrt(κ) yields density zero.
This is a constraint on the exact composite-block route, not on
the full spectrum. A target-dependent source must either allow
growing profiles beyond this regime or alter the reading maps.
It is the analogue of the rail example: a fixed little stock
of boundary types can be extended many times, yet still not
generate enough different absolute placements.
