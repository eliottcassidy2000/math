# Seven-return testing and a sharp 81-layer reverse decoder

**PROVED scoped algebra / FINITE-EXACT / CONDITIONAL / OPEN.** The uniform
seven-return inequality and original-source reset repayment remain open.
The full horizon-74 exclusion has native small controls and an audited
horizon-55 positive control; it has not received a second full horizon-74
implementation audit. All five return branches remain enabled.

## Inherited target and the new finite bound

Use the five strictly growing cheap-ternary return words from
[INTEGER_RETURNS](../ternary_integer_returns_20261029/INTEGER_RETURNS.md).
Each native odd step is n→(3n+1)/2^a, with a the actual 2-adic valuation.
The auxiliary queue starts at one and changes by −1+3[a divisible by 3].
Its termination is already proved; repeated growing complete returns are
a separate question. Let E_k be the positive odd sources with k consecutive
growing complete returns and s_k=min E_k. The exact source-lift recurrence
and its monotone pruning are inherited from
[LEAST_SOURCES](../ternary_integer_returns_20261029/LEAST_SOURCES.md).

The candidate is

    s_(k+7)+1 >= 7(s_k+1) for every k>=0.                 (1)

The earlier [ATOMIC packet](../atomic_tail_reset_20261029/ATOMIC.md) has
exact s_60=1090945607561904127 and the certified lower bound
s_67>=7636619252933328895. The present exhaustive cylinder search visits
763878259 prefixes and excludes every positive source strictly below

    B=49(s_60+1)−1=53456334770533302271

from E_74. Thus **s_74>=B**. This tests a two-stage consequence of (1).
Since s_67 is not known exactly, it does **not** verify the actual ratio
(s_74+1)/(s_67+1)>=7. Neither lower bound is asserted to be an exact minimum.

For shifted source m=(n+1)/2 use the original-source atom
w(m)=1/(2·4^floor(log2 m)); these weights sum to one. The inherited occupied
cell counts have initial values 1,2,3,4,6,9,13,18,24,31,40 and recurrence
N_L=N_(L−1)+N_(L−6)+3N_(L−11). Since floor(log2((B+1)/2))=64, the previously
proved all-source bound gives

    nu(E_74) <= N_74/4^74 + (1/2) sum_(j=64)^73 N_(j+1)/4^j
             = 641060333941629/2^148 < 1.797e−30.          (2)

This includes all positive odd sources, including those above B. It is
a finite tail bound for this failure language, not a proof of Collatz
convergence or of A_k→0. Every positive odd source below B has a local
non-growing return within 74 attempts; that return need not cross below
the original source.

## U128 exclusion engine and its arithmetic range

`window74.cpp` is the inherited endpoint-lift algorithm extended to an
exclusive source bound at most 2^66 and a horizon at most 74. It preserves
the original representative r, endpoint h, total odd clock J, and division
clock D. Appending branch (j,d,c,r_v) uses

    t = ((r_v−h)/2) 3^(−J) mod 2^d,
    r' = r+2^(D+1)t,
    h' = (3^j(h+2·3^J t)+c)/2^d.

If a nonzero t would place r' at or above the incumbent, it is pruned before
forming the large shift. A possible nonzero lift has D<=65 and J<=D, so
the multiplication by 3^J is in range. At depth d0, every retained actual
endpoint obeys h+1 <= 2^66(3/2)^d0. The tentative numerator at the next
step is at most 3^7·2^66(3/2)^73+16177 < 2^122. Stored powers through 3^74
also fit U128. Inverse powers are needed only through J<=7·73=511; the
table has 520 entries. The search enumerates all admissible next cylinders,
with no two-letter filtering and no sampling of sources.

`engine_checks.json` records an independent scan of all odd n<150000
through 16 growing returns, both inclusive and exclusive threshold controls,
and a native replay of the independently audited least source at horizon55.
The latter uses branch4 and checks the endpoint and both integer clocks.
The producer's source hash is pinned in that receipt. These controls do not
replace an independent full search at horizon74.

## An exact reverse decoder on the stable divisibility layer

In m=(n+1)/2 coordinates the five return maps are

| i | map B_i(m) | guard on m | J,D |
|---|---|---|---|
| 0 | 3m/2 | 0 mod 2 | 1,1 |
| 1 | 3^4(m+1)/2^6 | 63 mod 64 | 4,6 |
| 2 | 3^6(3m+11)/2^11 | 679 mod 2048 | 7,11 |
| 3 | 3^5(9m+25)/2^11 | 1135 mod 2048 | 7,11 |
| 4 | 3^4(27m+59)/2^11 | 1439 mod 2048 | 7,11 |

All statements use these actual guards. The positive layer D=81 Z_(>0)
is invariant under each legal map. At a nonzero letter i the 3-adic
valuation resets to e_i=(4,6,5,4); a zero letter adds one to it.

For positive m put v=v_3(m) and define the valuation-adjusted unit

    chi(m)=2^v (m/3^v) mod 81.                            (3)

Then chi(B_0(m))=chi(m) whenever B_0 is legal. For x in D,

    chi(B_1(x))=61, chi(B_2(x))=13,
    chi(B_3(x))=70, chi(B_4(x))=34.                       (4)

Proof: write B_i(x)=3^e(a x+b)/2^d. The bracket is a 3-adic unit, and
x=0 mod81, so (3) equals 2^(e−d)b modulo81. The four residues are distinct.

Consequently an endpoint m of a suffix B_i followed by z legal zero
letters, beginning at a source x in D, determines i and z uniquely:

    i from chi(m),       z=v_3(m)−e_i.

Undo the zeros to y=2^z m/3^z and then recover

    x=(2^d y−3^e b)/(3^e a).                             (5)

Reject a candidate unless z>=0, x is a positive multiple of81, and the
original guard holds. Conversely every accepted candidate replays to m
by the specified native maps. The recovered clocks are exactly J_i+z
and D_i+z. Thus the infinite collection of nonzero-letter/zero-run images
of D is pairwise disjoint and has an explicit partial inverse.

This is a suffix decoder; it is not a finite-modulus predictor of forward
behaviour. It uses the unbounded integer v_3(m), in addition to a unit
residue. Every indefinitely growing positive path would eventually enter
D: an infinite all-zero path would require a positive m divisible by
arbitrarily high powers of2, and after the first nonzero letter its output
is divisible by81. The decoder applies to subsequent nonzero letters.

The modulus and layer have a concrete sharp boundary. At precisions
3,9,27,81 the four constants in (4) are respectively

    (1,1,1,1), (7,4,7,7), (7,13,16,7), (61,13,70,34).

Moreover replacing D by 27 Z_(>0) produces a genuine endpoint collision:

    B_1(4671)=B_4(5535)=5913,
    v_3(4671)=v_3(5535)=3.                               (6)

So uniqueness is false at the immediately weaker divisibility layer,
not merely unproved by the chosen residue test.

## Why reverse descent has not paid the reset debt

Every accepted suffix inverse decreases its positive integer argument,
because it reverses strictly growing maps. This gives finite backward
decoding of any supplied finite endpoint. It does not prohibit an infinite
forward sequence of larger endpoints.

The original stopping source is essential even inside D:

    6399 --B_1--> 8100 --B_0--> 12150.                    (7)

If 8100 is the marked source, blindly applying the last-nonzero decoder
to 12150 jumps past it to 6399. A correct controller first checks the
marked source against the intermediate zero-run endpoints and stops at
the supplied clock. It cannot claim the extra backwards step as payment
earned by that source. In odd n-coordinates these are 12797,16199,24299.

This agrees with the stronger pre-existing collision analysis, THM-4555 /
[uniform switches and trailing-ones deletion](../../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md):
a smaller partner and a collision are useful only with the legal guard,
correct clocks, and an independently supplied certificate for the partner.

If (1) were proved uniformly, the inherited occupied-cell exponent s<1/2
would give source-weighted finite first-excursion growth moments:
7^(3/14)>3/2 beats the per-return growth envelope. This is a sufficient
**conditional** improvement. It still would not prove the terminal fall
pays all the way back to the original source.

The next precise helper problem is to combine (5) with a source-preserving
forward rank or a legal smaller-source switch. Any proposal must survive
(6), (7), all five letters, and arbitrarily long zero runs. A source-height
bound proved only after replacing the source by its decoded predecessor
does not meet that obligation.
