# Minus growth, complementary products, and Machin's two integer centers

The starting point is the [return-gate and all-plus construction](../chanpak_power_three_return_gates_20261118/README.md).
There are two exact repeatable minus moves here. One genuinely triples
the odd part but carries a fixed factor 17; the other stays on the
2,3-smooth centers and squares instead of tripling. Neither supplies
the missing pure 3^j minus tripling.

## 1. A genuine repeatable minus move with exact mode cancellation

Take deletion core D and apex outmask q with outgoing bit rows

    D=(30,156,8,128,12,671,959,277,575,159), q=31.

The only arcs from apex out-contacts to in-contacts are 1→7 and 3→7.
Thus a path returning to a clone must spend the common vertex 7;
across disjoint slots this can happen only once. Insert k ordered
unmarked twins immediately before the unchanged apex mark. With
h=H(deletion), κ=H(core including the mark), and W the two-slot
count whose first slot ends at an in-contact or is empty, the first
reads are

| k | h | W |
|---:|---:|---:|
| 0 | 91 | 1224 |
| 1 | 181 | 3672 |
| 2 | 361 | 11016 |
| 3 | 721 | 33048 |

There are at most two clone runs in h and at most three in W. The
first reads therefore determine both exactly in the onto-run basis
as in the return-gate theorem. They give

    h_k = 90·2^k+1,   W_k = 1224·3^k,   κ_k−2h_k=−1.

With the spine of width m≥3 tagged `11` followed by zeroes, the same
physical splice as before gives the strong family

    H = 17·3^(k+3)·2^(m+1)−1.

Every core vertex can reach the apex, and the apex reaches all via the
left side, proving strength. At k=0,m=3 the count is 7343=7·1049,
a genuine prime quotient (trial division through 32 suffices).
Relative to the previous minus example, both lower modes in W
have vanished, rather than become small. The odd 17 remains in the
leading coefficient and cannot be divided out by a spine edit of
this construction: the center has form aW2^t with positive integer a.

The ordinary positive odd shortcut does the exact next step:
for C=17·3^(k+3)·2^(m+1), v₂(3(C−1)+1)=1 and
the result is (3C/2)−1. This is the tournament with one extra clone
and one fewer trailing spine zero, as long as the resulting width
is at least 3. Thus choosing an initial width at least r+3 gives r
successive *ordinary* Collatz steps with literal tournament
witnesses from this one construction. It neither controls the
terminal low width nor proves anything on global convergence.

For genuine nonzero spine growth, a scoped obstruction on the positive seed explains why trying
a different fixed one-run flank is insufficient. Its reads are
W=36·3^k, h=12·2^k−1, R=W−h. Let a be a fixed left prefix and b a
fixed suffix of r bits, with G=2^r and c the reverse complement of b.
For a growing zero-run, F of that tag word has leading coefficient
a(1+c) in 2^m and its value has coefficient aG. For a growing
one-run these coefficients are (a+1)c and (a+1)G, respectively.
The spine formula has variable parts hF(u)+Ru (if the last tag
wins, there is an extra constant parity correction). Thus killing
the 2^k term requires c=G−1 in the zero case or c=G in the one case.
The latter is impossible. The former forces an all-zero suffix,
and the intercept is then 1 plus a fixed integer times h. It
cannot be −1 for all k as h grows. This is a fixed-prefix/suffix,
one-growing-run, no-right-flank obstruction, not a global exclusion.

## 2. Complementary products remove unwanted modes on smooth centers

If both C−1 and C+1 are known counts, ordinal gluing gives the
count (C−1)(C+1)=C²−1. It is a literal signed correction: the two
linear terms cancel exactly. Take C=3^j·2^m with j,m≥3 from one of
the known paired families. The universal plus construction provides
C²+1 at the new center, so this step repeats without end. It proves
both signs at 3^(2^r j)·2^(2^r m) for every r≥0, for example from
j=3,4,5,6,7 and every starting m≥3. There is no extra prime in
these centers. It does not increment j by one or preserve the
same width, as the true signed Collatz step would.

A polynomial rigidity explains the limitation. In a product of
factors of the forms X^d−1 and X^e+1 equaling X^n−1, exactly one
minus factor is allowed, by the simple root at 1. Its degree d
divides n. No plus factor vanishes at a root of unity of odd
order, so the largest odd divisor of n must divide d. Thus n/d
is a power of two. In particular a product consisting only of
a realized X−1 and plus-binomial factors cannot produce a
universal cubing law X³−1. It needs the odd cyclotomic factor
X²+X+1, the same factor that underlies the earlier
[Nagell–Ljunggren comparison](../chanpak_odd_center27_20261118/README.md).
This identifies a separate possible next target: a physical
realization of that cofactor for the dyadic/ternary centers.
The old Session 27 repunit audit already recorded Φ₃(2)=7 and
Φ₃(4)=21, the two forbidden counts. The next cofactor Φ₃(8)=73
is realized, and Φ₃(18)=343 gives the Nagell–Ljunggren cube.
This now explains a precise failure of the *factor-based cubing
method*: at C=4 it asks for 21 even though the whole target
4³−1=63 is realized by a different tournament. It does not explain
the global exclusions 7 or 21. The cubing step is also not the
same as multiplying the odd center by 3 while peeling a single
power of 2.

## 3. Exact landings at 57122 and 114244

Machin's identity follows from
(5+i)^4=4(119+120i)=2(1+i)(239+i). In the quadratic unit orbit,

    (1+√2)^7=239+169√2,     239²+1=57122=2·13⁴,
    2·57122=114244=338².

These are not just nearby magnitudes: the two centers both have
physical unit repairs on each side:

| Center | Lower count | Upper count |
|---:|---|---|
| 57122 | 57121 = 239² | 57123 = 3²·11·577 |
| 114244 | 114243 = 3·113·337 | 114245 = 5·73·313 |

Each factor is a literal path count in the replay. Their ordinal
products verify all four absolute landings, without claiming that
the two witnesses are related by an arc flip. The two centers
are even and cannot themselves be Hamiltonian tournament counts.

In the Pell system, x+y√2 moves by the unit 3+2√2 and preserves
x²−2y²=−1. Odd-power real parts start 1,7,41,239. In particular 7
is allowed in this system though forbidden as a tournament count.
No multiple of 3 can be such a real part (reduce the Pell equation
modulo 3), so 21 and 63 are *both* excluded there, while 63 is
allowed for tournaments. This is a countercheck on identifying the
two spectra. The exact overlap of labels is still striking:
119=(239−1)/2=7·17, so the odd prime left over in the minus compiler
also occurs in one Machin leg. That identity gives no way to cancel
17 in the tournament count. For why the quartic equation has only
its two positive solutions, see
[Steiner–Tzanakis on Ljunggren's theorem](https://users.math.uoc.gr/~tzanakis/Papers/LjunggrenEq.pdf).

## 4. Circulant Hadamard: exact zero correlations, not small ones

The [published Lean formalization](https://github.com/openai/math/blob/main/lean/docs/179.md)
classifies positive orders of real circulant Hadamard matrices as
1 and 4. At 1 the witness is [1]. At 4 the witness is the
circulant J−2I (first row −+++): every nontrivial cyclic shift
has two agreements and two disagreements. An explicit bridge to
a tournament is available **at 4**. Flip columns 2 through 4
and then all four rows, using indices 1..4. The first row is all
positive, the first column below the corner all negative, and
the bottom-right block has negative diagonal and positive
off-diagonal. Apply a 3-cycle to its three columns (including
their top entries). The bottom-right block now has positive
diagonal and one negative in each reflected pair: it is I+A
for a directed triangle. The matrix is precisely a skew Hadamard
border of that tournament, up to relabeling.

The shared proof lesson is a condition on exact lower modes:
Hadamard needs every nonzero autocorrelation to vanish; Machin
needs a specific unit determinant to obtain its little
correction; the tournament clone needs both lower exponential
coefficients of W to vanish and the signed defect to stay −1.
None of these conditions follows merely from their dominant
growth. The identity or classification in one problem does not
supply the missing realization in the others.

Run `python3 verify.py` for literal count factors, the minus-17 read,
the Gaussian-Pell arithmetic, and the order-four matrix conversion.
