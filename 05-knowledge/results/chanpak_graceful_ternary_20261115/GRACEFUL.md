# Graceful matching buffers and the exact obstruction on a single apex

**Proved:** an exact residual decoder from bounded-spacing switches; an
explicit warning that gracefulness without an ordering is insufficient;
a literal all-even-gap two-reader construction; and a complete
interaction test for two independent interval controls on one apex.
In particular, a one-apex two-control additive rectangle cannot cover
four consecutive odds. This is a scoped structural obstruction, not
a claim about arbitrary tournament gadgets or global Chan–Pak.

## 1. A graceful ruler with resource conflicts

The [earlier graceful and square-sum exchange note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/graceful_square_exchange_20261019.md)
distinguishes local midpoint exchanges, a two-end splice, and the
three-edge threshold for fixed-endpoint graceful exchanges. Its
Walecki path on `{0,...,m}` alternates
`0,m,1,m-1,2,m-2,...`. Its edge lengths are exactly
`m,m-1,...,1` in path order, not merely in some permutation.

**Bounded-spacing selector.** For integers m,k≥1, every integer in

    0 ≤ R ≤ C(m,k) = r m − k r(r−1)/2,    r = ceil(m/k),

is a sum of a subset of `{1,...,m}` whose chosen labels differ by
at least k. For k=2 the labels are the edge lengths of a matching
in the Walecki path, and C(m,2)=floor((m+1)^2/4).

**Proof.** Use all labels at most m−1, or take m and then use labels
at most m−k. The maximum is obtained by choosing m,m−k,m−2k,...
since each decreasing sequence with spacing k is coordinatewise
bounded by that greedy list. In the recursion the first interval
ends at or after m−1 because its single largest label is
available, and the second starts at m; hence the intervals meet
without a gap. Take C(t,k)=0 at t≤0. For explicit decoding, omit m
if R≤C(m−1,k); otherwise use m, subtract it, and recurse at m−k.
All steps operate on the actual target rather than enumerating witnesses.

Gracefulness alone does not give that conclusion for arbitrary
ordering: the graceful path (1,3,2,5,0,4) has lengths (2,1,3,5,4),
and its matching sums are `{0,1,2,3,4,5,6,7,9}`, missing 8.
Thus the bare numerical invariant is not the necessary
continuation data; the position of the conflicts is.

**Physical transfer criterion.** Suppose a literal host has switches
indexed by the Walecki lengths, and for every k-spaced choice S the
simultaneous switch is a legal tournament whose count is precisely
H₀+2 Σ_(i in S) i. Then its quadratic-scale repair interval is
guaranteed by the decoder above, even if other choices cannot be
combined. For a count *below* the baseline the same statement
works with all moves reversed or negative weights. This is a sufficient
condition, not a declaration that a host meeting it is available.
In particular it does not replace proving a physical absolute
baseline H₀ close enough to the target.

## 2. The exact one-apex seam is a product plus one parity term

Use a transitive spine with vertices 1,...,m in forward order and
one apex. A bit x_i is 1 when the apex beats i. Put h(w) for the
ordinary Hamilton-path count. The positive odd-cycle-packing
expansion gives

    h(w) = 1 + 2 ∑_(i<j) c(j−i) x_i(1−x_j),
    c(1)=1, c(d)=2^(d−2) for d≥2.                      (1)

Indeed no directed cycle avoids the apex and no packing uses two
cycles. For any endpoints i<j the spine segment order is fixed
and its intermediate vertices are freely chosen with even
cardinality. There are c(j−i) such choices. For the inherited
odd-cycle expansion see [Grinberg–Stanley](https://arxiv.org/abs/2307.05569)
and its use in the [old single-seam construction](../chanpak_target_address_20261115/seams/SINGLE_SEAM.md).

Split a word into U,V. Encode U by the ordinary binary integer a
(leftmost bit is highest); encode V by the reverse-complement
integer b (leftmost complemented bit is lowest). Let h(a) denote
h on any binary presentation of a; zero prefix padding has no effect
by (1). Converse plus reversing the spine identifies the count of V
with h(b). Then

    h(UV)=h(a)+h(b)−1 + ab + (a mod 2)(b mod 2).        (2)

The rank-one cross coefficient in h is exactly
2^(|U|−i+j−1) for the contribution of the i-th bit of U and j-th
bit of V (indices start at 1). This applies to all pairs except
the touching pair, whose weight is 2 instead of 1; its extra
coefficient is the last displayed product. This parity guard
cannot be silently dropped when asking for equality.

## 3. Full interaction test, including touching intervals

Each of two controls replaces a fixed-length block by one of two
distinct bit patterns, and all other bits remain fixed. In the
mixed difference compare settings as

    M = h_11−h_10−h_01+h_00.

If the changed sites of the first control all precede those of the
second and their nearest changed sites have distance d≥2, let u be the
last changed left site and v the first changed right site. Write δ_i
and η_j for signed bit changes. Formula (1) gives

    M = −2^(d−1) (Σ_i 2^(u−i) δ_i)(Σ_j 2^(j−v) η_j).

Both parentheses are odd because each facing extremal bit actually
changes by ±1. Thus M is **never zero over the integers**, even
though it vanishes at every dyadic precision below its valuation.
This is the literal path-count counterpart of the
[spaced checksum obstruction](../chanpak_bott_cayley_periods_20261115/SEPARATION.md).

For touching changed sites, split the full word at the boundary,
absorbing arbitrary fixed prefixes and suffixes into the two sides.
Relabel each setting so a₁>a₀ and b₁>b₀. The two differences are
positive integers because the presentations have fixed lengths.
By (2),

    M=(a₁−a₀)(b₁−b₀)+(a₁ mod2−a₀ mod2)(b₁ mod2−b₀ mod2),

where each parity term means the difference of the two residues.
Hence M=0 precisely when both numeric differences are 1 and
a₀,b₀ have opposite parity. This characterizes the otherwise
exceptional exact cancellation at distance one.

**It still cannot give four consecutive odds.** Assume without
loss that a=a₀ is even and b=b₀ is odd; (2) is symmetric. Then
the two steps in an additive rectangle are

    Δ_a = b+1−2 ceil(a/4),        Δ_b=a+h(b+1)−h(b).

The first identity follows because a→a+1 removes precisely the
cycles ending at the last site of the left word. Their weights
sum to ceil(a/4). For b odd let r≥1 count the trailing 1s of b.
The change to b+1 contributes 2^r from the new 1 and the r newly
zero lower digits. The coefficient of an old 1 at the next
position r+1 increases by 2(2^r−2); positions j≥r+2 have
strictly positive increase 2 c(j−r)(2^(r+1)−3). All other
contributions cancel. Thus h(b+1)−h(b)≥2, and equality holds
only for r=1 with no 1 in positions ≥3, namely b=1 or 5.
Also h(4)−h(3)=4 directly from (1).

Four consecutive odd counts in an additive rectangle would need
the two steps to have absolute values 2 and 4. If a≥4, Δ_b≥6.
If a=0, then to make Δ_a equal 2 or 4 requires b=1 or 3; Δ_b
is the same step in each case, not the other one. If a=2,
Δ_b≤4 requires b=1 or 5, so Δ_a is respectively 0 or 4 and
Δ_b=4. This exhausts all cases. In particular, even an
arbitrarily far nontrivial pair of fixed exterior blocks cannot
save this rectangle.

More broadly, three spatially disjoint interval switches already
have a nonzero mixed difference for the first and third. Those
controls are normally a legal matching in the graceful selector,
so this simple physical implementation of the quadratic buffer
fails. The theorem leaves open interleaved controls, nonlinear
selection, and marked multi-apex interfaces. The three-packet
modular decoder escapes precisely because it **uses**, rather
than discards, the interaction; that alone does not make it an
absolute decoder.

## 4. A positive marked reader and its missing absolute coordinate

There is a much simpler exact *difference* selector. Prefix an
arbitrary word w by two bits, and compare the prefixes 00 and 10.
Formula (2) with a=0 or a=2 gives exactly

    h(00w)=h(w),
    h(10w)−h(00w)=2+2D,
    D=Σ_(i=1)^m 2^(i−1)(1−x_i).

Thus an arbitrary positive even gap g uniquely specifies its
nontrivial spine word up to neutral padding, and the two
literal marked observations read every bit of D without any
approximation. The word length is O(log(g+1)). If
h(w)=1+2F(D), (1) gives 0≤F(D)≤D/2 by summing all
available contributions into each set bit (bit j≥1 contributes
at most 2^(j−1), and bit 0 none). Consequently the **upper** count
is at most (3/2)g, sharp within this reader for D=2^j. The lower
count is 1 exactly when w=0*1*, i.e. D=2^j−1. For gaps that are
not powers of two both counts have strong witnesses: remove the
zero prefix and the one suffix, which are ordinary transitive
source/sink pads. The remaining marked word starts with 1 and
ends with 0, giving a directed spanning cycle via the apex.
For power-of-two gaps this construction instead has a count-1
lower witness; there is no claim about other constructions for
those gaps. This is a concrete elementary two-reader mechanism,
not an external-priority claim for realizing even gaps.

What fails is substituting the gap for the absolute target.
The actual upper count is (3+2D+2F(D)); choosing D to make
the gap g also fixes the lower endpoint. For example this upper
lane misses the admissible count 11: D+F(D)=4 requires D≤4
and the first five values of D+F(D) are 0,1,3,3,6. The count
11 is attained by the cyclic join of two transitive 2-vertex
blocks and a singleton. Also D=2 and D=3 both give upper
count 9, but they have different exact gaps. This is an explicit
Berggren-style firewall: a marked pair has an exact inverse, but
the scalar projection loses its unique descent and its coverage.
Square-sum endpoint restrictions and graceful matching constraints
are useful here as tests of a lift, not as proofs of a physical
absolute anchor.
