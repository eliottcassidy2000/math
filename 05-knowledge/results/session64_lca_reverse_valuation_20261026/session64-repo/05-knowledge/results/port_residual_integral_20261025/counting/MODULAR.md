# A constant-tail quotient at fixed binary precision

2026-10-25. **PROVED:** modulo2^K, the integral moment transfer preserves
a constant-tail module requiring at most min(S,K) coordinates. An adaptive
sparse/moment implementation uses min(h,S,K) scalar weight coordinates
per nonterminal defect fibre, retaining the semantic sparse support.
**FINITE-EXACT:** `modular.py` checks counts, arbitrary constant-tail inputs,
and the precision/lattice hostiles below. This is absolute K-bit precision;
it does not assert the same budget for a quotient by a power of two.

## Inheritance and the exact observation

Use the integral gauge and notation in [README.md](README.md):

    Y_r=sum_m x_m*2^[r*(m-b)],   x_m=O(u)*w_m in Z,
    b=max_(d>=3) ceil(D(d)/(d-1)),   m>=b.

The incoming
`session63_integral_observers_root_arrangements_20261025/cancellation/RESULT.md`,
Theorem2, states the denominator-precision law: dividing an integral
numerator by a denominator of binary valuation e requires e additional
binary digits. Here every remaining denominator is odd, so e=0.

Fix K>=1 and work in R=Z/(2^K). Whenever r>=K, every term with m>b
vanishes, while the m=b term remains. Therefore

    Y_r = x_b = Y_K  in R,       r>=K.              (M1)

If S<K, retain all S channels. If S>=K, retain Y_1,...,Y_K, with Y_K
representing every higher channel. Terminal S=0 still has its one separate
coordinate Y_0. No higher-precision reconstruction is asserted.

This is not a statement about the entire integral product lattice from
the preceding theorem. Actual moment banks lie in the constant-tail
module after reduction. That module is generally larger than the image
of actual moment banks: their lower channels can satisfy further
congruences. We claim neither a minimal module presentation nor a rank
theorem over the nonfield R.

## Closure under every allowed transfer

The previous integral transition reads

    Y'_r += T * 2^(J0+r*delta) * Y_(ell*k+r),
    J0=A_ell(k)-v2(k!)+k*(ell*b-D(ell))>=0,
    delta=k+b-b'>=0.                               (M2)

T is the integer odd-factorial transfer including the signed beta(k,e).
Suppose the source is any member of the constant-tail module, not
necessarily a realizable moment bank. For any destination r>=K, the
source index ell*k+r is also at least K, so its value is the source tail.
If delta=0, the whole expression is independent of r. If delta>=1, its
power of two is at least K and it is zero in R. Thus every batch preserves
the constant-tail condition. Addition preserves it as well. This proves
closure for signed excess, coordinate deletion, skip steps, and merges.

For source indices above the stored range, read the tail coordinate.
Compute only r<=min(S',K). Those operations reproduce the complete moment
transfer modulo2^K exactly, without reconstructing any omitted channel.

Sparse conversion is also exact: reduce the same sums in(M1), keeping
only r<=min(S,K). No sign or positivity assumption is needed for these
module identities.

## Sparse support is semantic, not the set of nonzero residues

Use a sparse bank while h<=min(S,K), then convert to the compressed dense
module. Retain a sparse key even when its coefficient is zero modulo2^K.
Otherwise modular cancellation or divisibility could make the stored h
smaller than the actual support and invalidate the inherited monotone
support argument.

In the underlying exact counting transfer, each fixed u,m has the sign
(-1)^(u-m). Nonzero incoming terms at the same m never cancel as integers.
Every legal batch shifts m injectively by m->m+k. Thus hidden support
persists even when all its residues vanish. Once h>min(S,K), a destination
still has at least h>min(S,K)>=min(S',K), and dense storage never needs
sparse reconstruction. Sparse-only mergers convert at precisely that
threshold. The coefficient storage is min(h,S,K) per nonterminal fibre;
keys and semantic support labels are separate storage costs.

Ordinary A000568 still needs only its one actual-degree channel. The
constant-tail reduction improves the clone-excess vector representation;
it does not improve an already scalar ordinary state.

## Exact modular arithmetic and the denominator boundary

All quantities in

    T=beta(k,e)*O(u')/[O(u)*ell^k*O(k)]

other than beta are odd. Its integral value modulo2^K can therefore be
computed directly by modular inverses of the denominator's odd factors.
There is no need to form a large exact quotient first. The modular odd
factorials can be built from O(j)=O(j-1)*oddpart(j). The final division by
O(N) is also an odd-unit inverse. Thus the full absolute K-bit computation
needs exactly K-bit residue arithmetic; O(N) costs no additional digits.

By contrast, if the requested result is a(N)/2^B modulo2^K, reducing a(N)
modulo2^K first generally loses the answer. This recurrence can compute
a(N) modulo2^(B+K) and then divide the representative by2^B, provided the
proved baseline applies. The separate bounded-cycle precision algorithm
achieves normalized K-bit precision by cancelling each denominator and
baseline before reduction. These are distinct observation maps.

## Hostile boundaries and reproducible scope

* The arbitrary product-lattice vector (Y_1,Y_2)=(0,1), with K=1, is not
  a constant-tail input. An identity transfer preserves its unequal
  channels. Compressing this vector by pretending Y_2=Y_1 is invalid.
* The actual bank x_b=1 has every Y_r=1. The tail is retained, not set
  to zero. In particular a two-singleton completion of the empty prefix
  has a nonzero individual contribution modulo2, even though other signed
  completions can cancel it in the total order-two core count.
* A sparse bank whose integer coefficient is2^K is zero as a residue but
  still has a semantic cycle-count key. The implementation keeps that key
  through nonterminal transport; modular zero is not an absence proof.
* At n=50 the normalized ordinary count a(50)/2^25 is2 modulo4. Its
  absolute residue modulo4 is0. Dividing that already truncated zero by
  2^25 would give a false answer; absolute27-bit precision retains the
  required two normalized bits.

Run `python modular.py` and `python -O modular.py`. The checker compares
ordinary and core residues with the exact integer engine for n=0..40,
K=1..6, with an additional order44 core control. It checks arbitrary
constant-tail modules, sparse/dense commutation, zero-residue support,
and the normalized order50 hostile. No full tournament or Hamiltonian-
path enumeration and no timing comparison is performed.
