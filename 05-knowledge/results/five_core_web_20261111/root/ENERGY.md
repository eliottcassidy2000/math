# Binary path energy has an exact local form and a six-periodic radical

**PROVED IN THIS NOTE:** the integral banded form, determinant and radical
below. **FINITE-EXACT:** the accompanying controls. **OPEN:** integer target
coverage. This is an operation-bearing connection to the old homometry and
native-move work, not a claim that six-periodicity explains the five cores.

## 1. Recover the actual energy

The pair energy is inherited, not newly discovered here. See
[THM-5093, endpoint-half-bound-equality-and-binary-pair-energy](../../../../01-canon/theorems/THM-5093-endpoint-half-bound-equality-and-binary-pair-energy.md)
and [Session220's extrema proof](../../session220_binary_envelopes_and_selectors_20261107/extrema/PROOF.md).
For the binary digits `x_0,...,x_(m-1)` of n, least significant first,

    f(n) = sum_(i<j) w_(j-i) x_i x_j,
    w_1=2,   w_d=2^(d-1) for d>=2,
    g(n)=n+1-x_0-f(n).

Thus f is a weighted autocorrelation of the selected bit positions. June's
[signed-LRC homometry reflection](../../../../07-reflections/signed-lrc-as-cyclic-homometry-s705.md)
suggests preserving the whole difference multiset before specializing its
weights. Only that operation is inherited here; no LRC theorem transfers.

Reversal of exactly m bits preserves f. Hence, with r the reversed integer,

    g(r)-g(n)=r-n+x_0-x_(m-1).

Also `g(2^m-1-n)=g(r)`. To see the latter, the sum of weights incident to
bit i is `2^i+2^(m-1-i)-[i=0]-[i=m-1]`, while the full all-one energy is
`2^m-2`. Expanding the complement gives the identity. At m=1 it holds
directly. These identities retain a value displacement; energy equality
alone is not equality of a coloured-core response.

## 2. An integral change of coordinates makes the energy local

Let W_m have zero diagonal and off-diagonal entries w_(|i-j|), so
`f(x)=x^T W_m x/2`, now defined on arbitrary real or integer digit vectors.
Let D_m be lower bidiagonal, with diagonal1 and subdiagonal−2. Its
determinant is1. Direct multiplication gives

    D_m W_m D_m^T = -2 C_m,                              (1)

where C_m has the following nonzero entries (indices start at0):

    C_00=0,       C_ii=4 for i>=1;
    C_01=C_10=-1;
    C_(i,i+1)=C_(i+1,i)=-3 for i>=1;
    C_(i,i+2)=C_(i+2,i)=1.

Terms at distance at least3 cancel by the geometric weight recurrence.
The exceptional adjacent weight2 is exactly what leaves the second band.
This proves (1) at every m, including the truncated end rows.

Put `x=D_m^T y`. Then

    f(x)=-y^T C_m y,
    y_i=sum_(j>=i) 2^(j-i) x_j=floor(n/2^i)             (2)

for a binary vector. A complete, lossless local description is therefore

    y_i-2y_(i+1) in {0,1},      y_m=0,
    f=-y^T C_m y.

The energy has range two, but the integer registers y are unbounded.
This is compatible with the inherited finite-state selector obstruction:
locality of equations does not imply a bounded state alphabet.

## 3. Closed determinant and exact null direction

For m>=1,

    det(W_m)=(-2)^m a_m,

where

|m modulo6|a_m|
|---|---|
|0|1−m|
|1|0|
|2|m−3|
|3|m−5|
|4|−4|
|5|−m−1|

In particular, for m>=2 the form is singular **exactly when m=1 modulo6**.
Its kernel then has dimension one and is spanned by the length-m prefix

    v=(1,0,−2,−3,−2,0, 1,0,−2,−3,−2,0, ...,1).       (3)

Here is a finite recurrence proof, rather than extrapolation from a list.
For C_m y=0, choose y_1=a,y_2=b in one-based notation. The first two rows
give y_3=b,y_4=a−b. Every remaining recurrence is

    y_(j+2)=3y_(j+1)−4y_j+3y_(j−1)−y_(j−2).           (4)

The right boundary is `y_(m+1)=y_(m+2)=0`. Write y_j=p_j a+q_j b.
The initial six values are

    p: 1,0,0,1,2,2;
    q: 0,1,1,−1,−4,−6.

The characteristic polynomial of (4) is
`(z−1)^2(z^2−z+1)`. Substitution in (4) shows for every j that

    p_(j+6)=p_j,          q_(j+6)=q_j−6.

Eliminating y_3,...,y_m using the unit second-superdiagonal coefficients
gives

    det(C_m)=p_(m+1)q_(m+2)−p_(m+2)q_(m+1).

The sign is positive: the eliminated columns move past the first two
columns, making `2(m−2)` exchanges. The remaining two equations are the
right-boundary equations; replacing the penultimate residual by its
negative and the last residual by `3y_(m+1)−y_(m+2)` has determinant1.
The formula is also immediate for m=2. Inserting the six residue classes
gives the displayed table. At m=6k+1>=7, the two p coefficients vanish
and the two q coefficients equal `1−6k`, so b=0 and a is free. Transforming
the periodic p vector by D^T gives exactly (3). The m=1 case is direct.

## 4. Why a promising formal move is not yet a tournament operation

For m=6k+1, the vector v satisfies

    sum_i 2^i v_i=1,     v_0=1,     W_m v=0.

Consequently the polynomial extensions of both f and g are unchanged by
`x -> x+t v`, while the encoded integer changes by t. This looks like a
perfect count-neutral move in a continuous or signed-digit relaxation.

It has **no nontrivial binary realization** for m>=7: a difference of
binary vectors has every entry in {−1,0,1}, whereas v has entries−2 and−3.
Since the kernel is one-dimensional and its first entry is1, any integral
kernel difference is an integer multiple of v. Only the zero multiple is
a legal binary difference. The continuous flat direction is real, but
rounding it is the unresolved operation, not a harmless normalization.

This is a concrete transfer of the October
[native-move lattice lesson](../../session46_native_packet_decomposition_20261024/move-lattices/THEOREM.md):
an integer relation or a neutral quadratic direction does not by itself
give a path inside a capacity-constrained fibre. Here the exact missing
constraint is binary digit legality, rather than unnamed geometric data.

For a fixed coloured-core compiler, holding the opposite flank B fixed,
the formal move changes its polynomial extension by

    t*(LB+i+s),

because g is invariant, A changes by t and the low-digit coordinate by t.
Thus even before digit normalization, retaining the full core response
changes what “neutral” means. No all-target or even legal binary switch
is asserted from this identity.

## 5. Connection contract and next cheap test

Source: weighted binary autocorrelation, suggested by old cyclic homometry.
Target: the exact path-count compiler's nonlinear coordinate. Map: selected
positions to W_m energy, then the unimodular prefix transform (1).
Preserved: exact integer energy, address and all binary legality constraints.
Lost if y is bounded or relaxed: actual digits and the target count after
normalization. Sidecar: `y_i−2y_(i+1) in {0,1}` plus the marked core cut.

The new useful test is not another spectral-ratio fit. It is whether a
bounded core change can absorb the explicitly computable energy defect of
normalizing one of these signed moves. A successful rewrite must specify
actual new binary addresses and keep the seven-context response, rather
than merely remain in the nullspace of W. This remains OPEN.

The determinant's sixth-root factor and the golden overlap graph's
fifth-root geometry are different operators on different carriers. The
present result explicitly keeps them separate while making their lost
coordinates available for a future actual map.

Run `python3 -B 05-knowledge/results/five_core_web_20261111/root/energy_checks.py`.
The receipt records exact matrix controls, independent digit recurrences,
the explicit singular vector and a binary-legality hostile. The proof above
is the all-length justification.
