# A sparse physical spectrum can pass every previous local test

**Proved.** There is a power-sparse, multiplication-closed subset of the
actual tournament spectrum that preserves the previous literal +2 windows
on both sides, growing congruence precision, the entire fixed-modulus
family, and any prescribed finite set of achieved counts. It can still
omit almost every prime. This is not a claim that the unrestricted
spectrum omits any additional value.

The older [genus reflection](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/07-reflections/the-genus-is-the-local-global-gap-and-the-one-master-dichotomy.md)
asked what is not determined by local observations. The
[planar inverse work](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/kervaire_superstring_jacobian_20261007.md)
makes a related distinction between preserving finite observations and
having one literal global lift. Here the loss is measurable.

## The exact information cost on a height scale

Suppose every odd N in a scale [X,2X], and a comparable O(R)-enlargement,
has actual three-member windows with spacings +2,+2 on both sides,
within C R, with a designated member in its residue class modulo M.
Assume 1≤M≤R=o(X). Divide the scale into boxes of width between 2R and
2R+3. For every box and residue class with odd representatives choose
one representative and keep its *whole* pair of window triples.
For any target use the previous box's left window and the next box's
right window in its residue. These boxes each contain an odd representative
because their width is at least 2M. The full windows are strictly on
the required sides and within O(R). A distinguished center j chosen
by the residue mod 3 remains the same if 3 divides M.

Only O(XM/R) pieces are needed. Conversely one count can serve
O(R/M) odd targets by its modulus and proximity, so any set supporting
this test (even on one side) needs Ω(XM/R) distinct counts. This is
the sharp amount of data conveyed by the local theorem, not Θ(X)
as exact placement would require.

## Preserving both known kinds of local tests

Use the [polynomial clock](../chanpak_bott_cayley_periods_20261115/POLYNOMIAL_PLACEMENT.md)
with R=X^(3/8) and M comparable to X^(1/160), taking an allowed
3 times a power of two at that scale. The same M is valid in an
enlargement by reducing the exponent by an absolute constant,
since all target sizes are comparable there. The resulting witness
bank has O(Y^(101/160)) counts ≤Y. This exponent is 1−3/8+1/160.
The previously proved radius and precision up to constants survive.

To preserve the [fixed-CRT bounds](../chanpak_collatz_norm_clock_20261115/FIXED_CRT.md)
simultaneously, for integer e≥1 use M_e=6 lcm(1,…,e) and R=X^(1/3).
Start this lane only above a dyadic scale index J_e≥e, large enough
that the fixed-M_e theorem applies throughout the enlarged scale,
its retained counts are comparable to X, and M_e≤X^(1/e) and
M_e≤X^(1/3). Its error constant affects the distance but not the
number of boxes. At a given scale the sum of active M_e is O_ε(X^ε)
for every ε>0: finitely many e<2/ε contribute constants; each other
one is at most X^(ε/2) and there are at most log₂ X active indices.
Summing the banks over dyadic scales yields a union of size
O_ε(Y^(2/3+ε)) below Y. Each fixed m is eventually covered by
some M_e divisible by m, at radius O_m(N^(1/3)).

## Ordinal closure still need not create the missing values

Adjoin any finite verified census E, omit 1 as a generator, and close
multiplicatively. If A is the set of retained generators, then for
every s>2/3 the convergent sum Σ_(a∈A) a^(−s) gives

    Σ_(h∈monoid) h^(−s) ≤ Π_(a∈A) (1−a^(−s))^(−1) < ∞.

The right side overcounts factorizations, which is harmless. Ordinal
sum shows every h in that monoid is a physical tournament count.
The bound gives O_s(Y^s) counts ≤Y; the polynomial lane alone permits
s>101/160. Since prime counts below Y are ≫Y/log Y, most primes
can still be absent from the monoid. This is a logical countermodel
to the proposed sufficiency of the local tests, not to Chan–Pak.

## The missing datum can be stated as a precise lifting task

Choose two literal windows with designated counts on opposite sides
of N in a fixed residue M. If they lie in the same physical component
whose moves preserve that residue and change the designated count by
at most M, then the integer quotient (H−N)/M changes in steps of at
most one. The first crossing is necessarily an exact hit. The theorem
does not construct this component; it distinguishes the operation
needed from the existence of two independently chosen near-hits.
