# Seidel energy and dyadic quotient noncancellation

For every tournament quotient Q of order m>=2, the formal cover profile

    g_k=(-2)^(k-1)

is a substitution fixed point. Although all scalar cover coordinates forget
Q at this point, its full two-port response does not. The shared cubic bank
has the exact positive energy

    sum_(i!=j) S_ij [a_i a_j c_i] Delta
        = 2^(6-m) sum_(nonempty even I) ||S_I^(-1)1||^2 > 0.

The complete cover contrast is `(delta,2*delta,0,...)`, giving an odd ambient
multiplier. Set L=1+v2(energy), computed by exact rational inversion. Then
arbitrary unequal transitive child sizes that are positive multiples of2^L
are separated from quotient-converse targets in every fixed ambient. For a
transitive source, target child interiors may be arbitrary on the same sets,
after pairing queries by converse-complement duality. This is source-relative,
not all-pairs reconstruction.

For every odd carousel order n>=3, an exact closed form gives

    [a_0 a_((n-1)/2) c_0] Delta = -8(n-4).

Hence modulus16 suffices: child orders divisible by4 are already enough.
The proof uses an explicit independent-sign bijection, not finite-order
extrapolation. The admissible normalized-profile class is closed under every
quotient substitution, so it also contains recursive constructions.

`PROOF.md` gives all model and scope conditions. `verify.py` produces exact
complete block-uniform source/converse banks, including all zero and quartic
coordinates. Run normally and with Python optimization; the two receipts
must agree. Independent audit packets retain distinct literal named physical
banks and direct run-word controls.

General arbitrary child profiles remain outside this theorem. The geometric
profile is a rational/2-adic functional, not a claimed finite tournament.
