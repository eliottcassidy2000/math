# Independent audit: three-module clock identity

Status: **PASS, no correction required.** The result is a scalar Hamilton-count identity for three nonempty whole modules, not a full output bin-profile identity.

Frozen files in `session197-odd-cycle-path-polynomials/`:

- `THREE_MODULE_CLOCK.md`: SHA256 `25ac0883cb93cd55b1109160bdb253794413e3aa3b56682bccf28f1dd6542941`.
- `verify_three_module_clock.py`: `c1ffec079a9ae7596e6017f651f18c9469f7674e4298f249b3ca678b63540320`.
- `three_module_normal.stdout` and `three_module_optimized.stdout`: `7591c0f37a6287d18a641e13fdf1bcc2c572892e5a725d2ca861bd0e80a4bf8e`.

## Proof checks

The three cross-graph covers are literal and uniform. Their covering vertices4 and6 are fixed singleton blocks. The other enlarged modules remain isolated from the relevant cross edges, even if their internal tournaments are not transitive. In the scalar-clone argument the variable module is first regarded as one root, and the other two blocks are fixed. Thus the respective run bounds3,2,2 hold at every positive choice of the other sizes, not just on the interpolation grid.

The argument passing from separate bounded exponential spans to their tensor product is valid: successive Vandermonde interpolation expresses each coefficient function as a fixed linear combination of evaluations and preserves the other separate span bounds. The twelve actual graph counts therefore determine the entire positive-integer clone family. No zero-size interpolation value is silently used.

For fixed other actual blocks, the inherited whole-module fragment kernel is linear on the complete bin polynomial. The monomials x^m for m>=1 span the zero-constant subspace containing every nonempty tournament profile. Extending successively in T, U and V proves exactly the stated trilinear scalar identity. The nonempty-block restriction is necessary to this particular interpolation argument and is stated throughout.

The inverse rule chooses positive sizes b-2,a-2,1 and has exactly a+b+3 vertices. Strongness of the quotient lifts to every nonempty substituted block. Equal sizes of the paired modules when a is odd give the displayed physical anti-involution; no unjustified exclusion for even a is inferred. The modular divisibility condition and the explicit composite diagonal counterexample are correct. The new 7*4999 graph uses a changed binary module and does not contradict the fixed-amplitude72 family's separate obstruction.

## Independent physical controls

The source normal and optimized runs reproduce the frozen receipt byte for byte.

The separate `check_three_module.cpp` constructs vertices as explicit block/local-index pairs and counts every output by odd-cycle collections, not the source endpoint DP. It reproduces all twelve interpolation graphs and all64 actual input triples, including nontransitive blocks. It also checks the exact cross-edge lists after collapsing each variable module to a root for all192 combinations with arbitrary other blocks. This directly controls the uniform-cover claim outside the transitive interpolation grid.

Five lifted anti-involutions and their exact graph counts pass. The independent14-vertex witness rows agree with the source:

    [1022,1020,1016,1008,992,2048,416,288,2080,5600,10751,735,11775,3071].

Its independent odd-cycle count is34993. Direct division by every integer from2 through70 certifies4999 prime. Optimized and unoptimized library-checking independent builds give identical output.

Independent checker SHA256: `5f38429813ca4a4bceb3ddbd4a46968871640032c592e6f1653756347d24ade6`.
Independent receipt `three-module-independent.txt` SHA256: `9ecc53992bb2c18a86faeee0aeddad1d047d13e0edf85d02eb7ad8632da39ef4`.

The inherited whole-module linearity and run framework are appropriately credited. The all-order three-variable covers, trilinear identity and inverse family are the new consequences reviewed here; no literature-priority or unrestricted-spectrum claim is made.
