# Least missing counts, local invertibility, and a hidden winding coordinate

Chan–Pak remains open. The session proves a concrete local-to-global
failure, identifies which naive induction steps cannot fix it, and gives
two sharp distinctions between arithmetic and geometric symmetry.

## A mathematical failure at infinity inside one tournament grammar

[The modular escape theorem](MODULAR_ESCAPE.md) constructs strong
tournaments consisting of one apex and a transitive spine that realize
**every odd-compatible pair (path count, contact address) modulo every
modulus**, with even address where parity is specified. The proof first
sets binary digits. On the odd part it uses an exact affine commutator
of physical append words, with central step 3/4. A separate physical
reset bypasses the apparent obstruction at 3 without disturbing the
binary residue. These are simultaneous composite-modulus certificates.

The same grammar misses the exact count 13, although an explicit
five-vertex tournament outside it realizes 13. Any modular lift of
13 in this grammar must grow in active size as the modulus grows.
Even the full two-coordinate congruence state therefore fails to
capture absolute realizability. To turn a modular descent into
a proof one needs control of the actual witness height or a way
to carry the exact count, not additional local residue tests. The
[exact width envelope](HEIGHT_ENVELOPE.md) shows 13 sits
strictly between the width-four maximum and width-five minimum.

## What a smallest-counterexample proof must use

[Black-box substitution and the symmetry budget](
BLACK_BOX_AND_SYMMETRY.md) starts from the inherited reduction to
prime and 7p targets. A fixed whole-module construction that knows
only a witness's scalar H can only multiply copies of H: every
revisited module exposes higher ordered cover counts, which grow
under H-preserving padding. Even a finite bank on fixed input
counts and all neutral pads has only polylogarithmic range.

Every witness for a prime p≥5 is rigid. At 7p with p≥11 prime,
its automorphism group has order 1 or 7. There really are symmetric
7p witnesses (7·593 and 7·2297); the cofactor p is never geometric
symmetry in either prime lane. This does not rule out missing p.

[The winding lock](WINDING_LOCK.md) is a different kind of hidden
symmetry: two actual strong six-vertex cores share their scores,
automorphism group, and full cover polynomial. Reversing just the
second triangle reverses relative winding. It changes a painted
one-apex chart that reaches all three residues modulo 3 into one
that stays 3-divisible under every ordinal left/right continuation.
Thus even all ordinary cover moments do not detect the local
symmetry relevant to constructing a prime witness. For several
apices, [the sparse contact test](SPARSE_CONTACT_TEST.md) proves
that masks of size at most two per apex determine every exact
or modular cut response in a fixed base. This identifies a small
complete set of diagnostics, not a coverage certificate.

## The Jacobian comparison, with the exact hypotheses

The [Jacobian bridge note](JACOBIAN_BRIDGE.md) follows
[Gao's current primary account](https://arxiv.org/html/2608.00222v2#S3):
the complex constant-determinant counterexample in dimension three
extends to larger dimensions. The two-dimensional complex problem
remains open. The weaker everywhere-nonzero real condition had
already failed in dimension two, so is a different comparison.

An affine-fiber sweep with constant nonzero determinant is rigid in
any dimension; the new complex examples instead use a monomial twist
of the tangent sweep and fail properness. The useful proposed analogy
is the loss of **coherence at the boundary**, not a shared equation:
a local inverse or unit switch can look identical on every inspected
neighborhood while its lifts run off at unbounded size, or while a
missing orientation prevents the charts from gluing. For tournaments
the exact modular compiler and the relative winding example separately
witness those two mechanisms. Neither decides the planar Jacobian
problem or Chan–Pak.

Run `python mod_replay.py`, `python height_replay.py`,
`python replay_symmetry.py`, `python replay_winding.py`, and
`python replay_sparse_contacts.py` from this directory. Tests use literal
tournament subset counts for the stated witnesses and exact modular
arithmetic for the compiler. The all-modulus and all-order claims
are proved, rather than extrapolated from those tests.
