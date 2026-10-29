# Wedge coordinates, fixed-port recurrences, and integer descent

Session `wedge_descent_168_20261029`; source revision
`96e3f3fe347ac3060a09c88988809e8e5aacfaad`.

**Status:** PROVED elementary propositions with same-agent exact checks;
CITED exceptional group isomorphism; inherited independently audited
tournament classification. No new canon ID or independent-referee claim.
Collatz, ABC, LRC(14), and unrestricted tournament enumeration remain OPEN.

The new consequence of Session103 is a bounded-degree composition law.
For an n-vertex digraph with a marked Hamilton cycle whose auxiliary digraph
is acyclic, normalize its ordered-k-path observer by the canonical port
matching. The resulting native integer matrix B satisfies

    B = I+N,                 N^(k(n-k)+1) = 0,
    B^d = sum_{j=0}^{k(n-k)} binom(d,j) N^j,
    B^-1 = sum_{j=0}^{k(n-k)} (-N)^j.

Thus repeating the same port-matched observer has polynomial, rather than
exponential, dependence on composition depth. This concerns multiplication
of the declared observer tables; arbitrary tournament substitution is not
identified with that multiplication. The order-6 `ca` example with k=3
attains the degree bound 9. The full proof is [W1](THEORY.md#w1-fixed-port-nilpotence-and-exact-recurrences).

A closure criterion says exactly when the mechanism survives a menu of
operators: for nonnegative matrices with diagonal one, **every finite product
is unipotent if and only if the union of their off-diagonal supports is
acyclic**. A directed cycle constructs a failing product of length at most
the number of port states. This is [W2](THEORY.md#w2-a-decisive-test-for-closure-under-a-menu-of-port-operations).

The native hostile already occurs at four vertices. The seed's two-port
matrix has dimension12. Compose its normalized matrix with the same matrix
re-marked by the vertex permutation `(1 3)`: the product has trace22.
It has exponential repeated-composition growth. Both factors separately
have determinant1 and are unipotent. Keeping only these invariants cannot
support the recurrence under re-marking.

Other exact connections recovered and checked:

- **K4 to 168:** the even-weight subspace of F2^4 has dimension3. Its seven
  nonzero vectors are the six edges plus `1111`. Its full linear group has
  order168; fixing `1111` gives precisely S4 acting on the four vertices.
  Through the classical GL(3,2)=PSL(2,7) isomorphism this connects to the
  source-adjoined Paley7 switching group. The six-state language machine,
  with named ports and its outputs retained, instead has trivial symmetry.
- **Two different integral six-edge lattices:** the oriented cut/cycle
  split has quotient (Z/4)^2, order16; the inherited unsigned isotypic split
  has quotient (Z/2)^2 plus Z/6, order24. Orientation and integral gluing
  cannot be discarded because both presentations have six coordinates.
- **The user's two arithmetic rows:** reflection about c=2k(k+1) gives
  both balances. The square total is
  `k(k+1)(2k+1)(12k(k+1)+1)/6`. The same adjacent-block pattern has no
  positive-integer solution for cubes: its only positive center lies
  strictly between `3k(k+1)` and `3k(k+1)+1`.
- **169 and clocks:** in F169, norm and Frobenius ratio give 84 compatible
  pairs from168 nonzero elements; exactly one sign bit is missing. This
  field is different from the split quadratic ring O/13 in the older
  g58/Zagier calculation. Their identical state count169 must not identify
  their operations.
- **ABC square lift:** `(b-a)^2+4ab=(a+b)^2`, normalized to be primitive,
  has an exact radical-growth formula. It generally raises height; `1+3=4`
  is a fixed point. It provides a precise multiplication/addition bridge,
  and an explicit obstruction to using this lift itself as descent.

The descent lesson is constructive: a finite common acyclic support gives
an integer rank, explicit inverse, and recurrences. For Collatz the analogous
all-source guarded-return coverage remains missing. We retain exact first
descent words for999 positive odd sources, the nonroot 5x+1 cycle, and the
81/27 tournament pair. Finite verification is not universal coverage.

## Reproduce

From the repository root, with SymPy and python-flint:

```text
python -X utf8 -B 05-knowledge/results/wedge_descent_168_20261029/verify.py
python -X utf8 -B -O 05-knowledge/results/wedge_descent_168_20261029/verify.py --output <temporary-output.json>
```

Compare the second JSON byte-for-byte with [certificate.json](certificate.json).
The declared universe consists of35 native ordered-path matrices of orders4--7,
all512 binary 3x3 matrices, all4096 pairs of three-state off-diagonal
binary supports, all24 K4 vertex permutations, all720 permutations
of the six machine states, all168 nonzero F169 elements, reflection rows
1<=k<=300,125 offset triples, coprime pairs1<=a<b<=60, and odd sources3..1999.
Every truth gate stays active under Python optimization. No floating-point
arithmetic is used. The script independently compares matrix nilpotence with
longest paths in the support graph. The inherited 81/27 audit reconstructs
both tournaments and their groups, and checks their saved native response
matrices modulo3.

[THEORY.md](THEORY.md) gives proofs, quantifiers, and failure boundaries.
[RECOVERY.md](RECOVERY.md) gives the portfolio, source/map/loss ledger, and
primary literature. [manifest.json](manifest.json) pins source and output
bytes. [HANDOFF.md](HANDOFF.md) gives the next decisive tests.
