# Actual command words detect a full Boolean source law

**PROVED from the finite base certificate below; FINITE-EXACT base checks.**
The source family uses genuine strong atoms. Its ambient language retains
the full tournament and the selected Hamilton path after every operation.

## 1. Language and the eight-vertex certificate

In X_n the current path is0→1→…→n-1. Free coordinates are(i,j), j>i+1,
in lexicographic order, with bit1 meaning i→j. Define:

- C: reverse every free arc, leaving the current path arcs fixed;
- R_(i,j): remove the vertex at position i and insert it at position j,
  |i-j|>=2, provided the resulting word is a directed Hamilton path; then
  use that word as the new coordinate order.

All relocations and their inverses are named. Illegal commands are observable.
A length-d legality observation consists of the executable words of length
at most d. It does **not** include H unless explicitly stated.

Let A=X_8(509), B=X_8(253). Their outgoing-neighbor bitsets are

```text
A: 246,60,8,17,36,76,158,62
B: 246,28,8,17,36,78,158,62.
```

They differ only at the arc between1 and5: A has1→5 and B has5→1.
Both are strong and have trivial automorphism group. Exact counts are

```text
H(A)=171=9*19,       H(B)=173.
```

Of these paths,38 in A and40 in B use the changed edge consecutively;
the other133 paths are common. Thus a one-bit interpolation, with z=1
selecting A, is H(z)=173-2z. Permutations and subset DP give independent
counting implementations. Exhaustive40320-permutation automorphism checks
verify asymmetry; reachability verifies strongness.

The following complete guard table is the small proof certificate. Each
row lists all legal relocations from the indicated state, and the list is
the same for A and B. C is always legal.

| Already executed | Legal relocations |
|---|---|
| empty | (0,3),(1,7),(7,1),(2,7) |
| C | (0,2),(4,0),(5,0),(6,0),(0,7),(7,0),(3,1) |
| R_(0,3) | (0,3),(3,0),(0,7),(7,0),(4,1),(1,7),(4,7),(7,4) |
| R_(1,7) | (0,2),(6,1),(7,1) |
| R_(7,1) | (0,4),(1,7),(7,1),(2,7),(3,7) |
| R_(2,7) | (0,2),(2,0),(6,1),(2,7),(7,2) |

This proves equality of all depth-two legality observations. The word

```text
w = C, R_(3,1), R_(5,2)                         (chronological order)
```

is executable on A and not on B, so depth three is sharp for this pair.
The first two commands are legal on both. After them the rows are

```text
A: 6,244,200,243,229,197,129,1
B: 6,244,232,243,229,193,129,1.
```

The last relocation has vertex word(0,1,5,2,3,4,6,7). Its new arcs are
1→5,5→2,4→6: the middle arc is present only on the A branch.
All moved positions and all altered adjacency tests lie strictly inside
the eight-vertex block. No query uses an external or reset oracle.

The incoming `two_edge_transport_20261102/ESCAPES_AND_CARRIERS.md` has an
exhaustive depth-two injectivity certificate through n7. That is inherited
finite evidence, not an all-size theorem. The new pair refutes its natural
depth-two extrapolation; eventual all-size reconstruction remains OPEN.

## 2. All-size contextual extension

Write L→T→R for ordinal composition: every earlier block beats every later
one. Empty contexts are allowed. For any fixed L,R, replace T by A or B and
keep the concatenated selected paths.

**Theorem.** The two sources have identical depth-two legality observations
in the full ambient command language. Their Hamilton counts are respectively
171H(L)H(R) and173H(L)H(R), with H(empty)=1.

**Proof.** A Hamilton path in an ordinal sum must visit the blocks in order
and cannot return to an earlier block. Thus every legal relocation stays
within one block. Its extra boundary arcs, if any, are automatically legal.
Two successive relocations therefore reduce either to the base certificate
in one block or to independent unchanged guards in different blocks.
Words R,C and C,C cause no issue because C is total.

Only C,R needs an additional check: after C the old ordinal blocks need
not be modules. A guard can differ only if it tests the changed chord1--5
inside the eight-vertex block. This chord has length four, so it cannot be
the distance-two gap left by a removed vertex. Exactly four relocation
guards mention it. In local coordinates their new arcs are

```text
R_(1,4): 0→2,4→1,1→5
R_(1,5): 0→2,5→1,1→6
R_(5,1): 0→5,5→1,4→6
R_(5,2): 1→5,5→2,4→6.
```

Every one is internal to positions0..6; no context can change its truth
after C. The base certificate shows equality there. All other guards test
only equal arcs. This proves the observation assertion. The product formula
for H follows from the same ordered-block path decomposition. QED.

In particular the obstruction exists at **every order n>=8**, by appending
transitive singleton blocks. This does not assert that arbitrary padding
preserves arbitrary complement observations; the four internal guards are
the required sidecar.

## 3. A Boolean cube of nonisomorphic sources

For r>=1 and z in{0,1}^r, let

```text
T_z = T_(z_0) → ... → T_(z_(r-1)),   T_1=A, T_0=B.
```

Each block is strong, so these are exactly the SCCs in their unique linear
order. A and B have different H and are nonisomorphic. Therefore all2^r
words give different unlabelled tournaments. Since the cores are asymmetric
and the SCC order cannot be permuted, all T_z are asymmetric as well.

The contextual theorem, changing one coordinate at a time, shows that all
2^r sources have the same depth-two legality observation. The shifted word

```text
w_i = C, R_(8i+3,8i+1), R_(8i+5,8i+2)
```

is executable iff z_i=1. The moves preserve the block endpoints; every new
guard arc is internal, even after global C. Thus each bit is observable
at depth three with its source position retained.

Any deterministic exact observer for this language must distinguish all
2^r sources: two different words disagree at some i, and w_i distinguishes
them. This is a lower bound on an ambient closed carrier, not a claim that
these2^r sources are closed under C and re-marking. It amounts to at least
r bits on this family, not full quadratic reconstruction at order8r.

The scalar count retains only Hamming weight:

```text
H(T_z)=171^|z| 173^(r-|z|).
```

These are r+1 distinct values. Even after adjoining H to the depth-two
observer, one fibre contains binomial(r,floor(r/2)) inequivalent sources.
The order of the strong atoms is real information, despite identical atom
multisets, zero automorphism ambiguity, and equal total H.

## 4. Actual reset words and the all-modulus Boolean inverse

Each successful relocation is undone by its named inverse, and C is an
involution. Write w_i^{-1} for the three commands in reverse inverse order.
For S={i_1,...,i_k}, k>=1, use the literal word

```text
w_(i_1) w_(i_1)^(-1) ... w_(i_(k-1)) w_(i_(k-1))^(-1) w_(i_k).
```

It has length6k-3 and is executable exactly when z_i=1 for every i in S.
On a successful prefix the inverse returns the **full original source**,
not just its H value or isomorphism class. A failed guard stops the word.
The empty conjunction is the empty word.

Let a_Z be arbitrary integer weights on sources, with Z={i:z_i=1}. The
response to the S-test is

```text
m_S = sum_(Z containing S) a_Z.
```

This is the Boolean zeta transform, with inverse

```text
a_Z = sum_(S containing Z) (-1)^|S\Z| m_S.
```

Ordered by subset size its matrix is triangular with diagonal1, equivalently
a tensor power of[[1,1],[0,1]]. Hence it is unimodular and works over every
commutative coefficient ring, in particular modulo every positive integer.
No division or precision sidecar is needed for this particular inverse.
The underlying Boolean transform is classical; the present application is
its realization by these actual complement/re-mark/reset words.

All2^r conjunction rows are linearly independent, so a linear carrier for
arbitrary integer source weights preserving all these responses needs rank
at least2^r. If total mass is fixed and supplied separately, the affine
dimension is2^r-1. This is a linear-information claim, not a prohibition on
encoding many coordinates into a single arbitrary-precision integer.

For probability laws, the inverse must additionally satisfy a_Z>=0 and
sum a_Z=1. Unimodularity gives the algebraic inverse; nonnegativity and the
specified source family are the realization guards.

## 5. Higher correlations with a constant Hamilton-mean discrepancy

Give equal probability to each even-weight z, respectively each odd-weight
z. For every proper subset of the r coordinates the two marginals are
identical: an omitted coordinate pairs even and odd completions.

Nevertheless, using H(T_z)=product_i(173-2z_i),

```text
E_even H = 172^r+1,       E_odd H = 172^r-1.
```

Indeed the unnormalized total sum is(173+171)^r=344^r and the alternating
sum is(173-171)^r=2^r. Each parity class has2^(r-1) sources. Their mean
difference is exactly2 at every r. These are genuine positive source laws,
but **means are not individual Hamilton counts**. In integer multiplicity
coordinates the unnormalized difference is2^r; normalization is over Q,
and is not silently performed modulo2.

There is also an exact **two-state aggregate recurrence**. If E_r,O_r are
the unnormalized H-sums on even and odd words, then

```text
[E_(r+1)]   [173 171] [E_r]       [E_0]   [1]
[O_(r+1)] = [171 173] [O_r],      [O_0] = [0].
```

The eigenvalues344 and2 give the displayed formulas. This compression is
closed under appending an independent A/B choice and observing parity sums.
It is **not** closed under the individual selected-block tests w_i. Thus a
short scalar recurrence and an exponentially large joint-response space
coexist on the same actual sources without contradiction.

This gives an arbitrary-order version of the previous three-edge repair
cube, with actual strong-atom composition and actual later commands. The
classical parity/marginal mechanism is credited; neither the proof of
spectrum completeness nor unbounded required *command depth* follows.
The present family requires unbounded state memory with depth-three bit
tests, while conjunction depth grows with the number of specified bits.

## 6. Consequence for the remaining spectrum proof

The all-size7/21 exclusions and the two-prime-lane reduction remain inherited
from THM-1370 / `h-spectrum-omits-7-21-all-n` and THM-4094 /
`hamiltonian-matching-deficit-and-two-prime-lane-completeness`.
To finish sufficiency one still needs uniform actual sources for odd prime
p!=7 and7p for odd prime p!=3, or an equivalent uniform source-entry theorem.

This packet supplies a correct Boolean inverse **of joint response laws**.
It also rules out identifying those laws by depth-two legality or by total H
and an unordered atom inventory. It does not construct missing prime lanes,
extend the finite coverage ceiling, or prove a Jacobian conjecture. The
useful Jacobian comparison is precise: here the full response map has
integer determinant1, while a scalar projection loses source coordinates.
