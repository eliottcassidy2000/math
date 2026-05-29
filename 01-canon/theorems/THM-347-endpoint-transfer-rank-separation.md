---
id: THM-347
name: endpoint-transfer-rank-separation
status: PROVED
date: 2026-05-29
session: codex-2026-05-29-main-unification
depends_on:
  - THM-266
  - HYP-1772
  - HYP-1773
  - HYP-1774
---

# THM-347: Endpoint Transfer Separates Fiber Boundary from Rank

## Terminology

Let `X_n` be a fixed-path tiling cube and let

```text
E_n : X_n x {0,1}^{n-1} -> X_{n+1}
```

be endpoint extension. Let `q_n : X_n -> C_n` be any finite quotient of the
tilings and write

```text
F_n(c) = |q_n^{-1}(c)|
A_n(c,d) = #{(x,s) : q_n(x)=c and q_{n+1}(E_n(x,s))=d}.
```

Call

```text
partial_2(q_{n+1})(d) = F_{n+1}(d) mod 2
```

the **fiber boundary** of the child quotient. Call `q` **endpoint
parity-injective at n** if `A_n`, reduced over `F_2`, has full row rank
`|C_n|`.

## Statement

For every such quotient tower:

```text
sum_d A_n(c,d) = 2^(n-1) F_n(c)
sum_c A_n(c,d) = F_{n+1}(d)
```

and therefore, over `F_2` for `n >= 2`,

```text
sum_c A_n(c,-) = partial_2(q_{n+1}).
```

Thus endpoint transfer always exposes the odd child fibers, but it does not
automatically preserve all parent parity information. Full row rank is a
separate property of the quotient tower.

In particular:

1. For the unmerged tournament-isomorphism quotient, every child fiber is odd,
   because `F(C)=H(C)/|Aut(C)|`, Redei gives `H(C)` odd, and tournament
   automorphism groups have odd order. Hence the fiber boundary is every child
   class.
2. For the complement-merged tournament quotient, the fiber boundary is exactly
   the self-complementary child nodes: a fixed complement orbit has one odd
   class, while a non-fixed orbit is a pair of odd classes.
3. For the even-graph quotient, the child fiber is the labeled orbit size
   `(n+1)!/|Aut(G)|`. Hence the fiber boundary is exactly the set of graph
   orbits whose automorphism group contains the full 2-adic valuation of
   `(n+1)!`.

## Proof

The row sum counts endpoint extensions from a parent quotient class. Each
tiling in `q_n^{-1}(c)` has exactly `2^(n-1)` endpoint signatures, so

```text
sum_d A_n(c,d) = 2^(n-1) F_n(c).
```

The column sum counts child tilings by deleting the endpoint. Every child
tiling has a unique endpoint-deletion parent in `X_n`, so partitioning the
child fiber over parent quotient classes gives

```text
sum_c A_n(c,d) = F_{n+1}(d).
```

For `n >= 2`, every row has even total. Reducing over `F_2`, the coordinate
`d` of the sum of all rows is the column sum modulo 2, namely
`F_{n+1}(d) mod 2`. This proves the fiber-boundary identity.

The three examples are direct evaluations of `F_{n+1}(d) mod 2` in the three
quotients:

- tournament quotient: odd Hamiltonian-path fiber divided by odd automorphism
  group;
- complement quotient: one odd summand for fixed orbits, two odd summands for
  paired orbits;
- even-graph quotient: labeled orbit-size parity.

Finally, the boundary law alone cannot imply full row rank. Over `F_2`, the
matrix

```text
1 1 0 0
0 0 1 1
1 1 0 0
```

has even row sums and has a well-defined nonzero fiber boundary, but its row
rank is `2 < 3`. Therefore endpoint parity-injectivity is an additional
structural property, not a consequence of the transfer marginals.

## Consequence

The recent tournament/even-graph split should be read as a rank phenomenon, not
as a counting phenomenon. The double-counting transfer law survives every
quotient compatible with endpoint deletion. What changes is where odd fibers
live and whether the quotient retains enough parity information for the parent
rows to remain independent.

This isolates the current conjectural frontier:

```text
tournament and complement-merged quotients appear endpoint parity-injective;
the even-graph quotient is not.
```

The obstruction in the even-graph case is not endpoint recursion itself, but
the migration of parity from Hamiltonian-path fibers to automorphism 2-adics.
