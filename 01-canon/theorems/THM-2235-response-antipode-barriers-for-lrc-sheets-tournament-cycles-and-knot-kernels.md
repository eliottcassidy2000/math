---
id: THM-2235
title: "Response-antipode barriers for LRC sheets, tournament cycles, and knot kernels"
status: >
  PROVED. A measure-preserving fixed-XOR response can flip only balanced
  Boolean coordinates, and every such fixed-XOR response of an odd-order
  group is trivial.
  Consequently the raw four-checkpoint three-comb word has no nonzero
  measure-preserving XOR response, while 13-power root-sheet groups cannot
  access its top Walsh 2-torsion or negate a nonzero characteristic-zero
  lift vector. In the pinned tournament response, a permutation sigma has
  exterior cost at most N M(|C|-o(sigma)), where o(sigma) counts odd cycles;
  for a nonempty exterior ledger, exact antipodes are precisely even-cycle
  permutations supported on their alternating affine code. In the knot
  continuation-kernel semigroup, the only min-plus unit is the unknot. These
  are exact no-go/structure results,
  not closures of the remaining scalar LRC residue, tournament core transport,
  knot catalysis, or LRC(14).
source: codex-2026-07-25-response-antipode-barriers
depends_on:
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
related:
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
  - THM-2224-transfer-owner-word-temporal-union-bound
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-2226-three-checkpoint-bellman-sieve-and-eight-profile-residue
  - THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants
---

# THM-2235 -- exact barriers to response antipodes

THM-2225 succeeds because its cyclic response image contains the required
antipode. The same question has exact negative answers on three current
frontiers. The common lesson is not that response methods fail, but that
the surviving response is noninvertible or carries a cycle defect.

## 1. Marginal obstruction to a fixed-XOR response

Let a group `Gamma` act by measure-preserving bijections on a probability
space `(X,mu)`. Let

```text
Y=(Y_1,...,Y_p):X -> F_2^p
```

be measurable. Suppose that for every `g in Gamma` there is a constant
`alpha(g) in F_2^p` such that

```text
Y(gx)=Y(x)+alpha(g)                         a.e. x.   (1)
```

Then

```text
alpha(gh)=alpha(g)+alpha(h).                         (2)
```

Indeed, apply (1) successively to `h` and `g`; equality of two constant
translations of the same Boolean vector forces equality of the constants.
Thus `alpha` is a homomorphism.

If `alpha_i(g)=1`, measure preservation gives

```text
mu(Y_i=1)
 =mu(Y_i composed g=1)
 =mu(Y_i=0)
 =1-mu(Y_i=1).                                      (3)
```

Hence only coordinates with marginal `1/2` can occur in the support of
`im(alpha)`. In particular,

```text
mu(Y_i=1)!=1/2 for every i       implies alpha=0.    (4)
```

There is a second, independent obstruction. If `Gamma` is finite of odd
order, then its homomorphic image in the `2`-group `F_2^p` has order
dividing both an odd number and a power of two. Therefore

```text
|Gamma| odd                         implies alpha=0. (5)
```

The same conclusion holds elementwise for a group in which every element
has finite odd order.

## 2. The raw LRC checkpoint word has no response antipode

Let `d_1,d_2,d_3` be nonzero integers; in THM-2222's application they are
positive. Put

```text
D_a={x in R/Z:||ax||<1/14},
U=D_(d_1) union D_(d_2) union D_(d_3),
Sx=169x,

Y_r(x)=1_U(S^r x),                         0<=r<=3.  (6)
```

Haar invariance under every nonzero integer circle endomorphism gives

```text
P(Y_r=1)=measure(U)<=3/7<1/2                 for all r. (7)
```

Apply (4). For **every** measure-preserving group action on the circle
satisfying a fixed-XOR law (1) for the word `(Y_0,...,Y_3)`,

```text
alpha=0.                                               (8)
```

Thus no such action translates the raw checkpoint word by `1111`, by an
even-parity nonzero word, or by any other nonzero Boolean response. This is
stronger than saying that the visible `13`-sheet action lacks an involution.

For the actual root/deck groups, (5) supplies an algebraic version. Every
group `(F_13)^k` has odd order, so it has no nontrivial homomorphism to
`F_2^4`. More generally, if `g` has finite odd order and acts linearly on a
characteristic-zero vector space, then

```text
g v=-v                     implies                  v=0. (9)
```

If the order is `m`, then `v=g^m v=(-1)^m v=-v`.
Consequently a nonzero centered lift vector in THM-2218 cannot be sent to
its antipode by a `13`-power sheet translation.

THM-2222's identity

```text
L R=-R
```

does not contradict (8)--(9). The operator `L` sums over thirteen
preimages; it is not the Koopman operator of an invertible
measure-preserving action. After normalization its eigenvalue is `-1/13`,
not a sign character.

### The lawful replacement is a directed top correlation

Put

```text
f=1-2 1_U,
W_4=integral product_(r=0)^3 f(S^r x) dx.             (10)
```

If

```text
C_<4
 =sum_(A proper subset {0,1,2,3})
    (-1)^|A| integral product_(r in A)f(S^r x) dx,    (11)
```

then Boolean Fourier inversion gives the exact identity

```text
measure(intersection_(r=0)^3 S^(-r)U)
 =2^(-4)(C_<4+W_4).                                  (12)
```

The proper-character packet in (11) is the information retained by moments
through degree three; `W_4` is the missing top character isolated by
THM-2237. Equations (8)--(9) show that it cannot be manufactured by a group
antipode. It must be bounded as a directed correlation, or bypassed by a
noninvertible clause/owner transfer such as THM-2224 and THM-2226.

## 3. Tournament odd-cycle antipode deficit

Use THM-2221's pinned equal-block setting. Let `C` be the core index set,
let `mu(w)` be the exterior incidence-word multiplicities, and put

```text
M=sum_(w in {0,1}^C) mu(w).
```

For `sigma in Sym(C)`, let `o(sigma)` be the number of odd cycles of
`sigma`, with fixed points included. THM-2221 gives the whole-block exterior
response

```text
E_mu(sigma)
 =N sum_w mu(w)d_H(w,w composed sigma).              (13)
```

On a permutation cycle of length `ell`, the cyclic binary transition count
is even. Its maximum is `ell` when `ell` is even and `ell-1` when `ell` is
odd. Summing over cycles and then over words proves

```text
E_mu(sigma)<=N M(|C|-o(sigma)).                      (14)
```

Equality holds if and only if every word with positive multiplicity:

- alternates on every even cycle; and
- has exactly one equal cyclic adjacency on every odd cycle.

Define the alternating affine code

```text
A_sigma={w:w composed sigma=1-w}.                    (15)
```

Assume `M>0` for the antipode classification, so the supported ledger is
nonempty. (When `M=0`, the exterior response is vacuous.)
It is nonempty if and only if every cycle of `sigma` is even. Therefore the
exact antipode condition on the entire supported context ledger is

```text
w composed sigma=1-w for every w in supp(mu)

iff

every sigma-cycle is even and supp(mu) subset A_sigma. (16)
```

In that case `E_mu(sigma)=NM|C|`. A Gram variance detects that supported
words move; it does not test the affine-code membership in (16). The next
finite tournament computation should therefore solve these alternating
constraints on the low-`G(X)` transport permutations before using a
variance surrogate.

## 4. Knot continuation kernels have no nontrivial units

For knots, THM-2220 defines the exact min-plus continuation kernel

```text
P_K(A,B)=d_G(K#A,B)
```

and proves

```text
P_K tensor P_J=P_(K#J).                              (17)
```

Suppose `P_J` were an inverse of `P_K`, so that

```text
P_K tensor P_J=P_U,                                  (18)
```

where `U` is the unknot. Evaluate (17)--(18) at `(U,U)`:

```text
u(K#J)=P_(K#J)(U,U)=P_U(U,U)=0.                     (19)
```

Hence `K#J` is the unknot. Schubert prime decomposition makes the knot
monoid conical, so

```text
K=J=U.                                               (20)
```

Thus the continuation-kernel semigroup has no nontrivial units. In
particular, the mirror of a nontrivial knot is not an inverse kernel. The
faithful algebra is noninvertible min-plus composition; this result does
not forbid catalytic contraction or geodesic bypass.

## 5. Failure boundaries

1. Equation (8) does not rule out nonlinear or state-dependent word
   permutations, orbit packing with stabilizers, or Bellman
   supermartingales.
2. Equation (12) is an exact reformulation, not a bound on `W_4`.
3. Equation (14) concerns the pinned whole-block exterior term. It neither
   computes THM-2221's core kernel `G(X)` nor controls unpinned exchange.
4. Equation (20) rules out inverse kernels, not nonadditivity of unknotting
   number or positive catalysis.

The common survivor is a directed, defect-bearing response rather than an
antipodal group symmetry. QED.
