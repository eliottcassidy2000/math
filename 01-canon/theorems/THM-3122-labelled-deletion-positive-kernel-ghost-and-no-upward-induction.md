---
id: THM-3122
title: "Labelled-deletion derangement kernel ghost and no-upward-induction boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In every
  degree n>=2, factorial-normalized labelled deletion kills an explicit hook
  current whose Young-gap operator is exactly the Laplacian of the
  derangement Cayley graph of Sym_n.  That operator is positive semidefinite
  and nonzero.  Hence labelled deletion is not injective on coefficient
  currents and is not order-reflecting for carrier positivity: positivity
  after deletion cannot imply positivity before deletion.  This proves the need
  for a transverse pole/kernel selector; it is not an all-degree positivity
  or Gaussian-moment theorem.
source: root/gmc3000-audit-2026-08-02
audit: >
  Two independent hostile audits rederived the hook-kernel recurrence,
  averaged subgroup coefficient count, derangement-Laplacian identity,
  positivity factorization, generation/equality cases, and the repaired
  ambient-current/PSD-carrier no-induction scope, with no theorem defect.
  The exact companion checks the kernel identity, support coefficients,
  elementary transposition factorization, and every irreducible scalar
  through n=10.  Independent normal and optimized replays both match the
  stored output and declared LF hashes exactly.
depends_on:
  - THM-3119-factorial-normalized-labelled-deletion-and-young-carrier-order
related:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-2551-horizontal-transfer-transverse-projector-bicomplex-boundary
script: 04-computation/gmc_labelled_deletion_derangement_kernel_thm3122.py
output: 05-knowledge/results/gmc_labelled_deletion_derangement_kernel_thm3122.out
script_sha256: 835453ad3dab93cb0df7e8cf5fa823a3ba6b1268ce322417c365273d4eec1dc5
output_sha256: 67529918aeedbb1bc7346fbd19c5f6e2fd000ea0985b73034571fef9d27a13c6
hash_basis: LF-normalized bytes
---

# THM-3122 -- labelled-deletion derangement kernel ghost

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3119 proves that factorial conjugation turns algebraic block deletion
into genuine same-label stochastic deletion and that this deletion carries a
positive refinement boundary downward.  Downward preservation is not upward
reflection.  The missing kernel is already nonzero in degree two and has a
uniform positive description in every degree: it is the derangement Cayley
Laplacian.

## 1. The universal hook current

Use THM-3119's labelled-deletion operator

```text
A e_lambda=sum_(k>=1) k r_k(lambda)e_(lambda down k).          (1)
```

For `n>=2`, put

```text
lambda_j=(n-j,1^j),                         0<=j<=n-2,

kappa_n=sum_(j=0)^(n-2) (-1)^j binom(n,j)e_(lambda_j)
        +(-1)^(n-1)(n-1)e_(1^n).                            (2)
```

Then

```text
A kappa_n=0.                                                   (3)
```

Indeed, for `0<=j<=n-3`, the coefficient of
`e_(n-1-j,1^j)` is

```text
(n-j)(-1)^j binom(n,j)
 +(j+1)(-1)^(j+1)binom(n,j+1)=0.                              (4)
```

At the all-singleton target, the last two contributions are

```text
2(-1)^(n-2)binom(n,2)+n(-1)^(n-1)(n-1)=0.                    (5)
```

This also covers `n=2` directly.

## 2. The factorial carrier is the derangement Laplacian

Let

```text
w_lambda=product_i lambda_i!,
Kbar_lambda=w_lambda Lbar_lambda                              (6)
```

be the factorial Young carrier from THM-3119.  Since
`Lbar_(1^n)=0`, the operator carried by `(2)` is

```text
O_n=sum_(j=0)^(n-2) (-1)^j [n!/j!] Lbar_(n-j,1^j).            (7)
```

Let `Der_n` be the set of derangements in `Sym_n` and `D_n=|Der_n|`.
Then, as an exact central group-algebra identity,

```text
O_n=D_n identity-sum_(sigma in Der_n) sigma.                   (8)
```

To prove `(8)`, write `k=n-j`.  A set partition of type `(k,1^(n-k))`
is the choice of a `k`-set `B`, and its Young subgroup is `Sym(B)`.  If a
permutation `sigma` moves exactly `s` labels, its coefficient in the averaged
subgroup projection is

```text
binom(n-s,k-s)/[binom(n,k)k!]
 =(n-s)!/[n!(k-s)!]                                           (9)
```

when `k>=s`, and zero otherwise.  For a nonidentity permutation, `s>=2`.
Its coefficient in `(7)` is therefore

```text
-(n-s)! sum_(j=0)^(n-s) (-1)^j/[j!(n-s-j)!]
 =-(1-1)^(n-s).                                               (10)
```

This is zero unless `s=n`, when it is `-1`.  Each `Lbar` kills the trivial
representation, so `(7)` has augmentation zero.  Its identity coefficient
is consequently `D_n`, proving `(8)`.

## 3. Positivity and equality

The derangements are closed under inversion.  In every unitary
representation `rho`, `(8)` becomes

```text
rho(O_n)
 =1/2 sum_(sigma in Der_n)
      (I-rho(sigma))(I-rho(sigma))^* >=0.                     (11)
```

It is nonzero because the coefficient of every derangement in `(8)` is
`-1`.  Equality on a vector holds precisely when that vector is fixed by
every derangement.

For `n=3`, the derangements are the two three-cycles and generate `A_3`;
accordingly both the trivial and sign representations lie in the zero space,
while the standard representation is positive.  For `n=2` and every
`n>=4`, the derangements generate `Sym_n`.  For the only nontrivial case,
take

```text
alpha=(1 3 2 4 5 ... n),       beta=alpha^(-1)(1 2).          (12)
```

Both are derangements and `alpha beta=(1 2)`.  Conjugating gives every
transposition.  For `beta`, a fixed point `x` would say
`alpha(x)=(1 2)(x)`, which the displayed cycle order excludes.  Thus for
`n=2` and `n>=4`, `(11)` vanishes only on the
trivial isotypic subspace and is strictly positive on every nontrivial
irreducible representation.

Here conjugation is lawful because `Der_n` is conjugacy-stable, so the
subgroup it generates is normal; conjugating the displayed factorization
therefore puts every transposition in that subgroup.

## 4. The sharp no-upward-induction boundary

Equations `(3)` and `(11)` exhibit a nonzero current with a PSD Young-gap
carrier erased by labelled deletion.  In degree three this is already

```text
kappa_3=e_3-3e_(2,1)+2e_(1,1,1),
O_3=6(Lbar_3-Lbar_(2,1))>=0,             A kappa_3=0.          (13)
```

The inequality in `(13)` is strict on the standard representation and an
equality on the sign representation.

The opposite current `-kappa_n` has the same zero deletion but carries the
negative nonzero operator `-O_n`.  Hence

```text
A c has a PSD Young-gap carrier
   does not imply
c has a PSD Young-gap carrier.                                (14)
```

No argument using only the deleted current, even with exact cone
preservation and exact commutation with another deletion-compatible map, can
recover the missing sign.  An all-degree induction from THM-3115 therefore
needs a transverse observable that separates the `ker A` component.  The
active pole/Newton flag is one candidate; this theorem proves why some such
sidecar is logically necessary.

Here "positive carrier" means that the represented Young-gap operator is
positive semidefinite.  The theorem does not claim that `kappa_n` itself is
the boundary of a nonnegative Hasse one-chain; that smaller cone is only a
sufficient certificate for operator positivity.  The counterexample to
order reflection uses `-kappa_n`: its deletion is the zero Hasse boundary,
whereas its upper carrier is a negative nonzero operator.

THM-2551 has an analogous kernel-preservation boundary in the LRC transfer
bicomplex, but no map between the two problems is asserted.  The transferable
lesson is only that a flat commuting square cannot manufacture a predicate
stored in a positive kernel ghost.

## 5. Exact controls and scope

Run

```text
python 04-computation/gmc_labelled_deletion_derangement_kernel_thm3122.py
python -O 04-computation/gmc_labelled_deletion_derangement_kernel_thm3122.py
```

and compare with the stored output.  The companion uses exact integers and
rational numbers.  It checks `(3)` and the support coefficients in `(8)` for
`2<=n<=10`, checks `(12)`, and independently evaluates `(7)` on every
irreducible representation using the hook-content/Kostka scalar.

This theorem does not say that the live product-Gamma current has a negative
kernel component, nor does it refute the pole-flag route.  It supplies no
all-degree positivity, Gaussian-moment, NC2, LRC, JC, or DC conclusion.  It
is the sharp faithfulness and order-reflection boundary for the universal
labelled-deletion carrier.
