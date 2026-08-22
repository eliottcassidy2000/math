# U-spine prime toggles: clocks, cubical monodromy, and the next frontier

**Research synthesis, 2026-08-12.**  The proof source developed in this
session is THM-3346.  This note records the concept board, quotient losses,
and next experiments; it is not an additional proof source.  No claim below
closes LRC(14), supplies a global Berggren transducer, or creates a tournament
orientation.

## Outcome first

Three previously adjacent views now fit one typed diagram:

```text
U-spine position t
  z_t=(t+1)+it, C_t=4T_t+1, inradius=t
        |
        | compare two positions r,s
        v
two exact Gaussian channels
  gcd(C_r,s-r)          gcd(C_r,r+s+1)
  same local root       opposite local root
        |
        | fix an admissible grade N
        v
CRT root cube R_N = {t mod N : N divides C_t}
        |
        | quotient global conjugation and fill commuting squares
        v
K_N = Q_k^(2)/antipode
  H_1=Z/2, weighted systole=log N, large free H_2
        |
        | apply source-dependent ancestry paths
        v
Berggren tree: every square and the remaining Z/2 class die
```

The arithmetic gain is an exact splitter: every common prime power of two
U-spine hypotenuses chooses exactly one of the difference and reflected-sum
channels.  The topological gain is just as sharp: after all commuting prime
squares are filled, the allocation fibre retains one conjugation bit and no
other one-dimensional current.  The Pell gain is an operation-level compiler:
one square-triangular row gives two fixed-hypotenuse parents whose two
primitive products land at the adjacent square-hypotenuse U-spine depths.

The key frontier lesson is negative but useful.  Prime toggles create a rich
labelled two-dimensional carrier; they do not create a source-independent
ancestry action.  A consumer must keep height, Gaussian section, owner, phase,
and labelled columns.

## Anchor / Niche / Wildcard portfolio

| Lane | Object | Exact gain | Honest boundary |
|---|---|---|---|
| Anchor -- LRC(14) | determinant/Kelvin gate with labelled 13-column planes | a prime-local clock and content metric that can be attached to a coefficient direction | no lawful map from a residual LRC deck into a fixed root cube; no row closes |
| Niche -- primitive triples | consecutive spinors `z_t` and Gaussian multiplication | two-channel content factorization, primitive radius invoices, and a lossless fixed-grade root atlas | primitive/unordered triples still lose section and leg gauge |
| Wildcard -- topology | prime-allocation cube with commuting squares | antipodal quotient has `pi_1=Z/2`, `H_1=Z/2`, and explicitly growing `H_2` | ancestry tree kills the class; real/stable one-homology sees nothing |

The niche overtook the anchor because it produced a theorem and a new
obstruction.  Its best service to LRC is presently diagnostic: it says which
coordinates a future operation must preserve before a determinant certificate
can be compared.

## Inheritance pass

| Required item | Controlled source |
|---|---|
| closest proved mechanism | THM-3336: primitive Gaussian content cocycle and folded fixed-grade weights |
| native U-spine operation | THM-3334: the complete consecutive-parameter Berggren ray and its Gaussian collision fibres |
| canonical hostile | THM-3345: one prime direction has source-dependent Berggren paths and costs |
| corrected near miss | modular roots record only the selected `N`-primary content, not the full gcd of larger U-spine norms |
| least-used sidecar | the commuting two-cells of the raw factor-allocation cube |
| transverse compiler | THM-3341: consecutive square-hypotenuse depths from a square-triangular Pell row |

The minimal modular/full-content hostile is compact enough to remember:
`N=5`, `r=6`, `s=23`.  Both positions are roots modulo five.  Their modular
channels are `(1,5)`, but their full Gaussian contents are `(17,5)`.  The
first failed implication is therefore

```text
same modular root grade  !=  complete shared-norm factorization.
```

The strongest survivor is exact: the grade-`N` prime powers still split
between the two channels with no overlap.

## Live concept board

| Object | Representation | Invariant | Operation | Lost coordinate / next test |
|---|---|---|---|---|
| U-spine node | `z_t=(t+1)+it`, `Phi(z_t)` | `C_t=4T_t+1`, inradius `t` | Berggren `U` | a norm alone forgets Gaussian angle |
| pair of nodes | `(r,s)` | `gcd(C_r,C_s)=d_-d_+` | product / conjugate product | unordered pair swaps same/opposite channels |
| modular root | `x=2t+1`, `x^2=-1 mod N` | one sign at each split prime power | CRT sign toggle | residue forgets integral height |
| Gaussian parent | `gcd_Z[i](N,z_t)` | norm `N` up to units | source-dependent factor replacement | parent quotient forgets conjugation |
| root cube | `F_2^k` | prime labels and exponents | commuting coordinate toggles | graph alone forgets which commutators are filled |
| cubical quotient | `Q_k^(2)/antipode` | deck class in `Z/2`, free `H_2` | fill prime squares | ancestry functor kills all one-dimensional classes |
| Pell bridge | `(M_-,M_+,H)` | `M_-M_+=H`, folded weight `{M_-,M_+}` | two Gaussian products | no claim the two parents exhaust the grade |
| LRC deck | labelled primitive columns | determinant gate, saturation, owner/phase | column operation | no root-cube incidence map yet |

## Why the two-channel splitter is more than another gcd identity

The U-spine polynomial has discriminant `-4`:

```text
C_t=2t^2+2t+1,        2C_t=(2t+1)^2+1.
```

Modulo a shared prime power, `C_s=0` has exactly the two simple branches

```text
s=r                     or                    s=-r-1.
```

Those branches are already the two Gaussian operations.  The first makes
`z_r conjugate(z_s)` nonprimitive; the second makes `z_r z_s` nonprimitive.
Thus the gcd law is a local-orientation detector.  It converts the static
sum-of-two-squares factorization into a clock on the inradius parameter.

For `p=1 mod 4`, the clock has two Hensel branches at every precision:

```text
density(v_p(C_t)>=e)=2/p^e,
density(v_p(C_t)=e)=2(p-1)/p^(e+1).
```

There is no clock for `p=3 mod 4`.  Distinct prime clocks are independent by
CRT.  This is exactly the sort of sparse, prime-local address that could be
useful in a difficult frontier, but only after it is attached to the native
consumer rather than treated as a universal time variable.

## The topology says what survives controlled forgetting

The raw allocation cube is contractible, and its two-skeleton is simply
connected.  Global conjugation is the antipode.  For at least three prime
directions the antipode is free on every cell of dimension at most two, so
the quotient is a genuine double cover with one deck class:

```text
pi_1(K_N)=H_1(K_N;Z)=Z/2,
H^1(K_N;F_2)=F_2,
H^1(K_N;Z)=0.
```

The distinction between the graph and filled square complex is load-bearing.
The graph has many cycles.  Filling every commuting square kills all of their
free one-dimensional part.  What remains is the obstruction to choosing one
member of every conjugate pair coherently.  It is torsion, so no stable real
current detects it.

At three directions the quotient one-skeleton is `K4` and its three square
boundaries are precisely the three affine-V4 `C4` channels.  Filled together
they form a cubical `RP^2`, not three unrelated matchings.  In higher rank the
complex is not a surface; its free second homology has rank

```text
2^(k-4)(k^2-5k+8)-1.
```

This suggests a different use of the prime cube.  Do not search for a new
one-current after the tree has killed it.  Search for a labelled two-chain,
curvature, or compatibility obstruction evaluated by an external consumer.

## Pell composition is a factor-to-height compiler

For a square-triangular row `T_n=R^2`, define

```text
M_-=2n+1-2R,       M_+=2n+1+2R,       H=4R^2+1.
```

The selector spinor `u=(n+1)+in` and companion `v=2R+i` both have norm `H`.
Their two products have contents `M_+` and `M_-`, and primitive reduction
lands at the adjacent square-depth spinors.  In the fixed-grade parent cube,

```text
folded weight([u],[v])={M_-,M_+},
weighted distance=log M_-.
```

So the same complementary factors have three meanings:

1. Gaussian cancellation contents;
2. geodesic invoices inside the parent fibre; and
3. square roots of the two output hypotenuse grades.

The construction also reaches unbounded prime-toggle rank.  THM-3341's
strong-divisibility argument chooses Pell roots with arbitrarily many
pairwise-coprime split-prime divisors; the bridge adds the coprime preceding
root.  Hence the exact Pell grades support arbitrarily high-dimensional root
cubes, persistent conjugation torsion, and unbounded free `H_2`.

## Frontier leverage and cheapest decisive tests

### 1. LRC root-clock incidence test

Source: a residual THM-2056 coefficient deck.  Target: one `R_N` root cube.
Map each primitive column `c_i=(m_i,n_i)` to the root of a chosen split norm
only when `m_i+n_i h` has a controlled content and the normalized image deck
remains saturated.  Preserved predicate: exact determinant ratios with the
full content vector.  Destroyed information: height, owner, phase, and column
labels if projected too early.  Needed sidecar: `(d_i)`, Kelvin norm, owner,
phase, and the inverse row-basis map.  Cheapest test: run this only on the
finite residual cones and require a lawful inverse map on the physical speed
row.  Anything weaker is another sufficient-gate statistic, not a closure.

### 2. Two-chain rather than one-current probe

Source: commuting prime-toggle squares in `K_N`.  Target: a consumer-specific
defect or curvature group.  Preserved predicate: square commutation.  Lost
data: Gaussian section and ancestry source.  Needed sidecar: orient each
square using a labelled owner or phase.  Cheapest test: at `N=1105`, evaluate
the three filled `C4` channels against the actual six endpoint ancestry paths.
The flat tree functor must give zero; a useful consumer must distinguish at
least one square without violating that control.

### 3. Pell bridge compiler test on frontier predicates

Source: `(u,v,M_-,M_+)` from a square-triangular row.  Target: the two adjacent
square-depth triples.  Preserved predicates: exact norm, content, primitive
output, and U-spine depth.  Destroyed information: the rest of the parent
fibre and a fixed Berggren word.  Cheapest test: compare any proposed
predicate at the two parents and two outputs for the first rows
`(M_-,M_+)=(5,29),(29,169),(169,985)`.  A genuine compiler must explain both
channels and the degenerate `(1,5)` boundary.

### 4. Prime-labelled spectral sidecar

The folded cube graph has a known character spectrum, but its unweighted
version loses which coordinate is which.  Weight coordinate `j` by
`log(p_j^e_j)` and study the weighted Laplacian together with the unique
mod-two cover class.  The decisive hostile is complement symmetry: any
claimed scalar reconstruction must distinguish a subset from its complement
only after a conjugation section is chosen.

## Stopping boundary

This session does not produce a new LRC certificate, Berggren automorphism,
tournament, or Jacobian flux.  It does identify a sharper architecture:

```text
arithmetic clock + cubical two-cells + consumer owner/height/phase.
```

The clock alone is too small, the cube graph alone has fake free cycles, and
ancestry alone kills the remaining class.  The next serious move is to attach
the labelled two-cells to a native consumer and test a predicate that lives in
dimension two.  That is the frontier opened here.
