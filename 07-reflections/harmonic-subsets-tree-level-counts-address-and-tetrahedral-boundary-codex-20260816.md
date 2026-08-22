# Harmonic subsets: tree counts decide convergence, addresses decide drift

**Status: ELEMENTARY PROOF + VERIFIED-EXACT BOUNDARY SIDECAR; NOT A CANON
THEOREM AND NOT INDEPENDENTLY AUDITED.**

## 1. The precise form of “every subset is a harmonic subset”

Every subset `S` of the positive integers selects the subseries

```text
sum_(m in S) 1/m.                                        (1)
```

Binary shortlex indexing gives a second, useful description.  The words of
length `n` are exactly the integer block

```text
I_n=[2^n,2^(n+1)-1].                                    (2)
```

Thus an arbitrary subset of the natural numbers is exactly an arbitrary
colouring of all nodes of the full binary tree.  If

```text
a_n=|S intersect I_n|,                                  (3)
```

then THM-3499 proves the sharp coarse principle

```text
sum_(m in S)1/m converges
  iff sum_(n>=0) a_n/2^n converges.                      (4)
```

So the level-count/Kraft projection decides convergence.  It does not decide
the logarithmic coefficient, because reciprocal weights vary across a level.

## 2. Equal counts, two different regular coefficients

For every `n>=1`, split (2) into its first and last halves:

```text
E_n=[2^n,3*2^(n-1)-1],
L_n=[3*2^(n-1),2^(n+1)-1].                              (5)
```

Both contain exactly `2^(n-1)` integers, hence both have Kraft mass `1/2`
on every positive level.  Their harmonic masses are

```text
e_n=sum_(m in E_n)1/m,
l_n=sum_(m in L_n)1/m.                                  (6)
```

Riemann-sum or harmonic-number asymptotics give

```text
e_n -> integral_1^(3/2) dx/x = log(3/2),
l_n -> integral_(3/2)^2 dx/x = log(4/3).                (7)
```

Consequently the two regular languages “first digit zero” and “first digit
one” have the same count on every level but logarithmic densities

```text
delta_E=log(3/2)/log 2,
delta_L=log(4/3)/log 2.                                 (8)
```

They sum to one, but are unequal.  This is THM-3499's ordered-address
phenomenon in its smallest possible form.  It also explains why the
unweighted basin probability `1/2` is insufficient.

## 3. Equal counts and no logarithmic density

The boundary can be made sharper.  Let stage `k` have length

```text
r_k=2^(2^k)                                             (9)
```

in the *level index*, and select `E_n` throughout odd stages and `L_n`
throughout even stages.  Call the resulting subset `S_*`.  Every positive
level still has

```text
|S_* intersect I_n|=2^(n-1).                            (10)
```

In particular its Kraft series diverges with the same term `1/2` as both
regular controls.  However,

```text
(r_1+...+r_(k-1))/r_k -> 0.                             (11)
```

At the end of stage `k`, the current stage therefore dominates all earlier
levels.  Since its first level also tends to infinity, (7) is uniform on all
but a vanishing initial fraction of that stage.  If `N_k` is the last integer
through the last level of stage `k`, then

```text
log N_k=(r_1+...+r_k)log 2+O(1).                        (12)
```

Equations (7), (11), and (12) give the two subsequential limits

```text
lim_(k odd)  [sum_(m<=N_k,m in S_*)1/m]/log N_k
  =log(3/2)/log 2,

lim_(k even) [sum_(m<=N_k,m in S_*)1/m]/log N_k
  =log(4/3)/log 2.                                      (13)
```

Thus `S_*` has no logarithmic density despite having the same exact level
counts as two regular languages that do.  The hierarchy is now precise:

```text
level counts       decide harmonic convergence,
ordered addresses decide logarithmic drift,
finite-state law   forces the drift limit to exist.     (14)
```

## 4. Why four and six are one tetrahedral object

The four binary words at depth two are

```text
00,01,10,11.                                            (15)
```

Use them as the vertices of `K4`.  There are six unordered pairwise
comparisons, exactly the six edges of `K4`; orienting those edges gives a
tournament **on four vertices**.  This resolves the apparent sizes `4` and
`6`:

```text
4 = number of competitors,
6 = number of pairwise comparison channels.             (16)
```

They are not tournaments of two different orders.  A tournament on six
vertices would have `binom(6,2)=15` comparison edges.

The three nonzero linear forms on `F_2^2` split (15) into three `2+2`
partitions.  The two within-fibre edges of each split form one of the three
perfect matchings of `K4`.  In particular, the first-bit form gives

```text
{00,01} | {10,11},                                      (17)
```

which is exactly the early/late split in (5).  Cardinality and the abstract
matching forget which fibre comes first; the harmonic address restores that
lost orientation.  This is the same representation debt isolated by the
tetrahedral Haar/XOR atlas.

If one promotes the six `K4` edges to vertices, their intrinsic adjacency
graph is the octahedral line graph `L(K4)`: disjoint complementary edges are
opposites.  It is not a complete graph and therefore supplies no canonical
order-six tournament either.

## 5. Connection contract

| field | exact content |
|---|---|
| source | an arbitrary subset of `N`, equivalently a binary-tree node colouring |
| target | its harmonic subseries and level-count/Kraft series |
| map | binary shortlex index |
| preserved by counts | convergence or divergence |
| destroyed by counts | within-level order and logarithmic coefficient |
| required sidecar | ordered binary address, or a finite automaton that generates it |
| cheapest hostile | first-half and last-half languages have equal counts but coefficients (8) |
| stronger hostile | the staged language has the same counts and no coefficient at all |

This bridge does not recover Berggren ancestry, create an LRC current, or
transport a Jacobian flux.  It does explain why a subset of the harmonic
series is naturally tree-valued before it is scalar-valued, and why a
tournament/matching summary must retain an address sidecar whenever weights
depend on order.

## 6. Reproduction

```text
python3 04-computation/harmonic_shortlex_same_counts_no_log_density_probe_20260816.py
python3 -O 04-computation/harmonic_shortlex_same_counts_no_log_density_probe_20260816.py
```

The exact companion checks the first twelve harmonic blocks as rational
numbers, the THM-3499 endpoint squeeze, constant Kraft mass, the dominant
stage ratios, the `4/6/15` carrier census, and the three Walsh matchings.  Its
semantic ledger is

```text
e135ce689fe04866568fd98961a45f3ebc9fabca7c03c9bf92c50be762f2cb45.
```

Script/output LF SHA-256 are respectively

```text
e311a4e65491d414a9cb5a6a21d9c3113b811c74047a66f1f3e32ff698e3d2f2
e3e075bcaa94871a815815d6c59f4cbe88d46618f96b119fb79bc68dedfc46bd.
```
