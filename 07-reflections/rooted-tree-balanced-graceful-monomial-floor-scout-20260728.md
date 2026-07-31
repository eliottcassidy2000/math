# The factor-three tree range is a triangular certificate, not the observed coefficient floor

**FINITE-EXACT SCOUT; NO IMPROVED UNIFORM THEOREM.**

THM-2765's rooted leading monomial gives every `n`-vertex tree a
distinct-edge-distance labeling in `0,...,3n-5`.  That monomial is chosen for
triangularity, not balance.  The full obstruction polynomial

```text
Phi_T(X)=product_(i<j)(X_j-X_i)
         product_(e<f)(d_e^2-d_f^2)
```

has total degree

```text
D_n=binom(n,2)+2binom(n-1,2)=(n-1)(3n-4)/2.             (1)
```

Therefore any full-degree monomial has maximum exponent at least

```text
R_n=ceil(D_n/n).                                         (2)
```

If `Phi_T` contains a nonzero monomial attaining `(2)`, coefficient-grid
interpolation gives a labeling in `0,...,R_n`, roughly half the current
`3n` range.

The exact bounded-support expansion now checks all unrooted tree shapes
through six vertices.  The minimum coefficient width is

```text
n=2: 1;
n=3: 2;
n=4: star 3, path 4 (the average floor is 3);
n=5: 5 for all three shapes, exactly R_5;
n=6: 6 for all six shapes, exactly R_6.                 (3)
```

Thus every shape at `n=5,6` reaches the information-theoretic exponent floor;
the four-vertex path is the unique small obstruction in this census.  This is
strong evidence for a threshold form such as

```text
for n>=5, every tree graceful-obstruction polynomial has
a nonzero full-degree monomial of width ceil(D_n/n),       (4)
```

but `(4)` is **OPEN**.

## Why the naive balanced term fails

The edge-square Vandermonde in formal variables `Y_e^2` contains every
permutation of exponents `0,2,...,2(n-2)`.  Reversing those exponents against
the vertex Vandermonde appears to balance the load near `2n`.  Substitution

```text
Y_e=X_child-X_parent
```

creates cancellation, however.  Already on the rooted three-vertex path,
the proposed `X_parent^2` term from the child edge cancels the corresponding
term from the parent edge.  So the edge-permutation term cannot be pulled
through the rooted incidence map coefficientwise.

The missing theorem is a signed coefficient-selection statement: find a
balanced exponent vector whose total contribution survives all incidence
cancellations.  Plausible carriers are an Alon--Tarsi orientation count, a
sign-reversing involution leaving one flow class, or a Hall-type load balance
on the factor-support hypergraph.  More brute-force expansion is the wrong
next move: the exact seven-vertex state space already grows past the cheap
certificate budget.

## Information ledger

```text
SOURCE:          THM-2765 obstruction polynomial.
REPRESENTATION:  coefficient support / Newton width.
PRESERVED:       vertex injection and unsigned edge-distance separation.
LOST:            a canonical coefficient sign after Y_e=X_child-X_parent.
CHEAP HOSTILE:   the four-vertex path misses the average exponent floor.
POSITIVE SIGNAL: every n=5,6 shape attains the floor exactly.
NEEDED SIDECAR:  cancellation-aware balanced-monomial selector.
```

Reproduce with

```bash
python 04-computation/rooted_tree_graceful_monomial_width_scout.py
python -O 04-computation/rooted_tree_graceful_monomial_width_scout.py
```

against `05-knowledge/results/rooted_tree_graceful_monomial_width_scout.out`.
The LF-normalized SHA256 hashes are
`7f03c04a51b266e72d87e54c0c004a6221d405269c18a84246060e1a0a9a82d2` and
`55aa45687584190978eb94d13a4358e6c3fac96c4d99ef551dafb8f9922da635`.
