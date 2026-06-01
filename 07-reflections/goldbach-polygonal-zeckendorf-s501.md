# Goldbach, Polygonal Numbers, and Zeckendorf S501

The clean way to connect these ideas is to stop asking first about the
particular atoms and ask what kind of representation object is returned.

All four systems can be written as summand hypergraphs:

```text
target N
atoms A
hyperedges = finite multisets of atoms summing to N
```

The difference is the entropy and the certificate shape.

## Goldbach and Hardy-Littlewood

Strong Goldbach asks whether every even `N>2` has at least one edge in the
prime-pair graph:

```text
p + q = N.
```

Hardy-Littlewood says this graph should not merely be nonempty.  It should have
many edges, with a count controlled by:

```text
archimedean density  *  local singular series.
```

The S501 script uses the standard binary Goldbach heuristic as a diagnostic.  At
small scale the constants wobble, but the moral is visible: prime-pair rows are
not canonical.  They are redundant.  A proof wants a lower bound.

This is a very different feeling from Zeckendorf.  Goldbach is not looking for
the right carry rule.  It is looking for enough independent chances that at
least one survives the sieve.

## Helfgott

Helfgott's theorem proves the ternary Goldbach conjecture: every odd integer
greater than `5` is a sum of three primes.

In the hypergraph language this is the move:

```text
binary prime edges -> ternary prime hyperedges.
```

Adding one summand raises the dimension of the search space.  The S501 finite
counts show this immediately:

```text
101  has    38 unordered prime triples
1001 has  1095 unordered prime triples
1999 has  3105 unordered prime triples
```

So the third prime is not cosmetic.  It turns a thin pair-intersection problem
into a thick hypergraph-coverage problem.  Analytically, this is exactly why the
circle method has more room in the ternary problem than in the binary one.

## Fermat Polygonal

Fermat's polygonal theorem says every positive integer is a sum of at most `s`
`s`-gonal numbers.

This is also a representation-hypergraph theorem, but with a bounded depth:

```text
triangular atoms -> at most 3
square atoms     -> at most 4
pentagonal atoms -> at most 5
...
```

S501 verified the exact depth pattern through `600` for sides `3..8`.

This lives between Goldbach and Zeckendorf.  It is not usually unique, but it is
not probabilistic in the prime-distribution sense either.  It is a bounded
residue absorber: the structured quadratic atom set has enough shapes that `s`
copies erase all local obstructions.

## Zeckendorf

Zeckendorf is the opposite pole from Hardy-Littlewood.  It says every positive
integer has exactly one representation as a sum of nonconsecutive Fibonacci
numbers.

The extra datum is not just the atom set.  It is the conflict graph:

```text
F_i conflicts with F_{i+1}.
```

A valid representation is an independent set in this infinite path.  The no
adjacent rule is a local carry law:

```text
F_i + F_{i+1} -> F_{i+2}.
```

That carry law kills entropy.  Goldbach wants many representations and then
uses analysis to prove at least one exists.  Zeckendorf has many naive Fibonacci
representations, then uses confluence to collapse them to one.

## The Ladder

The proposed structure, recorded as HYP-1962, is:

```text
Hardy-Littlewood / Goldbach
  abundant representations, local singular series, analytic lower bound

Helfgott
  same prime atoms, one higher hypergraph dimension, enough room to prove coverage

Fermat polygonal
  bounded-depth structured basis, local residues absorbed by a fixed number of atoms

Zeckendorf
  canonical path-independent support, local carry graph confluent, entropy = 0
```

So the hidden axis is:

```text
redundancy -> bounded basis -> canonical peel
```

This belongs with the repo's natural-operation graph work.  Addition as an
unlabeled graph collapses to a transitive order, so the information is in the
fibers:

```text
all pairs summing to N
all triples summing to N
all polygonal packets summing to N
all carry-reducible Fibonacci packets summing to N
```

The operation fiber is the real object.

## New Questions

Can every representation theorem be placed by three coordinates?

```text
entropy:        log(number of representations)
local rank:     number/shape of residue obstructions
carry width:    graph width of local rewrite conflicts
```

Goldbach: high entropy predicted, high analytic difficulty.

Helfgott: higher entropy after adding a variable, proof possible.

Fermat polygonal: bounded depth, moderate entropy, residue absorber.

Zeckendorf: zero entropy after canonicalization, path-width one carry graph.

The next computational target is a small representation-hypergraph TDA: for each
target `N`, build the hyperedges for prime pairs, prime triples, triangular
triples, square quadruples, and Fibonacci representations, then measure connected
components, private atoms, local residues, and peelability.

## Sources

- Harald Helfgott, [The ternary Goldbach conjecture is true](https://arxiv.org/abs/1312.7748).
- G. H. Hardy and J. E. Littlewood, [Some Problems of "Partitio Numerorum"; III](https://archive.ymsc.tsinghua.edu.cn/pacm_download/117/5327-11511_2006_Article_BF02403921.pdf).
- Melvyn B. Nathanson, [A Short Proof of Cauchy's Polygonal Number Theorem](https://www.theoryofnumbers.com/melnathanson/pdfs/nath1987-55.pdf).
- [Zeckendorf representation](https://mathworld.wolfram.com/ZeckendorfRepresentation.html).
