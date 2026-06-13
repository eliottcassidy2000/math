# Goldbach, Polygonal Numbers, and Zeckendorf as Additive-Basis Regimes

The requested connection becomes clean if all four topics are put into one
object:

```text
representation hypergraph R_B:
  vertices = atoms in a basis B
  hyperedges over n = atom selections whose sum is n
```

Goldbach chooses `B = primes`.  Fermat's polygonal theorem chooses one
polynomial basis at a time.  Zeckendorf chooses Fibonacci atoms with a local
carry law.  Hardy-Littlewood is the asymptotic theory of the prime
representation hypergraph; Helfgott is the proof that the ternary prime layer
has no missing large or small odd targets.

## The Three Currencies

There are three ways an additive basis can pay for coverage.

```text
1. smoothing / many branches
2. bounded arity / enough summands
3. normal form / confluent carries
```

Goldbach is the smoothing extreme.  Primes up to `N` have density about
`1/log N`, so an even target has many candidate pairs, but the candidates are
arithmetically noisy.  Hardy-Littlewood predicts that the representation count
is a smooth main term times a singular series recording local prime-residue
compatibility.  Strong Goldbach remains open because the binary layer has too
little independent averaging for current tools.

Helfgott's weak Goldbach theorem moves one step up: odd targets use three
primes.  The parity obstruction is the same local parity story, but the third
variable gives a smoothing dimension.  That is the conceptual reason ternary is
proved while binary remains open.

Fermat polygonal numbers are the bounded-arity middle.  The `k`-gonal numbers
are much sparser than primes, roughly `sqrt(N)` atoms up to `N`, but the theorem
allows up to `k` summands.  The proof currency is not multiplicity on the scale
of primes; it is a finite arity budget large enough to erase residue gaps.

Zeckendorf is the normal-form extreme.  Fibonacci atoms up to `N` are only
`O(log N)`, far too sparse for smoothing.  Coverage comes from recurrence, and
uniqueness comes from local confluence:

```text
adjacent Fibonacci digits -> carry to the next layer
```

The no-adjacent rule is a path-independent-set condition, exactly the repo's
old `I(P_m,1)` Zeckendorf/fugacity lens.

## What The S501 Audit Saw

The script `additive_basis_goldbach_zeckendorf_s501.py` keeps the scope finite
and diagnostic.

Up to `5000`, binary Goldbach and ternary Goldbach had no missing targets in
the expected parity classes.  Ordered binary counts are already shaped by the
Hardy-Littlewood local factor: highly composite even targets have visibly more
representations than powers of two of comparable size.  Ordered ternary counts
are stable after normalization by `n^2/log(n)^3`, which is the smoothing signal.

For polygonals, the finite audit recovered the Fermat arity pattern for
`3`- through `8`-gonal numbers: the maximum needed number of summands was exactly
the polygon side count in each case.  For residues modulo `2..18`, triangular
triples, square quadruples, and pentagonal quintuples covered every residue.

For Zeckendorf, only `18` Fibonacci atoms cover all integers up to `5000`, with
maximum digit count `9`; every support satisfied the non-adjacent condition.
Residue coverage modulo `2..18` was complete, but the meaning is different from
the prime and polygonal cases: residues are covered by a deterministic automaton,
not by many independent representations.

## The Deeper Pattern

The same local-to-global question is being answered three ways:

```text
Hardy-Littlewood:
  local residue factors are multiplied into a density correction.

Fermat polygonal:
  local residue obstructions disappear after a bounded number of summands.

Zeckendorf:
  local carry conflicts are forbidden until a unique normal form remains.
```

This is the useful bridge back into the repo.  Tournament OCF already lives in
the same language:

```text
Zeckendorf: independent sets in a path at x=1.
OCF:        independent sets in a conflict graph at x=2.
Goldbach:   representation hypergraph over prime atoms.
Fermat:     representation hypergraph over polynomial atoms.
```

So the phrase "additive basis" should not be treated as one technique.  It is a
choice of proof economy.  If a problem supplies many local branches, use
smoothing.  If it supplies a finite packet of summands, look for bounded-arity
coverage.  If it supplies a recurrence with local conflicts, look for a
Zeckendorf-style normal form.

## Implication For LRC And Tournament Analysis

This helps classify proof attempts:

- LRC pressure-DAG methods are normal-form attempts: find a peel order or a
  carry order.
- Bridge-fiber dual certificates are bounded-arity attempts: finite rows should
  pay a finite invoice.
- Random perturbation or multi-gauge Tournament Analysis is a smoothing attempt:
  produce many relation-level witnesses and show that all local obstructions
  average away.

The current n=14/n=18 LRC work seems to sit between bounded arity and normal
form.  The exact local gate invoice is Fermat-like; the exported endpoint debt
looks Zeckendorf/Ostrowski-like.  A successful proof may need both currencies:
finite row weights to force export, then a carry normal form to show exported
debt cannot vanish.

## External Anchors

- Helfgott, "The ternary Goldbach conjecture", arXiv:1312.7748.
- Hardy and Littlewood, "Some problems of Partitio numerorum; III", Acta
  Mathematica 44 (1923), Springer DOI 10.1007/BF02403921.
- Zeckendorf, "Representation des nombres naturels par une somme de nombres de
  Fibonacci ou de nombres de Lucas", Bull. Soc. Roy. Sci. Liege 41 (1972).
