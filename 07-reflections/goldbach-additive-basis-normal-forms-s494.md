# Goldbach, polygonal numbers, and Zeckendorf as one additive machine (S494)

The useful move was to stop treating these as five separate famous theorems.
They all ask how a target integer is covered by additive atoms.

```text
target n
  atom set A
  legal hyperedges a_1+...+a_r=n
  local restrictions
  normal form or abundance
```

Goldbach is the sharp sparse version: `A = primes`, `r = 2`, even targets.  The
strong conjecture says every even target has at least one prime edge.
Hardy-Littlewood says much more: the number of edges should be archimedean
area times a p-adic singular product.  In repo language, this is a gap/debt
product: continuous volume multiplied by local residue debt.

Helfgott's theorem is the smoothing version.  Add one more prime.  The binary
edge problem becomes a ternary surface.  This is not just "three is easier
than two"; it is a dimension jump.  The circle method can average across the
extra degree of freedom and push the local residue structure into positive
global mass.

Fermat polygonal number theorem sits on the bounded-cover side.  For `k`-gonal
atoms, at most `k` summands cover every positive integer.  This is Waring-like:
structured polynomial atoms plus enough summands erase the obstruction.  The
S494 dynamic program only checks small `n`, but it mirrors the theorem cleanly:
triangulars hit max `3`, squares max `4`, pentagonals max `5`, and so on up to
the checked range.

Zeckendorf is the opposite pole.  It does not say "many representations"; it
says "one canonical representation."  Its obstruction is local and graph-like:
no adjacent Fibonacci indices.  That is exactly a path-independence rule, which
the repo already relates to OCF through independence polynomials.

The fresh idea is to put Zeckendorf under Goldbach, not beside it.  Given a
Goldbach pair `p+q=n`, overlay the Zeckendorf supports of `p` and `q`, then
compare the raw digit multiset to the canonical Zeckendorf support of `n`.
This gives a carry debt.  Some Goldbach pairs have debt zero:

```text
50   = 13 + 37
100  = 11 + 89
5000 = 673 + 4327
```

In those cases the prime edge is already compatible with the Fibonacci normal
form of the target.  Balanced pairs can have much higher carry debt.  So
"balanced", "prime-rich", "locally favored", and "normal-form compatible" are
different axes.

That feels like the real pattern:

```text
Goldbach/Hardy-Littlewood: abundance with p-adic local weights.
Helfgott: add dimension until positivity is provable.
Fermat polygonal: bounded arity gives universal cover.
Zeckendorf: local carry confluence gives canonical uniqueness.
```

For future repo work, the reusable object is an additive-basis profile:

```text
(atom set, arity, representation count, local product, carry debt, normal form)
```

This profile should transfer back into LRC and Tournament Analysis.  In LRC,
possible endpoint repairs are hyperedges; endpoint debt is the local product;
and a Zeckendorf/Ostrowski-style carry score may measure whether a repair is
compatible with the canonical boundary normal form.
