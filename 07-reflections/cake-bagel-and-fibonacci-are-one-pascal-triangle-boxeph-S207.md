# Cake, bagel, Moser, and Fibonacci as binomial-array readings

> **SCOPE CORRECTION — MISTAKE-222 (2026-07-21).** The exact identities in
> this note survive. “One Pascal triangle” means that each sequence is a fixed
> linear or diagonal reading of binomial coefficients. Bagel and Moser are not
> ordinary truncated-row sums. The equality `bagel-cake=T_n-1` does not by
> itself identify a torus handle with the `-1` in a g-bonacci shadow lattice,
> and no map to JC(2) or LRC(14) has been proved. The original polygonal-skip
> computation had an indexing error and has been removed.

*boxeph-2026-07-21-S207. Owner: relate the repository's figurate work to
Fibonacci and the cake/bagel cutting sequences. Builds on opus-S317,
klein-S313, and mac-mini-S137.*

## Exact binomial atlas

The displayed binomial formulas hold for `n>=0`, with a separately defined
bagel value at `n=0`. The shifted-row description of Moser starts at `n=1`.

| object | exact binomial reading |
|---|---|
| lazy caterer | `C(n,0)+C(n,1)+C(n,2)` |
| cake | `C(n,0)+C(n,1)+C(n,2)+C(n,3)` |
| Moser circle | `C(n,0)+C(n,2)+C(n,4)`; for `n>=1`, also `sum_(k=0)^4 C(n-1,k)` |
| bagel, `n>=1` | `C(n,3)+2C(n,2)+2C(n,1)=cake(n+1)-2` |
| Fibonacci `F_(n+1)` | `sum_k C(n-k,k)` |

The first two rows are prefixes of Pascal row `n`. Moser selects even columns
through depth four, equivalently the first five entries of row `n-1`. Bagel is
a weighted depth-three functional, equivalently the next cake number minus two.
Fibonacci is a shallow diagonal rather than a row functional. These formulas agree with the primary
sequence records [A000125](https://oeis.org/A000125),
[A000127](https://oeis.org/A000127), and
[A003600](https://oeis.org/A003600).

## The exact cake–bagel difference

For `n>=1`, subtracting the two binomial formulas gives

```text
bagel(n)-cake(n)
  = C(n,2)+C(n,1)-C(n,0)
  = n(n+1)/2-1
  = T_n-1.
```

This is an all-`n` identity, not merely a prefix match. It is the most useful
new relation in the session: it says exactly how the torus cutting count differs
from the ball cutting count. Equivalently, `bagel(n)=cake(n+1)-2`, which is the
strongest simple operation law exposed by the comparison.

Klein-S313 also contains a boundary term described as “deficit one.” Equality
of offsets is a prompt, not a mechanism. A genuine bridge would need a map from
the torus region increment to a shadow-lattice cell set, together with a proof
that the missing cell corresponds to the same boundary operation. No such map
is currently recorded.

## Fibonacci and g-bonacci

The standard diagonal identity

```text
F_(n+1) = sum_{0<=k<=n/2} C(n-k,k)
```

places Fibonacci on a diagonal of the same ambient binomial array. Separately,
the sequences with generating function `1/(1-x-x^(g+1))` satisfy
`a_n=a_(n-1)+a_(n-g-1)`; `g=1` is Fibonacci. With klein-S313's explicit
full-rank indexing,
`D_(infinity,g)(d)=sum_(k>=0) T_infinity(k,d-(g+1)k+1)`, these are exact
gap-diagonal readings. Its finite-rank shadows agree with the full kernel only
until their first deficit; the session script computes the full-kernel
recurrences, not the finite-rank shadows.

## What the shared array does and does not transfer

The binomial atlas is useful because row, weighted-row, and diagonal operations
can be compared in one coordinate system. It suggests concrete questions:

- Which finite-difference operators take cake to bagel or Moser?
- Does the `T_n-1` difference have a region-by-region bijective proof?
- Which boundary-cell or valuation operation explains the finite-rank first
  deficit, and can the bagel's two-region shift be realized by the same
  operation without discarding its source geometry?
- Which support or Dirichlet profiles are preserved by these readings?

It does **not** transfer an LRC or Jacobian predicate. The mac-mini-S137 golden-
corner and the LRC anti-golden comparison remain analogies until a source,
target, map, preserved predicate, and loss ledger are supplied.

## Reproduction

Run:

```bash
python3 04-computation/cake_bagel_figurate_fibonacci_boxeph_S207.py
```

The script checks the displayed finite prefixes and recurrence outputs. The
all-`n` cake–bagel identity is the algebra above; the computation is a referee,
not its proof.

Related: HYP-8820, MISTAKE-222,
`the-vandermonde-truncation-law-polygonal-vs-polyhedral-opus-S317`, and
`the-hurwitz-principle-jc2-golden-corner-lrc-ghost-convergents-gbonacci-macmini-S137`.
