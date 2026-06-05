---
id: HYP-2238
status: OPEN method hypothesis with S662 finite representation audit
source: codex-2026-06-05-S662
related:
  - HYP-2229
  - HYP-2237
  - HYP-2236
  - HYP-2235
  - HYP-2231
  - HYP-2230
  - HYP-2186
  - HYP-2184
  - THM-401
---

# HYP-2238: Pi Sum/Product/Fraction Trinity Carrier

## Claim

When `pi` appears as a carrier value, keep a trinity record instead of a scalar
record:

```text
infinite sum      = additive packet/moment face
infinite product  = local factor/norm face
infinite fraction = recursive boundary/owner face
```

The point is not that one face is "the best" definition of `pi`.  The useful
claim is methodological: the three faces preserve different proof payloads, and
collapsing all of them to decimal digits loses the side channel that usually
matters in this repo.

## Classical Faces

S662 audits five standard representations:

```text
Leibniz sum:      pi/4 = sum_{k>=0} (-1)^k/(2k+1)
Machin sum:       pi/4 = 4 arctan(1/5) - arctan(1/239)
Basel sum:        pi^2/6 = sum_{n>=1} 1/n^2
Wallis product:   pi/2 = product_{n>=1} (2n)^2/((2n-1)(2n+1))
Brouncker frac:   4/pi = 1 + 1^2/(2 + 3^2/(2 + 5^2/(2 + ...)))
Cot pole frac:    pi*cot(pi*z)=1/z+sum_{n>=1}(1/(z+n)+1/(z-n))
```

The convergence table is not the main theorem, but it gives a useful sanity
check.  With the S662 cutoffs:

```text
Leibniz N=10000          error ~ 1.0000e-4
Basel sqrt N=10000       error ~ 9.5490e-5
Wallis N=10000           error ~ 7.8535e-5
Machin 20 arctan terms   error ~ 8.2663e-30
Brouncker depth=100      error ~ 9.9007e-3
```

The numerical lesson is familiar: term choice and acceleration matter.  The
repo lesson is sharper: speed is not the same as retained information.

## Carrier Payloads

Sums expose additive packets, moments, and cancellation order.  This is the
face closest to LRC pair sums, Goldbach/Lemoine pair projections, Basel power
sums, and OCF coefficient sums.

Products expose local factors, zeros, norm constraints, and sieve channels.
This is the face closest to HYP-2229's sine product, Euler products, Wallis
factors, von Staudt denominators, local obstruction products, and CM/unit-norm
unit-distance packets.

Continued fractions expose convergents, continuants, and recursive boundary
state.  This is the missing face in HYP-2229's Basel note: it looks less like a
moment identity and more like a finite owner ledger.  Each convergent remembers
which boundary state was last retained, making it a natural analogue of S661's
desired paired carry/owner/deletion derivative.

The cotangent partial-fraction form is the same warning in pole language.  It
is a scalar identity only after one forgets which pole supplied which boundary
contribution; for LRC, that forgotten assignment is exactly the kind of
owner/carry side channel HYP-2237 is trying to retain.

Raw decimal digits are useful as checksums, but they are not proof carriers.

## S662 Tournament Analysis

`04-computation/pi_sum_product_fraction_s662.py` uses representation faces as
the tournament vertices:

```text
sum, product, fraction, raw_decimal.
```

The pairwise observable is:

```text
which side channel is preserved better?
```

The switch is a majority gauge over three coordinates:

```text
additive_packet, local_factor, recursive_boundary
```

with metric triples:

```text
sum         (3,1,2)
product     (2,3,1)
fraction    (1,2,3)
raw_decimal (0,0,0)
```

This produces a genuine trinity cycle:

```text
sum -> product -> fraction -> sum
```

and all three beat `raw_decimal`.  The fingerprints are:

```text
outscores={sum:2, product:2, fraction:2, raw_decimal:0}
score_hist={0:1, 2:3}
directed_3cycles=1, namely (sum, product, fraction)
sccs=[{sum, product, fraction}, {raw_decimal}]
hamiltonian_paths=3
```

That is the desired shape: the three representations are mutually correcting
faces, not a transitive ranking.

## Repo Transfer

HYP-2229 already gives the product/sum bridge:

```text
sin(pi*x)/(pi*x) = product_{n>=1} (1 - x^2/n^2)
```

The product face carries elementary factors; the log derivative turns them into
power sums `zeta(2k)`.  S662 adds the third face:

```text
continued fraction = recursive boundary state.
```

This suggests a practical pattern for future scalar-valued invariants:

```text
triune_value = (sum packets, product factors, fraction boundary state)
```

Use the sum face to see additive cancellation, the product face to see local
obstructions, and the fraction face to see recursive ownership.

For LRC `n=14`, this complements S661.  The no-leak theorem should not only
keep odd-wall sums and `C=27` product/gcd shell data; it should also keep the
continued-fraction-like boundary state of the carry recursion:

```text
Res_27/certificate card + attached carry/owner derivative.
```

So there are two compatible readings:

```text
global representation audit: sum -> product -> fraction -> sum
LRC14 carry-seam task:       fraction/carry first, product/gcd second, sum third
```

The first says no face globally dominates.  The second says that the current
`n=14` proof obstruction is specifically a branch/carry reconstruction problem.

For OCF/tournament work, the analogous warning is: `H(T)` is a scalar sum over
packets, the independence polynomial is the retained packet object, and a
continued-fraction/continuant encoding may be the right way to remember a
Hamiltonian path boundary through deletion or substitution.

## Assumption Challenge

Candidate Tournament Analysis vertices considered:

```text
values, formulas, terms, factors, convergents, proof lenses,
side-channel types, LRC residues, OCF packets, and owner states.
```

S662 uses representation faces because the session predicate is not numerical
equality.  The predicate is which information survives scalar collapse.

This destroys the fine internal structure of each individual formula.  That is
acceptable for this audit, but the next serious transfer should refine the
fraction face into actual convergent states or owner/carry derivatives.

## Next Experiments

- Build a small `triune_value` template for repo constants: scalar, sum packet
  list, product/local factor list, continued-fraction boundary data.
- Apply the template to LRC `n=14` floor atoms: odd-wall sums, `C=27` gcd
  products, and carry-owner continuants.
- HYP-2239/S663 now owns the concrete application lane.  Its first LRC script
  should use columns `(sum shadow, product card, fraction branch)` and test
  whether identical scalar/wall/product data but different carry-owner
  derivatives can disagree on floor-vs-strict status.
- Test an OCF continuant analogue: can a deletion sequence or substitution
  macro-word encode `H(T)` with a fraction-style boundary state rather than a
  raw scalar?
- Revisit pi/e trace-norm work (HYP-2211/HYP-2212) with the same trinity:
  trace sums, norm products, and branch/fraction sheets.

**See:** `04-computation/pi_sum_product_fraction_s662.py`;
`05-knowledge/results/pi_sum_product_fraction_s662.out`;
`07-reflections/pi-sum-product-fraction-trinity-s662.md`; HYP-2229, HYP-2237,
HYP-2236, HYP-2235, HYP-2231, HYP-2230, HYP-2212, HYP-2211, HYP-2186,
HYP-2184, THM-401.
