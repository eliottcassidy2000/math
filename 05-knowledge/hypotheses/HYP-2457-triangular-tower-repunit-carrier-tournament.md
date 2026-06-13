# HYP-2457 - Triangular-tower repunit carrier tournament

**Status:** OPEN synthesis; exact finite atlas.
**Source:** codex-2026-06-12.
**Companions:** HYP-2453, HYP-2450, HYP-2448, OPEN-Q-079.
**Computation:** `04-computation/triangular_tower_repunit_tournament_codex.py`;
stored output `05-knowledge/results/triangular_tower_repunit_tournament_codex.out`.

## Statement

The triangular-tower overlap families from HYP-2453 naturally carry
irreducible-polynomial data.

If an overlap block is the consecutive interval

```text
[a, a+L-1],
```

then its coefficient-row shadow is

```text
x^a R_L(x),  where  R_L(x)=1+x+...+x^(L-1).
```

The shift `x^a` is irrelevant for irreducibility, so the live datum is the
length `L`.

This turns the overlap atlas into an exact repunit/cyclotomic carrier:

1. `L` prime  =>  `R_L(x)=Phi_L(x)` is irreducible over `Z`.
2. `R_L(2)=2^L-1` prime  =>  Cohn certifies irreducibility of `R_L(x)`.
3. The fixed-path coefficient-row lens from HYP-2450 is literal here: the
   overlap family is a consecutive-support row, not merely an analogy.

The main qualitative claim is not that the unique exact hinge is the best
irreducibility carrier.  It is the opposite:

```text
exact overlap and irreducibility supply are different carriers.
```

The overlap-family tournament is therefore allowed to be nontransitive.

## Exact Finite Evidence

The stored run combines the whole-equation Pell family with the eight
endpoint-aligned side families from S713.

### 1. The exact hinge is reducible

The unique whole-side hinge is still

```text
B_3.L = A_4.R = [21,24].
```

Its length is `L=4`, so the carrier polynomial is

```text
R_4(x)=1+x+x^2+x^3=(x+1)(x^2+1),
```

which is reducible.  So the structurally most exact overlap is NOT the best
irreducibility lane.

### 2. The whole-equation Pell family is the richest prime-length lane

For the side-aligned whole-equation family `T_n=2T_m`, the block length is

```text
L = 2m+1.
```

Within the exact scan window `n<=10^6`, the observed lengths are

```text
5, 29, 169, 985, 5741, 33461, 195025, 1136689,
```

and the prime lengths among them are

```text
5, 29, 5741, 33461.
```

So this one family already contributes four exact cyclotomic-irreducible rows
in the stored window.  The first two also satisfy the base-2 Cohn condition:

```text
R_5(2)=31  prime,
R_29(2)=2^29-1=536870911  prime.
```

### 3. The side families contribute only the short prime lengths 2, 3, 5 in the stored window

Observed prime lengths by family:

- `BL prefix AL`: `5`
- `BL suffix AL`: `3`
- `BR prefix AL`: `3`
- `BR prefix AR`: `2`
- `BL prefix AR`, `BL suffix AR`, `BR suffix AL`, `BR suffix AR`: none in the stored window

So the side families excel at exactness/support, but their prime-length
irreducibility supply is much thinner in the scanned range.

## Tournament Analysis

Vertices are the family carriers themselves:

```text
whole_equation
BL prefix/suffix AL/AR
BR prefix/suffix AL/AR
```

Pairwise observable:

```text
majority(window_hits, prime_length_hits, mersenne_hits,
         exact_hits, support_ratio).
```

Tie path:

```text
larger support ratio, then larger prime-length count, then label.
```

Stored fingerprint:

```text
score_hist = {0:1, 1:1, 3:1, 4:3, 6:1, 7:2}
directed_3cycles = 6
SCC sizes = [7,1,1]
Hamiltonian paths = 33
edge_flips_vs_support_exactness_only = 8
ranking =
  BL suffix AL,
  BL suffix AR,
  whole_equation,
  BL prefix AL,
  BL prefix AR,
  BR prefix AR,
  BR prefix AL,
  BR suffix AR,
  BR suffix AL
```

The important point is the edge-flip count.  Prime-length information changes
`8` pairwise decisions relative to the support/exactness-only ranking.  So the
irreducible-polynomial carrier is not decorative; it genuinely reshuffles the
family geometry.

## Interpretation

The overlap atlas splits into three different notions of “best”.

1. **Exactness winner:** the `A_R` hinge families (`BL prefix AR`, `BL suffix AR`)
   contain the unique exact side equality.
2. **Prime-length winner:** the whole-equation Pell family supplies the richest
   exact cyclotomic lane in the stored window.
3. **Tournament winners:** `BL suffix AL` and `BL suffix AR` sit at the top of
   the mixed observable because they balance density, support retention, and
   exactness differently.

So the right slogan is:

```text
the overlap families are not ranked by one scalar;
they form a nontransitive certificate tournament.
```

## Assumption Challenge

Alternate vertices considered: raw intervals, lengths only, Pell equations,
coefficient rows, primes, Mersenne exponents, fixed-path tilings, and proof
obligations.

Chosen vertex set: overlap-family carriers.

Preserved: containment type, support-retention ratio, exact-equality events,
prime-length cyclotomic irreducibility, and base-2 Cohn certificates.

Destroyed: exact shift `x^a`, asymptotic prime production, and any support
channel not visible from the block length.

Challenged assumption: the most exact structural overlap should also be the
best irreducibility carrier.  The exact hinge length `4` shows that this fails
immediately.

## Next Moves

1. Attack OPEN-Q-079: do any of these Pell-family length sequences contain
   infinitely many primes?
2. Add exact base-`b` repunit certificates beyond base `2`, not just the
   Mersenne lane.
3. Push the same carrier into HYP-2450's fixed-path coefficient quotient and
   ask which family lengths persist as irreducible across richer support fibers.
4. Compare the prime-length carrier with HYP-2453's higher-moment fractional
   addresses: does a prime `L` correlate with a cleaner moment-lift ledger?
