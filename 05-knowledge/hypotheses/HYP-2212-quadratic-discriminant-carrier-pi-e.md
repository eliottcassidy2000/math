# HYP-2212: The pi/e Quadratic Carrier Has a Discriminant Sheet, So At Least Two of `e+pi`, `e*pi`, and `e-pi` Are Transcendental

**Status:** PROVED sharpening plus open proof routes.  S636 adds the
finite-height audit and carrier interpretation.

## Claim

Let

```text
S = e + pi,
P = e*pi,
D = e - pi.
```

Then no two of `S`, `P`, and `D` can both be algebraic.  Consequently at least
two of the three numbers are transcendental.

This sharpens HYP-2211's two-shadow statement.  The standard trace/norm pair
`(S,P)` has a hidden square-root sheet:

```text
S^2 - 4P = (e-pi)^2 = D^2.
```

The quadratic carrier is therefore not just trace and norm; it is
trace/norm/discriminant.

## Proof

Each coordinate pair reconstructs `e` and `pi` algebraically:

```text
S,D:  e=(S+D)/2,  pi=(S-D)/2.

S,P:  e and pi are roots of T^2 - S*T + P.

D,P:  e and -pi are roots of T^2 - D*T - P.
```

If any pair among `S,P,D` were algebraic, then the displayed formulas or
quadratic equations would make `e` and `pi` algebraic over `Qbar`, hence
algebraic over `Q`, contradicting the known transcendence of `e` and `pi`.

Thus at most one of `S,P,D` is algebraic, so at least two are transcendental.

## S636 Evidence

S636 adds `04-computation/quadratic_pi_e_carrier_s636.py` and stores
`05-knowledge/results/quadratic_pi_e_carrier_s636.out`.

The PSLQ height sieve is explicitly non-proof evidence.  At 100 digits it finds
no scalar relation of degree `<=8` and coefficient bound `<=1e8` for `S`, `P`,
or `D`; no pair relation of total degree `<=4` and coefficient bound `<=1e7`
for `(S,P)`, `(S,D)`, or `(P,D)`; and it correctly detects the true carrier
relation:

```text
S^2 - D^2 - 4P = 0.
```

The point of the sieve is not to prove irrationality.  It verifies that the
script can see the intended quadratic sheet while finding no small accidental
relations around it.

## Conditional Completion

If Schanuel's conjecture holds, then `e` and `pi` are algebraically independent.
Under that conditional, every nonconstant polynomial `f(e,pi)` with algebraic
coefficients is transcendental: if `f(e,pi)=alpha` with `alpha` algebraic, then
`f(X,Y)-alpha` would be a nonzero algebraic-coefficient polynomial relation.

So Schanuel would imply:

```text
e+pi, e*pi, e-pi are all transcendental,
and (S,P), (S,D), (P,D) are algebraically independent pairs.
```

This explains why a direct proof of the rationality questions is so hard.  It
is not enough to know individual transcendence.  One needs to rule out
specific algebraic curves containing `(e,pi)`, and the expected theorem is
algebraic independence.

## Exponential Lift Dead Ends

Two tempting lifts do not currently close:

1. If `S` were algebraic, then

   ```text
   exp(S) = exp(e)*exp(pi).
   ```

   Since `exp(pi)` is transcendental and `exp(S)` is transcendental for
   nonzero algebraic `S`, this gives no contradiction without a theorem about
   algebraic independence among `exp(e)`, `exp(pi)`, and `exp(S)`.

2. If `P` were algebraic, then

   ```text
   pi = P/e,
   exp(pi) = exp(P/e).
   ```

   Standard Lindemann-Weierstrass/Gelfond-Schneider tools do not apply because
   the exponent `P/e` is not algebraic.

This is the practical disproof map: every obvious route asks for an algebraic
independence theorem, not a one-variable transcendence theorem.

## LRC Carrier Reading

For LRC, `S` behaves like a trace/reset clock, `P` like a norm or
multiplicative shell, and `D` like the discriminant sheet identifying which
branch/observer owns which coordinate.  Forgetting `D` is the same
observer-blind collapse that HYP-2210 warns against: a quotient may remember a
trace and norm while losing the branch identity needed to certify the predicate.

The transferable schema is:

```text
For a generically finite quotient, record enough branch/discriminant data
to know which coordinate pair forces descent of the hidden source.
```

## Assumption Challenge

Alternate vertices considered: constants `e/pi`, shadows `S/P/D`, polynomial
relations, exponential lifts, log commensurability, PSLQ height sieves, LRC
clocks, relation lattices, and proof routes.

S636 chooses proof routes as Tournament Analysis vertices.

Preserved predicate: whether a route can rule out algebraic descent.

Destroyed data: exact branch identity, individual constants, and which shadow
carries transcendence.

## Next Tests

1. Generalize the three-shadow lemma to arbitrary generically finite maps:
   identify coordinate pairs that force algebraic descent of the hidden source.
2. Add discriminant-sheet labels to LRC reset-period and relation-lattice
   ledgers.
3. Search the repo for other triples where any two shadows reconstruct a known
   forbidden algebraic source.

## See Also

`04-computation/quadratic_pi_e_carrier_s636.py`;
`05-knowledge/results/quadratic_pi_e_carrier_s636.out`;
`07-reflections/quadratic-pi-e-carrier-s636.md`;
HYP-2211; HYP-2210; HYP-2154; HYP-2155;
`04-computation/transcendental_bases_s116h.py`;
`04-computation/logarithm_deep_s91g.py`.
