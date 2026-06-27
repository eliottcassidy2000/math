# LRC14 Conjecture 7.1 Refutation and Normalized-Arc Repair

The incoming polynomial-method bridge was directionally right but one step too
literal.

The paper's method proves `k=10,11,12`. Its Proposition 4.1 works over the
field `Z/(k+1)Z`, so for `k=13` the modulus `14=2*7` is exactly the composite
wall. The paper's fallback language of repeated lifts through prime factors
matches the project's descent: dyadic work for `c=2`, and THM-573 as the
level-7 lift sieve. HYP-2072's two-tier CRT sieve is also the right shape for
the named `I(k,p,1)` bottleneck.

The false synthesis was identifying the paper's Conjecture 7.1 with a direct
uniform largest interval in time.

THM-566 already warned against bounded denominators. S255 makes the warning
decisive: for

```text
S_B = {1,2,...,11,13,84*lcm(1,...,B)}
```

every denominator `d<=B` is killed by the last speed, but

```text
t = 1/12 + 1/(2*84*lcm(1,...,B))
```

is a strict witness once `B>=6`. Thus the rows are non-tight and primitive, so
Conjecture 7.1 is false for `k=13` as written. The direct largest safe component
also collapses with the apex: for `B=6` the largest direct component is
`1/5880`.

This failure improves the proof map. It says the project should not try to
prove a raw absolute-denominator theorem. The surviving target is a normalized
one:

```text
large-apex tuple
  -> split S = P union {V-e}
  -> prove a uniform slow-coordinate set G(P,E) with c>0 and bounded components
  -> use THM-565 finite-ruler sampling
  -> close the finite complement.
```

That is exactly where THM-530, THM-565, THM-573, and HYP-2072 meet. The
divisor-loaded counterexample sits in the `<=6` multiples-of-7 residual, so it
is not an annoyance off to the side; it is the stress test that forces the
right coordinate.

## Remaining proof obligations

1. Prove the normalized arc floor for the THM-573 residual, not a direct time
   arc floor. This means a uniform lower bound on `meas(G(P,E))` plus a
   uniform component-count/arc-count bound in the slow coordinate.
2. Integrate the HYP-2072 `I(k,p,1)` sieve with that normalized certificate.
   The sieve may use `c=2,7` residue data, but it must preserve the
   slow-coordinate witness object instead of pruning by small absolute
   denominator.

Tournament Analysis in the scout uses proof obligations as vertices rather
than runners. The normalized ruler arc beats the raw Conjecture 7.1 denominator
floor because it preserves the LRC predicate through divisor loading and still
interfaces with the paper's `c=2,7` lift bottleneck.

