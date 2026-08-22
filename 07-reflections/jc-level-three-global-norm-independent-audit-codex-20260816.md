# Independent audit of the fixed-map level-three norm divisor

**Verdict: PASS.**  Proved
[THM-3495](../01-canon/theorems/THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness.md)
survives every requested hostile gate.  In primitive integral normalization,

```text
N(H)=J/(2^35 L^7),
S_(F^3)=V(LHJ),
[Delta_3]=[-2J].
```

The older residue sidecar used `J_res=L^7N(H)=J/2^35`.  Its local statement
`[-J_res]` is therefore `[-2J]`, not `[-J]` for primitive `J`.  The older
reflection and the maintained frontier now state this conversion explicitly.

## Independent evidence

The audit used two disjoint companions in addition to replaying all three
producer scripts.

1. On the previously unused target slice `(b,c)=(1,1)`, exact
   regular-representation arithmetic gives

   ```text
   denominator N(H)=2^35 [a(27a-2)]^7.
   ```

   The corresponding primitive `J(a,1,1)` has degree 86, content one, is
   irreducible over `Q`, and is coprime to the specialized `L` and `H`.  Its
   coefficient-ledger SHA-256 is
   `f0b51842adf88383230a056f9bcdd19c434e17551940c23aea2fad2964c57447`.
   This independently fixes the sign and the full `2^35` normalization.

2. At the two rational intersections of that slice with `L=0`, namely
   `a=2/27` and `a=0`, the complete norm agrees exactly with

   ```text
   J_res|_L=(11^6 13^4/(2^17 3^15)) lambda^8 P(tau,lambda)
   ```

   at `(tau,lambda)=(1,1)` and `(1,3)`.  Thus neither the root-product sign
   nor the rational scalar is reversed.

3. The 527-term residual `P` was not trusted to the producer's bivariate
   factor call.  Its coefficients as a polynomial in `lambda` have gcd one
   in `Q[tau]`, excluding a `tau`-only factor.  The specialization
   `P(-1,lambda) mod 449` retains degree 84 and passes Rabin's irreducibility
   criterion: Frobenius closure at 84 and the three gcd tests associated to
   the prime divisors `2,3,7` of 84.  Hence `P` is irreducible over `Q`.
   Also `P(tau,0)` is nonzero of degree six, while `P(1,2) != 0`; the explicit
   `lambda^8` intersection is real, and the excluded seam
   `lambda tau=2` is not a hidden factor.

4. The global and tower normalizations meet in a new exact identity.  Since
   the leading coefficient of the degree-27 product is `N^2(L)`, the two norm
   laws imply

   ```text
   LC(P_3)=N^2(L)=J/(2^47 L^6 H).
   ```

   At `(1,1,1)`, this value equals the incoming primitive leading coefficient
   divided by its interpolation denominator, exactly.

5. A second coordinate-separation certificate avoids the characteristic-zero
   interpolation entirely.  Over `F_101`, target `(93,28,83)` has the split
   inverse fibres

   ```text
   3 first points,
   3+3+3 second points.
   ```

   Every one of the nine next cubic cores keeps degree three.  Grouped by the
   first point, their three degree-nine products are squarefree and pairwise
   coprime; the direct nine-cubic product is squarefree of degree 27.  All
   inverse denominators are units modulo 101.  This good-reduction witness
   independently proves that generic `x`-coordinate separation and block
   coprimality are nonempty in characteristic zero.

## Logical audit

The global proof has four correctly separated layers.

- THM-2473 makes `F` finite etale of degree three away from `V(L)`, so the
  norm is regular in `Q[a,b,c,L^{-1}]`; this controls denominator support.
- The generic-DVR Newton calculation gives two valuations `-7/2` and one
  valuation zero, so `v_L(N(H))=-7`; this controls the exact denominator and
  rules out cancellation.
- The dense part of irreducible `V(H)` has one irreducible finite image
  divisor.  Thus `L^7N(H)=uP^e`.  Either the inherited `(1,2)` slice or the new
  `(1,1)` slice rules out `e>1`, proving generic image degree one and absolute
  irreducibility of the image equation.  This argument, not computational
  squarefreeness, proves `J` prime.
- The full-degree squarefree 27-sheet specialization licenses THM-2582.
  Substitution of the primitive norm gives
  `[-LN(H)]=[-J/(2^35L^6)]=[-2J]`.

There is no silently lost exceptional hypersurface.  The exact global
resultant cancels the chart factor `S`; finiteness excludes every denominator
prime except `L`; the DVR residue excludes cancellation of `L`; and the
normalization factors `lambda` and `lambda tau-2` have been typed separately.

Finally, THM-2576's composition law gives the reduced closed set
`S_(F^3)`.  Density and continuity replace the two constructible images by
their closures `V(H)` and `V(J)`.  Since `L,H,J` are pairwise distinct
irreducibles, the set has exactly three components.  This does not determine
the scheme multiplicities of a nonproperness ideal.

## Scope

The verdict is only for this fixed three-dimensional map and depth three.
It proves neither an all-depth newest-prime law nor a classification of
Keller maps, and it has no implication for `JC(2)` or `DC(2)`.  The full
degree-27 discriminant was not expanded; only its square class and generic
nonvanishing were proved.

Reproduce with

```bash
python3 04-computation/jc_level3_global_norm_independent_audit_20260816.py
python3 -O 04-computation/jc_level3_global_norm_independent_audit_20260816.py
python3 04-computation/jc_level3_degree27_split_finite_field_audit_20260816.py
python3 -O 04-computation/jc_level3_degree27_split_finite_field_audit_20260816.py
```

Normal, optimized, and stored outputs are byte-identical.  Script/output LF
SHA-256 pairs are

```text
cb429497fbfbdf4cb538967bd472ab10051bd6fac33be10d14081b09e1543215
62fef84710ef14c3a021b196bb4ba3589ad672a3c275a2cfb3c24a25836aea54

9f8b6b4aad75b17eed8a6fd618cb661f3e61b85361b44ba676462148a2c5afc5
af2be3f27d524937836b93aacd5fa0700d5aa930c798dcf142b10206619fc520
```
