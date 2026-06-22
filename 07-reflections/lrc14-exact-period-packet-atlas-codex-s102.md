# LRC14 exact-period packet atlas

The finite-basis idea is now useful only after weakening it to an atlas.

The exact packet predicate is simple.  For denominator `D`, keep only units
`a mod D`; this packet is a level-`1/14` witness for `S` iff

```text
14 * min(s*a mod D, D - (s*a mod D)) >= D
```

for every speed `s`.  The capacity is `phi(D)`, but the LRC predicate is not
multiplicative in `D`, because the safe band is archimedean.

The S102 audit extends HYP-2865/HYP-2876.  Fixed denominator bases are killed
by divisor loading: `divload_B90={1,...,11,13,84*lcm(1,...,90)}` has
`N(21)=N(41)=N(53)=N(83)=N(89)=0`, and first opens at `D=97`.  So there is no
universal finite basis.  What survives is the exact-period packet chart.

The useful new evidence is that safe packets keep mod-7 structure before
scalarization.  In the tested mixed cases, quotient lenses explain packet
safety variance in the order

```text
mod14 > mod7 > chi_7 x affine_pair > affine_pair > chi_7 > parity.
```

`mod14` winning is not the point; it just reflects the `2*7` band.  The proof
relevant part is that `mod7`, `chi_7`, and the affine-pair quotient retain
signal exactly where HYP-2632 and HYP-2883 place the signed-current kernel.

The strong-component analogy becomes concrete here.  Tournament `H` is
multiplicative over strong components.  Exact-period capacity `phi` is
multiplicative over CRT factors.  But LRC safety has a multiplicativity defect:
for example `cover_84` has `rate(7*13)=1/36` even though the product of the
separate rates is `0`.  This defect should be retained as an atom, not treated
as noise.

Best next theorem shape:

```text
delete denominators killed by divisibility;
lift the HYP-2883 local signed-current balance on exact-period residue fibers;
retain CRT defects as labelled strong atoms;
route incoherent high denominators to the spectrum/L2/Part-A witness floor.
```

This does not prove LRC14.  It gives OPEN-Q-108 a cleaner discrete front end
and prevents another false "finite certificate basis" closure.
