# LRC Totient Curvature and Tournament Recursion Modes

The user's warning about the three recursion modes is the key point.  The same
letters do not mean the same sizes in the three formulas.

Full:

```text
A+B+C-D-E-F+G
A,B,C: n-1
D,E,F: n-2
G:     n-3
```

Even half:

```text
A+B-C
A,B: n-1
C:   n-2
```

Odd half:

```text
A+B-C+D-E-F+G
A,B: n-1
C,D: n-2
E,F: n-3
G:   n-4
```

This turns out to be exactly the right guardrail for the totient analogy.
`phi(n)` and `rho(n)=phi(n)/n` are multiplicative exact-period packet counts;
the `A..G` formulas are additive address recurrences for cell carriers.  They
are not supposed to agree.

The useful object is the residual.

At the LRC14 boundary:

```text
even_half n=14:
  A,B are size 13 prime carriers
  C is size 12 = 2^2*3
  actual n is 14 = 2*7
```

So the residual is exactly the Euler-factor curvature that says "you cannot
derive the composite 7-seam from the prime 13 and old dyadic/triadic carrier
without retaining labels."  The script records:

```text
even_half n=14 residual rho = -296/273
curvature = {2: 3, 3: 1, 7: 1, 13: -2}.
```

Odd size `15` is also instructive.  The odd half recurrence uses size `14` in
the `A,B` slots:

```text
odd_half n=15:
  A,B size 14
  C,D size 13
  E,F size 12
  G size 11
```

Thus LRC14 is an input carrier for the odd half recursion, not an output of it.
That corrects the tempting story where the odd recurrence somehow explains the
14-runner boundary directly.

The repo already had the pieces:

- HYP-2628/HYP-2629: `n=sum_{d|n}phi(d)` is exact-period packet bookkeeping,
  not a decorative totient identity.
- HYP-2630/HYP-2632/HYP-2883: copy mass is uniform over `F_7^*`; the useful
  residual is the `chi_7`/affine/signed-current layer.
- THM-442/THM-550/THM-554: the tournament recurrences are cell-address and
  score-partition recurrences, not Hamiltonian-path or OCF recurrences.
- HYP-2898: scalar additive energy fails; labelled Fejer/residual control is
  the robust LRC14 object.

This session's added synthesis is the name for the bridge: **Euler curvature**.
Apply the signed slot recurrence to a multiplicative coprime-density packet.
The nonzero signed remainder is not an error term to crush.  It is the finite
CRT/prime-exponent address that must be carried into the LRC support-six
residual.

The incoming KPS S31q and mac-mini S44 threads make this sharper.  KPS S31q
reads the sign words themselves as characters: Mobius for `+++---+`, Legendre
`chi_7` for `++-+--+`, and Eisenstein/`chi_3` for the `++-` omega mode.
mac-mini S44 reads resonance killing as totient-weighted Farey killing:
divisibility by `b` kills all `phi(b)` primitive denominator-`b` points, with
`Phi(14)=64`.  So the live bridge has two coordinates, not one:

```text
character sign channel  x  exact slot size / prime-exponent curvature.
```

KPS S31r adds a third coordinate: parity.  The even mode `A+B-C` is the
Eisenstein/pronic fold and the odd seven-term mode is the Legendre/square
channel.  In that language LRC14 is literally

```text
14 = 2*7 = Eisenstein(even fold 14 -> 7) o Legendre(odd apex 7)
```

on top of the Mobius floor.  The exact residual `-296/273` for even-half
`rho` at `n=14` is the local arithmetic scar of that fold: the operator has
exposed the apex `7`, but the slot recurrence still samples size-13 and size-12
carriers unless the exact-period labels are retained.

That is why the owner's warning about different subtournament sizes matters.
The same sign character sampled on the wrong slot sizes gives the wrong
Euler-factor curvature.  S31q names the characters; S44 names the
totient-weighted killing; S31r names the parity composition; S113 pins the
exact slot sizes and residuals at the LRC14 boundary.

The proof target I would carry forward:

```text
exact-period phi packets
  -> Mobius / chi_7 / chi_3 character decomposition
  -> parity composition 2q -> q
  -> Euler-curvature residual of the full/even/odd slot atlas
  -> mod-7/coimage/chi-labelled signed-current packets
  -> HYP-2890 Fejer residual leak
  -> scalar cap/floor only at the last step.
```

That is also the best answer to the coprime-density question.  Coprime density
is multiplicative, but LRC loneliness is archimedean.  The multiplicative side
supplies packet capacity; the archimedean side supplies the Fejer/cap geometry;
the curvature residual is the interface where neither side may be forgotten.
