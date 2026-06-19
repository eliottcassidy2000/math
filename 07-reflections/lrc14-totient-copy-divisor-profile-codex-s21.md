# LRC14 Totient-Copy Divisor Profile - S21

The copy-rule prompt turned out to be a very clean way to say what the
squarefree-profile route was trying to say.

The equation

```text
sum_{d|n} c(d) = n
```

is `1*c=id`, so `c=mu*id=phi`.  The important part is not just the identity.
It is the model: on a denominator `q` grid, the residues split by exact reduced
denominator `d|q`, and the packet size is `phi(d)`.  So the divisor profile is
not a static list of factors.  It is an exact-period census.

That changes the LRC14 reading of the last few sessions:

```text
raw denominator -> exact-period phi packets -> squarefree masks -> coimage tail
```

This is better than raw divisor counts and better than jumping straight to the
radical.  Raw divisor counts give every divisor weight `1`, which ignores unit
density.  Straight radical projection loses prime powers.  The raw `K_14`
denominator `1260=2^2*3^2*5*7` keeps both:

```text
selected p^a layer contributes p^a-1.
```

So relative to squarefree `210`, the raw profile is amplified by `3` on the
dyadic coordinate and by `4` on the triadic coordinate.  That is exactly what a
pair of `6` blocks should be doing if the Hill row `7,6,6,5` is a denominator
carrier and not just a numerological echo.

The nicest new resonance is:

```text
phi(2520) = 576
full-mask mass inside the raw 1260 profile = 576.
```

The half-clash denominator `2520` is therefore not arbitrary.  It counts the
all-four-prime exact-period mass of the raw `1260` profile, and then the
`tau <-> 1-tau` symmetry halves the denominator back to `1260`.

The more proof-flavored surprise came after weighting actual safe centers.  If
we count strict-safe residues `a/Q` for subsets of `{1,...,13}`, then `Q=210`
catches only `11` of the `13` one-drop AP cores.  The misses are exactly drops
`6` and `12`.  But `Q=1260` catches all `13` one-drop AP cores, while the full
AP13 still has zero strict-safe residues.

That is a tiny model of the whole proof problem: the tight object remains tight,
but every one-runner deletion has a raw-denominator witness.  It also explains
why radical `210` is not enough even though it is the right coimage address.
`210` is the address; `1260` is the exact-period packet resolution needed to
see the local transversality.

This does not prove LRC(14), but it makes the proof route less blurry.  HYP-2626
should probably be re-indexed by exact-period packets before coimage projection,
not by raw residues alone.  The unit seam `(Z/14Z)^* -> F_7^*` is the smallest
example: `phi(14)=6`, and the quotient is literally a six-unit exact-period
packet becoming the six nonzero mod-7 residues.

The challenged assumption this session was that squarefree masks are the atoms.
They are not quite the atoms.  Exact-period `phi` packets are the atoms; masks
are the compressed address after the zeta/Mobius transfer has recorded the raw
denominator information.  That distinction is exactly where the factor `2`
between `1260` and `2520`, and the lost dyadic layer in `315`, live.
