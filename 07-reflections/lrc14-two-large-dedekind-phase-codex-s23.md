# LRC14 Two-Large Dedekind Phase

**Source:** codex-2026-06-19-S23, HYP-2632 / T880.

HYP-2630 said the repeated tail was asking for a two-large
cotangent/Dedekind bound that did not erase the quadratic character.  The new
thing is that the phrase has become computationally concrete.

The right finite object is not just the two large residue pair `(a,b)`.  It is
the additive frequency `m` used to enforce the relation hyperplane
`sum a_i r_i = 0 mod 7`.  Expanding that indicator gives a finite
Dedekind/cotangent packet:

```text
S_d(a) = (1/7) sum_m sum_T (-1)^|T| chat(0,T)^(d-6)
         prod_i D_T(m a_i).
```

Once this `m`-address is retained, the numbers snap into place.  In units
`U=147/(16*pi^6)`, the dominant `4+2` row satisfies

```text
a=2,4:     S = -25U
a=3,5,6:   S = -18U
a=1:       S = 0
```

Equivalently, on the punctured row,

```text
2S/U = -43 - 7 chi_7(a).
```

That is the QR/NQR split from HYP-2630 with the mechanism still visible.  Copy
mass was uniform; the phase was hiding in the additive-frequency shell
distribution.

The `4+1+1` packet is also compressed: unordered unit pairs have signed mass
only in `{0,U,8U}`.  The exact zero rows `(3,6)` and `(4,5)` are the warning
label.  Their blind two-large residue absolute matrices are large, so an
absolute envelope cannot see the cancellation.  A proof must sum by frequency
first.  The incoming HYP-2632 refinement makes the finite address sharper:
the zeros are the unit-domain part of the affine lane `a+b=2 mod 7`, and off
that lane the high/low split is classified by
`Q(a,b)=ab*(1+3(a+b))-1` through `chi_7(Q)`.

The next proof move is therefore narrower than "prove a Dedekind bound":

```text
finite wall deletion
-> retain two-large packet address
-> retain affine lane + Q selector on 4+1+1
-> split reciprocal tail by additive frequencies/conjugate shells
-> apply signed summation-by-parts there
-> only then take absolute values.
```

This is still not LRC(14).  It is the shape of the analytic closer.  The
finite coefficient layer now says the target constants should be `25U` for
`4+2`, `8U` for unit `4+1+1`, and `4U` for the small zero-cusp halo, with
the `4+1+1` constants organized by the affine lane and `Q` selector rather
than a generic 159-class coimage envelope.

Assumption challenged: the tournament vertices are not runners, raw arcs, or
even just residue classes.  The useful vertices at this step are frequency
shells, two-large residue pairs, coimage classes, finite wall ledgers, and lift
obligations.  The quotient preserves signed support-six coimage mass and
destroys raw witness times; that is the right loss.
