# Independent hostile audit: the desheeted alternate diagonal is a delta-cell lift

**Status: FINITE-EXACT INDEPENDENT AUDIT ACCEPTS THE SCOPED DESHEETED
COMMON-RESIDUAL PACKAGE.  LRC(14) remains open.**  The candidate was not
imported.  The audit started from the hash-pinned THM-3514 endpoint engine and
the old point-diagonal API conventions, then rebuilt the residual sweep with
a different bitmask algorithm.  Candidate artifacts were read only after the
clean-room semantic surface had been fixed.  This reflection now incorporates
the exact owner-geometry correction recorded as MISTAKE-417; it withdraws the
original mixed-cell interpretation while preserving the independently checked
alternate diagonal and residue bridge.

## Verdict

The coordinate

```text
91t = 7a+r,       y = r/7 = 13t-a,       t=(y+a)/13
```

does provide a character-independent common endpoint base.  On that base,
the old point diagonal is exactly the same-sheet sector, while allowing all
ordered sheet pairs gives a larger common-residual coupling.  Same-sheet plus
cross-sheet equals the full coupling before every Fourier transform.

The full residual bridge is nonzero but differs from the original Cartesian
endpoint bridge.  However, the endpoint owner factor forces every exact
character row into cell zero.  Each apparent `7 x 13` table is therefore
`delta_0(ell)R(t)` and has matrix rank one.  Its formally complete Fourier
support and 72 centered mixed modes are a separable outer-product artifact,
not genuine cell/residue mixing.  Thus the package constructs a natural
alternate endpoint diagonal but no septimal spectral closure.

## Why both the word and phase descend

Write an endpoint point in sheet `a` as `t=(y+a)/13`.  Then

```text
169t = 13y+13a = 13y mod 1,
```

so the delayed factor is the same `Q(13y)` on every sheet.  The target
frequency has the load-bearing factorization

```text
742586 = 13^2 * 4394.
```

Consequently its sheet phase is integral, and the residual frequency is
`742586/13 = 57122`.  The audit checked both identities on every transformed
endpoint, not only symbolically or at sample points.

For exact integer arithmetic, if the original endpoint grid has length `T`,
the audit uses

```text
u = 169(x-aT/13),          0 <= u < 13T.
```

Then `y=u/(13T)`, the word coordinate is `u mod T`, and the endpoint phase is
`root^(742586 u)`.  This simultaneously eliminates fractions and exposes the
sheet cancellation.

## The owner hostile and exact cell collapse

The load-bearing condition was already present in every endpoint interval:

```text
OWNER:  ||13t|| < 1/14.
```

On the desheeted branch `t=(y+a)/13`, integer sheet translation disappears:

```text
||13t|| = ||y+a|| = ||y|| < 1/14.
```

The right-hand side is exactly `cell_0`.  The interval sweep now requires
`ell=0` on every active segment before evaluating a character or reducing
modulo the split prime.  Consequently all 2,197 character rows have geometric
occupancy `(2197,0,0,0,0,0,0)`, rows 1 through 6 of every inverse table vanish
in characteristic zero, and the three table ranks are exactly one.  The
earlier wording "zero only mod p" was false: the zeros are forced over the
rational interval model.

## Guard restoration and the multiplicity mechanism

For each `(alpha,beta)`, the guard was removed before atomization, exactly as
in THM-3514.  The audit independently derived the guard value on every one of
the `39 x 13` atom/twist states from the literal inequality.  It then compared
the selected unguarded atom union with the directly built fully guarded set
for all 78 triples over the six inherited control pairs and all thirteen
guard shifts.  Every interval endpoint agrees.

After desheeting, let `n(y)` be the number of guarded endpoint sheets active
at a residual point.  The three integrands are then

```text
same:   Q(13y) n(y),
cross:  Q(13y) n(y)(n(y)-1),
full:   Q(13y) n(y)^2.
```

The clean-room implementation sweeps residual interval endpoints while
carrying the active sheets as a bitmask.  This realizes the identity
`n+n(n-1)=n^2` pointwise, so same plus cross equals full before character
summation, residue inversion, or spectral testing.

There is an intrinsic directed-fibre interpretation: `n` counts loops,
`n(n-1)` counts ordered nonloop arcs, and `n^2` counts the complete directed
fibre with loops.  This is a valid tournament-adjacent decomposition because
the binary observable is literal equality of residual base points.  It does
not supply an orientation, flux, or ancestry current.

## Exact recovery and spectra

Summing cells in the same-sheet character bank reproduces the historical
point-diagonal digest exactly:

```text
771545a5cb1f0f03459b8d351de668ad950ece5fcb985fa61d599d643de3303f.
```

Its two inverse values and bridge are

```text
q_H  =   633668780131603861,
q_q5 =405160484437854840264,
bridge=167726070588785644466.
```

The independent full and cross character digests, all table digests, and all
spectrum digests agree exactly with commit `6d50f8b6a`.  The full coupling has

```text
q_H  = 95783417057771114126,
q_q5 =124341028951542154618,
bridge=543695274352737840377,
```

and its cross-sheet bridge is `375969203763952195911`.  The Cartesian bridge
is `389266878372286537904`, so equality fails exactly as the typing predicts.

For each of same, cross, and full, the formal two-dimensional spectrum still
has census

```text
(total, DC, F_7 axis, F_13 axis, mixed)=(91,1,6,12,72).
```

This does not show two-coordinate interaction.  If `R(t)` is the cell-zero
row, then exactly

```text
table(ell,t) = delta_0(ell) R(t),
ANOVA(ell,t) = (delta_0(ell)-1/7)(R(t)-mean(R)).
```

Both matrices have rank one.  The seven Fourier transforms of `delta_0` are
all `1`, so every septimal frequency simply repeats the same nonzero
one-dimensional residue coefficient.  The audit checks this equality at all
91 frequencies and verifies that the residue profile alone has 13/13
nonzero modes.  At `(1,0,6)` the repeated full-sector value
`289814661037836286866` is therefore a useful residue nonvanishing witness,
but it contains no independent cell information.

## Strict boundary

The desheeted base uses actual endpoint interval integrands, but it has not
been identified with the canonical THM-2471 collision base.  In particular it
contains no collision root, source packet, word/source sheet, horizon, or
chronological arrival.  It also performs only the refined residue pushforward
and does not isolate `C(a;X,m)`.

Therefore none of the following follows:

- equality with the original Cartesian endpoint current;
- a physical LRC current or THM-2512 bridge;
- a grouped exact-address coefficient or row exclusion; or
- LRC(14).

Any useful next map must first introduce a refiner not implied by `OWNER` and
pass two cheap hostiles: at least two cells must be occupied in characteristic
zero, and the cell-by-residue matrix must have rank at least two.  Only then
would lifting through a typed collision record, chronology, and exact address
be worth interpreting as a mixed current.

## Reproducibility

Normal and optimized transcripts are byte-identical.  The repaired semantic
digest is `9d2070f27bcac8cf576a75bc542222be1e31679bc12e989c2f9bc276e8dd872c`.
The final script and output digests are recorded in the results index.
