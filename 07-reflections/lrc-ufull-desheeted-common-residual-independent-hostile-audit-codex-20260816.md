# Independent hostile audit: desheeting constructs an alternate diagonal, not the Cartesian bridge

**Status: FINITE-EXACT INDEPENDENT AUDIT ACCEPTS THE SCOPED DESHEETED
COMMON-RESIDUAL PACKAGE.  LRC(14) remains open.**  The candidate was not
imported.  The audit started from the hash-pinned THM-3514 endpoint engine and
the old point-diagonal API conventions, then rebuilt the residual sweep with
a different bitmask algorithm.  Candidate artifacts were read only after the
clean-room semantic surface had been fixed.

## Verdict

The coordinate

```text
91t = 7a+r,       y = r/7 = 13t-a,       t=(y+a)/13
```

does provide a character-independent common endpoint base.  On that base,
the old point diagonal is exactly the same-sheet sector, while allowing all
ordered sheet pairs gives a larger common-residual coupling.  Same-sheet plus
cross-sheet equals the full coupling before every Fourier transform.

All three resulting cell-by-residue tables have complete `7 x 13` spectrum,
and their output interactions retain all 72 mixed modes.  The full residual
bridge is nonzero but differs from the original Cartesian endpoint bridge.
Thus the package constructs a natural alternate endpoint diagonal; it does
not identify that diagonal with the Cartesian current or with LRC ancestry.

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

For each of same, cross, and full, the two-dimensional spectrum has census

```text
(total, DC, F_7 axis, F_13 axis, mixed)=(91,1,6,12,72).
```

Double-centering each full output table gives `(72,0,0,0,72)`, a stronger
check than the candidate's required full-table control.  At the fixed residue
class `(1,0,6)`, all seven septimal Fourier modes are nonzero in all three
sectors.  The full-sector value is the repeated nonzero reduction
`289814661037836286866`.

Rows 1 through 6 of each cell table reduce to zero in the certified split
field.  The audit records these only as mod-`p` zero reductions.  A single
prime cannot certify that the corresponding characteristic-zero coordinates
are zero; only the seven nonzero Fourier reductions lift to nonvanishing.

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

The useful next map would have to lift this residual coordinate through a
typed collision record and compare it with the already audited source-cell
tensor without discarding either chronology or exact address.

## Reproducibility

Normal and optimized transcripts are byte-identical.  The semantic digest is
`dac1a968808aaf3bf5c1f2208f62fd1c68e55b4e17af2e12aa65d8a9809a969e`,
and the pinned script LF digest is
`ede780b135f4032be364b43a9543ede259686d1b298256a258a666cdcae083f2`.
