# LRC14 Exposure-Poset Creative Pass

This pass reframes the proof search as an exposure problem.  A strict LRC14
counterexample would not merely be a row with low scalar value.  It would be a
source packet that remains invisible to every channel the repo has built:
qdiv, AP/GW boundary, Haar open mass, endpoint bridges, Fejer/Toeplitz duals,
K33/C27 state labels, twist/Ramanujan witnesses, and moment barriers.

The new computation audits certificate channels rather than runner sets.  It
checks `12015` bounded AP-neighborhood and hard-frontier rows with exact
safe-component data and Fejer-vector PSD tests.  The output is simple:

```text
zero-safe rows                2: AP, GW 12->24
positive-safe rows            12013
positive rows with Fejer dual 12013
unexposed rows                0
```

The hard q>=14 front is familiar: `P10+GW`, `drop(6)->63`, K33 `12->36`,
`P10+K33`, `drop(6)->86`, and nearby two-swap hard faces.  This is useful
because it turns the "try another proof" request into a concrete theorem
target: prove no hidden exposure kernel exists after labels are retained.

## What Changed

Earlier route-atlas and source-packet synthesis said that the missing object
is not a scalar.  HYP-2988 makes the statement operational:

```text
row -> labelled source packet -> exposure channel set.
```

If the channel set is empty, we have a genuine new F7 sector.  If it is not
empty, the proof obligation is to formalize the channel.  In this pass, every
audited positive row has both `OPEN_HAAR_BRIDGE` and `FEJER_PSD_DUAL`; the
remaining issue is not discovery of a witness but familywise certification.

This shifts the summit target from:

```text
find an invariant that separates all rows
```

to:

```text
prove exposure-channel completeness for labelled packets.
```

## Tournament Analysis

Vertices are exposure channels.  I explicitly considered runners, gaps, fixed
sections, wall crossings, residues, Fourier modes, endpoint events, packet
families, and proof obligations.  The chosen quotient preserves q-threshold,
open/boundary status, safe-component size, Fejer degree and margin, packet
family, and state labels, while destroying full endpoint-owner incidence.

The channel tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
```

The leading channels by this finite exposure gauge are open-Haar bridge and
Fejer PSD dual.  AP/GW and K33 are low-volume channels but remain load-bearing:
they are not weak because they score low; they score low because they are
small exceptional exits.

## No-Hidden-Kernel Lemma

The clean theorem target is:

```text
Every primitive q>=14 non-AP/GW source packet either has positive open Haar
exposure with a familywise Fejer interval certificate, or carries a named
C27/K33 state-lift label.
```

HYP-2981 says the numerical interval work is not huge: even the smallest
selected Fejer margin only asks for a finite atom-L1 precision budget.  HYP-2988
says where to attach those certificates: not to an unstructured row list, but
to the exposure kernel of the source-packet sheaf.

After the rebase, the exposure kernel should be read as the common residual of
the newer analytic and deformation lanes too.  HYP-2982 warns that
squarefree/inverse-unit analytic weights forget prime-power and endpoint
labels; HYP-2983 says smoothing choices must carry Kaczynski boundary approach
classes; HYP-2984's kernel-homotopy lane requires preserved packet
certificates or boundary-defect atoms; HYP-2984's Farey scheduler routes
`e=14p-q` and product-collapse `p` into AP/GW, C27 petal, or K33/state-lift
families; HYP-2985's admissible-smoothing dispatcher names which smoothing
clock may certify each packet family; and the incoming HYP-2986
tope-wall lane plus the HYP-2987 certificate-handoff atlas add endpoint-cell
and zipper-arrow versions of the same obstruction.  HYP-2988 is therefore the audit layer asking
whether any packet survives all of those labelled exits simultaneously.

## Direct Falsifier

The first serious new object would be a row with:

```text
qdiv > 14,
safe_mu = 0,
not AP/GW,
no K33/C27/petal label,
no Toeplitz/Fejer negativity,
no moment/twist/Ramanujan handoff.
```

No such row appeared in this bounded pass.  If one appears later, it should be
promoted immediately to the F7 Johnson-harmonic/source-spectrum sector rather
than treated as another scalar anomaly.
