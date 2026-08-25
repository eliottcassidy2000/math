# LRC(14) `d=2` two-phase all-height probe

**Status:** `PROMOTED TO THM-4049 + FINITE-EXACT + VERIFIED-EXACT`. This is
the discovery report for
[THM-4049](../../01-canon/theorems/THM-4049-lrc14-d2-two-phase-residue-firewall.md).
The theorem closes an explicit residue-defined family, not the physical
`THM-3818` `d=2` branch and not LRC(14).

## Inheritance pass

- **Closest proved mechanism:**
  [THM-4041, `d=2` affine defect-edge boundary](../../01-canon/theorems/THM-4041-lrc14-d2-affine-defect-edge-boundary.md),
  especially its exact condition `G(H) subset Sigma_(alpha,beta)` and its
  retention of the two common lift labels.  The common-label ancestor is
  [THM-4004, three-detuned divisor-comb profile](../../01-canon/theorems/THM-4004-lrc14-three-detuned-divisor-comb-profile.md).
- **Canonical hostile:** THM-4041's typed row with
  `H={1,...,10,12}`, exceptions `(1,11)`, and pack phase `1/11`.  That
  particular safe pack phase is fully spoiled although the row is lonely.
  It forbids choosing an arbitrary cited pack witness.
- **Corrected near miss:** THM-4004/4041 correct the formerly advertised
  minimal two-lift hostile `(3,7)` to `(3,5)`.  More importantly for this
  probe, the finite bank through exception height `79` cannot be extrapolated
  merely because no two-endpoint hostile occurred there.
- **Least-used relevant sidecar:** the exact control-bank audit
  `.scratch/lrc14_d234_pack_safe_intersection_20260824/REPORT.md` found the two
  closed pack-safe phases `{1/14,5/56}` for `d=2`, but used them only on 759
  certificate-positive odd pairs through height `79`.  The present probe
  attacks its all-height extension directly.

The meta-patterns used were **Attack a proposed bound before extending it**,
**Respect symmetries by searching orbit representatives**, and **Verify the
consequence, not the model's own assumptions**.  The height extension is a
complete residue problem modulo `112`, and the companion prints full-row
clearances at the actual lift times.

## Anchor / Niche / Wildcard and concept board

- **Anchor:** decide the physical `d=2` inclusion
  `G(H) subset Sigma_(alpha,beta)` without forgetting pack phase.
- **Niche:** dualize from describing the whole spoiled set to finding a small
  closed hitting bank `T subset G(H)` that no odd exception pair can spoil.
- **Wildcard:** quotient the four lift times by their common denominator and
  ask whether the height-`79` observation is actually a finite mask identity.

| live concept | object / representation | predicate | operation / invariant | loss / obstruction | cheapest test |
|---|---|---|---|---|---|
| affine defect edge | open set `Sigma_(alpha,beta)` | can both labels be spoiled? | odd pair, half-shift | forgets `G(H)` | intersect with named phases |
| divided pack | closed set `G(H)` | is a named pack phase safe? | evaluate residues mod `56` | named phases need not be pack-safe | exact residue filter |
| dual phase bank | four lift times over two `y` values | does some lift survive? | complement/hitting dual | forgets all other phases | enumerate all odd masks |
| periodic mask carrier | odd residues mod `112` | can one pair cover both banks? | quotient by common denominator | forgets speed height, harmless here | all 1,596 residue pairs |
| physical producer | filtered tuple `(s,t,p,q,u)` | does its projected `H` enter the residue class? | projection to `(H,E)` | the proved `91^12` box is not tractably enumerated | execute projection or sharpen the bound |

The new mask identity strengthens the dual-bank and periodic-carrier items.
It does not solve the physical-producer item: a physical divided pack need
not lie in the residue class below.

## Theorem signal: a four-time residue firewall

Put

```text
A = (Z/56Z) minus {0,11,14,22,23,28,33,34,42,45}.     (1)
```

Let `H` be any finite set of integers whose residues modulo `56` lie in `A`,
and let `alpha,beta` be odd integers.  Then one of

```text
x in {1/28, 15/28, 5/112, 61/112}                    (2)
```

satisfies

```text
||2h x|| >= 1/14  for every h in H,
||alpha x|| >= 1/14,  ||beta x|| >= 1/14.             (3)
```

Consequently every distinct-positive-speed family

```text
2H union {alpha,beta}                                 (4)
```

of this type is `1/14`-lonely.  In particular this holds for the eleven-pack

```text
H_0={1,2,...,10,12}                                   (5)
```

and **all** odd exception pairs, with no exception-height bound and no affine
certificate hypothesis.

## Proof

Write

```text
y_1=1/14,   y_2=5/56,
x_(i,j)=(y_i+j)/2,   j=0,1.                           (6)
```

For every `h`, both labels preserve the pack phase:

```text
||2h x_(i,j)||=||h(y_i+j)||=||h y_i||.                (7)
```

At `y_1`, pack safety is equivalent to `14` not dividing `h`.  At `y_2`, it
is equivalent to the circular residue of `5h mod 56` having magnitude at
least `4`.  Combining these two conditions gives exactly the ten forbidden
residues in `(1)`.  Thus every `h in H` is safe at all four times.

An odd speed sees the two labels in one bank at phases differing by `1/2`, so
it cannot be dangerous at both.  At the first bank,

```text
alpha dangerous at x_(1,0)  iff alpha == +/-1  (mod 28),
alpha dangerous at x_(1,1)  iff alpha == +/-15 (mod 28). (8)
```

If the pair does not cover both first-bank labels, a first-bank lift proves
`(3)`.  If it does cover them, one exception, call it `w`, must satisfy
`w == +/-15 (mod 28)`.  At the second bank the lift numerators over denominator
`112` are `5` and `61`.  Since `61 == 5 (mod 28)`,

```text
5w == 61w == +/-19 (mod 28).                          (9)
```

Danger at either second-bank lift would require the numerator modulo `112`
to have circular magnitude less than `8`; because it is odd, its reduction
modulo `28` would lie in

```text
{+/-1,+/-3,+/-5,+/-7},                               (10)
```

contradicting `(9)`.  Therefore `w` is safe at both second-bank labels.  The
other odd exception cannot cover both labels by itself, so one second-bank
label survives.  Equation `(7)` then proves `(3)`.

This argument is all-height.  The residue audit is a redundant exact check,
not the source of the universal quantifier.

## Boundary, hostiles, and scope

- The pack inequality is closed (`>=1/14`) and exception danger is open
  (`<1/14`).  The pack `H_0` reaches equality through `h=12` at `y_2`; this is
  a positive boundary control, not a discarded wall.
- Deleting the `y_2` pack condition is genuinely unsafe for this four-time
  bank: `H=(11)`, exceptions `(1,15)` make all four displayed times fail.
- Deleting the `y_1` pack condition is likewise unsafe: `H=(14)`, exceptions
  `(1,11)` make all four displayed times fail.
- Those two rows are hostiles to the **fixed four-time certificate only**;
  they are not Lonely Runner counterexamples.
- Oddness is load-bearing.  It makes the two exception labels antipodal and
  is exactly the physical parity supplied by the `d=2,c_2=9` boundary.  With
  the valid canonical pack `H_0`, the distinct even exceptions `(22,28)`
  defeat all four displayed times; this is a parity hostile to the fixed
  certificate, not an LRC counterexample.
- The result does not prove that an arbitrary physical divided eleven-pack
  satisfies `(1)`.  It therefore does not close THM-4041's whole physical
  image, any other `11+2` branch, or LRC(14).

## Connection contract

```text
source:    a d=2 row 2H union {alpha,beta}
target:    two labelled danger masks on each of y_1,y_2
map:       evaluate the four lift times and reduce speeds mod 112
preserves: existence of a displayed 1/14-safe full-row phase
destroys:  every phase outside the four-time bank and the physical tuple origin
sidecar:   H mod 56 lies in A; alpha,beta are odd
test:      all 1,596 odd residue pairs with repetition, by Fraction and
           independent integer-residue masks
```

The next direct test is not to enlarge the odd-pair height again.  It is to
project a genuinely physical THM-3818 `d=2` producer onto the divided pack's
modulo-`56` residues and measure how often `(1)` holds. THM-3818 already
supplies a finite `91^12` box; the missing work is an executed projection
census or a sharper bound that makes that census tractable.

## Exact artifact

Run from the repository root:

```text
python3 -B 04-computation/lrc14_d2_two_phase_all_height_probe_20260824.py
python3 -B -O 04-computation/lrc14_d2_two_phase_all_height_probe_20260824.py
```

Both modes byte-match
`05-knowledge/results/lrc14_d2_two_phase_all_height_probe_20260824.out`.
The frozen transcript retains its pre-promotion `NOT CANON` provenance line;
the audited mathematical statement is now canonical in THM-4049.
The companion checks all `56` odd residues modulo `112`, all `1,596`
unordered pairs with repetition, direct Fraction masks against an independent
integer-residue implementation, the exact pack residue firewall, five
positive controls (including an exception beyond height `79`), the two
pack-hypothesis hostiles, and the even-exception parity hostile.  It prints
the actual full-row witness and clearance.

```text
script SHA-256: c24612429463a4d0caf88ccaf5d301a186b0d2346bd284141920aeb52e0699dd
output SHA-256: 10d8f230a1212e4b6a6e061fb7866b6cd12335f5c530f21989b70e9365a56adf
semantic digest: 364f63bbeec172c18250e856905f376a47264a43569d0206972684c0e6590536
hash basis: raw LF bytes
```
