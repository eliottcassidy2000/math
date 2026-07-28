# The odometer digit is a handoff coordinate

THM-2684's nilpotence is exact, but it is nilpotence of the scalar dilation
`D(x)={13x}`.  It should not be read as nilpotence of every chronology that
projects to dilation.

THM-2657 supplies the missing coordinate.  Its lawful carry/root lifts are
translations `k/13^6`, with `k` a unit modulo thirteen.  Compose one after
dilation:

```text
T_k(x)={13x+k/13^6}.
```

Near the central tooth, write `x=1/2+u`.  For a fixed lift `k`, with
`tau=k/13^6`, the first two displacements are

```text
u_1=13u+tau,
u_2=169u+14tau.
```

Their zeroes differ by `|tau|/169`.  Thus every nonzero fixed lift opens a
positive three-event clock-flip window that scalar `D` did not have.  But the
next zero lies on the same side as the second: a fixed affine lift crosses the
central clock wall at most once.  It cannot support two consecutive clock
flips.

The decisive operation is therefore not “add a translation”; it is “retain
the translation digit on every edge.”  The two valid lifts

```text
k_0=-14,                 k_1=14
```

give the exact cycle

```text
1/2+1/13^6  --T_-14-->  1/2-1/13^6
1/2-1/13^6  --T_+14-->  1/2+1/13^6.
```

The nearest-seven clocks alternate `4,3`; the predecessor carries alternate
`7,5`; both future digits are `6`.  Relative to THM-2657's quotient map, the
two lift labels are root steps `11,2`.  A perturbation grows by a factor of
thirteen per edge, so the exact cycle is repelling, but for every finite
horizon `H` the initial interval of radius

```text
1/(3*13^6*13^H)
```

remains inside the central tooth and keeps every consecutive clock distinct.
This is positive support at every prescribed finite depth, not merely one
exceptional point.

The new object is a controlled affine cocycle, not a single circle map:

```text
state       = (x, odometer lift k),
transition  = x -> {13x+k/13^6},
observable  = nearest-seven clock plus base-thirteen carry/future digit.
```

Projecting away `k` recovers THM-2684's nilpotent scalar picture and destroys
the alternating repair.  This is the same kind of loss as forgetting gain on
a graph edge or orientation on a knot band: the projected vertices remain,
but composition changes.

The remaining hostile is physical and sharp.  The three-tooth envelope is a
union over labelled rails; it does not say that one source/owner/deep packet
survives along the alternating affine fibre product.  Present factors are not
translation-covariant, and primitive units have not been transported.  The
next decisive test is therefore the actual `(-14,+14)` labelled rail product,
beginning with the THM-2640 `h=10,c=0,r=2` atom and retaining every source,
carry, half-edge, and global-content label.
