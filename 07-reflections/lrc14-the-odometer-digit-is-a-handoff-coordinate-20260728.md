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
positive flip window for nearest-clock values *along the affine orbit*.  But
this is not yet the stored event edge: the intrinsic shallow and owner are
`c_7(Dx)` and `c_7(D^2x)`, and a handoff must identify the owner at `x` with
the shallow at `T_kx`.

This distinction kills the first attractive cycle.  The valid lifts

```text
k_0=-14,                 k_1=14
```

exchange `1/2+1/13^6` and `1/2-1/13^6`, but their stored edges are `4->4`
and `3->3`, and neither owner equals the following shallow clock.  The full
rail bank is positive at both points, while the THM-2640 present packet is
already zero there because the speed-66 safe factor excludes a neighbourhood
of `1/2`.  This is a useful hostile, not an escape.

The arithmetic also gives the repair.  Put `S=13^5` and use

```text
k_0=-(S+1),             k_1=S+1,
a=(S+1)/(14*13^6).
```

Then

```text
1/2+a  --T_(-(S+1))-->  1/2-a,
1/2-a  --T_(S+1)----->  1/2+a.
```

Since

```text
13a =1/14+1/(14S),
169a=13/14+13/(14S),
```

the intrinsic stored edges are exactly `4->3` and `3->4`.  Each owner is the
next event's shallow clock, so this uses the existing nonconstant edge grammar
rather than inventing a new observable.  The predecessor carries again
alternate `7,5`, both future digits are `6`, and the THM-2657 quotient labels
are root steps `11,2`.

A perturbation grows by a factor of thirteen per edge, so the exact cycle is
repelling, but for every finite horizon `H` the initial interval of radius

```text
1/(3*13^6*13^H)
```

remains inside the central tooth, keeps each stored edge nonconstant, and
preserves every owner-to-next-shallow identification.  This is positive
clock/state support at every prescribed finite depth, not merely one
exceptional point.

The new object is a controlled affine cocycle, not a single circle map:

```text
state       = (x, odometer lift k),
transition  = x -> {13x+k/13^6},
interface   = intrinsic shallow/owner clocks plus carry/future digit.
```

Projecting away `k` recovers THM-2684's nilpotent scalar picture and destroys
the alternating repair.  This is the same kind of loss as forgetting gain on
a graph edge or orientation on a knot band: the projected vertices remain,
but composition changes.

The remaining hostile is physical and sharp.  The three-tooth envelope is a
union over labelled rails; it does not say that one source/owner/deep packet
survives along the alternating affine fibre product.  Present factors are not
translation-covariant, and primitive units have not been transported.  The
next decisive test is the actual `(-(13^5+1),13^5+1)` labelled rail product,
retaining every source, carry, half-edge, present factor, and global-content
label.  The small-lift hostile shows why checking only the aggregate envelope
or the orbit clock is insufficient.
