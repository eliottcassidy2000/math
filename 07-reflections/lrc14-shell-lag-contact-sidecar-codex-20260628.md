# Shell magic sees where the long gaps sit

HYP-3245 gave a clean slogan:

```text
AP = triangular lag law
traps = outward lag transport
```

That is true, but only as a first projection.  HYP-3246 says the shell packet
still remembers something ordinary lag transport forgets.

The bounded-bank collision is exact and hard to unsee:

```text
E_A = (0,1,2,3,4,12,13,14)
E_B = (0,1,2,10,11,12,13,14)
```

They have the same ordinary support autocorrelation.  They also have the same
residue histogram and even the same residue word mod `7`,
`(0,1,2,3,4,5,6,0)`.  So the naive story would be that shell magic should not
distinguish them.

But it does:

```text
magic(E_A) = 211/98
magic(E_B) = 71394/35035
```

The only visible combinatorial difference is where the long gap `8` sits in
the anchored gap word:

```text
E_A gaps = (1,1,1,1,8,1,1)
E_B gaps = (1,1,8,1,1,1,1)
```

That is the whole lesson.  Shell magic does not only see how much pair mass
moved to large lags; it sees where the contact defect sits in the ordered
endpoint word.

This turns HYP-3245's open warning into an exact controlled-forgetting rule.
Ordinary support autocorrelation is not a legal terminal quotient for
HYP-3228.  Residue histogram mod `7` is not a legal terminal quotient either.
The missing finite repair is an ordered contact-support sidecar.

The bounded counts make the point sharper:

```text
ordinary lag profile              -> 1677 mixed shell fibers
lag profile + residue histogram   -> 62 mixed shell fibers
lag profile + residue histogram + contact support -> 0 mixed shell fibers
```

The striking part is what does **not** help.  Gap multiset information alone
does nothing to those last `62` fibers.  So the hidden coordinate is not the
size data of the long gaps.  It is their placement.

That connects three threads that looked merely adjacent before:

- HYP-3204's ordered-tail exchange pricing
- HYP-3243's circle endpoint arrangement cells
- HYP-3244's lift/compress descent sidecars

All three were really saying the same thing.  A scalar summary of a row must
remember an ordered contact coordinate if it wants to talk to the shell packet.

So the better formulation is no longer

```text
shell magic is a lag-space scalar.
```

It is

```text
shell magic = lag transport + ordered contact sidecar.
```

The next useful step is to compare that contact sidecar directly against the
finite endpoint cells / trap chambers from HYP-3243 and the tiling/half-tiling
descent sidecars from HYP-3244, instead of asking lag space to do a job it is
not allowed to do alone.
