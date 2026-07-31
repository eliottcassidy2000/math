# THM-2868 to THM-2861: the recombined frequency edge

Status: `FINITE-EXACT SCRATCH RESULT`; not canon and no LRC(14) row
exclusion.

THM-2868 supplies, on one fixed signed q3 endpoint selector,

```text
U_(r+1)=omega^3 U_r,        V_(r+1)=V_r.
```

Recombine the two Prony branches:

```text
S_r=U_r+V_r.
```

Projectively put `t_r=U_r/V_r`, so `S_r/V_r=1+t_r` and
`t_(r+1)=omega^3 t_r`.  The oriented Hermitian edge

```text
E_r=(1+t_(r+1))(1+t_r^(-1))
```

then satisfies, identically,

```text
E_r=omega^3 conjugate(E_r).
```

Writing the THM-2861 scalar as

```text
c_r=A-B omega^(3r),
```

the exact conductor normalization gives

```text
t_0=-B/A=xi^955,       1+t_r=c_r/A.
```

Moreover `V_0=P A` and `S_r=P c_r`.  Since
`A conjugate(A)=1`, the projectively normalized edge is exactly

```text
c_(r+1) conjugate(c_r),
```

not merely another edge of the same phase type.  The full signed-current
edge is its fixed nonzero multiple `P conjugate(P)`.  This reconstructs the
THM-2861 edge from THM-2868's 26 lawful raw multiplier samples and local
Prony transport.  Its cyclic Fourier support is exactly `{0,3,10}`, and all
thirteen values are nonzero and distinct.

There is a small additional simplification.  The phase line is already
known and `1+omega^3` is nonzero, so

```text
E_r = omega^3/(1+omega^3) (E_r+conjugate(E_r)).
```

Thus the coefficient-field symmetric trace determines the oriented edge
after retaining the separately proved forward frequency orientation.  The
hostile is exact: reversal has the same trace.  No physical ordinary-union
or polarization measurement is constructed, so this identity does not
bypass THM-2380's realization gate or itself manufacture the orientation.

What this repairs is the *frequency-dual coefficient* adjacency: the two
Prony summands are reconstructed from the same signed q3 selector and
recombine into the full three-channel detector.  What it does not repair is
the type change from multiplier section to ancestry section.  The variable
local charts are not one positive packet or one marked triangle, q11
cancels, q7 has no E3 support, and no physical ancestry co-shift follows.

Reproduce:

```text
python3 .scratch/lrc_frequency_hermitian_bridge_20260728/probe.py
python3 -O .scratch/lrc_frequency_hermitian_bridge_20260728/probe.py
```
