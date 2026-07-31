# Fragment decode: the depth-two tree is `3+1+3=7`, not 9

> **CORRECTION (root, 2026-07-28).**  The earlier note preferred an
> Arithmetic-Kakeya interpretation and described a `mu_3`-fixed collision.
> Exact term matching and reconstruction identify the fragment with
> THM-2473's sporadic Keller-map inverse tree.  No finite points collide:
> two sheets escape to infinity.  The fixing symmetry is an involution and
> the generic monodromy is `S3`, not `mu_3`.  The former AK transfer is
> withdrawn.

The truncated owner fragment was:

> The depth-2 tree is 3+3+1 = 7, not 9: the mu3-fixed branch P0 has a
> degenerate fiber — two p...

Its exact host is the sporadic Keller map of
`01-canon/theorems/THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy.md`:

```text
u=1+xy
F=(u^3 z+y^2 u(4+3xy),
   y+3xu^2z+3xy^2(4+3xy),
   2x-3x^2y-x^3z),
det J_F=-2.
```

At `v*=(-1/4,0,0)`, the first fibre is

```text
P-=(-1,3/2,13/2),  P0=(0,0,-1/4),  P+=(1,-3/2,13/2).
```

The three next fibres have sizes

```text
P-: 21119x^3-404x-208                 3 distinct points
P0: z=0, y=0, x=-1/8                  1 finite point
P+: 20929x^3+532x-208                 3 distinct points
```

so the positional count is

```text
3+1+3=7.
```

The core eliminant is

```text
E(x)=Lx^3+(4-3bc)x-2c.
```

At `P0`, its leading coefficient vanishes and the polynomial becomes
`4x+1/2`.  The missing two sheets escape along the Jelonek nonproperness
surface.  Since `det J_F=-2`, the map is etale and a finite
ramification/collision explanation is impossible.

The exact symmetry is

```text
sigma=diag(-1,-1,1):
sigma(P0)=P0,  sigma(P-)=P+,  sigma(P+)=P-.
```

Thus the branch stabilizer visible here is `C2`; THM-2473 proves generic
monodromy `S3`.  A `mu_3` label does not preserve the actual fibre data.

## Consequence for Arithmetic Kakeya

This fragment supplies no identification-gluing certificate for the
Arithmetic-Kakeya forcing game.  The AK idea remains an independently
testable construction principle, but the Keller degeneration destroys
finite support by escape to infinity rather than merging two vertices.
Any future transfer must specify a source/target map and show that forcing,
legality, and certificate score survive it.

Exact reproduction:

```bash
python3 .scratch/unknotting_mu3_20260728/keller_depth2_exact.py
```

The script passes ordinary and optimized Python with byte-identical output.
