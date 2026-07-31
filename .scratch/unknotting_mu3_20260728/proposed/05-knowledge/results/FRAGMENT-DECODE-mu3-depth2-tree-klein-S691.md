# Fragment decode: depth two is 3+1+3=7, not 9

**Status: SUPERSEDED / EXACTLY RESOLVED by
`THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy.md`.**

The fragment refers to the depth-two inverse-image tree of the sporadic
Keller map

```text
u=1+xy,
F=(u^3z+y^2u(4+3xy),
   y+3xu^2z+3xy^2(4+3xy),
   2x-3x^2y-x^3z),
det J_F=-2.
```

For `v*=(-1/4,0,0)`,

```text
F^-1(v*) = {
  P-=(-1,3/2,13/2),
  P0=(0,0,-1/4),
  P+=(1,-3/2,13/2)
}.
```

At the next level:

```text
|F^-1(P-)|=3,
F^-1(P0)={(-1/8,0,0)},
|F^-1(P+)|=3,
```

so the exact finite count is `3+1+3=7`.

The missing two points do **not** collide.  Since `det J_F=-2` everywhere,
finite ramification is impossible.  The core eliminant

```text
E(x)=Lx^3+(4-3bc)x-2c
```

has `L(P0)=0` and specializes to `4x+1/2`; two sheets escape to infinity
on the Jelonek nonproperness surface while the unique finite survivor is
`x=-1/8`.

The symmetry fixing `P0` and swapping `P-` with `P+` is the order-two
involution `sigma=diag(-1,-1,1)`.  Generic monodromy is `S3`, not `mu3`.
Thus “mu3-fixed” is, at best, informal ternary imagery and must not be used
as the exact group action.

Exact reproduction:

```bash
python3 .scratch/unknotting_mu3_20260728/keller_depth2_exact.py
```

The former preferred Arithmetic-Kakeya/Ward-recurrence readings and the
claim that “two preimages coincide” are retracted.  Any independent
identification-gluing idea in the AK workbench may remain a hypothesis,
but this fragment supplies no evidence for it.
