# THM-2788 modular odometer/Heisenberg hostile audit

**Verdict: PASS.  THE INCOMING UPSTREAM PROMOTION REPAIRS THE SCRATCH
EVIDENCE PACKAGE AND SHARPENS THE SCOPE CORRECTLY.**

No tracked file was edited.

## 1. Reproduction

The discovery companion

```text
.scratch/modular_odometer_heisenberg_thm2788/probe.py
```

replays under ordinary and optimized Python to the same LF SHA-256 claimed
in its report:

```text
script  bd64ea35ea093519c23b12484e1cf3857a74361d8633b541868916de62f3995b
output  2792a695832dc2c05227813bbf2e3182dce33c251856711146320f564f9ea35e
```

The scratch directory does not contain `probe.out`; its output hash is
reproducible but its stored transcript is missing.  The later upstream
canonization repaired that packaging under

```text
04-computation/lrc14_modular_odometer_heisenberg_bockstein_thm2788.py
05-knowledge/results/lrc14_modular_odometer_heisenberg_bockstein_thm2788.out
```

in `origin/main`.  In a detached worktree, ordinary and optimized Python
both byte-match the stored canonical output:

```text
canonical script
  d414bf2afb6aa3e40de9378ae20f03db1cb7bff75f59f13a60ac96e56cb95a89
canonical output
  99ad33904617d45d76a285de5467b96408dc164839cb4168905c7fe678db8f66
```

The independent audit does not import the discovery script.  Its ordinary
and optimized transcripts both byte-match:

```text
.scratch/thm2788_hostile_audit/audit.py
  64dbdd8c32821aced592bf01eb777a856cd475f91c91952bc752b1c65809ad9d
.scratch/thm2788_hostile_audit/audit.out
  0bbdb5046651c65da3d0a55ffda7060c277706c179f1f896084dd63ae617bc1c
```

Both independent scripts have zero `assert` AST nodes.

## 2. Presentations and cocycle orientation

Fix the convention

```text
[g,h]=g h g^(-1) h^(-1),
s(a,b)=Y^b X^a.
```

Write a central normal form as

```text
(a,b,c)=Z^c Y^b X^a.
```

The independently reconstructed laws are

```text
H:
(a,b,c)(a',b',c')
 =(a+a', b+b', c+c'+a b'),

M:
(a,b,c)(a',b',c')
 =(a+a', b+b', c+c'+a b'
   +floor((b+b')/p)),
```

with all coordinates reduced modulo `p`.  Direct permutation composition
with

```text
X(n)=(1+p)n,
Y_phys(n)=n+1,
Y_0(v,w)=(v+1,w),
Z(v,w)=(v,w+1)
```

matches these laws.  Therefore the report's positive cocycle signs are
correct:

```text
c_H((a,b),(a',b'))=a b',
c_M=c_H+floor((b+b')/p).
```

The carry term is the standard positive Bockstein for the digit section
`{0,...,p-1}`.  It is symmetric, so both antisymmetrizations are

```text
a b'-a' b.
```

For both groups,

```text
Z(G)=[G,G]=<Z>,
[X,Y]=Z.
```

The modular presentation is therefore

```text
M_p = C_(p^2) semidirect_(1+p) C_p,
Y_phys^p=Z,
X Y_phys X^(-1)=Y_phys^(1+p).
```

No cocycle-orientation correction is needed.

## 3. Odd-prime power and order boundary

For odd `p`, exhaustive normal-form powers give

```text
(Y_0^b X^a)^p=1,
(Y_phys^b X^a)^p=Z^b.
```

Thus

```text
H_p:  1^1, p^(p^3-1),

M_p:  1^1,
      p^(p^2-1),
      (p^2)^(p^2(p-1)).
```

The exact tested spectra are:

```text
p=3:  M 1^1 3^8 9^18;
p=5:  M 1^1 5^24 25^100;
p=7:  M 1^1 7^48 49^294;
p=11: M 1^1 11^120 121^1210;
p=13: M 1^1 13^168 169^2028.
```

The corresponding `H_p` spectra are `1^1 p^(p^3-1)`.  Hence the odd-prime
groups are nonisomorphic.  The report's count `p^2(p-1)` of order-`p^2`
elements is correct.

## 4. The binary boundary

At `p=2`, the two permutation subgroups of `S_4` are literally equal and

```text
Y_phys=(X Y_0)^(-1).
```

With rotation `r=Y_phys` and reflection `s=X`,

```text
r^4=s^2=1,       srs=r^(-1),
```

and the spectrum is

```text
1^1 2^5 4^2.
```

Thus the common group is `D_8`, not `Q_8`.

One qualification should be stated: the two displayed cocycles are not
equal in the fixed quotient basis at `p=2`.  The equality of action groups
uses the generator—and hence quotient-basis—change above.  Calling this
“the same `D_8` after generator change” is exact.

## 5. Minimal faithful degree

The lower bound proof is valid for both groups.  Every index-`p` subgroup is
normal and is the inverse image of a hyperplane in

```text
G/[G,G] isomorphic F_p^2.
```

There are `p+1` such kernels, each of size `p^2`, and every one contains
`Z=[G,G]`.  Therefore an action whose orbits all have size `1` or `p`
annihilates `Z` on every orbit and cannot be faithful.  A faithful action
must have an orbit of size at least `p^2`.

Both displayed actions on `Z/p^2` are transitive and faithful, with
point-stabilizer size `p` and trivial action kernel.  Hence

```text
mu_perm(H_p)=mu_perm(M_p)=p^2.
```

The independent script enumerates all `p+1` index-`p` kernels and the two
faithful actions for every tested prime.  No correction is needed.

## 6. Exact THM-2779/2782 scope

THM-2779 already proves:

- the carry-suppressed `H_p` action;
- `mu_perm(H_p)=p^2`;
- the `p=2` `D_8` quadratic-refinement boundary;
- the two-digit `156+13` carry-wall law at `p=13`.

THM-2788 should cite rather than re-promote those as new.  Its genuinely new
payload is:

```text
the physical group M_p and its presentation;
the explicit carry/Bockstein cocycle c_M-c_H;
the odd-prime p-power map and complete order-spectrum separation;
the p=2 equality M_2=H_2 after generator change;
mu_perm(M_p)=p^2.
```

THM-2782 constructs a physical three-point `+13` central segment, not a
physical `+1` successor cycle.  Since both extension classes share `Z`, that
segment is compatible with both and cannot detect the carry cocycle.  The
new obstruction applies only to a **full faithful successor action** that
tries to identify `Y_phys` with `Y_0` while preserving `X`, `Z`, and the
group law.  It does not rule out a nonfaithful quotient, a graded lift that
retains carry, or an unrelated physical current.

Nothing here proves that an LRC current realizes `M_13`, attaches endpoint
origins to the physical address, supplies a root-deck intertwiner, excludes
a row, or proves LRC(14).

## 7. Incoming intersection, wreath, and tower addenda

The upstream promotion adds three claims not present in the original scratch
report.  Independent controls support all three.

For every tested odd prime,

```text
H_p intersect M_p=<X,Z>,                 order p^2.
```

The exact intersection sizes are `9,25,49,121,169` for
`p=3,5,7,11,13`; at `p=2` the intersection is the whole common `D_8`.

Putting `C=O Y^(-1)`, the conjugates

```text
C_r=Y^r C Y^(-r)
```

increment the high digit only on the column `v=r`.  They commute, have order
`p`, and satisfy

```text
product_r C_r=Z,             product_r C_r^r=X.
```

Their disjoint supports make them an independent `F_p^p` column basis, so

```text
<H_p,M_p>=C_p wreath C_p,                 order p^(p+1).
```

Finally, for `p=13` and depths `k=2,...,6`, the audit checks

```text
|<O,X> on Z/13^k|=13^(2k-1),
[X,O^(13^r)]=O^(13^(r+1)),
Z(G)=<O^(13^(k-1))>,
|ker(G_(13,k)->G_(13,2))|=13^(2k-4).
```

At depth six this gives the exact invoices

```text
|G|=13^11=1792160394037,
|kernel|=13^8=815730721.
```

Therefore the upstream distinction is correct: reduction modulo `169`
retains the physical modular carry class, while replacing `O` by the
carry-suppressed `Y` is a second, coarser shadow operation.  THM-2782's
failure after `+169` is typed by the next address commutator
`[X,O^13]=O^169`; this identifies a group-coordinate source for the failure
but still does not prove factor/current covariance under that commutator.

## 8. Recommended canonical wording

Suggested promoted title:

```text
Physical modular odometer, carry-Bockstein cocycle,
and odd-prime Heisenberg separation
```

Suggested theorem digest:

> Let `p` be prime and let `X(n)=(1+p)n`, `Y_phys(n)=n+1` on
> `Z/p^2`.  Then `M_p=<X,Y_phys>` is
> `C_(p^2) semidirect_(1+p) C_p` of order `p^3`, with
> `Z=Y_phys^p` and `[X,Y_phys]=Z`.  Relative to the section
> `s(a,b)=Y^bX^a`, its extension cocycle is
> `a b'+floor((b+b')/p)`, whereas the carry-suppressed Heisenberg action of
> THM-2779 has cocycle `a b'`.  The two commutator forms agree.  For odd
> `p`, the quotient power maps are respectively `(a,b)->bZ` and zero, and
> the groups have spectra
> `1^1 p^(p^2-1) (p^2)^(p^2(p-1))` and
> `1^1 p^(p^3-1)`.  At `p=2` the two action groups are the same `D_8` after
> the generator change `Y_phys=(XY_0)^(-1)`.  The modular action also has
> minimal faithful permutation degree `p^2`.
>
> For `p=13`, the physical successor and carry-suppressed successor agree
> off the 13-point low-digit wall and differ there by `Z`.  This is an
> abstract group/action boundary.  It constructs no physical LRC current or
> endpoint/root intertwiner.

Suggested canonical evidence paths:

```text
04-computation/lrc14_modular_odometer_heisenberg_bockstein_thm2788.py
05-knowledge/results/lrc14_modular_odometer_heisenberg_bockstein_thm2788.out
```

Recommended dependencies:

```text
depends_on:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
```

The core modular-group proof is elementary and self-contained.  THM-2779 is
needed for the comparison and to avoid presenting inherited results as new;
the expanded tower/semantic interpretation uses THM-2657 and THM-2782
directly.
