# Primitive Gaussian content and the spinor-surface frontier

**Research synthesis, 2026-08-12.**  Status labels in this note are local and
literal.  THM-3336 is `PROVED + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED`; the congruence-surface family below is `FINITE-EXACT` at five
primes and remains a theorem candidate.  Nothing here proves LRC(14).

## Portfolio and inheritance

The session kept three lanes live.

| Lane | Starting object | Closest proved mechanism | Load-bearing loss |
|---|---|---|---|
| Anchor | LRC(14) determinant gate | THM-2056 Kelvin/Farey certificate and THM-2596 Gram-owner cocycle | phase, owner switches, saturation, and global exit |
| Niche | multiplication of primitive Pythagorean spinors | THM-3333 Gaussian square, THM-3334 fixed-norm fibre, and THM-3341's state-dependent self-square transplant | Gaussian content and odd/odd parity content |
| Wildcard | reduction of the Lorentz/Farey carrier modulo odd primes | THM-2626's 84-state frame bundle and THM-3088's retained radial fibre | integral height, Kelvin denominator, and physical column incidence |

The canonical hostile inherited from THM-3333 is that primitive unordered
triples collapse distinct spinor gauges.  The corrected near miss is
THM-2596's endpoint-scalar recursion: it becomes exact only after restoring
norm/Gram and owner data.  The least-used sidecar is THM-3088's
`F_p^*/{+-1}` coordinate inside a projective line; it turns out to be exactly
the radial frame needed by the wildcard surface.

The live concept board was:

```text
Gaussian multiplication; primitive content; Farey faces; Lorentz/radius
polarization; fixed-norm Boolean fibres; LRC determinant gates; modular
frame surfaces.
```

The strongest pull came from applying an operation, rather than looking for
another scalar identity.  Raw multiplication preserves too much and is only
a similarity.  Primitive reduction is where new arithmetic appears.

## 1. The operation lane: from sums of squares to a content-curved face

For primitive `s=a+ib`, left multiplication is the integral similarity

```text
G_s=[a -b],             det G_s=a^2+b^2=N.
    [b  a]
```

After dividing `G_su` by its coordinate gcd `d_s(u)`, determinants become

```text
det(mu_s(u),mu_s(v))
  =N det(u,v)/(d_s(u)d_s(v)).
```

The key insight is that `d_s` is not noise.  A Bezout reduction produces a
root `h^2=-1 (mod N)` and the exact formula

```text
d_s(m,n)=gcd(m+hn,N).
```

So a representation of `N` as a sum of two squares selects one modular cusp
linear form.  This is the operation-level version of THM-3334's factor-choice
fibre: conjugation changes `h` to `-h`; coprime multiplication glues roots by
CRT; a shared split prime with opposite choices collapses into content.

On a Farey triangle `(u,v,u+v)`, the three contents are pairwise coprime.
Writing `kappa=N/(d_ud_vd_w)`, primitive normalization gives a weighted
mediant and opposite determinant labels

```text
d_w W=d_u U+d_v V,
(ell_u,ell_v,ell_w)=kappa(d_u,d_v,d_w).
```

The apparently lost multiplier norm is then reconstructed from one image
face:

```text
kappa=gcd(ell_u,ell_v,ell_w),
N=lcm(ell_u,ell_v,ell_w).
```

This changes the interpretation of the Farey graph.  Gaussian multiplication
is not a graph endomorphism except for units.  It is a rational
commensurator overlay whose edge labels are primitive cusp-width/intersection
invoices.  The abstract triangulation forgets the overlay; the labels recover
it.

### Radius packaging

The same content is visible without abandoning the triangular-number/radius
lane.  Package the inradius and three exradii into

```text
H(A,B,C)=1/2 [ A+B-C   A-B+C]
                 [-A+B+C  A+B+C].
```

Then

```text
det H=-<X,X>_L/2,
H(Phi(G_su))=G_s H(Phi(u)) G_s.
```

The absent transpose is real: this is a two-sided Gaussian action, not a
congruence.  Pairwise determinant polarization of the `H` matrices recovers
the three edge labels, hence `N`.  Thus triangular in/exradii do not merely
decorate a right triangle; as a four-coordinate matrix they can reconstruct
the index of a Gaussian-transformed Farey face.

### Generalization boundary

For an arbitrary primitive-entry integer matrix `A`, Smith normal form still
gives the content cocycle, pairwise-coprime face contents, weighted mediant,
and determinant/index reconstruction with `|det A|`.  It does not give a norm
similarity.  The counterexample

```text
A=diag(1,2), u=(1,1):  ||Au||^2=5 != 4=det(A)||u||^2
```

separates the general cyclic-lattice mechanism from the specifically
Gaussian sum-of-two-squares structure.  Brahmagupta composition, the root of
minus one, the radius covariance, and Kelvin-gate transport belong only to
the latter.

## 2. Exact leverage and exact obstruction for LRC(14)

For a primitive direction `u` and labelled columns `c_i`, primitive Gaussian
reduction transforms the THM-2056 ratio as

```text
|det(mu_s(u),mu_s(c_i))|/||mu_s(u)||^2
 =(d_s(u)/d_s(c_i)) |det(u,c_i)|/||u||^2.
```

This is an exact new coordinate for deck operations: columns whose modular
cusp forms absorb different parts of `N` receive different weights.  The
companion verifies both directions of certificate reversal on positive,
saturated, distinct 13-column decks, first at `N=2` and then at `N=101`.
At `N=2`, the norm changes `130 -> 65` in one direction and `65 -> 130` in
the other while `D=1`; at `N=101`, it changes `101 -> 1` and back.

The surprising hostile is that the exact lonely-runner margins move opposite
the certificate in all four controls.  For example, the `N=2` pass-to-fail
deck has exact margin `1/4 -> 1/2`.  This does not contradict THM-2056: its
gate is sufficient, not necessary.  It sharply rules out interpreting the
content-weighted gate as a monotone safety statistic.

There is a second typing boundary.  Multiplying every column by `G_s` without
normalizing is only a change of basis of the same rational plane; parameters
transform contragrediently.  Dividing columns by different contents generally
creates a different plane.  The lawful controls therefore compare different
saturated decks, not two presentations of one LRC instance.

The useful next experiment is correspondingly narrow:

1. start with each finite residual THM-2056 fan cone;
2. enumerate Gaussian multipliers whose content pattern preserves saturation;
3. retain the full vector `(d_1,...,d_13)`, the Gram owner, and the clock
   packet rather than only the new determinant ratio;
4. test whether any operation maps a failed cone into an already discharged
   labelled cone with a lawful inverse map on the physical speed row.

Until step 4 succeeds, content curvature is an exact reparameterization and
obstruction, not an LRC closure.

## 3. Wildcard: a finite spinor surface hiding in the 84-state carrier

The companion
`04-computation/gaussian_spinor_congruence_surface_probe_20260812.py` makes a
separate finite-exact observation.  For odd prime `p`, take

```text
V_p=(F_p^2-{0})/{+-1},
[u]~[v] iff det(u,v)^2=1,
```

and fill every graph triangle.  Exact enumeration gives:

| `p` | `V,E,F` | genus | lifts of base projective triangles |
|---:|---:|---:|---:|
| 3 | `4,6,4` | 0 | `1:4` |
| 5 | `12,30,20` | 0 | `0:10, 2:10` |
| 7 | `24,84,56` | 3 | `1:56` |
| 11 | `60,330,220` | 26 | `1:220` |
| 13 | `84,546,364` | 50 | `0:182, 2:182` |

For all five primes the probe verifies degree `p`, two faces at each edge,
each vertex link a single `C_p`, connectedness, orientability, and the
Gaussian-square Lorentz identity.  Projection to `P^1(F_p)` is a
`(p-1)/2`-sheet graph cover of `K_(p+1)`, not a surface/flag cover.

At `p=13`, the vertex set is not new: it is exactly THM-2626's
`14 x 6` projective frame bundle.  The new finite-exact structure is the edge
and face incidence.  In its affine frame coordinates,

```text
(infinity,eta)~(x,theta) iff eta=theta,
(x,eta)~(y,theta) iff (x-y)^2=eta theta.
```

The carrier partitions into six vertex-disjoint wheel disks `W_14`; their 78
faces are exactly the 78 directed Paley arcs.  Forgetting the radial frame
projects to `K_14` and creates 182 fake triangles.  For primes `1 mod 4`, a
quadratic-character two-graph selects the base triangles with two lifts; for
primes `3 mod 4`, the character is oriented and every base triangle has one
lift.  These are different predicates and must not be fused into one
tournament statement.

### Why this does not yet transfer to LRC

Three exact hostiles stop the tempting modular shortcut:

```text
(1,0) and (1,1183) are the same mod-13 state, but the one-tail
THM-2056 gate respectively fails and passes;

(5,1) has norm 26=0 mod 13, so the Kelvin denominator is singular;

the standard one-tail column deck reduces to eleven repeated projective
points, one exceptional point, and one zero column, not a 13-cycle link.
```

The surface may still organize finite currents.  At `p=13`, its graph cycle
rank is `463`, compact surface homology has rank `100`, and puncturing all 84
vertices separates `83` cusp currents from `100` compact currents.  But a
physical map would need an owner/height-coloured lift; residue alone cannot
carry the Kelvin predicate.

The cheapest next test is therefore not another genus computation.  On live
residual LRC rows, require all 13 relation columns to be nonzero modulo 13 and
to occupy the 13 distinct neighbor slopes of one framed owner; only then
compute the induced link current, together with integral norm and owner
sidecars.  Failure of this incidence predicate is a stopping reason, not a
negative LRC verdict.

## 4. Updated concept board

The session changed every lane:

- **Gaussian multiplication** is now an exact content-labelled operation,
  not a presumed primitive-triple monoid.
- **Farey topology** survives as an abstract triangulation but acquires index
  labels under primitive cyclic embeddings.
- **Triangular radii** become a rank-one matrix representation that reads
  those labels by determinant polarization.
- **Fixed-norm Boolean fibres** select roots of minus one and predict where
  shared primes cancel under composition.
- **LRC gates** are provably content-sensitive and nonmonotone under changes
  of primitive deck.
- **Modular frames** support a genuine finite surface, but only after retaining
  the radial sheet that projectivization erases.

The frontier lesson is operational: sums of two squares become useful when
their representations are treated as maps with contents, roots, and covering
indices.  Norm-only coincidences are too compressed.  The most promising
hard-frontier move is to combine the content vector with the already-required
Gram owner and clock/phase packet, and demand an exact predicate-preserving
map on a live residual row before drawing a safety conclusion.

## Reproduction

```bash
python3 04-computation/primitive_gaussian_content_curvature_thm3336.py
python3 -O 04-computation/primitive_gaussian_content_curvature_thm3336.py
python3 04-computation/gaussian_spinor_congruence_surface_probe_20260812.py
python3 -O 04-computation/gaussian_spinor_congruence_surface_probe_20260812.py
```

The normal and optimized transcripts must respectively byte-match:

```text
05-knowledge/results/primitive_gaussian_content_curvature_thm3336.out
05-knowledge/results/gaussian_spinor_congruence_surface_probe_20260812.out
```
