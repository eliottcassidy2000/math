---
id: THM-2271
title: "Expiration support forces a weighted owner-absorber cut"
status: >
  PROVED + VERIFIED-EXACT. In every one of the 150 strict first-depth-one
  scalar profiles, THM-2263 supplies a labelled exclusive-owner stratum
  whose post-expiration support has Haar mass at least
  15041431/70270200. After adjoining the guard complement as an explicit
  absorber label, the exact global cover makes source and target service
  eligibility total. A finite measurable root-branch section then realizes
  THM-2267's singleton owner cut with support-Haar transition energy at
  least 5002831/70270200. Thus nonzero transition energy is genuinely
  forced, not merely proposed. The energy can be paid entirely by guard or
  unit absorbers, is not THM-2266's signed shallow response, and gives no
  profile exclusion. The 15 repeated-first profiles remain below one-comb
  capacity and are outside this conclusion. The exact missing coordinate
  is positive ancestry-aware return to the blocker residual.
source: codex-2026-07-25-expiration-owner-absorber-cut
depends_on:
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
related:
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2268-two-shell-private-owner-trident-and-raw-carry-cocycle-no-go
script: 04-computation/lrc14_expiration_owner_absorber_cut_thm2271.py
output: 05-knowledge/results/lrc14_expiration_owner_absorber_cut_thm2271.out
script_sha256: d238f51a61ac81f9786730f0d6440deafe2a714b9c5badb3eeb206139240db6e
output_sha256: d7510603b3d9cbbc397403fdbf90fc77a051c13e564aa42659c6660a15af16aa
hash_basis: working-tree bytes (LF)
---

# THM-2271 -- a strict owner expiration pays a positive absorber cut

THM-2263 raises the strict-row labelled expiration-image floor above one
danger-comb capacity.  THM-2267 supplies the exact transition cut kernel.
The two results compose, but only after making every target state legal:
the guard complement must be retained as an absorber label, and the
post-expiration image must be joined back to its marked source by an exact
root-branch section.

This produces a substantial positive transition energy on every strict row.
It does not yet produce a blocker-to-blocker handoff or a contradiction.

## 1. Scalar cover and the marked expiration image

Use the scalar notation

```text
D_s={x in R/Z:||sx||<1/14},
C_H={x in R/Z:||Hx||>1/7},

A_0=C_H minus union_(i=1)^5 D_(q_i).                 (1)
```

Assume the scalar cover

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(k=1)^3 D_(c_k)             (2)
```

almost everywhere, and a strict first-depth-one profile

```text
c_k=13^(lambda_k)u_k,          13 does not divide u_k,

(lambda_1,lambda_2,lambda_3)=(1,b,c),
2<=b<c,                        5<=c<=19.             (3)
```

For each blocker label, let

```text
E_j
 =A_0 intersection D_(c_j)
       minus union_(k!=j)D_(c_k).                   (4)
```

Thus, modulo null endpoints, every point of `E_j` is serviced by exactly
one of the nine labels used below, namely `c_j`.

Put

```text
T(x)=13x mod 1,

f_j=T^(lambda_j+1),
B_j=f_j(E_j).                                       (5)
```

THM-2263 proves that some labelled `j` satisfies

```text
mu(E_j)>=15041431/593783190                         (6)
```

and, using the exact expiration expansion,

```text
mu(B_j)>=15041431/70270200.                         (7)
```

Every danger comb has Haar mass `1/7`.  Therefore

```text
Y_j=B_j minus D_(c_j)

satisfies

mu(Y_j)
 >=15041431/70270200-1/7
 =5002831/70270200
 >0.                                                (8)
```

Equation (8) is stronger than saying that some owner switch occurs at a
torsion witness.  It is a positive Haar floor on the target support where
the marked source owner is unavailable.

## 2. Total target eligibility requires an absorber

Let the service-label set be

```text
L={g,q_1,...,q_5,c_1,c_2,c_3},                      (9)
```

where `g` is the guard-complement absorber.  At a circle point `z`, define
the complete eligibility set

```text
S(z)
 ={g:z notin C_H}
  union {q_i:z in D_(q_i)}
  union {c_k:z in D_(c_k)}.                         (10)
```

The global cover (2) makes `S(z)` nonempty almost everywhere:

- outside `C_H`, the label `g` is available;
- inside `C_H`, at least one of the eight danger labels is available.

At a source point `x in E_j`, equations (1) and (4) give

```text
S(x)={c_j}.                                         (11)
```

At a target point `y in Y_j`, equation (8) gives

```text
c_j notin S(y),                 S(y) nonempty.       (12)
```

Without `g`, a target outside the guard could have empty eligibility and
would not define a legal transition state.  Adding the absorber is therefore
not cosmetic: it is the sidecar that makes the owner-transition problem
total while recording exactly why a point left the blocker residual.

## 3. A finite measurable branch section

Write

```text
N=13^(lambda_j+1).
```

The map `f_j(x)=Nx mod 1` has the labelled inverse branches

```text
r_y=(y+r)/N,                         0<=r<N.          (13)
```

For every `y in B_j`, at least one branch lies in `E_j`.  Choose the least
such branch index:

```text
r(y)=min{r:r_y in E_j},
s(y)=r(y)_y.                                         (14)
```

This is a measurable section of `f_j` over `B_j`.  In the present setting it
has a more elementary finite description.  Every set in (1), (2), and (4)
is a finite Boolean combination of rational intervals.  For each fixed
branch `r`, the condition `r_y in E_j` is therefore a finite Boolean
combination of intervals in `y`.  Taking the least legal branch and then
refining by all memberships in (10) partitions `Y_j`, modulo finitely many
endpoints, into finitely many positive-measure cells

```text
Y_j=disjoint union_alpha C_alpha.                   (15)
```

On each cell:

- the selected root branch `r(y)` is constant;
- the source eligibility is the singleton `{c_j}`;
- the complete target eligibility is one fixed nonempty subset of
  `L minus {c_j}`.

For every `C_alpha`, make two transition vertices, one for its selected
source cell and one for its target cell, and join them by an edge of weight

```text
w_alpha=mu(C_alpha).                                 (16)
```

This is a finite weighted eligibility graph of the exact type in THM-2267.
The weights measure target **support**, not preimage multiplicity.

## 4. The singleton owner cut

Apply THM-2267 with the owner cut

```text
Q={c_j}.                                             (17)
```

Equation (11) forces every source vertex to the `Q` side, while equation
(12) forces every target vertex to the complementary side.  Each edge in
(16) must cross the binary cut.  Consequently its binary cut energy is

```text
kappa_Q
 =sum_alpha w_alpha
 =mu(Y_j)
 >=5002831/70270200.                                (18)
```

THM-2267's cut inequality gives the same lower bound for the full
owner-switch energy:

```text
omega_G>=kappa_Q>=5002831/70270200.                 (19)
```

Thus all `150` strict profiles force positive weighted transition energy.
The conclusion is invariant under which of the three blocker labels is
selected by THM-2263.

The terminology matters.  The energy in (19) is an eligibility-switch
energy on an exact source-to-target section.  It is not the raw relative
integer-carry increment of THM-2268, whose locally nonzero jumps form an
exact cocycle and may have zero loop holonomy.

## 5. Why positive energy is not a contradiction

The cut in (17) has only two sides:

```text
the marked blocker c_j
versus
every other blocker, every unit mask, and the guard absorber. (20)
```

All of the mass in (18) could, in principle, be paid by `g` or by the five
unit labels.  Neither LRC nor the existing Bellman inequalities assign a
forbidden cost to that handoff.  In particular:

1. THM-2261 proves that raw guard-plus-owner expiration has the whole circle
   as image.  No smaller named target follows without the global service
   sidecars.
2. THM-2266's signed shallow response

   ```text
   s=integral R 1_(D_(u_1))
   ```

   is not the cut energy (19).  Its bounded-relation alternative does not
   say where the marked expiration support lands.
3. THM-2268 forces a global cyclic owner word, but does not glue its private
   torsion witnesses into the branch section (14) with positive Haar weight.

Accordingly (19) excludes no profile.

## 6. The exact missing return coordinate

The blocker residual at a target point is `A_0`.  Define the
ancestry-aware return set

```text
H_j
 =E_j intersection
   f_j^(-1)(A_0 minus D_(c_j)).                      (21)
```

Equivalently,

```text
f_j(H_j)=B_j intersection A_0 minus D_(c_j).         (22)
```

At a target in (22), neither the guard absorber nor a unit label is
available.  The global cover therefore forces one of the other two blockers.
Hence a positive lower bound on (21) would upgrade (19) from a generic
owner-to-absorber cut to a genuine blocker-to-blocker handoff.

No existing result proves even

```text
mu(H_j)>0
```

for the THM-2263-selected label.  This is the precise gluing coordinate
left after composing THM-2261, THM-2263, THM-2266, and THM-2267.

### Cheapest decisive computation

The two sharp incoming spectra point to the same first control:

```text
profile (1,3,5),
all normalized blocker cores equal,
reduced shallow relation K_0=1.                     (23)
```

The profile `(1,3,5)` is the unique worst strict overlap ledger in THM-2263
and the weakest interior Bellman row in THM-2266.  The relation `K_0=1`
also contains THM-2261's local hostile geometry, so it attacks the proposed
bridge at its hardest branch.

The decisive finite problem is:

```text
retain an exact marked E_j branch through f_j;
enforce the pointwise global cover at every descendant;
retain the exact related-pair root law for K_0=1;
minimize mu(H_j), separately for j=1,2,3,
subject to the THM-2263 labelled source/image floors.       (24)
```

An exact zero optimum should be frozen as an ancestry/gluing stopping
witness.  A positive optimum would be the first lawful input to a
blocker-only cut and should then be tested over THM-2266's `1,822` small
primitive relation shapes.  Starting with all shapes would be premature:
the local membership and relation data are known to admit zero-future
hostile behavior before the global ancestry constraints are installed.

## 7. Repeated-first boundary

For the `15` profiles `(1,1,c)`, THM-2263 gives only

```text
mu(B_j)>=5229541/70270200
          <1/7.                                     (25)
```

Subtracting one-comb capacity yields

```text
5229541/70270200-1/7
 =-4809059/70270200<0.                              (26)
```

Thus the capacity argument in (8) forces no positive singleton-owner cut
in the repeated-first branch.  Its crossing of the smaller `1/14` threshold
needs a different named half-comb or two-owner sidecar.

## 8. Exact arithmetic audit

The companion verifies the `150/15` profile split, both labelled owner
floors, the `169/20` expiration factor, the strict cut floor in (18), and
the repeated-row failure in (26).  Reproduce with

```bash
python3 04-computation/lrc14_expiration_owner_absorber_cut_thm2271.py
python3 -O 04-computation/lrc14_expiration_owner_absorber_cut_thm2271.py
```

Normal and optimized transcripts are byte-identical to the stored output.
QED.
