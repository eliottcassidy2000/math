# The missing two-index carrier is spectral, not a packet label

Status: **FINITE-EXACT typing obstruction plus positive physical-macro
control.**  This note does not construct an endpoint current, decrement a
ledger, exclude a scalar row, or prove LRC(14).

## Inheritance pass

The closest proved mechanism is the exact relation-residue pushforward of
`THM-2334-relation-residue-current-and-character-twist-pushforward.md`, with
the factorwise target covariance of
`THM-2350-owner-pivot-dual-dipole-normal-form.md`.  The canonical hostile is
the full-`X` annihilation in
`THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary.md`:
refining an endpoint pair is useless if its entire pair-to-target line sum is
still zero.  `THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary.md`
is the corrected near miss: its moving-endpoint support is one-sided and its
status explicitly leaves the left residue missing.  The least-used positive
sidecar is `THM-2625-canonical-endpoint-current-full-transvection-sector-survival.md`,
which really does split a marked coefficient into separately descended left
and right spectral factors before forming a joint endpoint current.

The recent physical-carrier theorems
`THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary.md`,
`THM-2701-literal-singleton-word-one-step-dilation-nilpotence.md`,
`THM-2707-full-physical-lift-fibre-common-simplex-and-packet-scc.md`, and
`THM-2710-central-half-phase-literal-word-nilpotence-and-prescribed-clock-invisibility.md`
make a different kind of object: packets, support overlaps, and physical
translations.  The central question was whether that object could supply the
two indices missing from THM-2563.  It cannot do so merely by relabelling its
vertices.

## The type dictionary

Six nearby alphabets must remain separate.

```text
x                         circle/event variable
w_i and named roles       runner speeds and semantic factor roles
n in Z/(13^6)             physical deck address
n mod 13                  residue class of that physical address
carry, root, owner        event-dependent packet labels
s,t in the target dual    independently chosen character-twist parameters
u,beta,v in Z^d           global Fourier multiindices in THM-2334
r=u+13^6 beta+m e_c-v     exact relation address
eta.(u-v) mod 13          target residue after the chosen allocation
```

In particular, `u` and `v` are summation indices of Fourier coefficients.
They are not pointwise observables of `x`; they are not runner, packet, carry,
root, or physical-lift labels.  Thus an expression such as

```text
1_(left relation index at x = u) 1_(right relation index at x = v)
```

is ill-typed.  An exact relation address also does not recover its two
endpoints: THM-2334 groups a whole free transitive gauge orbit of
decompositions over each address.

The faithful two-index object is instead a **spectral cospan**.  At a fixed
marked frequency triangle `(X,m,Y)`, independently twist the left and right
endpoint factors, retain the two arrays `P_s` and `Q_t`, and only then take
their two-dimensional inverse finite Fourier transform.  In the notation of
THM-2625 this produces

```text
Lstar(l) = sum_s P_s chi_s(-l),
Rstar(r) = sum_t Q_t chi_t(+r),
J(l,r)   = C Lstar(l) Rstar(r).
```

The target current is the line sum over `l-r=q`.  This order matters:
THM-2350 covariance must hold factorwise, and a chosen pivot/allocation is
part of the data.  A packet orbit can be the physical ancestry on which such
factors are built, but it cannot replace them.

## An all-depth transversality obstruction

Let `P=13`.  At depth `r`, a physical translation `x in C_(P^r)` moves a
blocker/graft factor pair in the direction

```text
x |-> (b x, g x),
```

whereas the required target dipole moves it in the opposite direction

```text
y |-> (-y,y).
```

Equality is the kernel equation for

```text
Phi_(b,g)(x,y)=(b x+y,g x-y).
```

Its determinant is `-(b+g)`.  In every canonical pair checked here,

```text
b in {2197,742586} = 0 mod 13,
g in {27,40,53,66,14} = 1 mod 13,
```

so `b+g` is a 13-adic unit.  Therefore `Phi_(b,g)` is an automorphism of
`C_(13^r)^2` for every `r>=1`; the two cyclic directions intersect only at
the identity at every depth.  Restricting `y` to the order-thirteen target
subgroup cannot create a kernel.

At the physical depth `R=13^6`, the same proof is the elementary pair of
congruences

```text
b k/R =  s/13  (mod 1),
g k/R = -s/13  (mod 1).
```

Adding gives `R | (b+g)k`.  Since `gcd(b+g,R)=1`, this forces `k=0`, and then
`s=0`.  The neutral coordinates independently give the same obstruction:
the guard `H=1` already forces `R|k`; even if the guard is erased, neutral
source `c1=13` forces `13^5|k`, after which the blocker phase is integral and
again `s=0`.

The load-bearing hostile control is positive after deleting one factor.  If
the oppositely moving graft coordinate is forgotten, the blocker equation
alone has `gcd(b,R)` physical solutions for every target colour: `2197` for
the first dipole and `371293` for the second.  It is precisely the paired
opposite direction—not a lack of physical lifts—that kills the target
action.

This is the useful finite-torus meaning of the “Kakeya needles” picture:
the physical diagonal-translation needle and the target anti-diagonal needle
are transverse in the two-factor phase torus, and the toothpick refinement to
deeper powers of thirteen remains transverse.  No Kakeya measure theorem is
being invoked.

## What the full packet SCC does give

An independent exact rescan of all `13^6=4,826,809` addresses recovers the
THM-2707 bank of `3346` packets, split among eleven residue classes as

```text
0:304, 1:305, 2:304, 3:305, 4:304, 5:305,
6:304, 9:301, 10:304, 11:305, 12:305.
```

Exactly the `304` residue-zero packets retain the following atom.  The check
is on the whole inherited open interval, not only its centre.  It verifies
both components of that atom:

1. the base support at `q_n={13x}+7n/R`; and
2. the delayed prefix centred at `{Rz}T`, with the entire pulled-back radius
   `13 R T I_radius`.

Every one of those `304` endpoints has a physical two-step loop through each
of the `3042` packets in the other ten residue parts, hence

```text
304 * 3042 = 924768
```

rooted two-step loops.  More generally the complete directed multipartite
law gives a two-step path between every ordered pair of residue-zero
endpoints.  This is a strong positive physical macro: the frozen following
atom survives and the square of the lift relation is complete on that bank.

It is nevertheless not the missing current.  Each nontrivial macro edge is
a global deck translation, hence lies outside both nonzero target-dipole
subgroups by the all-depth lemma.  The two edge translations of a closed
macro sum to a target-neutral physical displacement; this does not define
the separate Fourier residues `eta.u` and `eta.v`.

## Two exact non-identifiability controls

The eight actual half-cycle support rows provide two cheap hostile tests.

First, let `S` be one such right-shift support.  Write
`ubar=eta.u`, `vbar=eta.v` for abstract projected residues.  The two abstract
joint tables

```text
J_diag(ubar,vbar)  = 1_S(vbar) 1_(ubar=vbar),
J_fixed(ubar,vbar) = 1_S(vbar) 1_(ubar=0)
```

have the same marginal after forgetting `ubar`.  The first is supported only
at `q=ubar-vbar=eta.(u-v)=0`; the second has nonzero-`q` support.  This is an
information-theoretic control only: it does not identify a physical shift
with a Fourier residue.

Second, use the dual data the calculation really knows: the fixed-left row
`K(0,t)`.  Put that observed profile in row zero of two `13 x 13` arrays and
set every unseen entry to zero in the first.  In the second, additionally
fill the unseen diagonal entries `(s,s)`, `s!=0`.  The observed row is
identical, yet their common-diagonal primitive target-character verdicts are
opposite on every one of the eight rows.  If `0` is present, the first
diagonal is a nonconstant singleton profile while the filled diagonal is
constant; if `0` is absent, the first diagonal vanishes while the filled
off-zero diagonal has every primitive character nonzero.

Consequently, no amount of sharpening the one-sided right-shift table alone
can determine the joint endpoint target current.

## Minimal sidecar and next decisive test

The smallest object that could close this gap has four parts:

1. a selected owner-aligned pivot/allocation;
2. a factorwise THM-2350-covariant pair of endpoint spectral orbits built on
   one common physical ancestry;
3. a fixed marked triangle `(X,m,Y)`, or an equally explicit physical
   normal/jet/orientation that supplies its weight; and
4. a proof that the resulting bi-Fourier coefficient survives the
   pair-to-target line sum.

The fourth item is essential.  Summing all `X` recreates the pointwise
danger-to-safe product annihilated by THM-2568.  Full packet support, more
recurrence, or a larger SCC does not repair that cancellation.

THM-2625 is the correct positive prototype: it separately descends the two
endpoint factors and exhibits a full joint spectral array at one marked
triangle.  Its limitation is equally exact: it is a canonical non-cover
control with no elementary chronology.  The high-leverage experiment is
therefore not another packet census, but a transplant test: on one of the
`304` frozen-following packet ancestries, can one define separately covariant
endpoint factors at a fixed marked triangle and prove one target line sum
nonzero?  Failure should be recorded at the first of covariance, common
ancestry, semantic endpoint typing, or line-sum survival.

## Reproduction

Run

```bash
python3 04-computation/lrc14_physical_lift_vs_target_dipole_intersection_probe_20260728.py
python3 -O 04-computation/lrc14_physical_lift_vs_target_dipole_intersection_probe_20260728.py
```

Both modes reproduce
`05-knowledge/results/lrc14_physical_lift_vs_target_dipole_intersection_probe_20260728.out`
byte for byte.  SHA-256:

```text
eb7bb7b1dd116bf276bfc774e001bd14505f079d8da52353dba369cbd50e1610  script
2ee95333f109550555342c57849d79f32ff51fe362942b97208d33fac2fe403c  output
```
