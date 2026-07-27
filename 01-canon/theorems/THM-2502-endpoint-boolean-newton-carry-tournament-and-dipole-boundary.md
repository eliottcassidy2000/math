---
id: THM-2502
title: "Endpoint Boolean Newton top cell, carry tournament, and dipole boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The 128 complete local
  endpoint masks of THM-2452 are the seven-cube. For each fixed
  blocker word, THM-2445's six first-failure labels have
  1,16,8,4,2,1 complete extensions; after union over the four blocker
  words the six layers contain 4,64,32,16,8,4 atoms. The all-safe
  ghost is the unique atom in the all-unit-safe, all-blocker-safe
  cell. Orienting two masks by their
  first differing role gives only the transitive lexicographic
  tournament; its role-h edge counts are 2^(13-h), while its successor
  chain is coloured by the 2-adic ruler h=7-nu_2(k+1). In Boolean
  Gregory--Newton coordinates the ghost indicator is the degree-seven
  top monomial: every lower mixed difference at zero and every
  coordinate-zero face misses it, while the full sevenfold difference
  is one. The lawful target quotient of THM-2350 has disjoint dipole
  representatives and no nonzero singleton-role character. Splitting
  the two roles of one active dipole therefore supplies only one
  universal four-cell zero-sum equation, with a three-dimensional
  ambiguity; moreover
  danger and safe complements have opposite nontrivial Fourier
  coefficients and identical jump nodes. Thus uncoloured tournament
  transitivity, lower Newton faces, singleton complement transport,
  and node-only Prony separation cannot recover semantic root
  orientation. THM-2452 bypasses all four losses by copying a complete
  atom and releasing the old exact frequency. No scalar row is removed.
source: mac-mini-2026-07-26-endpoint-newton-carry-boundary
depends_on:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2448-right-endpoint-cospan-transition-atlas
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
related:
  - THM-362-natural-operation-graph-shadows
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
  - THM-2438-poisson-newton-ternary-half-and-harmonic-divisor-incidence
  - THM-2456-two-root-replica-uniform-offset-boundary
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
  - THM-2462-mixed-radix-root-phase-orbit-universality
  - THM-2501-switching-fourth-moment-signed-c4-and-gram-energy
script: 04-computation/lrc14_endpoint_newton_carry_boundary_thm2502.py
output: 05-knowledge/results/lrc14_endpoint_newton_carry_boundary_thm2502.out
script_sha256: 7f6bea32a306237f942a8da70231bc5d0cfc64e7937546142e980dc14aa71bf1
output_sha256: e5901d57d1e81c83b4d1073447bedfc4a6e0e2fc9c28daaa049e8df960f243a2
hash_basis: working-tree bytes (LF)
---

# THM-2502 -- the endpoint cube is transitive only after losing its carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2452 restores the missing endpoint by refining to all `128`
complete Boolean masks, copying one mask to the other leg by
idempotence, and selecting a fresh exact frequency only afterwards.
Three tempting reformulations of this bank are exact but too coarse:

```text
first-difference tournament,
lower Gregory--Newton face data,
one-role complement derivative.                               (1)
```

The point of this theorem is to locate their common loss.  The
unlabelled relations are transitive; the discarded coordinates are
the carry depth, top Boolean atom, dipole partner, complex amplitude,
terminal word, and absolute root orientation.

## 1. The 24 cells refine six first-failure layers of the seven-cube

Use THM-2445/2452 notation.  The five guard/unit bits come first,
followed by the two nondeep blocker bits.  Write `1` for safe and `0`
for danger.  Thus a complete local mask is a word

```text
x=(x_1,...,x_7) in {0,1}^7.                           (2)
```

For `i=1,...,5`, the first-failure label `i` fixes

```text
x_1=...=x_(i-1)=1,       x_i=0,                      (3)
```

and leaves `5-i` guard/unit tail bits free.  For each fixed two-bit
blocker word `sigma`, it therefore has

```text
2^(5-i)
 =16,8,4,2,1                                             (4)
```

complete extensions.  The all-unit-safe label `i=0` has one complete
guard/unit extension for each fixed `sigma`.  Thus the exact
per-`sigma` extension counts are `N_i=1,16,8,4,2,1`.

The `24` THM-2445 cells are the pairs `(i,sigma)`.  Unioning the four
blocker words at each fixed `i` gives six coarser layers with sizes

```text
4,64,32,16,8,4,       with total 128.                  (5)
```

The word

```text
mathbf 1=(1,1,1,1,1,1,1)                              (6)
```

is the unique all-safe ghost, lying in the single `(i=0,
sigma=all-safe)` cell.  The other `127` atoms retain a literal first
danger, blocker danger, or both.  This is a refinement identity, not a
semantic THM-2305 word identification.

## 2. The first-difference tournament is exactly lexicographic order

For distinct words `x,y`, let `h(x,y)` be their first differing
coordinate.  Orient the pair from the word with `0` there to the word
with `1` there and colour the edge by `h`.

The orientation is exactly lexicographic order, hence a transitive
tournament.  For a fixed colour `h`, choose the common prefix, then
the two ordered differing bits, then the two independent suffixes.
The number of unordered edges is

```text
2^(h-1) * 2^(7-h) * 2^(7-h)
 =2^(13-h).                                            (7)
```

Thus the colour census is

```text
h:       1    2    3    4    5   6  7
edges: 4096 2048 1024  512  256 128 64,               (8)
```

whose sum is `8128=binom(128,2)`.

Identify a word with its seven-bit integer `k`.  The transitive
reduction of the total order is the successor chain

```text
0 -> 1 -> ... -> 127.                                  (9)
```

For `0<=k<127`, incrementing `k` resets exactly
`nu_2(k+1)` trailing ones.  Therefore the first-difference colour of
the edge `k -> k+1`, with coordinate one most significant, is

```text
h(k)=7-nu_2(k+1).                                     (10)
```

The simple shadow (9) looks like the uniform successor shadow of
addition by one.  Its edge colour is instead the binary ruler/carry
position.  THM-2448's tail recursion

```text
T(0)=1,       T(n)=2+2T(n-1)=3*2^n-2                 (11)
```

is the corresponding first-mismatch carry tree: two off-diagonal
branches stop, while two diagonal branches retain the next tail bit.
Forgetting edge colour and the remaining tail turns this entire
structure into an information-poor transitive tournament.

## 3. The ghost is the top Boolean Gregory--Newton coefficient

Let `X_i` denote the safe bit and define the ghost indicator

```text
G(X_1,...,X_7)=product_(i=1)^7 X_i.                  (12)
```

For `S subseteq {1,...,7}`, let `Delta_S` be the mixed Boolean forward
difference in the coordinates of `S`, evaluated at the origin.  The
multilinear Gregory--Newton expansion gives

```text
Delta_S G(0)=0              if |S|<7,
Delta_{1,...,7}G(0)=1.                                 (13)
```

Equivalently, the only nonzero multilinear coefficient of `G` is the
degree-seven monomial itself.  Also `G` restricts to zero on every
coordinate-zero face `X_i=0`.

Consequently any proposed endpoint detector made only from lower
mixed differences at zero, or only from coordinate-zero face
restrictions, has `G` in its kernel.  The all-safe ghost is not an
irregular exception to Newton calculus; it is precisely the top cell
that lower Newton data are designed to forget.  THM-2452 avoids this
kernel by retaining a complete atom and using `P_x^2=P_x`, not by
reconstructing the top coefficient from its faces.

There is an exact operation-graph formulation.  The seven-bit function
algebra is

```text
B=Q[X_1,...,X_7]/(X_i^2-X_i)
  isomorphic to Q^128.                                  (14)
```

Its complete-mask indicators `P_x` are the primitive orthogonal
idempotents:

```text
P_x P_y=0              for x!=y,
P_x^2=P_x.                                             (15)
```

Thus the fully labelled multiplication table has input vertices
`{P_x}`, output vertices `{0} union {P_x}`, and sends `(P_x,P_y)` to
`P_x` on the diagonal and to `0` off it.  It consists of idempotent
loops together with typed zero-products; multiplication on all of `B`
is not being called diagonal.  THM-2452 is precisely multiplication in this
semisimple algebra before Fourier extraction.  THM-2445's
first-failure cells are prefix-code sums of the primitive idempotents,
and Gregory--Newton is the Möbius change of basis between monomials and
these atoms.  More static polynomial identities in `B` cannot restore
clock, ancestry, root section, or terminal-word data discarded when
the physical packet was mapped into `B`; those coordinates belong to
the external action/embedding sidecar.

## 4. Lawful target characters are dipoles, never singleton roles

THM-2350 proves that every lawful target character has the unique
omitted-unit gauge representative

```text
ell(s,t)
 =s(e_a-e_(k_a))+t(e_b-e_(k_b)),       (s,t) in F_13^2,  (16)
```

where the two dipoles have disjoint support.  If exactly one of `s,t`
is nonzero, (16) moves two roles; if both are nonzero, it moves four.
Adding a nonzero gauge vector turns on every ungrafted unit
coordinate.  Hence no nonzero lawful target character moves exactly
one local role.

Choose one active dipole axis in (16): `(a,k_a)` when `s!=0`, or
`(b,k_b)` when `t!=0`.  Split the two roles of that dipole into the
four cells

```text
gg, gd, dg, dd.
```

At a nontrivial character on the chosen dipole axis, summing its four
statuses removes dependence on that axis.  The only universal linear
relation supplied by complementarity alone is

```text
J_gg+J_gd+J_dg+J_dd=0.                                (17)
```

The kernel of the row `(1,1,1,1)` is three-dimensional.  In
particular, fixing one transition such as `J_gd` does not determine a
matched diagonal such as `J_gg`.  For example,

```text
(-1,1,0,0),        (0,1,0,-1)                        (18)
```

both satisfy (17) and have `J_gd=1`, but have different `J_gg`.
These are coefficient-space hostile vectors proving non-identifiability;
they are not asserted to be physical LRC packets.
This statement is deliberately about the two roles of one active
dipole; an arbitrary cross-pair of moving roles need not have a
zero-sum parent.  Thus the one-role identity `J_g=-J_d` cannot be imported into the
lawful target quotient.  A coefficient-first repair would need a
dipole-polarized reference or two-copy square retaining the other
three cells, not a singleton Newton derivative.

## 5. Complements also defeat node-only separation

For a danger indicator `d` and its safe complement `g=1-d`, every
nontrivial cyclic Fourier coefficient satisfies

```text
g_hat(k)=-d_hat(k),             k!=0.                (19)
```

The two step functions have the same jump locations, with opposite
jump coefficients.  Therefore a Prony or resultant test that retains
only the endpoint-node polynomial sees identical node support on the
two sides.  A mismatch label alone supplies no node transversality.
This is the complement boundary already isolated in THM-2408, now
seen inside each Boolean endpoint role.

## 6. Exact consequence and stopping boundary

The four views fit into one loss diagram:

```text
complete endpoint atom
  --forget edge colour/tail--> transitive tournament
  --discard top difference--> lower Newton faces
  --project target action--> dipole four-cell equation
  --discard jump signs--> common Prony node set.                (20)
```

None of these quotients can recover a canonical root orientation or a
semantic owner/repair word.  The intrinsic atom relation is
THM-2457's directed co-support graph, not a tournament.  These
quotients also do not obstruct THM-2452:
aggregate idempotence keeps the complete atom, sums over the old exact
frequency, and then reselects a fresh exact address.  That operation
deliberately bypasses every quotient in (20).

Subsequent THM-2457 proves that literal truth bits do not determine the
semantic word and supplies common-root service from a directed atom
edge.  THM-2459 retains aggregate drift and such service together in a
Boolean union of at most four atoms.  THM-2460 proves that an already
retained word copies with no further loss and identifies maximizing-entry
root-image inclusion as the exact sufficient sidecar; it does not make
that inclusion automatic.  THM-2462 further shows that fixed-speed
phase compatibility alone cannot exclude the finite uniform-offset
atlas.  The remaining LRC(14) obligation is to use the full scalar-cover
coupling to exclude the physical THM-2456 mixtures or force the retained
drift/service packet into a canonical same-word root image.  No scalar
row is removed, and LRC(14) is not proved.

## 7. Exact companion

Run

```text
python3 04-computation/lrc14_endpoint_newton_carry_boundary_thm2502.py
python3 -O 04-computation/lrc14_endpoint_newton_carry_boundary_thm2502.py
```

The dependency-free exact companion:

- enumerates all `128` masks, all `8128` unordered pairs, and the
  seven first-difference colour counts;
- verifies transitivity and the `127` successor ruler colours;
- reconstructs all `24` `(i,sigma)` cell sizes, their six aggregate
  layers, and the unique ghost;
- evaluates all Boolean mixed differences of the top monomial;
- verifies the dipole support dichotomy over all `168` nonzero
  characters of `F_13^2`;
- checks the four-cell rank-three ambiguity and the exact cospan tail
  recurrence through all five omitted-tail depths; and
- verifies the complement polynomial identity behind (19).

Normal and optimized runs must byte-match

```text
05-knowledge/results/lrc14_endpoint_newton_carry_boundary_thm2502.out
```

QED.
