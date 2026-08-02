---
id: THM-3141
title: "Quartic V4 modular congruence shadow and Gamma(3) sidecar boundary"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/codex-thm3088-push-2026-08-02
audit: >
  An independent hostile audit reconstructed the pointed A4 action, the
  PSL2(F3)/PGL2(F3) identifications, the first length-six congruence relation,
  all Farey and THM-2056 gate controls, and the rooted-C3/equal-factor scope.
  It caught and repaired noncanonical-pointing, physical-safety,
  cyclic-readout typing, and transcript-wording overclaims.  Fresh immutable
  normal and optimized runs agree after LF normalization with the stored
  eleven-line transcript and both declared hashes.
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
  - THM-3121-path-cover-walk-content-substitution-kernel
  - THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform
script: 04-computation/modular_v4_psl2f3_congruence_shadow_thm3141.py
output: 05-knowledge/results/modular_v4_psl2f3_congruence_shadow_thm3141.out
hash_basis: LF-normalized bytes
script_sha256: 1730c911852c4680bc85ec7eeeee17861dfd4a58265e221e46c194b5e925b1ea
output_sha256: d87ef15554abb466da6cc63988f3bd1e4f044dc08b8ed9b4aed3ce5a9d935a48
---

# THM-3141 -- quartic `V4` modular congruence shadow and `Gamma(3)` sidecar boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The exact structure behind the binary/ternary analogy is presently a modular
**congruence shadow**, not a faithful modular-group action.  THM-2950 supplies
an unbased `V4` torsor, with no canonical order-three automorphism.  After a
base half-system, a nonzero translation, and a cyclic orientation are chosen,
the four half-systems admit the displayed `PSL2(F3)=A4` action.  Its missing
coordinate is precisely the level-three kernel detected by `T^3`.

## 1. The pointed `V4` torsor

Write the half-system torsor as

```text
Omega = F2^3 / <(1,1,1)> = {o,x,y,z}.                    (1)
```

Choose the origin `o`, the nonzero translation `x`, and an oriented labeling
of the three conjugate-pair directions.  Let `S` translate by `x`, and let
`R` cyclically permute the three nonzero pair-flip directions:

```text
S=(o x)(y z),             R=(x y z).                     (2)
```

Then

```text
S^2=R^3=(SR)^3=1,
<S,R> = V4 semidirect C3 = A4.                           (3)
```

Indeed the conjugates of the nonzero translation `S` under `R` generate all
of `V4`, and `R` supplies its order-three automorphism.  This is the
orientation-preserving half of THM-2950's full affine symmetry
`V4 semidirect S3=S4`.

An explicit projective identification is

```text
o -> [1:1],  x -> [1:2],  y -> [1:0],  z -> [0:1]
                                                        in P1(F3),          (4)
```

under which `(2)` is induced by

```text
S = [ 0 -1 ]             R = [ 0 -1 ].                  (5)
    [ 1  0 ]                 [ 1  1 ]
```

Consequently, as abstract groups,

```text
A4 = PSL2(F3),             S4 = PGL2(F3),
S4/V4 = S3 ~= PSL2(F2).                                   (6)
```

The `3`-statement is generator-compatible with the displayed modular
matrices: the pointed four-sheet frame is the mod-`3` projective line.  The
`2`-statement is only an exceptional **abstract** isomorphism here.  The
selected `S` in `(2)` lies in `V4` and therefore becomes trivial in
`S4/V4`, whereas the standard modular `S` modulo `2` is nontrivial.  Thus the
cubic resolvent is the `S4/V4` quotient, abstractly isomorphic to
`PSL2(F2)`, but it is not a compatible mod-`2` reduction of the displayed
`S,R` action without a different lift.

## 2. The first missing modular coordinate

Over the integers, the matrix product in `(5)` is projectively

```text
SR = [1 1] = T,
     [0 1]
```

up to the central sign.  Therefore `(SR)^3=T^3` is nontrivial in
`PSL2(Z)`, but becomes the identity modulo `3`.  Exact reduced-normal-form
enumeration in `C2*C3` finds no nonempty identity of syllable length below
six; at length six it finds `(SR)^3`, `(SR^2)^3`, and their cyclic starts.
Thus `(SR)^3=1` is the first extra relation on the torsor, not merely one
relation among many shorter accidents.

Equivalently, the intrinsic action factors through

```text
PSL2(Z) -> PSL2(F3)=A4,                                  (7)
```

and forgets the `Gamma(3)` coordinate.  The smallest concrete witness is

```text
T^3(0,1)=(3,1)!=(0,1),                                  (8)
```

although the two slopes define the same point of `P1(F3)`.

## 3. Farey flanks survive only one way

For primitive integer vectors `u,v`, a Farey edge has
`|det(u,v)|=1`.  Reduction modulo `3` leaves this determinant nonzero, so

```text
Farey neighbors -> distinct points of P1(F3).             (9)
```

The converse is false: `(1,0)` and `(1,4)` have determinant `4` but distinct
projective residues.  Moreover the fixed THM-2056 sufficient polar
certificate is not preserved
by the modular generators.  In its one-tail chart put

```text
D(a,b)=max(13|b|,|a-12b|),       d=(-100,-7).             (10)
```

Then

```text
91D(d)=8281 <= ||d||^2=10049,
91D(Sd)=118300 > ||Sd||^2=10049,
91D(Rd)=126581 > ||Rd||^2=11498.                         (11)
```

Thus exact unimodular adjacency survives, but the fixed THM-2056 determinant
certificate does not.  This witness makes no claim that `Sd` or `Rd` is
physically unsafe: outside the polar polygon means only that this sufficient
certificate is silent.  Transporting the polar polygon and quadratic metric
along with the slope would make the *certificate* equivariant, not prove a
currently available physical LRC move.

## 4. The ternary root is separately forgotten

THM-3121/3134 supply a cyclic-substitution transform whose **unrooted
numerical profile/endpoint-jet output** does not retain a chosen block root.
Let `S1,S2,S3` be the labelled input tournaments, let `Psi` denote the actual
substitution, and let `J` denote its complete numerical endpoint-jet readout.
For unequal labelled factors, cyclic rotation of the input triple gives an
isomorphic substituted tournament; for equal or canonically identified
factors, THM-3134 does give the substituted tournament a cyclic-wreath `C3`
automorphism.  In either case put

```text
rho(S1,S2,S3)=(S2,S3,S1).                                (12)
```

in the precise sense

```text
J(Psi(rho(S1,S2,S3))) = J(Psi(S1,S2,S3)).                (13)
```

For unequal labelled factors the two substituted objects underlying `(13)` are
identified only after the cyclic relabeling/isomorphism is forgotten; `(13)`
does not assert equality of rooted labelled tournaments.  The complete
numerical jet retains every path-cover count but still forgets
which input block was designated as root.  A rooted quotient block is
therefore independent of the endpoint-jet sidecar, even when an equal-factor
output tournament has its genuine cyclic automorphism.  On the naive product
of the four-point torsor with a labelled jet triple, the diagonal action still
obeys `(SR)^3=1`; if `R` is put only on the jet factor, it commutes with `S`
and degenerates even further to `C2 x C3`.

## 5. Exact evidence and scope

The companion performs only finite integer arithmetic.  It checks the two
generator orders, all `12` generated permutations, all reduced normal forms
through the first length-six collision, the projective matrix realization,
`744` bounded Farey-neighbor controls, the determinant-`4` converse hostile,
and all three exact gate values in `(11)`.  Normal and optimized raw
transcripts agree, and their LF-normalized bytes match the stored output.

The theorem target is deliberately structural.  It proves neither a faithful
`PSL2(Z)` action nor a degree-four Keller exclusion.  It supplies no common
physical atom and no LRC gate.  A genuine combined carrier still needs:

1. a `Gamma(3)`/Farey lift on which `T^3` remains visible;
2. a rooted `C3` block in addition to the endpoint jet; and
3. a common physical coupling to the pointed `V4` lift.

Without all three, the pointed four-sheet action remains a mod-`3`
congruence shadow, while the binary/resolvent and ternary/substitution towers
remain sidecar grammars rather than two generator-compatible faithful free
factors acting on one object.

THM-2951 supplies the sharp reconstruction boundary behind item 3: a full
fifth compound can reconstruct multiplication, but no natural signed-pair
linear contraction canonically selects the pointed balanced-third-compound
phase needed for the common physical coupling.
