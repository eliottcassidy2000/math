---
id: THM-2894
title: "Unmarked residual semilattice order and group-clutch no-go"
status: >
  PROVED + LEAN-VERIFIED.  Literal deletion by danger labels factors through
  the free commutative idempotent semilattice of finite label sets.  It is
  therefore blind to order and repetition, and remains order-blind after
  every later unmarked continuation.  Any union-multiplicative invariant
  from this carrier to a group is trivial.  Consequently THM-2888/2893's
  unmarked union/residual algebra cannot itself retain THM-2889's ordered
  quaternionic central sign.  This is a bridge guardrail, not a row
  exclusion or an LRC(14) proof.
source: root/lrc-residual-semilattice-no-go-2026-07-29
depends_on:
  - THM-2887-quaternionic-arf-lift-of-the-semantic-v4-and-global-carry-no-go
  - THM-2888-eight-body-first-apex-global-pair-cap-atlas
  - THM-2889-dicyclic-reverse-action-joint-carrier-and-skew-lift-separation
  - THM-2893-complement-cap-finite-core-flag-lemma
lean_module: 04-computation/lean/TournamentH7/TournamentH7/LRCResidualSemilatticeNoGo.lean
lean_sha256: 4bf7150a11c41b90e6fbd249e76acb12c0593df0b249c6494e06d41b19a6cac9
hash_basis: LF-normalized bytes
---

# THM-2894 -- unmarked residual semilattice order and group-clutch no-go

**PROVED + LEAN-VERIFIED.**

The direct THM-2888/2893 route and the THM-2889 horn route use genuinely
different carriers.  This theorem identifies the exact information loss
between them.

## 1. The universal unmarked residual carrier

Let `V` be a set of allowed labels, let `C` be a carrier, and attach a danger
set `D_v` to each `v in V`.  For a finite word `w` in the labels, write
`supp(w)` for its finite support and define the literal residual

```text
R_C(w)=C minus union_(v in supp(w)) D_v.                 (1)
```

Thus `(1)` factors as

```text
free words on V
  --support--> P_fin(V)
  --literal subtraction--> residual sets.               (2)
```

The middle object is the free commutative idempotent semilattice:

```text
S union T = T union S,        S union S = S.             (3)
```

Consequently, for all words `u,v,z`,

```text
R_C(uv)=R_C(vu),
R_C(uu)=R_C(u),
R_C(zuv)=R_C(zvu).                                      (4)
```

The last equality is the persistence statement: no later unmarked deletion
can recover order once support has forgotten it.

More generally, if `F` is any function of the residual set alone, then

```text
F(R_C(uv))=F(R_C(vu)).                                  (5)
```

This does not require `F` to be additive, multiplicative, scalar-valued, or
continuous.

## 2. Every multiplicative group image is trivial

Let `G` be a group and suppose

```text
phi : P_fin(V) -> G
```

satisfies

```text
phi(S union T)=phi(S) phi(T).                            (6)
```

Idempotence gives

```text
phi(S)=phi(S union S)=phi(S)^2.
```

Left cancellation in `G` yields

```text
phi(S)=1                                                (7)
```

for every finite `S`.  No normalization assumption at the empty set is
needed.  Hence a nontrivial group clutch cannot be a
union-multiplicative invariant of unmarked residual support.

The root-imported Lean module
`LRCResidualSemilatticeNoGo` proves:

```text
insert_order_blind
insert_idempotent
continuation_order_blind
no_order_sensitive_lift
union_multiplicative_group_invariant_trivial.
```

All five theorems are sorry-free; their `#print axioms` reports contain only
standard Mathlib foundational axioms and no `sorryAx`.

## 3. Application to the two LRC carriers

THM-2888 and THM-2893 use exactly the carrier in `(2)`.  Their observables

```text
U_C(S)=measure(C intersect union_(v in S) D_v)
```

and their literal residuals depend only on the finite set `S`.  In
particular, the heavy-pair relation is intrinsically undirected.  The exact
verifiers sort label tuples, check deletion against simultaneous union, and
retain no ancestry, section, or order field.

That loss is harmless for their target predicate: whether a finite family of
danger sets covers `C` is itself support-based.  It is nevertheless fatal to
an attempted descent of the THM-2889 central clutch.  In the canonical
quaternionic section,

```text
QA*QB=-QAB_canonical,        QB*QA=+QAB_canonical.       (8)
```

The two values in `(8)` are distinct, whereas the two corresponding unmarked
residuals are equal by `(4)`.  Therefore no map determined by the residual
set can assign the two ordered lifts in `(8)`.  The stronger statement `(7)`
rules out trying to repair this with a multiplicative group-valued statistic
on finite support.

This does **not** say that an unordered pair of externally labelled `V4`
directions cannot be assigned the symmetric commutator-curvature bit
`det(QA,QB)=1`.  In characteristic two that bit is unchanged by swapping the
directions.  What the residual carrier cannot recover is which speed
occurrence is the `QA` or `QB` axis, which axis carries the selector versus
carry character, or the absolute ordered lift relative to the canonical
`QAB` section.  More precisely, finite support retains the speed identities,
but the residual algebra itself supplies no canonical or physical assignment
of those occurrences to the semantic axes.  The clutch in `(8)` is of this
last, order-sensitive type.

This also explains why orienting the heavy graph cosmetically is invalid.
There is no asymmetric pairwise observable in the direct carrier from which
such an orientation could be recovered.

## 4. The minimum lawful bridge

A sidecar capable of retaining the full physical THM-2889 clutch must contain
at least:

1. an ordered word or path field which does not factor through finite
   support;
2. the marked ancestry and quaternionic section gauge;
3. a coefficient action on one common physical packet; and
4. a proof that forgetting the sidecar still preserves the literal coverage
   predicate used by the direct route.

Once physical semantic edge labels are supplied, THM-2887 gives the exact
minimal algebraic lift.  Let `H=V4=F2^2` and use the cocycle displayed there
in the chosen canonical-`QAB` section:

```text
alpha(x,y)=x1*y1+x2*y2+x1*y2.                          (9)
```

Replace an unordered support by a strict oriented Hasse flag

```text
empty=S_0 subset S_1 subset ... subset S_m
```

whose edge `S_(i-1)->S_i` has a supplied semantic label `a_i in H`.  Starting
from `(h_0,e_0)=(0,0)`, accumulate

```text
h_i=h_(i-1)+a_i,
e_i=e_(i-1)+alpha(h_(i-1),a_i).                        (10)
```

Then

```text
(R_C(S_m),h_m,e_m)                                    (11)
```

is a `Q8`-lifted residual state, and forgetting the last two coordinates
preserves the literal residual exactly.  Swapping adjacent labels `a,b`
changes the central coordinate by

```text
alpha(a,b)+alpha(b,a)=det(a,b).                        (12)
```

For two paths carrying the same fixed semantic directions `a,b` in opposite
orders, `(12)` measures the change of lift.  In particular, for
`a=QA,b=QB`, the two histories have the same residual and `V4` endpoint but
opposite central lifts.  After the `V4` endpoint, fixed semantic edge labels,
ordered path, section/cocycle, and initial lift are supplied, the minimal
endpoint fibre retaining which central lift occurred is `C2`.  Formula `(10)`
computes the bit from the retained path; it becomes separately necessary only
when path history is compressed to its endpoint.  Once the labels are
supplied, `Q8` is the complete semantic extension state space, but it does not
construct the missing speed-to-edge assignment.  Retaining the full THM-2889
address clutch requires its `(C169^2) semidirect Q8` state.

THM-2887 rules out the specific normalized base-independent map
`ell:C13 -> V4` which would satisfy both the carry coboundary law and the
required horn values.  Hence the labels in `(10)` cannot come from that
residue-only ansatz; a lawful assignment must retain some additional
contextual coordinate, such as packet, ancestry, origin, or event data, or
abandon one of those carry/horn laws.  No claim is made here against every
residue-to-`V4` function outside THM-2887's hypotheses.  An order field by
itself is only decoration.

The cheapest decisive physical test remains THM-2889's typed edge

```text
e9=(-9,+9,QB)
```

on the same marked packet as `e8`, with via/direct clutch

```text
(13,13,-1).                                           (13)
```

If `(13)` is not realized physically, the added orientation carries no
theorem-bearing information.

## 5. Scope

This theorem proves a quotient-loss boundary.  It does not weaken the
THM-2888/2893 heavy-flag reduction, whose target is already unmarked.  It
does not construct the ordered sidecar listed above, close a heavy-triangle
residual, exclude an LRC row, or prove LRC(14).
