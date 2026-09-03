# September 3, 2026 intake: rank-two Poisson gauge and percolation gluing

## Status at a glance

The supplied arXiv identifiers concern different subjects and do not combine
into a theorem about planar Jacobian maps.

| source | status used here | imported content | not imported |
|---|---|---|---|
| Christopher D. Long, [arXiv:2608.23777v1](https://arxiv.org/abs/2608.23777) | **PREPRINT v1; formulas independently exact-checked** | explicit nonautomorphic rank-two Poisson map; four-variable symplectic Keller map; appendix construction in the fourth Weyl algebra | planar `JC(2)`, `DC(2)`, or a second independent counterexample mechanism |
| Kozma--Nitzan, [arXiv:2401.12397v1](https://arxiv.org/abs/2401.12397) | **PREPRINT v1; conjectures remain conjectures in that paper** | finite-graph gluing conjectures, proved small-interface cases, and `Conjecture 3 => theta(p_c)=0` | a proof of Conjectures 1--4 or unconditional critical continuity |
| Leder, [fixed Lean artifact](https://github.com/anthropics/formal-math/tree/795efb86f191735c5481675763537cfb4ff37e55/percolation) | **EXTERNAL FORMAL CLAIM / SOURCE-REPORTED MECHANICAL CHECKS / NOT INDEPENDENTLY REBUILT OR REFEREED HERE** | claimed additive gluing theorem, Conjecture 3, and `theta(p_c)=0` for bond percolation on `Z^d`, `d>=2` | Kozma--Nitzan Conjectures 1, 2, or 4; site percolation; other lattices; quantitative critical behavior |

## 1. Long's rank-two Poisson counterexample

Long works in `P_2=Q[x,q,p,z]` with `{p,x}={z,q}=1`.  The
preprint gives explicit `(R,T,D,S)` satisfying the six canonical Poisson
relations, with

```text
R=x(2-3xq),
```

and an exact reduced three-point fibre over `(0,1/8,0,0)`.  Thus its point map
is a nonautomorphic symplectic Keller map in four variables and the standard
two-pair Poisson conjecture is false.  Appendix B constructs a
nonautomorphism of `A_4`; its rank is four Weyl pairs, not `A_2`.

### Exact relation to THM-2044

This is not a second noninvertibility mechanism.  Let

```text
s=xq,
F=q^3(63s^2+318s+601)/840,
U_F:(x,q,p,z) |-> (x,q,p+F_x,z+F_q),
tau:(r,t,d,s) |-> (r,s,d,-t).
```

The exact certificate
[THM-4397](../../01-canon/theorems/THM-4397-rank-two-poisson-counterexample-symplectic-gauge-equivalence.md)
proves

```text
Phi_THM2044 = tau o Phi_Long o U_F.                    (1)
```

Both gauge maps are polynomial symplectomorphisms and lift to exact Weyl
automorphisms; the only mixed source relation is `F_xq=F_qx`.  Consequently
finite `A_2` quantizability or polynomial termination is invariant under this
gauge change.  The paper gauge is nevertheless much smaller:

```text
                  terms (R,T,D,S)       total degrees (R,T,D,S)
Long              (2,47,139,22)         (3,15,23,11)
THM-2044          (2,35,246,78)         (3,21,48,30)
```

It is therefore a useful implementation gauge for the open coupled-`D`
quantization calculation, without changing THM-2049's abstract Ore boundary
complex or its chosen nonterminating grade ladder.

The common first coordinate has no polynomial planar Jacobian mate by
THM-2045.  Rank two here means two canonical pairs and four commutative
variables.  Neither the preprint, `(1)`, nor its `A_4` appendix proves `DC(2)`
or planar `JC(2)`.

## 2. Kozma--Nitzan's gluing program

For independent bond percolation on a finite weighted graph, the paper's main
post-FKG conjecture is

```text
P(0 connected b)
 >= P(0 connected A) min_(a in A) P(a connected b).    (2)
```

It also states a stronger pre-FKG form.  Both are open in the paper.  Proved
content includes the pre-FKG form for `|A|=2`, several typed `|A|=3` cases,
and a recursive class where `0` is close to `A`.  Theorem 6 proves that the
weaker near-one gluing Conjecture 3 would imply `theta(p_c,Z^d)=0` for bond
percolation in every `d>=2`.

The proof architecture of the reduction is important.  It conditions on the
complete state of an interface, deletes closed interface edges, contracts
open interface components, and only then applies the local gluing statement
on the resulting graph.  This preserves boundary state instead of taking an
early scalar average.

### Source-level cautions

The v1 TeX needs several repairs before being treated as a polished proof
source:

- Theorem 2's `m_S` is displayed as an event but subsequently used as its
  probability.
- In Theorem 12 the printed equality `f(p_e)=-m` should generally be the
  sufficient inequality `f(p_e)<=-m`.
- The score `D_e` in Theorem 12 requires `0<p_e<1`; zero and one probabilities
  need deletion/contraction or a limiting argument.
- Several conditional probabilities require positive-conditioning or
  degenerate-case branches.
- Proofs of Theorems 7--9 and 11 are omitted, and no code or data accompanies
  the stated numerical evidence.

These issues do not by themselves refute the paper's reduction, but they bar
silently treating every displayed auxiliary theorem as established.

## 3. The later Lean claim

A fixed 2026 commit of `anthropics/formal-math` contains a large Lean project
claiming the unconditional result.  Its new finite-graph inequality is

```text
P(a connected b)>=1-t for every a in A
  => P(o connected b)>=P(o connected A)-t.             (3)
```

Taking `t=delta=epsilon/2` gives Kozma--Nitzan Conjecture 3.  The artifact says
that `(3)` follows from a conditioned slack hierarchy, then formalizes the
Kozma--Nitzan reduction and the required classical percolation inputs.  It
reports a pinned Lean/Mathlib build, comparator acceptance, no added axioms,
and no `sorry` outside two deliberate challenge-file placeholders.  It also
states explicitly that there has been no independent refereeing.

This intake inspected the statement surface and audit record but did not run
the roughly 13 GB build.  Therefore the repo status is deliberately narrower
than the artifact's self-description:

```text
EXTERNAL FORMAL CLAIM / SOURCE-REPORTED MECHANICAL CHECKS /
NOT INDEPENDENTLY REBUILT OR REFEREED HERE.
```

The claim is a material post-2024 status change and should be checked before
future text calls `theta(p_c)=0` open for bond percolation on `Z^3,...,Z^10`.
It does not settle the multiplicative inequality `(2)`.

## 4. Typed transfers to current frontiers

### Planar-JC row eleven and chart seams

| field | transfer contract |
|---|---|
| source | Kozma--Nitzan Lemma 10: condition on the full interface state, then delete/contract |
| target | the `Phi!=0` / `Phi=0` split after THM-4395 and the exceptional-quartic raw/admissible transgression split |
| map | interface state to pivot/seam stratum; deletion/contraction to localization, saturation, and the stage-preserving quotient |
| preserved predicate | each boundary stratum has a lawful continuation and overlap-compatible residual ideal |
| destroyed information | branch identity, conductor class, normal jet, and descent of a target representative |
| required sidecar | normalization/conductor labels, saturation ideal, pivot minors, transition maps, admissible-kernel image |
| cheapest test | compute row eleven separately at `Phi!=0`, `Phi=0`, and the first vanishing pivot minor; clear denominators and compare overlap ideals |

This is an architecture, not a probabilistic theorem inside algebraic
geometry.  It agrees with the incoming quartic calculation: unrestricted
representative motion spans the raw normal quotient, whereas the
solution-preserving kernel motion does not pay the retained seminormal class.

### Consumer-independent response

Kozma--Nitzan's three-terminal response system suggests looking for one exact
left inverse that cancels every row-eleven bracket/depth consumer.  Only its
linear algebra transfers.  Positivity and conditional correlation have no
automatic meaning over the characteristic-zero coefficient field.  A moving
minor or the `Phi=0` pivot loss is a seam, not a license to divide by `Phi`.

### LRC owner-first decomposition

The claimed formal proof localizes a union bound by assigning the first relay
before summing.  The precise LRC analogue is to retain endpoint owner, phase,
and arrival state before aggregating packet failure mass.  A Walsh/Bessel
bound transfers only on a genuine product cube and remains an average bound;
it cannot prove a universal lonely time without an exact lower quantum for any
nonzero failure mass.

That lower quantum actually fails on the full first open shell.  The exact
family

```text
w_m=(1,m,16m-1), c=C=(1,-16,1), m=5 (mod 6), m>=17
```

has primitive minimal relation norm `18`, but its live carrier component is
`3/[7(16m-1)]`, which tends to zero.  This is an algebraic obstruction to a
uniform-gap transfer, not an obstruction to `LRC(14)` itself.  The independent
audit and canonical status are handled separately from this source intake.

### Tournament firewall

The paper's directed analogue fails on four vertices.  An arc flip is not an
undirected open/closed edge, and tournament contraction generally leaves the
tournament category.  Any FKG-style tournament transfer must retain both
endpoint hypotheses, contraction orientation, and realization back in the
target class.  Otherwise the four-vertex directed gadget is a mandatory
hostile.

## 5. Next exact tests

1. Complete the weight-at-most-fourteen row-eleven residual on both `Phi`
   strata, then test the late weight-eighteen response separately; do not infer
   an all-weight obstruction from a capped boundary.
2. Use Long's smaller gauge for the coupled `D` Weyl relations, while treating
   THM-2049 termination as gauge-invariant.
3. For the exceptional quartic, construct the actual moving graph-family
   quotient before comparing its image with the raw fixed-`x` normal sidecar.
4. Replace the failed uniform-gap Walsh/Bessel route by an owner-conditioned
   inequality that can scale with the exact carrier component; any proposed
   bound must survive the family `w_m=(1,m,16m-1)` above.
