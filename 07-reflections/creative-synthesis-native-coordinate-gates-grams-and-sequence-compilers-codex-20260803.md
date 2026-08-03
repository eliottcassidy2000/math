# Creative synthesis: native coordinates, first gates, Gram sidecars, and sequence compilers

**Status:** current research synthesis across
[THM-3316](../01-canon/theorems/THM-3316-prime-right-boundary-interpolation-forces-scalar-rigidity.md)--[THM-3324](../01-canon/theorems/THM-3324-tournament-deletion-response-gram-ordered-join-compiler.md).
`LRC(14)`, `JC(2)`, `DC(2)`, `FC(3)`, and the global tournament problems
remain OPEN.  The point of this note is to identify what genuinely moved,
which new objects are reusable, and which next tests can decide something.

## 1. What the new block actually establishes

| lane | new proved object | frontier after the result |
|---|---|---|
| LRC prime cells | prime right-boundary interpolation is a scalar ramp; sparse blockers have a weighted unique-protector fragility test | composite `14`, physical lifting, and the global rung remain open |
| LRC projected atlas | fourth `z1=216` ruler prefix closes; full-dual arity is smaller than integral-template arity on two packets | ledger `373153`; `349` wall rows in `29` families; next complete family has `19` rows |
| planar Jacobian | the exceptional quadratic critical deck is a connected algebraic etale germ over both clutch slopes | no global component, rational section, mate, or inverse |
| divergence obstruction | a separate gradient-unimodular family has the exact annihilator ladder `(P^r),(P^(r-1))` | poles still block a polynomial mate; no carrier map to the critical deck |
| factorial conjecture | mixed triangle moments have one Hesse kernel; support `<=4` is excluded; the support-five degree-21 map has exact rank `10980` and quotient dimension `1670` | full support five is still the sole open chart; degree `29` is only count-eligible |
| tournaments | the switching second moment is `P` plus a deletion Gram; the full two-coordinate response Gram closes under joins and diagonalizes | third tensors, substitution response, and SCC-order sidecars remain open |
| sequences | Hesse coefficients have a one-sum/lattice recurrence; tournament walks and repeated-join Gram channels have product-power closed forms | bit complexity and broader operation interfaces are separate questions |

These are not shadows of one theorem.  They share a productive research
procedure, while their carriers and target predicates remain distinct.

## 2. Complete the native coordinate before reading complexity

Three lanes now delimit the rule sharply.

### LRC: add the affine constant

THM-3320 begins with an integral zero-constant majorant.  Two packets appear
to require three threshold layers.  Once the equality-row dual is written in
its full affine form,

```text
sum_(t in T) 1[c(P)>=t] <= a+sum_(i in P)w_i,               (1)
```

both have the smaller certificate

```text
T=(1,5),       a=1/2,       w=(1/2,1/2,1/2,3/2).           (2)
```

Exact proper-support tables prove that support two, not three, is intrinsic
to the full common-table problem.  The constant row was not decoration; it
was a missing coordinate in the dual.

### Jacobian: release both external slopes

At the THM-3306 fixed slice, the two linear-subresultant coefficients have a
unit Jacobian in the internal coordinates `(x,c)`.  THM-3319 releases both
external clutch parameters `(d,k)`.  The relative Jacobian criterion then
makes `V(a,b)` etale over the clutch plane at the exceptional degree-36
point.  The old fixed deck is the special fibre of a moving algebraic germ,
not an isolated accident.

### Factorial: reject a coordinate that is not a symmetry

The Hesse surface `uv=r^3` has a formal torus, but THM-3321 checks the moment
functional rather than the carrier alone.  The five nonzero pure terms of
`M_3` have three distinct torus weights.  The torus therefore cannot normalize
a second coefficient.  Completing a coordinate system is lawful only when
the added action preserves the target predicate.

The combined rule is:

```text
before declaring high arity or isolation,
  add the native affine/relative coordinate;
before quotienting by it,
  test the actual functional or predicate.                 (3)
```

## 3. Test the first obstruction gate, not the most sophisticated invariant

The Jacobian results now give a clean two-gate taxonomy.

```text
THM-3319: gradient ideal proper
          -> the divergence class mu(P) is not defined;

THM-3318: gradient ideal unit
          -> mu(P) is defined but is nonzero torsion.       (4)
```

This prevents two common wastes: computing a late invariant before checking
its domain, and treating all negative mate results as the same mechanism.
The carriers differ and no reduction between the families is known.

The same procedural gate appears elsewhere:

- in FC, projective boundary closure must precede an affine normalization;
- in LRC, exact marginal feasibility must precede naming circuit arity; and
- in tournaments, the operation must decide whether `P`, `s=(P,zN)`, scalar
  `D_T`, or the full response Gram `Gamma_T` is the first closed interface.

## 4. Gram sidecars measure what scalarization forgets

For tournament switching, THM-3322 proves

```text
E[N_d(z)N_d(w)]
 = P-determined terms - 2*D_T(z,w),                         (5)

D_T(z,w)=sum_i p_i(z)p_i(w),                               (6)
```

where `p_i` is the centered characteristic polynomial after deleting vertex
`i`.  At order seven, two tournaments have the same `P` but a rank-one
difference of deletion Grams.  Thus `D_T` is precisely the second-order
coordinate lost by spectrum.  THM-3324 then asks what the **next operation
consumes**.  Put

```text
r_v=(p_v,z*nu_v)^T,       Gamma_T=sum_v r_v(z)r_v(w)^T.    (7)
```

The answer is not scalar `D_T` but the full two-coordinate response Gram:

```text
Gamma_(X join Y)
 =H_Y(z)Gamma_X H_Y(w)^T+H_X(z)Gamma_Y H_X(w)^T.           (8)
```

Order-six masks `73,83` agree on `(P,N,D,E,E^T)` and differ only in the final
`F=sum z*nu(z)w*nu(w)` block; joining a singleton exposes that difference in
`D`.  This is a sharp example of a quotient that suffices for one consumer
(the switching second moment) but not for its next native operation (join).

This is structurally comparable, but not identical, to two other sidecars:

- LRC's common-table primal remembers whether separately valid tail demands
  coexist under one marginal system; and
- the Jacobian quadratic deck remembers the two conjugate common roots that
  scalar elimination cannot select individually.

In all three cases the right question is not “can the scalar be sharpened?”
but “what finite labelled Gram/table/deck restores the consumer?”

## 5. Two exact sequence compilers

The factorial lane now has the rational kernel

```text
C(a,b)=[X^aY^b](1-X^3-Y^3-3XY)^(-1),
mu(a,b)=2*a!*b!*C(a,b)/(a+b+2)!,                            (9)
```

with

```text
C(a,b)=C(a-3,b)+C(a,b-3)+3C(a-1,b-1)+[a=b=0].             (10)
```

This replaces a six-index barycentric expansion by either a one-sum formula
or constant work per lattice cell.

The tournament lane has a different compiler.  For a selected switch,

```text
G(2z)=N(z)/(P(z)-zN(z)),                                   (11)
```

so its total-walk sequence is constant-recursive.  Under order join,

```text
P_join +/- zN_join=(P_1 +/- zN_1)(P_2 +/- zN_2).          (12)
```

Equation `(12)` is an efficient closed form for repeated joins, but its
commutativity forgets factor order: `K1 triangle C3` and `C3 triangle K1`
have the same `(P,N)` and different source/sink structure.  A fast compiler
is not automatically an injective classifier.  THM-3324 extends the compiler
to every diagonal Gram channel.  If `u_a=P+a*zN`, then for
`J=T_1 triangle ... triangle T_k`,

```text
Gamma_hat_(J,ab)
 =sum_i Gamma_hat_(T_i,ab)
        product_(j!=i)u_(T_j,a)(z)u_(T_j,b)(w).           (13)
```

For `k` identical factors this is
`k*Gamma_hat_(T,ab)*(u_(T,a)(z)u_(T,b)(w))^(k-1)`.  The
closed form accelerates exact evaluation while making the order loss
transparent.

## 6. LRC: the next combined test

THM-3316 and THM-3317 should not be merged into THM-3320 by vocabulary alone.
The prime-cell theorems concern boundary interpolation/protectors; THM-3320
concerns a projected common-status table.  The cheapest lawful bridge is a
typed experiment on the next complete family:

1. batch the nineteen `gcd72/L7056` rows by exact capacity signature;
2. for each surviving signature, record the unique endpoint protector and
   the sum of gcd shift costs from THM-3317;
3. split prime local factors, where THM-3316 forces a scalar ramp, from the
   composite `14` interaction;
4. run the ordinary common-table screen unchanged; and
5. ask only whether the protector/ramp packet predicts a valid dual or a
   smaller residual universe.

Success means an explicit map from a cell/protector packet to a status-table
inequality.  Failure means a minimal row/signature on which the prime-local
data agree but common-table feasibility differs.  Either outcome identifies
the missing composite sidecar.  No heuristic score transfer is admissible.

## 7. Four decisive next probes

### LRC anchor

Close or sharply stop the full `gcd72/L7056` family using capacity-signature
reuse.  The deliverable is one of:

- nineteen empty exact screens and the new ledger;
- a finite list of residual signatures with exact feasible tables; or
- a demonstrated batching obstruction that identifies which row label is
  not captured by capacity data.

### JC local-to-global probe

Let `Z=V(a,b)` in the boundary-localized `(x,c,d,k)` space.  Saturate its
closure by the known owner and finite-clutch walls, then compute:

- the degree and component containing the THM-3319 etale germ;
- the first branch/discriminant locus of `Z -> A^2_(d,k)`; and
- the first parameter curve on which the quadratic critical deck ramifies or
  changes gcd degree.

The hostile stopping condition is a wall or extra component before any
global rational section appears.

### FC support-five probe

THM-3323 shows that degree `21` attains the universal ceiling exactly: rank
`10980`, quotient dimension `1670`, and no hidden early degeneracy.  Do not
build the full degree-29 rectangle yet.  Continue the sparse quotient degree
by degree from the frozen pivot basis and decompose its dual through the Hesse
surface `uv=r^3`, while imposing THM-3310's phase and modulus barriers.  The
decisive output is an exact saturation certificate, a structural recurrence
for the quotient, or its first departure from the formal count.  Degree `29`
remains only the first count-eligible fallback.

### Tournament wildcard

The join law and its diagonal closed form are now solved.  Derive the third
marked-response tensor under join and compare it with the third switching-cube
moment, whose degree-two Walsh characters first close in triangle-shaped
triples.  In parallel test one nontransitive substitution quotient while
retaining quotient owner and SCC order.  The hostile target is a minimal pair
agreeing on the proposed tensor but diverging after one more substitution or
in a triangle Walsh coefficient.

## 8. Connection boundary

The synthesis transfers procedures, not conclusions:

- the LRC affine constant is not a Jacobian clutch parameter;
- the Jacobian quadratic deck is not the tournament deletion deck;
- the Hesse kernel does not make the factorial ideal a tournament spectrum;
- prime-cell rigidity has no physical LRC lift yet; and
- none of these results proves `LRC(14)`, `JC(2)`, `DC(2)`, or `FC(3)`.

The durable advance is a more discriminating workflow: complete the native
coordinate, test the first gate, retain the smallest labelled sidecar, derive
the operation response, and freeze the first hostile where closure fails.
