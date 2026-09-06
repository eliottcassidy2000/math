# Consumer-complete quotients and natural addresses

**Synthesis status, 2026-09-03.** This note combines PROVED/FINITE-EXACT
THM-4361, THM-4363--4368, and THM-4386--4394 with the proved quotient and
tournament controls THM-840, THM-851, THM-2294, THM-4088, and THM-4095. It
proves only the elementary factorization and residue lemmas below and records
exact corollaries of those theorems. It asserts no LRC-to-JC-to-tournament
transfer, no LRC(14) decrement, and no planar-Jacobian consequence.

## 1. The right question about a natural-number code

Let `q:X->Y` be an observation and `c:X->C` a named consumer. Then

```text
c=f after q for some f:q(X)->C
iff
q(x)=q(y) implies c(x)=c(y).                              (1)
```

For a sidecar `s:X->S`, the repaired state `(q,s)` determines `c` exactly
when

```text
ker(q) intersect ker(s) subset ker(c).                    (2)
```

These are representative-independence in one line. They are the arbitrary-
consumer form of THM-840's kernel criterion.

If `e:Y->N` is injective, then `ker(e after q)=ker(q)`. Consequently replacing
the `n`th odd square by `n`, or Cantor-pairing a structured finite tuple, is
harmless when the code is injective and its decoder is retained. It is a
renaming, not a compression theorem. If overlaps are allowed, the code has
become a genuine quotient and must pass `(1)` for the next consumer.

For deterministic operations `T_a:X->X`, start with `E_0=ker(q)` and refine

```text
x E_(r+1) y
iff x E_r y and T_a(x) E_r T_a(y) for every a.            (3)
```

On finite `X`, this decreasing chain stabilizes at the largest common
operation congruence contained in `ker(q)`. Equality in `E_r` is exactly
equality of every observation `q(T_w x)` for words of length at most `r`.
Thus a static code should not be called a state for a process until its
operation closure has stabilized.

A parallel interface lemma handles gluing. For finite maps `f:A->I` and
`g:B->I`, with fibre sizes `alpha(i),beta(i)`,

```text
|A times_I B|=sum_(i in I) alpha(i)beta(i).               (4)
```

Marginal nonemptiness or cardinality does not determine compatibility;
interface support determines existence and interface multiplicity determines
the count. In a stacked linear system, the mixed left-null relations are the
corresponding unpaid gluing consumers.

## 2. Exact source/map/consumer ledger

| Source | Quotient/map | Named consumer | Destroyed datum and hostile |
|---|---|---|---|
| THM-4363's four `h=420` rows | 828 statuses plus all 282 completed physical chains | first missing-component exit | missing physical prefix; `P=761,1015` have the same local role packet but teeth `(761,21),(1015,28)` and different exits |
| THM-4365's odd tail | the same constant quotient `Q` | metric exit `E_x` | residue/address and reciprocal scale; `P=48679,95873` even share `rho=1` but have different exits |
| THM-4367's strict-active cells | primitive ray `(a,b)`, equivalently `(kappa,b)` | metric exit `E_x=kappa/(14b)` | scale `g`; the quotient is exact for the exit, but the four scales over `29/1050` have different physical binders and split after `P->P+2` |
| THM-4396's three-sheet comb | finite two-factor spectrum plus exact ordered-pair geometry | local failure-mass upper certificate | third-sheet incidence; `(1,5,11)` has positive slack on an open interval for every finite pair degree |
| THM-4409's graph refinement | third-sheet component identities, contact graph, and separate vertex capacities for `g,1-g` | sharper local failure-mass upper certificate | overlap endpoints, edge integrals, and coupling of the two flows; `(1,19,79)` is already a strict fibre |
| THM-4364/4368 diagonal packets | boundary pair `(u,v)` or its triangular address | complete fixed-order trace | source exponents inside the fibre `mu(N,n0)`; reducing `(u,v)` to its primitive ray instead destroys trace scale |
| THM-4366's `U=0` pullback | two fixed-bank rank-one hierarchy flags | joint source restriction | ambient independence does not prevent the two opposite-module consumers from becoming proportional by `-4/3` on the selected source graph |
| THM-4361's six sign-sheet points | `Y=P^2` | signed `H_C` | sign sheet; the two points over one `Y` have opposite values |
| THM-2294/4088/4095 metric data | sign/order tournament | margin or arithmetic type | gaps, scale, valuation, offsets, and honest tie provenance; `{1,2}` and `{1,3}` have the same order tournament but different margins |
| THM-851 coloured arrows | defect and endpoint-node marginals | factorization/composition count | defect-to-factor-node incidence; the masks `8,24` give the exact finite hostile |

The common fact is not a shared carrier. It is that each tempting quotient
erases an interface read by a specified later consumer.

## 3. LRC: finite phase plus an unbounded metric counter

On THM-4365's tail, write

```text
1303P=47194n+rho,                  -23597<rho<=23597,
j=(rho+23595)/2 in {0,...,23596}.                        (5)
```

Advancing through odd parameters has the exact skew-product update

```text
j(P+2)=j(P)+1303 mod 23597,
n(P+2)=n(P)+1_[j(P)>=22294].                             (6)
```

The proof is direct: add `2606=2*1303` to `(5)`; crossing the centered
modulus increments the Euclidean quotient. Since `gcd(1303,23597)=1`, all
23,597 finite phases occur. Exactly 3,370 are active, but phase alone is not
a complete metric state:

```text
P=48679 and P=95873 have Q equal and rho=1,
yet their first exits are 18817/681506 and 37059/1342222. (7)
```

The residue selects the symbolic branch and binder type; `P`, equivalently
the address together with the residue, evaluates the reciprocal correction.
This is the precise architecture

```text
periodic finite state + affine natural-number address + reciprocal metric.
```

At `|rho|=3371`, tie provenance matters: the `3371` and `P` teeth both bind
at the safe open endpoint even though neither is active there.

THM-4367 identifies the exact compression available only on the strict-active
part. Put

```text
c=3371-rho,       g=gcd(c,P),       c=ag,       P=bg.
```

Then `a` is even, `b,g` are odd, `gcd(a,b)=1`, and

```text
a+1303b=3371kappa,          gkappa=1 (mod 14),
E_x=kappa/(14b).
```

Thus `(kappa,b)`, equivalently the primitive ray `a/b`, is the exact quotient
for the metric exit. Its fibre is one bounded arithmetic progression of
scales `g`; the sharp maximum size is 241. Equality occurs exactly at

```text
a=2,   b=401+6742t,   kappa=155+2606t,
t>=1,  t mod 7 in {0,1,2,5}.
```

This is the positive counterpart to the quotient warnings: a gcd reduction
really can be consumer-complete. It is nevertheless static. The complete
first repeated class

```text
P=11625,12675,13725,14775,              E_x=29/1050
```

maps under `P->P+2` to four distinct exits. Therefore exact knowledge of the
present metric exit is not an operation state. Its missing coordinate is the
scale `g` (or `P`), which also selects the physical binder. Moreover the
reciprocal ray swap `(a,b)->(b,a)` leaves the admissible parity sheet: it
turns the required even/odd pair into odd/even. The Stern--Brocot reflection
is an ambient symmetry here, not an internal symmetry of this LRC family.

## 4. Jacobian depth: two valuations, one triangular address

For fixed diagonal intercept `ell`, put `s=ceil(ell/2)`. THM-4364/4368 show
that a normalized nonzero single diagonal packet factors through

```text
(N,n0),                 N=c+e,        n0=b+c+2e,          (8)
```

and its complete order-`q` row trace has generating function

```text
F_q(z)=(-1)^(n0-s) z^(n0-s)(1-z)^(N-q-1).                (9)
```

At `q=0`, the valuations at `z=0` and `z=1` recover `n0` and `N`. Hence
neither start nor extent is decorative: together they are the exact trace-
equivalence state for one normalized packet. Complete fixed-`q` streams are
invertible integral transforms of the raw diagonal, so changing simplex
order does not add information when the entire stream is retained.

The positive boundary coordinates and their natural address are

```text
u=n0-s+1,       v=N,       Addr(u,v)=T(u+v-2)+u.         (10)
```

This is a genuine bijection, not a quotient: triangular blocks recover the
ordered pair exactly. A pair is realized by a source monomial precisely when

```text
N>=ceil(ell/3),          n0>=N,          n0+N>=ell.      (11)
```

Even this exact trace address forgets which source monomial produced it. The
complete fibre is indexed by

```text
max(0,ell-2N)<=e<=min(N,n0-N),
mu(N,n0)=min(N,n0-N)-max(0,ell-2N)+1.                   (12)
```

Every member expands to the identical packet. This is a proved overlap in the
source-to-trace map, not an accidental collision caused by numbering.

Boundary swap is an ambient reciprocal involution. With `tau=u+v-1`,

```text
Addr(u,v)+Addr(v,u)=tau^2+1,
F_R(z)=(-1)^(v-u)F_0(1-z).                              (13)
```

The diagonal `u=v` gives the honest fixed-point ties in odd `tau` blocks.
But, conditional on the starting pair satisfying `(11)`, the reflected pair
is source-realizable only when

```text
u>=ceil(ell/3),                 u-v<=s-1.
```

At `ell=10`, `(u,v)=(3,6)` is valid and its reflection is not. Also `(4,5)`
and `(8,10)` lie on one primitive rational ray but have different traces,
while `(4,5)` and `(4,6)` have the same strict order bit but different
traces. Thus the gcd/ray compression that is exact for THM-4367's LRC metric
consumer is specifically unsafe for this Jacobian trace consumer.

The other triangular structure is a rank clock. For
`rho=ceil(ell/3)`, the THM-4364 hierarchy bank at `(ell,m,d)` has exact
ambient-functional rank

```text
max(0,rho-max(0,ell+d-m)).
```

It enters one order at a time, has ramp sizes `1,...,rho`, and therefore has
`T(rho)` marked ramp positions. Consecutive fixed-row orders have a unit
Pascal minor. This makes the next-row tactic precise: evaluate only the newly
entering order first, then test whether source restriction makes it redundant
with a consumer on another module.

THM-4366 supplies both outcomes as controls. On the `U=0` pullback, two
ambiently distinct rank-one fixed banks become proportional by `-4/3` after
source restriction, so the new order is an alternative certificate rather
than a new cut. Nevertheless the six strict row-ten source points are killed
by a row-eleven residual coprime to their cubic. Ambient rank, restricted
rank, and the mixed source ideal are three different consumers.

THM-4361 adds the same warning in a sign direction. On its six-point locus,

```text
H_C=-P d(Y)/D,                     Y=P^2.                 (14)
```

The signed value does not factor through `Y`, but `H_C^2` and the predicate
`H_C=0` do. A sidecar can be essential for one consumer and irrelevant for
another. Separately, each projected-depth block can be consistent while the
stacked system is not: the mixed cokernel, not two marginal yes/no flags, is
the correct gluing object.

## 5. Tournament lesson without a false transfer

THM-2294 makes zero determinant an honest tie and treats symmetric character
data as a colour, not an orientation. THM-4088/4095 prove that order
tournaments forget metric and arithmetic data. THM-851 shows that composition
reads defect-to-factor incidence rather than marginal colour volume.

Accordingly, a tournament is a valid natural-number state only relative to a
consumer invariant under:

```text
metric variation inside one orientation cell,
tie-gauge variation,
and rearrangement of erased interface incidence.          (15)
```

Failing any test in `(15)` requires the corresponding sidecar; it does not
license an arbitrary orientation completion.

The reduced rational `a/b` can be drawn as the directed edge `a->b` between
natural-number vertices. With the labelled endpoints retained, reciprocal
reflection reverses this edge and no information is lost. Keeping only its
orientation inside a larger tournament forgets both endpoint labels and
metric value; allowing unreduced endpoints also introduces a scale fibre.
THM-4367 proves that one particular LRC consumer is constant on that scale
fibre. THM-4368 proves that one particular Jacobian consumer is not. This is
the rigorous overlap behind the analogy; there is no cross-problem theorem.

### The raw LRC carrier and the one honest tournament bit

THM-4386 supplies a second exact addressing example.  For a primitive
three-speed comb, nearest-integer data `n` modulo the common translation
`n -> n+tw` have the chart-independent image

```text
C=w cross n in Lambda_w={C in Z^3:C dot w=0}.            (16)
```

The map `Z^3/Zw -> Lambda_w` is an integral isomorphism.  Thus `C`, not a
chosen relation chart, is the consumer-complete component address for the
local length `L_w(C)` and the raw-carrier sum.  If a primitive relation
`c dot w=0` is selected, a Bezout section rewrites one defect fibre as

```text
C=C_delta+k c.                                          (17)
```

Here `(delta,k)` is a chart-labelled affine coordinate: changing the Bezout
section translates `k`, and changing the relation changes the whole chart.
It is not an intrinsic group law on physical components.  THM-4393/4394's
multiple-relation rays make the warning literal: the same raw `C` can receive
different `(delta,k)` presentations, so sector totals must not be added.

There is nevertheless one intrinsic binary relation in the distinct-owner
gate.  Write

```text
o_i=-w_i^(-1)n_i mod 3,        {o_1,o_2,o_3}=F_3,
i -> j iff o_j-o_i=1 mod 3.                              (18)
```

This orients the three labelled speed vertices as a directed 3-cycle.  An
affine relabelling `o_i -> o_i+t` preserves it, while negation reverses it.
The two live raw-carrier cosets modulo three are therefore the two opposite
orientations.  This is an honest tournament observable because the pairwise
rule is intrinsic after the positive-speed gauge and owner labels are fixed.

The coefficient residues say exactly when a relation chart can still see
that orientation.  Since `w_i` is a unit modulo three, a primitive relation
has either all three `c_i` nonzero modulo three or exactly one zero modulo
three; two zeros would force the third.  In the all-nonzero case the three
`c_iw_i` are equal, so distinct owners force

```text
delta=c dot n=0 mod 3.                                  (19)
```

Each defect fibre then retains two live `k` classes, the two reversed owner
cycles.  The scalar defect forgets their orientation.  If exactly one
coefficient is zero modulo three, the other two weighted residues are
opposite; distinct owners instead force `delta!=0 mod 3`.  Its sign relative
to the zero-coordinate label retains the orientation, and only one `k` class
is live on each defect fibre.  This is a proved elementary residue dichotomy;
it does not say that either relation type bounds the comb by itself.

THM-4392 is the dual version of the same information loss.  Fourier transform
of the **unoriented union** of the two raw cosets gives character weights `6`
on equal weighted residues and `-3` on all-distinct residues.  The transform
retains the two-coset union exactly but no longer names which directed owner
cycle a primal component used.  Conversely, keeping only the tournament
orientation discards raw scale, the exact lattice point `C`, and the metric
length `L_w(C)`.  It therefore cannot determine failure measure.  There is no
LRC-to-tournament theorem here: the useful statement is a loss ledger linking
one genuine orientation bit to the raw and Fourier representations.

THM-4405 stress-tests the orientation/defect distinction across all 40
one-zero presentation shapes through norm 20.  Every new row maximum lies in
the all-unit `(1,1,2)` sector, yet the common equality comb `(1,5,11)` moves
its two masses between defect pairs `+/-1`, `+/-2`, `0`, and `+/-3` as the
relation chart changes.  What persists is the two-entry raw-carrier
dictionary; what moves is the scalar coordinate assigned by a chosen
relation.  In tournament language, the directed owner cycle survives while
its numerical chart label does not.  This is a finite structural compression,
not a claim that future coefficient norms inherit the same winners.

THM-4413 shows that the same orientation sidecar has a second, genuinely
metric consequence.  For a raw carrier `C`, the three roof margins satisfy

```text
p_i=3(w_j+w_k)-14|C_i| = |C_i| (mod 3).
```

The distinct-owner tournament is equivalent to `3` not dividing any `C_i`;
odd speeds make every `p_i` even.  Hence the physical owner lattice misses the
continuous roof boundary by at least two integer units.  The orientation has
not determined the metric length—it still loses height and the exact carrier—
but, together with parity, it excludes tangency and yields the sharp adaptive
floor `min(p_*,6w_2)/(14w_2w_3)`.  This is a useful pattern for other problems:
a finite label need not encode a magnitude to force transversality when the
boundary functional has a rigid residue alphabet.

THM-4396 gives a pointwise version of the same discipline.  After smoothing
two sheet indicators, write `X=f_i f_j-g_i g_j`.  Forgetting the exact third
indicator but retaining only `0<=f_k<=1` has the optimal envelope

```text
sup_(0<=z<=1) zX=X_+.                                   (20)
```

This “quotient-sharp” statement is relative to the declared information loss;
it is not global optimality of the resulting numerical bound.  The equality
comb `(1,5,11)` proves the distinction: an open region has `f_k=1` and
`X<0`, so replacing `f_kX` by `X_+` loses a strictly positive amount for every
finite pair of Fejer degrees.  A smaller bound must restore some third-sheet
incidence, couple pair views, or change the representation.  Merely raising
the cutoff cannot repair a quotient that has erased the sign-bearing
interface.

THM-4409 performs the next refinement step without restoring the full triple
partition.  On each sheet it retains the bipartite contact graph between pair
components and third-sheet components, together with vertex integrals of `g`
and `1-g`.  Actual overlap integrals form a feasible fractional flow, so two
maximum-flow capacities preserve the upper-bound consumer.  At degree zero
the formula collapses to a rational component-length capacity.  The sidecar
is exactly sufficient on `(1,5,11)` because every graph is a nested matching,
yet it is not consumer-complete: `(1,19,79)` has physical mass `108/10507`
and best graph bound `8/553`.  This is a useful middle quotient—strictly more
informative than sheet blindness and strictly cheaper than overlap geometry.

THM-4414 reveals that the max-flow language at degree zero is itself a
temporary representation.  Each danger-sheet family is six-separated, and
intersecting two such families preserves six-separation.  A sharp general
factor-two lemma therefore makes all edgewise minimum loads feasible at once:
the contact graph has no Hall competition at any height.  What remains is the
meet-envelope hinge

```text
min(|I|,|J|)-|I intersection J|
 =min(|I minus J|,|J minus I|).
```

This is the same operation seen in older residual-capacity and interval-Gram
work, now with a clean loss ledger: graph plus lengths forgets crossing
placement, while the other two pair-overlap facets restore the exact
box-spline roof.  The graph is a sparse forest, so tournament completion would
destroy contact zeros and hinge magnitudes without helping the mass consumer.
The open arithmetic target is consequently an explicit choice among three
raw-carrier projection sums, not a generic network optimization.

THM-4420 finds a lane where that choice disappears.  On
`w=(1,m,2m+sigma)`, the kernel relation and a roof width strictly below one
force every live carrier onto `k(-sigma,-2,1)`.  The exact natural address is
then just the ordered list of positive `k` with every multiple of three
deleted, plus reflection sign; counting this list bounds all three projection
sums simultaneously.  This is an example of a quotient becoming complete
because a transverse integer coordinate has width below one.  The extension
to `(a,m,2m+sigma a)` pinpoints the next obstruction: divisibility only forces
that coordinate into `aZ`, so an additional residue/phase sidecar may be
needed before ordinal counting is again lossless.

THM-4421 tests whether THM-4413's residue-forced transversality transfers to
the exceptional-quartic collision period.  In one integral compiler it does:
`F(h)=8O(h)+18E(h)`, so the nonzero odd-coefficient residue forces the sharp
gap `F(h) in 2Z minus 6Z`.  The hostile is the missing invariance sidecar.
The visible kernel direction `9x-4x^2` preserves the collision, actual
polynomial transgressions annihilate the first Smith class, and the lawful
target rescaling `C -> C/3` removes all three-torsion.  Unlike LRC, where the
speed/owner lattice and roof functional are canonical, this Jacobian lattice
is a coordinate choice.  The transferable principle is therefore stricter:
a residue alphabet forces geometric transversality only when its integral
lattice is intrinsic and stable under every allowed coordinate change.

## 6. Generated sharp tasks

1. **LRC operation closure.** THM-4367 resolves the static exit-collision
   classification. Start from `(kappa,b)`, adjoin the least sidecar that makes
   `P->P+2` deterministic, and run refinement `(3)`. The four-point
   `29/1050` fibre is the first mandatory hostile; parameters separated by
   47,194 test whether a proposed finite phase has silently dropped the
   affine counter.
2. **LRC physical state.** Classify the joint consumer `(metric exit, binder
   set, binder index)` on strict-active and boundary cells. Determine whether
   `(kappa,b,g)` is minimal after the two boundary flags are included, rather
   than merely sufficient.
3. **Jacobian rank-clock step.** For each source stratum, evaluate the single
   hierarchy order entering at the next row before rebuilding all older
   consumers. Then compute its class modulo the restricted old-bank span and
   the mixed left-null quotient; THM-4366's `-4/3` relation is the redundancy
   control, while its coprime row-eleven residual is the genuine-cut control.
4. **Boundary-orbit enumeration.** Count the `(u,v)` for which both a packet
   and its reciprocal reflection satisfy `(11)`, by triangular address block;
   retain the two source-fibre sizes `mu(u,v)` and `mu(v,u)`. Fixed points must
   remain honest ties rather than being oriented by convention.
5. **Linear-combination firewall.** THM-4368 classifies single normalized
   packets, not sums. Find the smallest two-packet cancellation sharing one
   boundary valuation, then determine what Hankel/valuation sidecar recovers
   the multiset or prove that no bounded sidecar can.
6. **Tournament incidence closure.** At order eight, refine endpoint colours
   by actual defect-to-factor incidence and compare operation closure with the
   full deck. Test whether triangular block reflection supplies a useful
   involutive hostile without treating the analogy as a carrier map.
7. **Foundational firewall.** Promote `(1)--(4)` and injective ordinal
   invariance only after checking their exact set-theoretic quantifiers. These
   lemmas separate renumbering, static quotienting, and composition-compatible
   state before another domain-specific code is introduced.
8. **Component-network height frontier.** Determine whether the minimum of
   THM-4409's three rational degree-zero capacities stays at most `6/77` at
   arbitrary height, or find the first exact hostile above `79`.  Separately,
   refine the strict fibre `(1,19,79)` by adding one coupled-flow or edge-mass
   coordinate and test whether it is the minimal consumer-complete repair.

The reusable rule is: name the consumer, test one quotient fibre, retain the
first interface it reads, and only then encode the surviving structured tuple
by a natural number.
