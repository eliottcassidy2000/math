# Five-core shells, reversed addresses, and the neglected Pollock five

**2026-11-10. PROVED IN THIS NOTE:** exact high-bit shell transport,
its complete jump multiplicities, and a joint dyadic rectangle law.
**FINITE-EXACT:** small literal-tournament controls and larger arithmetic
controls in `shell_checks.json`. **OPEN:** the target-occupancy conclusion
needed for Chan–Pak. No new all-height coverage is claimed.

The [continued target-selection investigation](continuation/PROOF.md)
derives the full internal-bit derivative and centered doubling inequality,
then studies3560 actual witnesses to the322 inherited small-parent holes.
It records the useful one-step joint signal, two failed completion ideas,
and the exact source-preservation gate still needed.

The historical search found a more useful additional five than another
five-object symmetry: the June Pollock **five tetrahedral summands** work.
It changes a coverage problem into an obstruction to simultaneous failure
at related addresses. Transferring that move to the five-core compiler
exposes an exact reversed-address coordinate. This coordinate is absent
from a marginal jump count and from a mere list of five core names.

## Inheritance and recovered side lanes

The closest mechanism is the inherited seven-coordinate positive compiler
in [POSITIVE_COMPILER.md](../../chanpak_five_core_bridges_20261110/algebra/POSITIVE_COMPILER.md).
The hostile is a prescribed target with no matching source: every even
jump can occur somewhere without a source at the required starting H.
The corrected near miss is the old golden projection: contraction in one
coordinate can hide expansion in the conjugate coordinate. The least-used
sidecar recovered here is the **joint defect pair**, rather than a scalar
defect or a marginal set of increments.

Historical sources actually read in this lane, beyond Session252 and the
latest five-core packet:

* March23 [heap/tournament/tableau notebook](../../../../07-reflections/heap-tournament-tableau-triangle.md):
  retain its suggestion to compare conflict pieces and admissible schedules;
  its stronger geometric identifications are not dependencies here.
* May16 [Eureka/Zeckendorf/simplex-cuboid notebook](../../../../07-reflections/eureka-zeckendorf-simplex-cuboid.md):
  separate the *representation law* from the number of objects it represents.
  This recovers the choice between bounded-arity and normal-form coverage.
* June13 [HYP-2491](../../../hypotheses/HYP-2491-pollock-tetrahedral-defect-pair-descent.md),
  [HYP-2497](../../../hypotheses/HYP-2497-pollock-sierpinski-carry-scout.md),
  and their [defect-pair reflection](../../../../07-reflections/pollock-tetrahedral-defect-pair-descent.md):
  the useful exact map is displayed in Section5 below. Their tail conjectures
  remain inherited OPEN statements, not inputs to our theorem.
* The later [Farey/Fibonacci/additive-basis synthesis](../../../../07-reflections/farey-fibonacci-additive-basis-carrier-s169.md):
  distinguishes smoothing, bounded arity, sparse-carry normal forms, and
  Farey addresses. It motivates asking which resource the core bank uses.
* October7 [frieze/Collatz reset note](../../frieze_collatz_reset2_20261007.md):
  its five-flip coordinate change preserves a marked word; a matrix ear
  identity alone need not preserve legal arithmetic endpoints. This directs
  the proposed five-core rewrites toward actual addresses and native cuts.
* October18 [golden visible/conjugate embedding](../../fivefold_boundary_20261018g.md):
  keep both coordinates when trying to decrease a height. Our binary
  reversal is an explicit discrete companion coordinate, not an asserted
  identification with Galois conjugation.

The live concepts are core/cut response, fixed-address shell, reversed
address, jump multiplicity, simultaneous target failure, and proper height.

## 1. The exact shell identity

Use the inherited function

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod2).

For `k>=1` and `0<=n<2^k`, let `r_k(n)` reverse the **k binary digits**
of n, retaining leading zeroes. Put `t_k(n)=floor(n/2^(k-1))`.

**Proposition1.**

    g(2^k+n)-g(n) = 2^k-r_k(n)-t_k(n).                 (1)

The top bit is necessary. For k=1, the two differences are2 and0.
For k>=2 write n=2q+b. If b=1, the difference reduces to the
`k-1` case at q. If b=0, it is the `k-1` difference plus `2^(k-1)`;
the added high prefix has even parity at this stage. Since

    r_k(2q+b)=b*2^(k-1)+r_(k-1)(q),
    t_k(2q+b)=t_(k-1)(q),

these are exactly the two cases of (1). This proves the identity at every
length. The recurrence is also valid at q=0, by its initial values.

The result is a genuinely two-coordinate shell law: a fixed high bit adds
`2^k`, but its H contribution depends on the low word read backwards.

## 2. Every even jump occurs, with exact multiplicities

Since the top bit of n is the low bit of r, (1) becomes

    D_k(n)=2^k-r-(r mod2).

**Corollary2.** On the complete k-bit address universe,

    {D_k(n):0<=n<2^k} = {0,2,4,...,2^k}.

Each interior value occurs twice;0 and2^k occur once. Indeed an interior
value2j has exactly the reverse-address preimages

    r=2^k-2j,       r=2^k-2j-1.                       (2)

Discard an out-of-range value at the two endpoints. Reversing these
addresses gives an explicit inverse for every prescribed even jump.
No search through all k-bit words is needed.

This is not a prescribed-source theorem. The baseline `g(n)` changes with
the inverse address in (2), and only the resulting correlated pair is
available. That limitation is the new proof target, rather than a reason
to discard the exact jump law.

## 3. Transport through any coloured tournament core

For the inherited positive signature `(h,d,L,i,o,s,e)`, write

    Phi(A,B)=h+d(g(A)+g(B)-2)+LAB+iA+oB
               +s(A mod2)+e(B mod2).

Here every coordinate is fixed by an actual core and apex cut; `d=H(C)`.
Fix B and put `alpha=LB+i`. The parity of n is unchanged by adding2^k.
Therefore (1) gives the all-height law

    Phi(2^k+n,B)-Phi(n,B)
      = alpha*2^k + d(2^k-r_k(n)-t_k(n)).              (3)

**Corollary3.** At fixed B the high-bit jump spectrum is exactly

    alpha*2^k + 2d*{0,1,...,2^(k-1)},                 (4)

with the same endpoint/interior multiplicities. There is a symmetric right
flank law, replacing alpha by `LA+o` and n by the low right address.

For the five cores, the grid spacings are respectively

    2,6,10,18,18.

The two H9 cores therefore have the same *marginal* jump spacing. Their
marked cuts, linear offsets and baseline correlations distinguish them.
This is an exact explanation of why a histogram of available increments
cannot replace the two different core charts.

With leading-zero padding retained, (3) is a literal one-vertex extension:
insert one new first vertex in the transitive left flank; it beats every
old source vertex and loses to the apex. The old flank has exactly k
vertices, even when n has fewer significant bits. Omitting neutral padding
preserves H but obscures that the operation added exactly one vertex.

## 4. The joint shell geometry has exact dyadic occupancy

There is useful structure before projecting to the jump. Consider

    P_k={(n/2^k,r_k(n)/2^k):0<=n<2^k}.

**Proposition4.** Every half-open dyadic rectangle inside `[0,1)^2` of side lengths
`2^-a,2^-b` with `a+b<=k` contains exactly `2^(k-a-b)` points of P_k.

Proof: its horizontal address fixes the a highest bits of n. Its vertical
address fixes the b highest bits of r, hence the b lowest bits of n.
These bit sets are disjoint, leaving precisely `k-a-b` free bits. QED.

This supplies a complete bounded-prefix/bounded-suffix joint law at every
scale. The full H baseline is not constant on those rectangles: it depends
on interior weighted bit interactions. Accordingly Proposition4 does not
prove target occupancy. It suggests a sharp experiment: which joint cells
remain entirely forbidden after **all five** core/cut responses are allowed?
A finite-state single-address picture is already obstructed by the
inherited selector theorem; the interior baseline remains a real arithmetic
coordinate, not a dispensable hidden state.

## 5. The Pollock five gives the right kind of next statement

Let `Te_k=k(k+1)(k+2)/6` and let D4 denote integers with no representation
as a sum of at most four tetrahedral numbers. The old Pollock map is

    N=Te_k+r = Te_(k-1)+(r+tri(k)),
    tri(k)=Te_k-Te_(k-1).

Thus two competing fifth-summand choices fail together only when
`r` and `r+tri(k)` are both in D4. It is a **joint defect-avoidance**
problem. No theorem on the number of marginal gaps alone settles it.

The transfer to Chan–Pak is now explicit:

|Connection component|Here|
|---|---|
|Source object|Pollock's two neighboring shell decompositions|
|Target object|A core/cut's old and high-bit-extended flank witnesses|
|Map|Replace the shell difference by the exact expression (3)|
|Preserved predicate|Both observations come from actual legal constructions|
|Lost information in a scalar quotient|Which baseline belongs to which jump|
|Necessary sidecar|`(core,cut,B,k,n,r_k(n),Phi(n,B))`|
|Cheapest decisive test|Enumerate simultaneous target misses on joint address cells, rather than extending a marginal jump histogram|

One precise open gate is the following. Given an odd N above the inherited
finite base, find an actual core/cut, B,k and n with `0<=n<2^k` such that

    Phi(n,B) = N-(LB+i)*2^k
                    -d*(2^k-r_k(n)-t_k(n)).             (5)

The right side is a proposed **baseline**, and (5) is the correlated
occupancy obligation. Corollary2 ensures any admissible marginal jump has
one or two explicit inverse addresses, so a target-driven search can work
from candidate jump values directly. Equation(5), by itself, remains an
equivalent reformulation of a high-bit witness, not a solved theorem.

A stronger useful next result would be a tail exclusion for simultaneous
failure of several chart/address choices linked by (1) and reversal. This
is the exact analogue of Pollock's no-long-defect-pair lemma, and weaker
than proving any one chart universally covers its own target band.

The golden/Platonic/Fermat/Ogg links already recovered by Session252 are
structural background. The additional five here contributes a specific
research operation: spend the extra chart to avoid a **paired failure**.
The present contribution is the exact shell transport needed to state
that operation without guessing its arithmetic gap.

## 6. Reproduction and finite scope

Run `python3 -B 05-knowledge/results/five_core_web_20261111/fivefold/shell_checks.py`.
The checks compare recurrence and a separately expanded binary-pair formula,
all shell identities through k14, dyadic rectangles through k9,44 inherited
core/cut signatures at B0..7 and k1..8, and literal small tournaments with
leading-zero padding. The finite checks support the explicitly proved
identities; they do not claim an infinite target cover. The output reports
the actual universe and the witness for the missing-top-bit hostile.
