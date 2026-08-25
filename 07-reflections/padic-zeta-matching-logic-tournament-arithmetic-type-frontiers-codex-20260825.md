# P-adic zeta, matching logic, and arithmetic type through tournament quotients

**Session date:** 2026-08-25
**Status:** research synthesis with four promoted exact theorems.  The external
`22` p-adic singleton claims remain **AUTHOR-CLAIMED / UNREFEREED**.  LRC(14),
planar `JC(2)`, Mahler `3/2`, the Rule-30 prizes, and the named open
irrationality/transcendence problems remain **OPEN**.

## 1. Outcome first

The incoming sources produced four durable results, each at a different
logical boundary.

1. **THM-4088 (order/type firewall).**  Every scalar strict-order tournament
   is transitive and its fibre realizes every prescribed finite assignment of
   rational, exact-degree algebraic, and transcendental labels.  Conversely,
   any positive-margin lonely-runner witness persists on an explicit open
   ball, so times of each arithmetic type are dense there.  The AP13 unit
   skeleton proves that equality-only witnesses need not have this property.
2. **THM-4089 (global margin obstruction).**  Both continuous variables of
   the pinned modular-template margin have unique global optima.  At those
   optima the next cells `(2,31),(3,13),(5,7),(7,5)` remain rigorously
   negative.  Progress beyond the claimed `22` cells therefore needs a
   changed template, saving, or cost—not better tuning of `xi,Y`.
3. **THM-4090 (minimal sort-flow obstruction).**  Chen--Rosu's three-sort
   matching-logic counterexample reduces to two sorts and one symbol.  One
   sort is complete in the cited fragment; two sorts already fail because a
   total semantic constraint crosses `b->a` while the conclusion lives at
   `b` and no feed path returns.
4. **THM-4093 (gauge and zeta tangent).**  Edge weights `d_i/d_j` are exactly
   diagonal similarity, so every closed-walk spectral invariant is unchanged.
   The first non-gauge coefficient of a tournament's adjacency zeta is its
   directed-triangle count, with sharp p-adic valuation
   `v_p(Z_T(x)-1)=3v_p(x)` iff `p` does not divide `c_3(T)`.

These are not four proofs of one grand analogy.  They identify four different
maps and, crucially, four different losses.

## 2. Inheritance pass and concept board

The bounded inheritance pass used the following current objects.

| lane | closest proved mechanism | hostile / repaired near miss | least-used sidecar |
|---|---|---|---|
| p-adic irrationality | Apéry integrality/nonvanishing/decay gate; published packet existential results | a finite gcd or positive last margin does not prove a singleton | full valuation vector, global source, Casoratian/height |
| rational edges | THM-4057 ordered coprime arcs and Stern depth gauge | coprimality graph is not a tournament; `[4]` completions change `c_3` | endpoint height and non-exact cycle decoration |
| Stern packets | THM-4059/4061/4068/4071 owner-resolved inverse parity | signed cancellation is not positive Jensen energy | full owner DFT and total variation |
| Berggren | THM-3756 outer/inner odd-square ordinals | outer square alone collides | inner ordinal `s` |
| tournament zeta | THM-1926 strong-core determinant factorization | MISTAKE-409 forbids generic “coboundary implies singular” | actual diagonal similarity and `c_3` |
| LRC(14) | THM-358/613 strict-margin interval mechanism | AP13 has rational equality points only | `(margin,max speed)` and fixed time scale |
| matching logic | cited backward core/double cover | forced tournament completion invents reachability | coordinate labels and sort-feed graph |

The live board was:

```text
global arithmetic source; valuation clock; owner-resolved Fourier rows;
rational-edge gauge; directed triangles; sort-feed reachability;
strict witness interval.                                      (1)
```

Every successful connection either retained one of these coordinates or
proved that forgetting it was fatal.

## 3. What the p-adic source actually changes

The pinned repository
[`b46a1770901551961710e155d775aae7c5ea39e7`](https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7)
contains a substantial proposed proof, not just a list of decimals.  Its
architecture is:

```text
raw Eichler rows / Dwork relation
 -> BGG algebraic connection and genuine positive global source
 -> small-prime Hasse kernel saving
 -> large-prime torsor product digits and no-backflow
 -> Bost/CDT height inequality plus continuation
 -> owner-counted modular Jensen energy
 -> positive real margin.                                  (2)
```

The exact certificate reproduces all `22` final positive margins.  That
replay is **FINITE-EXACT** for the numerical implication only.  It does not
check the arrows in (2).  The smallest margin, at `(5,5)`, is about `0.1318`,
so a missing factor or normalization can matter.

Deleting only the manuscript's Hasse saving flips exactly three displayed
witnesses—`(2,29),(5,5),(7,3)`—negative.  This isolates where that one new
ingredient is numerically indispensable; it does not prove the other
nineteen cells by older methods because all remaining arrows in (2) still
need audit.

THM-4089 then solves the next question globally.  Its objective separates as

```text
M_(p,s)(xi,Y)=G_(p,s)(Y)-(s+1)tau_(p,s)(xi).              (3)
```

The prime-window/Hasse side is globally convex with rational minimizer

```text
xi*=12(s(s+1)-1)/(12s^2+(s-1)(p+1)),                     (4)
```

and the Jensen side is globally concave with one critical point.  Exact
derivative brackets and tangent bounds make all four adjacent global maxima
negative.  The research consequence is procedural: stop spending compute on
`xi,Y` retuning for those cells.  Change one of the structural terms in (3).

## 4. Three lawful tournament carriers in the p-adic proof

### 4.1 Pivot orders: a transitive proof skeleton

If the Bost construction supplies distinct pivot orders
`n_1<...<n_r`, their scalar order tournament is canonical and transitive.
Its edge count gives exactly

```text
sum_i n_i >= 0+1+...+(r-1)=binom(r,2).                   (5)
```

This is a lawful use of a tournament: the intrinsic pairwise observable is
order.  It retains filtration rank and loses every local height, denominator,
and coefficient.  The adelic evaluation weights are the required sidecar.
THM-4088 proves that the bare order can carry any arithmetic-type labels.

### 4.2 Jensen bottom rows: coprime arcs, not a tournament

Primitive modular bottom rows `(c,d)` with `p|c` map to reduced arcs `d/c`.
The shell has `phi(c)` owners.  This preserves coprimality and exact
denominator while losing the full matrix/coset and collision height.  The
carrier has missing noncoprime pairs and reciprocal orientations; completing
it to a tournament adds information not present in the Jensen formula.

THM-4059 refines a shell by its owner and Stern sign.  Jensen takes the
unsigned quotient

```text
{(d,epsilon_c(d)):d in U_c}  ->  #U_c=phi(c).             (6)
```

A signed Jensen laboratory may reveal Fourier cancellation, but it cannot
replace a positive zero mass without an unsigned total-variation sidecar.

### 4.3 Vertex ratios: spectrally pure gauge

For any adjacency matrix `A` and nonzero vertex labels `d_i`, the proposed
rational edge weight gives

```text
W_ij=A_ij d_i/d_j,                 W=DAD^(-1).            (7)
```

Open-path products remember the endpoint ratio; every closed path cancels it.
Therefore traces, characteristic polynomial, determinant zeta, and all
closed-walk counts are unchanged.  Applying `v_p` turns (7) into the exact
additive coboundary `h_i-h_j`.

This answers one part of the user's numerator-to-denominator proposal
sharply.  The *orientation* of a reciprocal edge can encode a chosen Stern or
order selector.  The raw *ratio weight* on a fixed orientation cannot create
spectral arithmetic.  A useful decoration must have nontrivial cycle
holonomy—approximation error, continued-fraction digit, owner character,
height, or a recurrence transfer matrix.

## 5. Rational, algebraic, and transcendental patterns

There are three distinct statements, only one of which is a positive detector.

1. **Bare finite orientation is blind.**  By THM-4088, every transitive order
   tournament realizes any prescribed vertexwise arithmetic types.  A pure
   tournament is invariant under relabelling, so no unretained number type can
   be recovered from it.
2. **Finite rational-tree address is blind.**  Every nonempty open real arc
   contains a finite Stern--Brocot cell and representatives that are rational,
   quadratic irrational, exact-degree algebraic, and transcendental.  The
   cell must allow terminating completions; a strict cylinder of infinite
   continued-fraction words contains no rationals.  Termination, eventual
   periodicity, and an infinite growth/nonperiodicity certificate are three
   different tail predicates.
3. **A nonconstant rational function preserves transcendence forward.**  For
   a nontransitive tournament, `P_T(u)=det(I-uA)` is nonconstant.  If `alpha`
   is transcendental, then `1/P_T(alpha)` is transcendental; otherwise alpha
   would satisfy a polynomial over the algebraic numbers.  This is THM-4093's
   exact positive statement.  It uses a transcendental *input* and does not
   infer that input from the tournament.

The transitive boundary is sharp: `P_T=1`, so the same zeta collapses every
input to `1`.  For algebraic inputs, rational and algebraic-irrational cases
can also collide under a rational function.  A finite tournament therefore
does not classify the three types even though one nonconstant associated
function preserves transcendence in the forward direction.

## 6. The exact `p=13` shell and the character correction

Normalize Beukers' six reflection-orbit carriers by

```text
u_a=omega(a)^(-2) H_13(3,a,13),       u_(13-a)=u_a.      (8)
```

Then `U_13/{+-1}~=C_6`.  Taking the unique odd representative
`d in {1,3,5,7,9,11}`, setting `r=7,s=(d+1)/2`, and applying THM-3756 gives
the complete outer Berggren shell `B+C=13^2`:

| `d` | primitive triple `(A,B,C)` | word | Stern/Legendre sign |
|---:|---|:---|---:|
| `1` | `(13,84,85)` | `LLLLL` | `+1` |
| `3` | `(39,80,89)` | `ML` | `+1` |
| `5` | `(65,72,97)` | `RM` | `-1` |
| `7` | `(91,60,109)` | `LLR` | `-1` |
| `9` | `(117,44,125)` | `LRR` | `+1` |
| `11` | `(143,24,145)` | `RRRRR` | `-1` |

In generator-two order the odd representatives are

```text
1,11,9,5,3,7
```

and the signs alternate.  This is exactly the order-two (`k=3`) Fourier row
of `C_6`, a Paley colour/graph coordinate rather than a tournament.

The referee caught a crucial branch error.  In the normalized `u` coordinates,

```text
sum_a u_a = zeta_13^Delta(3)=L_13(3,omega^(-2))          (9)
```

is the trivial Fourier row.  Beukers' principal branch is instead

```text
L_13(3,1)=sum_a omega(a)^2 u_a,                           (10)
```

an order-six row.  The Stern/Legendre row is orthogonal to the trivial row and
isolates neither (9) nor (10).  The two-coordinate abstract hostiles make the
loss unavoidable:

```text
(u_0,u_1)=(sqrt(2),-sqrt(2)) : trivial sum rational,
                                      order-two row irrational;
(u_0,u_1)=(sqrt(2), sqrt(2)) : trivial sum irrational,
                                      order-two row rational.              (11)
```

Thus the shell is a lossless index bridge but not an irrationality selector.
The next lawful object is the full even-character DFT paired with simultaneous
Padé/Casoratian nonvanishing and height.

## 7. Matching logic as a test of global transport

The matching-logic paper's one-sort theorem localizes a global theory by
closing all hypotheses under backward coordinate paths.  Its double cover
works because the basic language cannot distinguish two copies of the outside
phantom.  Fixpoints, free set variables, and nominals break different parts of
that mechanism.

THM-4090 shows the many-sort obstruction at the smallest possible sort count:

```text
b --f--> a,
Gamma at sort a forces |M_b|=1,
phi at sort b says |M_b|=1,
but no proof rule transports Gamma against the feed arrow.                (12)
```

This resembles the p-adic source only at a typed meta-level.  A true global
semantic implication may fail to travel through a calculus whose local feed
graph lacks the return path; a true final real margin may fail to travel
backward into irrationality when the geometric/lattice interfaces are
unverified.  There is no direct theorem map between the subjects.

The native matching carrier is a labelled directed dependency graph with
loops, missing pairs, and possibly both directions—not a tournament.  On the
special case where its unlabelled direct dependency graph is a finite
tournament, backward-closed cores are suffixes of the transitive SCC
condensation.  This is a useful localization skeleton, but it forgets
coordinate labels, tuple witnesses, symbol interpretations, and formula
semantics.

## 8. Connections back to other repo frontiers

### LRC(14)

THM-4088 is a post-processing theorem, not an existence theorem.  When a row
has strict margin `eta`, every arithmetic type occurs in its witness interval.
At a tight boundary, THM-358 shows the opposite.  Hence rational clock
certificates are sufficient for existence when available, but their
rationality is not an invariant of the whole strict witness set.  THM-4092's
new parity-weighted comb compiler supplies positive intervals on additional
families; THM-4088 immediately enriches those intervals without closing the
remaining LRC(14) rows.

### Planar Jacobian

The p-adic draft and the JC response program both use Hermite--Padé-looking
matrices, but MISTAKE-222's firewall remains decisive.  JC needs polynomial
entry, divisor contact, a constant-Jacobian implication, and algebraization;
the p-adic argument needs denominator lattices, nonzero linear forms, and
adelic/archimedean decay.  Shared syntax supplies no transfer.  The useful
method import is narrower: draw the feed graph of proof obligations and test
whether each local source can reach the global target before optimizing its
numerical margin.

### Mahler `3/2`

THM-4082's long common binary prefixes are another finite-address quotient.
They do not decide ordinary termination, just as a finite Stern cell does not
decide arithmetic type.  The shared research move is to retain the infinite
tail predicate; there is no claimed numerical or dynamical conjugacy.

### Rule 30

Finite cylinders and Haar-random spatial laws do not automatically describe
the distinguished temporal orbit.  This is the same quantifier boundary as
finite p-adic digit balls and finite continued-fraction cells.  A transfer
needs a temporal current, not an order tournament on observed prefixes.

### Hamiltonian-path tournament spectrum

The concurrent THM-4094 turns the matching-logic warning into an exact
tournament statement without conflating the subjects. For vertex deletion,
the **full** two-sort path incidence retains the exact increment

```text
H(T)-H(T-v)=sum_P(a_v(P)-1)+O_v,
```

whereas selecting one extension per old path discards both excess insertions
and orphan paths. This is the same controlled-forgetting grammar as sort flow,
but here the source, target, loss map, and repaired sidecar are literal
tournament objects. Combined with SCC multiplicativity and the proved gaps
`7,21`, it reduces the open global `H`-spectrum conjecture to two strong prime
lanes, with first unforced construction targets `613` and `623`. THM-4088
prevents a cosmetic shortcut: assigning rational, algebraic, or transcendental
vertex labels cannot create those strong atoms unless the arithmetic
decoration controls the insertion deficit or cycle structure intrinsically.
The common minimal activation is `C_3`: a scalar-order tournament is
transitive, with `Z_T=1`, `H(T)=1`, and zero deletion deficit regardless of
the arithmetic types decorating its vertices; `C_3` instead has `c_3=1` and
`H(C_3)-H(K_2)=2`. Thus `H>1` certifies nonconstant tournament zeta, but
`H`, `c_3`, and `nu_p(H)` observe adjacency/SCC structure rather than the
arithmetic type of external labels.

### Exact field transport versus optimization

The incoming THM-4095 locates a positive result hidden behind THM-4088's
no-go. For fixed finite `S` and rational `beta`, the piecewise-affine observer
does not merely preserve coarse rational/algebraic/transcendental type:

```text
Q(F_S(t)-beta)=Q(t).
```

The active branch `(v,k,epsilon)` reconstructs `t`. The loss occurs only when
that branch and value are replaced by scalar order, or when the quantifier is
changed by optimizing over speed sets. Primitive two-speed optimized margins
then form a discrete spectrum with gap `(0,1/15)`. This cleanly separates
three levels: exact field transport for a fixed observer, topological density
inside a strict witness component, and a gapped spectrum after optimization.

THM-4096 gives the analogous boundary on the corrected regularization lane.
The normalized Cover14 carrier is a nonnegative rational affine combination
of the correctly typed twisted values `L_p(-1,omega_p^2)` only at `m=1`; its
required moment is negative for every `m>=2`. Because the construction places
cross-prime rational shadows in a common real affine line, it deliberately
forgets p-adic topology. The next lawful experiment must stay at one prime and
retain character plus valuation; otherwise it can only rediscover the affine
hostile.

## 9. Generated next tasks

### Anchor -- specialist audit of the `22`-cell draft

1. Reconstruct the `(5,5)` and `(7,3)` BGG rows and genuine multiplier in a
   CAS, reduce at several auxiliary primes, and verify Hasse degree/kernel on
   the actual global source.  Use `(2,3)` as a positive control and a doubled
   unipotent action as a hostile.
2. Reimplement the torsor product-digit no-backflow proposition over finite
   flat coefficient rings, including zero-divisor special fibres.  Delete one
   carry or one full-product digit to force the predicted failure.
3. Reconcile the modular Jensen formula with the prior `Gamma_0(2)` lune:
   owners, orbifold multiplicities, active cosets, and normalization must agree
   before using its energy globally.

### Niche -- character-resolved `p=13` isolation

4. Compute the complete six-row even-character DFT of every available Hurwitz
   linear-form template.  Record rank, Casoratian, integral Smith data,
   archimedean height, and p-adic decay for the desired row; the order-two row
   alone is now closed as insufficient.
5. Build an owner-resolved signed-Jensen diagnostic on small shells
   `c=3,5,9,15`, but carry unsigned total variation alongside it.  Seek a
   domination inequality, not a visual cancellation pattern.

### Wildcard -- non-exact tournament decorations

6. Classify minimal edge decorations on Stern tournaments modulo diagonal
   gauge.  The target coordinate is fundamental-cycle holonomy; test whether
   continued-fraction error, digit, or owner characters create a nontrivial
   p-adic zeta family with controlled denominator growth.
7. Correlate THM-4059 lower-star imbalance with directed triangles and strong
   cores.  `T_5,T_9` are mandatory hostiles: both have balanced apex but
   nontrivial, different zeta.

### Cross-lane proof engineering

8. Mechanize THM-4090's two-sort witness and generalize the rule induction to
   a precise feed-cut lemma.  Keep this scoped to the displayed basic Hilbert
   system; a cross-sort/global rule is an explicit escape.
9. Build an LCM minimal-denominator atlas for the external Eichler rows,
   retaining every prime valuation and Dwork residual.  Feed the vector—not a
   gcd scalar—into the repo's Apéry searches.
10. For each remaining LRC/JC/Rule-30 certificate, annotate the proof-obligation
    feed graph and identify one local statistic that cannot reach the global
    conclusion.  Replace that statistic before enlarging its finite census.

## 10. Closing status

The session made concrete progress without promoting a headline it did not
prove.  The p-adic draft is now a sharply pinned, auditable research claim;
its next four cells have a complete formula-level obstruction.  Tournament
structure has one exact arithmetic-positive role—the nonconstant zeta
preserves transcendence of a supplied input—and two exact firewalls: scalar
order forgets label type, and vertex ratios are diagonal gauge.  The most
promising frontier is therefore not another bare tournament encoding, but a
character- and height-resolved family with nontrivial cycle holonomy and a
global nonvanishing/denominator packet.
