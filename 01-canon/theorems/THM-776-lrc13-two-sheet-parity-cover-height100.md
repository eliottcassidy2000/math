---
id: THM-776
title: Finite-exact exclusion of the two-sheet parity-cover packet through height 100
status: PROVED (finite-exact rational atomization and exhaustive CNF; not a uniform theorem)
source: codex-2026-07-14-S8
depends_on:
  - THM-769   # exact folded two-sheet criterion
related:
  - THM-668   # pair-crossing exact maximin principle used in the direct audit
  - THM-772   # primitive divisor-complete quotient transfer
  - THM-774   # equivalent folded-diamond odd-pair obstruction
  - THM-775   # proved dyadic deletion descent toward a finite base
  - HYP-6820  # n=12 sporadic-branch audit
verification: 04-computation/lrc13_s2_parity_cover_h100_codex_S8.py
  (+ 05-knowledge/results/lrc13_s2_parity_cover_h100_codex_S8.out)
---

# THM-776 — The two-sheet parity-cover packet is empty through height 100

## Statement

Let `A` be a set of twelve distinct positive integer speeds with

```text
max(A)<=100
```

and suppose that exactly ten members of `A` are even and two are odd.  Then

```text
M(A)=max_t min_(v in A) ||vt|| > 1/13.                    (1)
```

Equivalently, no such set is tight at the twelve-speed threshold.  In
particular, the primitive two-exception deep packet of THM-769 is empty through
height 100.

The conclusion is stronger than the primitive deep-branch application: the
finite census assumes neither primitivity nor the divisor-completeness forced
by THM-772.  It is, however, only a bounded theorem.  It does not exclude a
two-sheet packet of height greater than 100.

## 1. The exact folded obstruction is a hitting set

Write uniquely

```text
A=2U union {x,y},       |U|=10,       x<y odd,
U subset {1,...,50},    x,y<=99.
```

Put

```text
G_U={tau in R/Z : ||u tau||>1/13 for every u in U}.        (2)
```

For an odd `w`, whenever `||w tau||<=2/13`, the nearest integer
`N_w(tau)` is unique.  Give `w` the sheet colour

```text
epsilon_w(tau)=N_w(tau) mod 2.
```

Define the closed opposite-colour eligibility locus

```text
H_(x,y)={tau : ||x tau||,||y tau||<=2/13
                and epsilon_x(tau)!=epsilon_y(tau)}       (3)
```

and put `B_(x,y)=(R/Z)\H_(x,y)`.

In THM-774's folded-diamond coordinates

```text
a=(x+y)/2,       b=(y-x)/2,
```

the same closed set is

```text
H_(x,y)={tau:||a tau||+||b tau||>=11/13}.                 (3a)
```

Thus the computation below is simultaneously an exact parity-colour census
and an exact finite census of the folded-diamond obstruction.

The two lifts of `tau` are

```text
t_0=tau/2,       t_1=(tau+1)/2.
```

At both lifts an even speed `2u` has clearance `||u tau||`.  An eligible odd
speed kills exactly the lift bearing its colour.  Hence the two odd speeds
kill both core-safe lifts exactly on `H_(x,y)`.  Consequently

```text
M(A)<=1/13
iff G_U subset H_(x,y)
iff B_(x,y) subset union_(u in U){tau:||u tau||<=1/13}.    (4)
```

The last condition is a finite hitting-set problem once its rational boundary
cells are exposed.

## 2. Lossless rational atomization

Use all points in `[0,1]` of the following two forms:

```text
||u tau||=1/13,       1<=u<=50,
||w tau||=2/13,       1<=w<=99 odd.                        (5)
```

After duplicate removal these give exactly

```text
6,876 endpoints and 6,875 open atoms.                     (6)
```

All predicates in (2)--(3) are constant on each open atom.  This midpoint
test is exact despite the endpoint asymmetry: `G_U` and `B_(x,y)` are open, so
any failure of (4) is open and contains a whole atom.  An isolated endpoint
can never be the only failure.

For a bad atom `I subset B_(x,y)`, define its coverer set

```text
C(I)={u in {1,...,50}: ||u tau||<=1/13 on I}.              (7)
```

Then (4) holds exactly when `U` intersects every `C(I)`.  If
`C(I) subset C(J)`, the obligation from `J` is redundant.  Removing all such
supersets leaves an inclusion-minimal monotone CNF without changing its
transversals.

There are

```text
C(50,2)=1,225
```

odd pairs.  Their compressed obligation hypergraphs have between 332 and 362
minimal clauses.  Across the full atlas the clause cardinalities are

| coverer-set size | number of minimal obligations |
|---:|---:|
| 4 | 5,995 |
| 5 | 6,205 |
| 6 | 135,811 |
| 7 | 257,611 |
| 8 | 11,259 |
| 9 | 9,459 |
| 10 | 1,178 |

## 3. Exact transversal result

For every odd pair, the compressed hypergraph is UNSAT with a cardinality
constraint of at most eleven selected quotient speeds.  The exact census is

```text
transversal-number histogram:  {12: 1,225}.               (8)
```

The lower bound in (8) was obtained with CaDiCaL 1.9.5 on the at-most-eleven
CNF.  A separately encoded at-most-ten CNF was also UNSAT for all 1,225 pairs
under Glucose 4.  The upper bound is elementary and checked atom by
atom: `{1,...,12}` covers the whole circle by closed `1/13`-danger teeth.

The canonical instance/atlas fingerprints are

```text
rational endpoint grid:
  6da3195bf93d5446d6eb66e6bcbf29d65e3069e0d12d3a13551db7addcf6d7ac

compressed constraint atlas:
  a2d415506978ad031c9f3b77ae7d751598a3658c378e832c25828e44ae28bd96
```

Since a packet core has only ten quotient speeds, (8) makes (4) impossible.
Thus some bad atom contains a `tau` safe for all of `U`; one of its two lifts
is then strictly `1/13`-safe for every speed in `A`.  This proves (1).

### Independent small-height audit

As a separate check on the folded encoding, all

```text
C(12,10) C(12,2)=4,356
```

ten-even/two-odd sets of height at most 24 were evaluated by the exact
pair-crossing/self-cusp formula for `M(A)`.  There are zero tight rows.  The
closest row is

```text
A={2,6,8,10,12,14,15,16,17,18,20,22},
M(A)=2/23 at t=6/23.                                      (9)
```

Thus the direct maximin engine, folded atom engine, and two Boolean solvers
agree on their common range.

An independent referee regenerated the full height-100 rational grid and the
entire folded-diamond obligation atlas, reproducing both fingerprints and all
clause histograms.  A separate exhaustive search through height 30 obtained
transversal number 12 for all 105 odd pairs, and a different totalizer
encoding proved the at-most-eleven CNF UNSAT for the representative pairs
`(1,3)`, `(49,51)`, and `(97,99)`.  No DRAT/LRAT traces are stored, so the
fingerprints are reproducibility anchors rather than independently checkable
proof certificates.

## 4. Interaction with the divisor-transfer theorem

THM-772 says that a hypothetical primitive packet in this branch has a
primitive core `U` containing a multiple of every modulus `2,...,12`, with
`x,y<=11 max(U)`.  The height-100 theorem does not need those restrictions.
An exact subset DP quantifies the saving:

```text
10-subsets U of {1,...,50} divisor-complete for 2,...,12:
    885,427,750

primitive such cores:
    884,640,190

primitive filtered core/odd-pair packets before parity:
    1,083,684,232,750.                                    (10)
```

The first divisor-complete ten-core already has maximum 12.  Directly
enumerating (10) is the wrong coordinate system.  The atom-incidence duality
reverses the quantifiers: fix one of only 1,225 odd pairs, then rule out every
core simultaneously by its transversal number.

THM-774 supplies a complementary exact slice.  Its widest-component argument
rules out every `U subset {1,...,19}` while allowing the two odd exceptions to
be arbitrarily large.  The present theorem allows `U subset {1,...,50}` but
restricts the odd exceptions to at most 99.  Neither finite theorem contains
the other; together they close a low-core/unbounded-odd rectangle and a larger
bounded-height square in the same folded incidence space.

## 5. Tournament and information-preservation audit

On every component of `H_(x,y)`, use the two odd runners as vertices and
orient

```text
x -> y  iff x owns sheet 0,
y -> x  iff y owns sheet 0.                               (11)
```

The switch/gauge `tau -> 1-tau` flips both nearest-integer parities and hence
reverses (11).  Each observed tournament is the transitive tournament on two
vertices:

```text
score histogram (0,1), directed cycles 0,
SCC sizes (1,1), Hamiltonian-path count 1.
```

The fixed order `x<y` is the tie Hamiltonian path on same-sheet cells.  Across
the complete height-100 atlas there are exactly 9,358 good components in each
orientation, and 27 odd pairs have no good ownership component at all.  The
average circle measures per odd pair are

| state | average measure |
|---|---:|
| `x -> y` | 0.0232141528654 |
| `y -> x` | 0.0232141528654 |
| same-sheet tie | 0.0481158158318 |
| exactly one eligible | 0.426296372259 |
| neither eligible | 0.479159506178 |

This challenges the default assumption that runners, arcs, or pairwise
orientations are the proof vertices.  The runner tournament records ownership
where ownership succeeds, but it destroys the simultaneous-cover condition
and cannot see the transversal number in (8).  Sheet vertices alone are also
too coarse.  The exact carrier has **bad atoms/proof obligations as vertices**
and quotient speeds as incident killers.  Passing from geometry to the
inclusion-minimal incidence hypergraph destroys atom positions and widths, but
preserves the LRC predicate (4) exactly.

## 6. Honest frontier

The bounded result identifies a robust two-tooth deficit: through quotient
height 50, every odd-pair failure locus needs twelve core danger families,
whereas the packet supplies ten.  Nothing here proves that the transversal
number stays above ten when larger quotient speeds are admitted.

The uniform two-sheet problem should now be stated as follows:

> Can a primitive divisor-complete ten-set `U`, together with odd
> `x,y<=11 max(U)`, hit every bad-atom obligation of `B_(x,y)`?

THM-772 supplies the arithmetic restrictions; THM-774 and THM-776 supply two
complementary finite-exact base regions and the predicate-preserving object.
THM-775 proves that every failure of hereditary primitivity is a dyadic seam,
but it does not yet force the terminal hereditarily primitive quotient into
either certified region.  A uniform proof still needs that quantitative
landing step, or an unbounded lower bound on the corresponding
obligation-hypergraph transversal number.
