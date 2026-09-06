# A finite hereditary gcd sieve through seven-speed cores

**Status: PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN RUNNERS +
FINITE-EXACT + INDEPENDENTLY AUDITED.** The analytic reduction is
unbounded in speed height. The finite computation enumerates complete
arithmetic signatures after a proved cutoff; it does not enumerate or
construct counterexamples to LRC(14).

Let `V` be a primitive set of thirteen distinct positive integer speeds,
and put `M(V)=max_x min_(v in V)||vx||`. If `M(V)<1/14`, then every subset
`P` of the indicated size satisfies:

| `|P|` | largest possible `gcd(P)` after this sieve | number of retained gcd values | number of sorted signatures, including gcd one |
|---:|---:|---:|---:|
| 12 | 1 | 1 | 1 |
| 11 | 2 | 2 | 2 |
| 10 | 4 | 4 | 5 |
| 9 | 9 | 7 | 19 |
| 8 | 32 | 16 | 110 |
| 7 | 96 | 43 | 1,217 |

In particular, a primitive row containing seven speeds whose gcd exceeds
96, eight whose gcd exceeds 32, or nine whose gcd exceeds nine is safe.
For an arbitrary row, replace these inequalities by
`gcd(P)>B_|P| gcd(V)`; common dilation preserves `M`.
The bounds are sharp for the declared arithmetic relaxation, not a sharp
classification of actual failures. Explicit safe rows attain its maximal
signatures. LRC(14) remains **OPEN**.

## 1. Inheritance, the retained object, and the new comparison

The preceding [effective-clock proof](lrc14_effective_clock_empty_core_sep06.md)
closed the first three rows. Incoming
[THM-4447 — composite-clock capacity and small-clock reduction](../../01-canon/theorems/THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction.md)
independently recovers that boundary and retains the exact count chambers
and divisor absorption. The common-fibre mechanism is older:
[THM-2064](../../01-canon/theorems/THM-2064-multitail-sheet-capacity-and-dyadic-seam.md)
credits concurrent THM-2060. Neither orbit capacity nor absorption is
claimed as a new operation here.

For overlap credit, the predecessors are
[THM-1217 — mixed-period beat masks and Hunter tree](../../01-canon/theorems/THM-1217-mixed-period-beat-mask-hunter-tree.md)
and [THM-860 — primitive Hamming-six ramification](../../01-canon/theorems/THM-860-primitive-hamming-six-finite-ramification-reduction.md).
The latter has all relative cuts and a fibre-Hunter strengthening on its
AP-centred six-colour carrier at threshold `1/13`. Its cap is not a theorem
about arbitrary seven-speed cores at threshold `1/14`.

Here the source is an actual subset `P` and its complementary tails. The
map keeps `c=gcd(P)` and `g_i=gcd(c,w_i)`, then absorbs divisible tails at
every divisor of c. It preserves the validity of a sufficient free-sheet
certificate. It loses tail units, phase addresses, and all divided-body
speed ratios. Balanced residue counts restore a phase-uniform part of the
pairwise intersection data before the scalar union bound is taken.
The new finite hierarchy combines that bound with an inductive pivot on
the number of tails.

The live board is: subset gcds; effective orbit orders; divisor absorption;
balanced CRT intersections; tree overlap credit; actual safe phases.
The cheap hostile is the [nine-body boundary](gcd_nine_audit_empty_core_next_sep06.md):
clock nine with tails `(1,5,6,7)` covers all nine lifts of `y=1/2`, while
an actual body is safe there and the full row is safe at another time.
The missing sidecar remains the whole body-to-tail phase map.

Targeted statement, constant, and synonym searches recovered the cited
ramification and capacity precedents, but no existing `9,32,96` hierarchy
for arbitrary subsets of a primitive thirteen-speed strict counterexample.
No external-priority claim is made. After consulting `CORE-PAPERS.md`, the
[current primary preprint](https://arxiv.org/html/2604.23906v2), Theorem1.3
and Definition2.1, were checked: the former supplies lower-runner LRC;
the latter's properness sieve retains leave-one-out gcds, and is not the
smaller-subset classification asserted here.

## 2. Padded capacities and an exact CRT pair bound

For a clock `a>=2`, a body-safe phase has a complete labelled fibre
`x_j=(y+j)/a`, `j mod a`. A tail with `h=gcd(a,w)` has effective order
`q=a/h`. Let

```text
k(q)=ceil(q/7),               B_a(w)=h k(q).
```

The actual bad orbit residues are a cyclic consecutive block of at most
`k(q)` residues after multiplication by a unit. This uses the open danger
arc of length `1/7`; when `7|q`, the capacity is exactly `q/7`.
Pad that block to exactly `k(q)` residues in the same unit order. The
padded labelled bad set contains the original set and has cardinality
`B_a(w)`. Padding is only an upper-bound device; it need not be realized
by the same physical phase or by every tail simultaneously.

Consider two padded sets of orders `q,r` dividing a. Write

```text
b=gcd(q,r),       k(q)=A b+u,       k(r)=B b+v,
0<=u,v<b.
```

Their intersection is at least

```text
L_a(q,r) = a/lcm(q,r) *
           [b A B + A v + B u + max(0,u+v-b)].            (1)
```

Proof. Reducing a consecutive unit-step block modulo b gives counts
`A` or `A+1`, with exactly u high classes; the other block similarly
has B or `B+1`, with v high classes. The high-class intersection is at
least `max(0,u+v-b)`. For each compatible residue pair, CRT supplies
exactly `a/lcm(q,r)` labelled preimages. Expanding the dot product of the
two residue-count vectors proves (1). Choosing unit step one and suitably
translated high-class intervals attains this pair minimum in the padded
relaxation. Simultaneous attainment for many pairs is not assumed.

For coprime q,r, (1) is simply
`a k(q)k(r)/(qr)`. At `a=288`, orders16 and9 force twelve common labels.
This is information the separate tail capacities discard.

## 3. A tree bound and every divisor of the body clock

Given tail gcds `h_i` at clock a, put `q_i=a/h_i`. Define

```text
U(a;h)=sum_i h_i k(q_i)
       - max_(spanning trees T) sum_(ij in T) L_a(q_i,q_j). (2)
```

The empty or singleton edge sum is zero. The actual union of all tail
danger sets has size at most U. Indeed it is contained in the union of
the padded sets. At any label belonging to `m>=1` padded sets, the induced
edges of a tree number at most `m-1`; summing gives the usual Hunter
upper bound with actual padded intersections. Replacing each intersection
by its lower bound (1) preserves the upper bound. Maximizing the tree
credit is legitimate even when its individual pair minima cannot all be
attained together. An upper bound `U<a` supplies a free label.

Now start with `c=gcd(P)`, complementary tails `w_1,...,w_d`, and
`g_i=gcd(c,w_i)`. For every `a|c`, `a>=2`, absorb all a-divisible tails:

```text
R_a=(c/a)(P/c) union {w_i/a:a|w_i},
D_a={w_i:a does not divide w_i}.
```

This is exactly the same physical row `V=aR_a union D_a`. Full primitivity
makes `D_a` nonempty, so `|R_a|<=12`. Cited lower-runner LRC provides a safe
body phase, with clearance strictly greater than `1/14`.
The surviving tail gcds are `h_i=gcd(a,g_i)`, for `a` not dividing `g_i`.
Consequently every strict counterexample must satisfy all the finite cuts

```text
U(a; sorted(gcd(a,g_i): a does not divide g_i)) >= a
                                  for every a|c, a>=2.  (3)
```

All child-subset cuts are already included: adjoining some tails changes
the exact body gcd to a divisor of c. Every divisor of that new gcd was
already tested in (3), with all its divisible tails absorbed. This is a
flat divisor test, not a claim that forgetting those tails preserves the
original chosen body phase. Each repackaged body receives its own cited
safe phase.

## 4. Why the complete arithmetic universe is finite

Let `S_d` denote the gcd values surviving the declared sieve for a body of
`13-d` speeds. Set `S_0={1}` by full primitivity. Suppose `d<=6` and the
preceding set `S_(d-1)` is known. Every complementary tail gives

```text
g_i=gcd(P union {w_i}) in S_(d-1).
```

At the original clock c, ignoring overlap and absorption is still a valid
necessary relaxation:

```text
sum_i beta(c/g_i)>=1,          beta(q)=ceil(q/7)/q.         (4)
```

Some q therefore has `beta(q)>=1/d`. Since
`ceil(q/7)<=(q+6)/7`, this implies

```text
q <= floor(6d/(7-d)).                                     (5)
```

Thus c belongs to the explicitly generated finite set

```text
{qg:g in S_(d-1), 1<=q<=floor(6d/(7-d)), d ceil(q/7)>=q}.  (6)
```

For each c, enumerate every sorted d-tuple from
`{g in S_(d-1):g|c}`, require `gcd(g_1,...,g_d)=1`, and apply every
cut (3). This induction proves completeness before computation. Every
actual counterexample signature must be present; the converse is not
claimed. The resulting gcd sets are

```text
S_1 = {1},
S_2 = {1,2},
S_3 = {1,2,3,4},
S_4 = {1,2,3,4,6,8,9},
S_5 = {1,2,3,4,5,6,8,9,10,12,15,16,18,24,30,32},
S_6 = {1,2,3,4,5,6,8,9,10,11,12,15,16,17,18,20,22,23,24,25,
       27,29,30,32,33,34,36,40,44,45,46,48,50,51,54,58,60,64,
       66,72,88,90,96}.
```

The primary candidate/proposal counts are
`(1,1),(2,2),(5,14),(11,161),(31,2721),(138,397280)` for d1 through6.
These are finite arithmetic words, with no speed-height or phase cutoff.
The [full retained signatures](lrc14_recursive_gcd_empty_core_next_sep06.json)
are part of the consequence certificate, not merely an intermediate count.

The sharp nine-body cap also has an elementary short proof in the
[independent nine-body note](gcd_nine_audit_empty_core_next_sep06.md).
Retaining all exact ten-body/absorption conditions removes two of the
20 nontrivial profiles left by just its scalar caps; 18 plus c1 agree here.

## 5. What overlap changes, and the precise stopping boundary

Without tree credit the same recursive scalar/divisor sieve gives maxima
`1,2,4,9,36,288`. The new credit first changes the eight-speed result.
For example, at `c=36` the gcd word `(1,4,4,6,9)` has capacities
`(6,8,8,6,9)` summing to37. Orders9 and4 force an intersection of two,
so its union is at most35, below36. At `c=288`, the word
`(1,1,4,4,18,32)` has total capacity290, while orders16 and9 force an
intersection of twelve. Thus that old scalar survivor cannot cover.
These examples explain the first failed implication of a scalar-only
sharpness claim; they do not assert unsafe rows before the improvement.

The largest retained clock32 word is `(1,1,2,4,4)`. At clock96 the only
largest words are `(1,1,4,4,6,12)` and `(1,3,4,4,6,12)`.
All these signatures, as well as the maximal clock9 word, are realized
by explicitly safe primitive thirteen-speed rows. Take body
`c{1,...,13-d}` and tails

```text
w_i=g_i[1+c(1+14(i+1))],             i=0,...,d-1.          (7)
```

The source checks distinctness, each gcd, and safety at `x=1/(14c)`.
For these named words, every g is at most12 and every tail phase is
`(g_i+1/q_i)/14`, between `1/14` and `13/14`; body phases are also safe.
These realizations prove arithmetic consistency and block interpreting
retained words as counterexamples.

At **seven tails**, this finite-clock mechanism stops. For every integer
`c>=2`, the word `g=(1,...,1)` of length seven survives (3). At each
`a|c`, `a>=2`, all seven residual orders equal a, total capacity is
`7ceil(a/7)>=a`, and every pair lower bound is zero because
`2ceil(a/7)<=a`. The last inequality is immediate at a2 and follows from
`ceil(a/7)<=(a+6)/7<=a/2` for a>=3. Thus every integer clock survives
this arithmetic relaxation. Formula (7) with six body speeds gives an
explicitly safe realization at every such clock. This is a proved boundary
of the retained information, not a suggested LRC counterexample. The
independent referee also checked this all-clock strengthening; the eleven
prime controls are a finite subset of it.
A further descent requires phase/unit constraints or another operation.

## 6. Reproduction and independent audit

The [primary source](../../04-computation/lrc14_recursive_gcd_empty_core_next_sep06.py)
uses only standard-library integer arithmetic, exact divisibility, and a
Kruskal maximum spanning tree. It imports no mathematical producer.
Its 3,760 literal pair-mask cases cover every order pair from1 through28,
both orientations and every relative common-residue shift, and verify
both the lower bound and attainment in the padded pair relaxation.
The source also checks scalar hostiles, absorption, all named maximum-clock
safe realizations and eleven prime controls for the seven-tail boundary.

```bash
python3 -B 04-computation/lrc14_recursive_gcd_empty_core_next_sep06.py --write-profiles 05-knowledge/results/lrc14_recursive_gcd_empty_core_next_sep06.json
python3 -B -O 04-computation/lrc14_recursive_gcd_empty_core_next_sep06.py
```

Normal and optimized stdout agree with the
[frozen output](lrc14_recursive_gcd_empty_core_next_sep06.out): 4,618 explicit
gates, plus the complete finite compiler. The profile JSON SHA256 is
`7af5425db35516037e03bb2d114c34bd1f2c8971cdee91866dc77b754b13f557`.
The independent audit uses every clock up to the coarse analytic bound,
a Prim tree implementation, and separate literal unit/block intersections.
The [independent audit](gcd_pair_hunter_audit_empty_core_next_sep06.md)
passes the full written proof and every retained signature, not only the
counts or maxima. It has 1,579,350 explicit gates and uses all clocks through
`1,2,8,32,135,1152`; all six profile arrays agree exactly. The root also
independently read the nine-body proof and pair/Hunter argument: **PASS**.

Frozen primary raw-byte hashes:

```text
source 30d9e767cc502ec3bb818edf07d8ce8adf1546a4224abe731aaa0ee37d287063
output 86ad4e7f2a61c4efaa5e538646132ae9305c4fab849e295fa20e9809fae958c5
```
