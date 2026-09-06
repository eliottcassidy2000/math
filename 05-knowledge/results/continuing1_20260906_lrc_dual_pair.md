# A dual pair-gcd inequality closes actual decoder classes without a unit

**Status: PROVED relative to the named proved and CITED suppliers;
FINITE-EXACT controls; INDEPENDENTLY AUDITED.** The [independent root audit](continuing1_20260906_lrc_audit.md)
accepts the full proof and reconstructs its high-minimum actual entry. This note
does not claim a universal decoder closure or LRC(14). No theorem ID or
external priority is claimed. Repository inheritance was checked at clean
`677bde8eb`; all new files were written outside the repository.

Let Q=91^6=567,869,252,041. Consider a primitive, positive, distinct
thirteen-speed physical row in the box `sum speeds<=Q^2`, with exactly two
**actual** THM-3818 decoder components and equality `W_(Q,3)=V_dec`. Write it
as

    tV union gU,  gcd(V)=gcd(U)=gcd(t,g)=1,
    |V|=a, |U|=b, a+b=13, 2<=a<=6<=b,
    u=min U, L=max U.

Neither primitive shape is required to contain one. Choose two distinct
labels A<B in V, and put

    D=gcd(A,B), p=A/D, q=B/D.

The main sufficient condition is

    7(b+1) D L <= aQ,
    7(b+1) B   >= a G_b,                         (1)

where the inherited hypothetical-failure ceilings are

    (G_7,G_8,G_9,G_10,G_11)=(90,30,9,4,2).

**Theorem.** Every stated actual entry satisfying (1) has a common time with
clearance at least 1/14. For the branch `g<=G_b` the argument supplies a time
with clearance strictly greater than 1/14. The `g>G_b` branch uses the
inherited weak-safety supplier and does not export strictness.

In particular, in the balanced 6+7 case it is enough to find a V-pair with

    B>=10,    gcd(A,B) max U <= floor(3Q/28)=60,843,134,147.   (2)

This applies to high-minimum six-components and to seven-components without
a unit. The product is a selected pair gcd times the **global primitive
maximum of the other component**. Replacing either factor by a collective
gcd or a selected local maximum is not justified.

## 1. Inheritance, concept board, and correction lineage

The closest mechanisms are the actual crossing prohibition and internal
pair-height bound in **THM-3818**, Sections6.4--6.5,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`;
the finite-box supplier **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`;
and the audited exact signed-box theorem
`05-knowledge/results/overnight12_20260906_lrc_gcd_semigroup.md`.
The full component/rank normalization is independently audited in
`05-knowledge/results/overnight14_20260906_lrc_general_decoder.md`.

The completed, audited note
`05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.md` supplies
G_b for a hypothetical primitive failure. Its ceilings are not conditions
on every safe row. The lower-dimensional phase input is **CITED**:
Sungkawichai--Trakulthongchai, *Eleven, twelve, and thirteen lonely runners*,
[arXiv:2604.23906v2, Theorem1.3, printed page2](https://arxiv.org/pdf/2604.23906v2).
Its primary abstract and theorem were reopened for this session. It supplies
clearance1/(k+1) for k<=12 nonzero speeds; its external computer-assisted
proof was not repeated here.

The previous
`05-knowledge/results/overnight15_20260906_lrc_larger_unit.md` uses a pair
(1,L) on U and an outside V-label to force t large. Here the pair is on V,
and an outside U-label forces t large by making the minimal outside
coefficient exceed Q. This reverses which quantity leaves the coefficient
box. It removes the unit assumption for the stated new subclasses.

The incoming `05-knowledge/results/creative_20260906_lrc_bridge.md` was
read in full. It computes attainable endpoint-walk cancellation factors
exactly, but does not force a successful endpoint or a new safe class.
The present theorem consumes a different selected pair, across the two
component roles, and reaches a safe time. No divisor-walk compiler is
reproved. The inherited connected-gcd descent lemma, rather than a new
walk classification, supplies the automatic table in Section4.

The concept board is: bounded three-label relations; minimal outside
coefficients; component-scale caps under failure; pair versus collective
gcd; global safe arcs; and primitive star factors. The source-to-target map
is a bounded relation obstruction to a lower bound on t, then to t lifts of
one V-safe phase. This preserves all V-clearances. Coprimality and the whole
U-safe arc are the sidecars preserving coherent U-clearance. It discards
the initially chosen phase, does not preserve a prescribed owner, and makes
no body-measure inference. The decisive tests are literal cleared equations,
actual atlas graphs, both mixed-triple orientations, and exact physical
phases.

Two provisional mistakes were caught and repaired before freezing:

* I proposed using `t<=31,950` for a whole six-component. Reopening
  `overnight12_20260906_lrc_decoder_descent.md` Section3 refuted that step:
  its six-subset must grow along an actual edge inside a component with at
  least seven vertices. A whole six-component cannot do so. No t cap is
  used below. Section5 gives an actual equality row with t>31,950 satisfying
  every genuine seven-through-twelve subset cap.
* The root caught an overly permissive provisional atlas evaluator. The
  actual atlas requires **every** prime of p+q to be2mod3, with exponent at
  most two, including the final prime after trial division. It is not
  enough to restrict exponents only at inert primes. The frozen program
  rebuilds all5,855 pairs independently by multiplicative inert sums. The
  three primary controls survived unchanged; one provisional residual
  control split and was replaced by the strictly checked row below.

These are correction lineage, not live assumptions. Targeted current-result
and canon searches for the selected pair-gcd product, reverse crossing,
nonunit cutoff, and coprime-star reduction found no matching closure in the
searched surfaces. Existing unit, endpoint and phase-gluing overlaps remain
credited; this is not an external priority claim.

## 2. Exact crossing and the stronger native criterion

Actual equality identifies V_dec with the direct sum of the two internal
weighted kernels. Every nonzero mixed relation on at most three labels lies
outside it: one component contains a single occupied label, so its partial
weighted sum cannot vanish. Equality excludes every such row of height<=Q.
The actual box gives q<=Q for every internal primitive pair.

Recall the signed box

    B_Q(p,q)={p r+q s : r,s integers, |r|,|s|<=Q},
    R=Q(p+q)-(p-1)(q-1).

For coprime `1<=p<q<=Q`, every integer in [-R,R] belongs to the box and
`R>Qq`. These are the inherited elementary signed-semigroup statements,
not approximations. For completeness, complement x to Q(p+q), choose the
coefficient of p in0,...,q-1 modulo q, and the other coefficient is between
0 and2Q. Subtracting both from Q gives the required signed coefficients.

Use the physical outside label gu and set

    delta=gcd(tD,gu), c=tD/delta, x=gu/delta.

If c<=Q and x<=R, a signed-box witness `x=p r+q s` gives the literal row

    c(gu)-r(tA)-s(tB)=0.                         (3)

All coefficients have magnitude<=Q, and c>0 distinguishes the nonzero
outside coordinate. Thus this crossing is forbidden by entry. In particular,

    x<=R  ==>  c>Q  ==>  t>Q delta/D>=Q/D.       (4)

The physical gcd delta cannot be erased from (3). It can be dropped in the
last inequality because it is a positive integer.

A stronger **native sufficient criterion**, with no failure assumption, is

    gu/gcd(tD,gu) <= R,
    7(b+1)LD <= aQ gcd(tD,gu).                    (5)

Under actual equality these imply strict safety by (4) and the next section.
They depend on the actual scales. The scale-free theorem (1) is a uniform
sufficient consequence of (5) in the hypothetical-failure branch.

## 3. Coherent phase lifts and proof of the main theorem

The cited lower-dimensional supplier gives phases eta,zeta for V,U with
clearances at least1/(a+1),1/(b+1), respectively. The U minimum-clearance
function is L-Lipschitz. Hence its closed arc around zeta of radius

    [1/(b+1)-1/14]/L = a/[14(b+1)L]

is weakly1/14-safe. The full arc length is

    ell=a/[7(b+1)L].                              (6)

The t physical times `(eta+j)/t`, `0<=j<t`, preserve every V-clearance.
In U's clock their images are `g eta/t+g j/t` modulo one. Since gcd(g,t)=1,
these are a complete translated t-grid. If t ell>1, one point lies in the
interior of the safe arc and every U-clearance is strictly greater than1/14.
All V-clearances are already strict at that target.

This proves (5): (4) gives `t>Q delta/D`, and its second line implies
`(Q delta/D) ell>=1`. Equality in the product threshold is allowed.

For (1), if g>G_b, the inherited subset-gcd supplier already yields weak
safety. Otherwise g<=G_b and

    x=gu/delta <= G_b u <= G_b L
      <= G_b aQ/[7(b+1)D]
      <= QB/D = Qq < R.                           (7)

The last non-strict inequality is exactly the second condition in (1).
Thus (4) gives t>Q/D, while the first line of (1) gives `(Q/D)ell>=1`.
The coherent grid produces a strictly safe time. This proves the theorem.

The original unit version survives with its extra gcd: when1inU, the
outside label is g, so delta=gcd(tD,g)=gcd(D,g). Because g<=90 and R>Qq,
entry forces `tD/gcd(D,g)>Q`. Consequently

    7(b+1)LD <= aQ gcd(D,g)

closes that branch. This was the initial dual proposal. The nonunit theorem
uses the retained physical gcd and target bound instead of pretending it
has that same simplified delta.

## 4. Automatic nonunit classes and a structured six-component family

For a connected primitive a-component, the inherited boundary-edge gcd
lemma implies that **every** pair gcd is at most355^(a-2). Start with the
two vertices and adjoin vertices along actual outgoing edges until the
whole component is reached. At each step the gcd loss divides an oriented
primitive edge coefficient, which is at most355. The final gcd is one.
There are a-2 additions. This argument takes place inside V and makes no
claim that V can grow into U.

For a=2,3,4 the B-threshold in (1) is one, so every pair qualifies. For a=5
it is three; choose a pair containing maxV, which is at least five because
the labels are positive and distinct. The theorem therefore gives:

| Actual split | All primitive U with maxU at most | Unit required? |
|---|---:|---|
| 2+11 | 13,520,696,477 | No |
| 3+10 | 62,323,312 | No |
| 4+9 | 257,485 | No |
| 5+8 | 1,007 | No |

These bounds are `floor(aQ/[7(b+1)355^(a-2)])`. The same crude comparison
for6+7 gives L<=3, incompatible with seven positive distinct U-labels.
Thus no automatic all-six conclusion is drawn from that comparison.
The1+12 case has no V-pair and is outside this theorem.

A useful six-component family retains more arithmetic. Let q_1,...,q_5>=2
be pairwise coprime, let p_i>0 with gcd(p_i,q_i)=1, and suppose p_i+q_i is
an actual admissible atlas sum<=356. Put

    P=product_i q_i,  V={P} union {p_i P/q_i : 1<=i<=5},

and require six distinct labels. The displayed star is an actual connected
atlas graph. Its collective gcd is one: for a prime dividing q_i, the ith
leaf has valuation zero. For distinct leaves i,j one has exactly

    gcd(p_i P/q_i,p_j P/q_j)
      = [P/(q_i q_j)] gcd(p_i,p_j).                 (8)

To see the last identity, remove the common factor P/(q_iq_j). Coprimality
of q_i,q_j and of each numerator with its own denominator shows that
`gcd(p_i q_j,p_j q_i)=gcd(p_i,p_j)` prime by prime.

If two leaf numerators are coprime, their pair gcd is a product of three
denominators, at most355^3. Each leaf is at least2^4=16, so B>=10 is automatic.
Therefore every stated actual6+7 entry of this star form closes for

    maxU <= floor(3Q/[28*355^3])=1,359,              (9)

without a unit in either primitive shape. No claim is made that every
six-component has this star arithmetic or two coprime leaf numerators.

For the sharper displayed star, take

    q=(179,181,183,185,187), p_i=356-q_i.
    V=(185370716505,189592176609,193905909065,
       198314972625,202822562745,205114343115).

Every displayed edge has sum356=2^2*89, and the shape is primitive. Its
selected first two labels have D=5,929,017 and reduced pair(31265,31977).
Thus **every stated actual entry with this V and arbitrary primitive
seven-component U of maxU<=10,261 closes**. Its minimum exceeds the previous
unit theorem's balanced cutoff60,843,134,147. Fixed Q and the physical box
give a finite, enormous class; this is not described as infinitely many
in-box rows.

## 5. Actual controls, hostiles, and open boundary

The program uses the strict all-scale atlas, all internal height gates and
both mixed-triple orientations. For six-plus-seven there are exactly231
mixed triples. It also has a separate no-crossing certificate: if
`t D_min>2Q g maxU`, any nonzero smaller-component partial sum on one or two
labels is too large to cancel the larger partial sum. Internal equality
plus no mixed crossing identifies W=V_dec by the inherited general decoder.

For the displayed V, two genuine equality controls are:

| U | t, g | Physical sum | Strict physical phase |
|---|---|---:|---|
| (2,3,4,6,9,12,18) | 3,448,015; 1 | 4,051,833,733,739,682,014 | 1379207/6896030 |
| (2,3,4,6,18,486,9234) | 1,768,827,685; 1 | 2,078,585,993,174,527,392,593 | 141506215/707531074 |

Both are inside Q^2, neither shape contains one, both actual graphs have
exactly components6+7, and all231 mixed tests are negative. Their exact
clearances are1379197/6896030 and70752184/353765537. These phases are
independently evaluated on all thirteen physical speeds; the proof does
not assume those particular phases are supplied in a general row.
Every subset-gcd maximum for sizes7,...,12 is respectively `(3,3,3,3,1,1)`,
within `(90,30,9,4,2,1)`. In particular, their t values refute extending the
eleven-core-derived31,950 cap to a whole six-component, even on genuine
actual equality rows obeying all those scalar caps.

For comparison, setting t=g=1 with the first U preserves the actual graph
and all scalar caps but produces231 mixed crossings. The first literal row
in the output has outside coefficient5,929,017 and inside coefficients
567,869,224,579 and−555,225,046,329. This is a genuine graph control showing
that a convenient partition plus scalar caps does not establish entry.

The selected condition is not necessary for safety. The primitive connected
six-shape

    V'=(68215446373,80770436503,83857204777,
        87189279139,94716411151,103665993307)

comes from the complementary products of `(127,139,151,157,163,193)`. Its
actual strict atlas is connected, minimum exceeds the old balanced cutoff,
and minimum pair gcd is418,499,671. With the second U above, t=25,059,565,
g=1, the physical row is an in-box actual equality entry, all genuine subset
caps pass, and all231 mixed tests are negative. But the new product gate
fails; its scalar L cutoff is only145. It is nevertheless strictly safe at
10023827/50119130 with clearance5007296/25059565. Thus the first failed
implication would be from failure of a sufficient product inequality to
unsafety. The surviving open obligation is to close additional actual
entries for which every eligible pair has a large gcd-times-global-maximum,
or to exploit more of the native delta and exact signed-box targets.

The arithmetic-only controls separately show why q<=Q matters, why the
B-threshold is used for automatic target containment, and why coprime scales
are needed for a complete phase grid. They are not asserted to be actual
decoder entries. A global six-component pair-gcd bound strong enough to
close arbitrary U is not proved here.

## 6. Reproduction and finite scope

The standalone source and transcript are
`04-computation/continuing1_20260906_lrc_dual_pair.py` and
the corresponding `.out`. The source uses no producer imports and only the
Python standard library. Every check is active under optimization.

```text
python -B 04-computation/continuing1_20260906_lrc_dual_pair.py
python -B -O 04-computation/continuing1_20260906_lrc_dual_pair.py
```

The finite controls are all5,855 strict atlas pairs; all481 coprime pairs
`1<=p<q<=B<=16` and all164,259 points across their full signed support plus
one point on each side;5,112 literal physical normalizations in the stated
small scale bank;2,511 exact coprime-grid controls;1,200 complete numerator
choices for the explicitly fixed five-denominator identity bank; and the
four explicitly stated physical rows. The last includes both a graph-only
hostile and a safe actual entry outside the sufficient gate. This is not a
height census of LRC inputs. The analytic proofs establish the arbitrary
parameters; the controls test normalization, strict endpoints and entry.

The final normal and optimized runs pass347,547 always-active gates with
identical LF output bytes. Source and output hashes are supplied in the
freeze message and independent audit; no upstream file was edited.
