# The minimum connecting-tree consumer closes the entire upper clock range

**Status: PROVED + FINITE-EXACT; [independent audit accepted](continuing8_20260906_lrc_minimum_tree_audit.md).**
The analytic result below converts an unknown connected actual decoder graph
into a minimum-tree overlap certificate on its seven positional sheet labels.
Its complete finite consumer excludes every inherited clock at least 12,000:
549 exclusions, leaving 7,646 necessary clocks with maximum 11,995. This is
a scoped entry reduction; LRC(14) remains open.

## 1. Inheritance and the missing quantifier

The proved grid/forest inequality is inherited from
[third_20260906_grid.md](third_20260906_grid.md), using **THM-2060**,
`01-canon/theorems/THM-2060-crt-tail-coset-saturation.md`, and **THM-2064**,
`01-canon/theorems/THM-2064-multitail-sheet-capacity-and-dyadic-seam.md`.
The strict ratio atlas is **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
The complete projected profiles and the prior necessary clock array are the
audited inputs of
[continuing7_20260906_lrc_residue_domains.md](continuing7_20260906_lrc_residue_domains.md).
The latter leaves 8,195 clocks, maximum 14,868.

The canonical hostile is a selected grid killed by a tail divisible by its
scale. The corrected near miss is independent minimization of located
intervals, which loses their shared translation. The least-used operation is
to minimize over **all possible connecting trees** after conditioning every
edge on its two actual sheet margins. The earlier forest inequality and
maximum-tree bounds for a known actual graph are not new claims here.

The live concepts are: all 126 complementary profiles; margin-conditioned
strict ratios; placed open intervals; unknown actual graph; positional
minimum trees; and native realization of surviving relaxations. No tournament
orientation is imposed on the symmetric intersection relation.

## 2. Analytic minimum-tree theorem

Let `t>=1`, let the seven actual positive speeds be `u_i`, and set

    d_i=gcd(t,u_i),
    X_alpha={(alpha+j)/t mod1:0<=j<t},
    D_u={x:||ux||<1/14},
    E_t(d)=sum_i d_i ceil(t/(7d_i))-t.

Assume an actual graph on these seven labels is connected and each of its
edges has primitive ratio `p<q` in the strict atlas: `gcd(p,q)=1`,
`p+q<=356`, and every prime factor of `p+q` is `2 mod3` with exponent at most
two. The atlas has 5,855 ratios.

For each unordered margin pair `(a,b)`, put `e=gcd(a,b)` and `n=t/e`. Consider
every strict ratio satisfying

    {e gcd(n,p), e gcd(n,q)} = {a,b}.                 (1)

For such a ratio, let `m(n;p,q)` be the **minimum over every real translation**
of the number of points of a translated `n`-grid in `D_p intersect D_q`.
Define

    w_t(a,b)=min_(compatible p,q) e m(n;p,q).          (2)

An empty candidate set means that edge is impossible. Retain seven separate
vertices even when some margins coincide. Let `M_t(d)` be the minimum
spanning-tree weight for these pair weights on the seven positions.

**Theorem.** If the graph of possible edges is disconnected, no actual
connected graph realizes these margins. Otherwise every translated grid has
at least

    M_t(d)-E_t(d)                                    (3)

weak-safe points. In particular a strictly positive value proves a common
grid point of clearance at least `1/14`.

**Proof of the sheet map.** For an actual edge write its two speeds as
`D p,D q`, with coprime `p,q`. Then

    gcd(d_i,d_j)=gcd(t,D)=e.

Writing `D=eD'` and `t=en` gives `gcd(D',n)=1`. Multiplication by `D`
sends the actual `t`-grid onto a translated `n`-grid, each image having
multiplicity `e`; the coprime cofactor `D'` permutes its points. Moreover
`gcd(t,Dp)=e gcd(n,p)`, so this edge occurs in (1). Its actual intersection
count is therefore at least the lower bound (2).

Choose any spanning tree of the actual connected graph. At a grid point
lying in `r` danger sets, the induced forest has at most `max(r-1,0)` edges.
Thus the sum of its six actual pair intersections is bounded above by the
total excess of danger multiplicity over union size. Each edge separately
has its uniform bound (2), giving

    multiplicity excess >= sum_(edge in actual tree) w_t(d_i,d_j)
                        >= M_t(d).

The inherited marginal ceiling bound subtracts `E_t(d)`, proving (3).

The independently minimizing ratios in (2) need not coexist in one physical
row, nor need their minimizing translations agree. They supply lower bounds
for every actual edge and every actual translation, which is the sufficient
direction used above. No converse realization is inferred. A **maximum** tree
over merely possible edges would invent edges that the actual graph might
not contain, and would invalidate the argument.

## 3. Exact strict-endpoint count

Intersect the literal open danger intervals for `p` and `q`. With common
denominator `L=14pq`, join the two origin pieces on the circle and write the
disjoint resulting intervals as `(a_j/L,b_j/L)`, allowing the origin
interval's left endpoint to be negative.

For one interval, write

    n(b-a)=hL+r,    0<=r<L,
    A=na mod L,    B=nb mod L.

If `r>0`, its grid count equals `h` plus the indicator of the **open**
projected arc from `A` to `B`. If `r=0`, the count equals `h` except at its
coincident endpoint, where it is `h-1`. Sum the bases `h`; sweep all projected
walls, retaining start, end, and coincident-point defects. At a wall, end
arcs and punctures are absent; immediately after it, start arcs are present.
The minimum over walls and intervening chambers is the exact `m(n;p,q)`.

The independent literal control evaluates all those walls and chamber
midpoints directly using

    # integers in (na/L-alpha, nb/L-alpha)
      = ceil(nb/L-alpha)-floor(na/L-alpha)-1.

The separate interval minimum
`sum_j(ceil(n(b_j-a_j)/L)-1)` is used only as a lower bound. If that bound is
already at least the current minimum in (2), skipping the more expensive
located count cannot change that minimum. No monotonicity in `t` of the
located count is assumed.

## 4. Complete finite universe and consequence

The input alphabet consists of the 42 proved sheet states, restricted to
those dividing the current clock. For each sorted seven-multiset, keep
repeated entries as separate positions, require total gcd one, and for every
nonempty proper subset `I` require

    (c, sort(gcd(c,d_j):j outside I)) in P_(7-|I|),
    c=gcd(d_i:i in I).                               (4)

These are all 126 inherited full profile conditions. The producer uses only
necessary partial-completion pruning derived from this same table; it does
not prune on a partial minimum-tree weight. Adding a vertex can reduce a
minimum tree, so the latter would be unsound. For small domains, a second
unpruned full-multiset enumeration checks completeness. The independent
referee uses unpruned enumeration for the entire domain bank.

The declared exclusion universe is **every** member `t>=12000` of the pinned
continuing7 array, not a selected list of owners or favorable residues. There
are 549 clocks, 122 distinct divisor domains, 1,083,014 distinct-domain valid
words, and 1,532,310 word/clock evaluations. All have positive (3).
The smallest margin is 93, at

    t=12672, d=(2,8,9,9,36,48,64), E=108, M=201.

The complete range therefore leaves 7,646 necessary clocks, maximum **11,995**.
The new array's canonical JSON SHA-256 is

    8ffc6d14b3883cf7e02c3ab02ddca5339d909a8051411def9096dee83b0aaed7.

**Physical consumer.** In a hypothetical primitive thirteen-speed failure,
choose six labels `A` and the complementary seven `B`, with the actual strict
graph on `B` connected. Let `t=gcd(A)` and `g=gcd(B)`. The inherited reduction
has `g<=90`, `gcd(t,g)=1`, the profiles (4), and the continuing7 clock array.
The same CITED proper-six supplier used in the inherited entry theorem gives
a safe phase for `A/t`. Its `t` lifts preserve those six clearances; in the
`B/g` clock they form exactly a translated `t`-grid. Formula (3) then supplies
a common weak-safe lift whenever `t` is one of the excluded clocks.

Consequently every such hypothetical failure has `t` in the new array. This
does not assume decoder equality, boxed speeds, or connectivity of the
six-label graph. It also does not supply an entry theorem for complements
whose actual strict graph is disconnected.

## 5. The unresolved lower-clock object

The same complete consumer was also tested at `t=7200`, as a deliberately
hostile control. Among all 76,814 valid full words, exactly 15 survive (3).
The worst is

    d=(4,9,9,16,18,24,32), E=74, M=0.

Its certificate retains an entire zero-weight minimum tree and each edge's
compatible native ratio. These are necessary projected survivors, not
unsafe speed rows. The one-pair route and this full-tree route therefore
still leave a real compatibility problem: which low-credit ratio choices
can coexist with a single global valuation assignment, the omitted physical
profiles, and a common phase? The bound intentionally discarded those data.

The source-target map retains each pair's margins and the existence of an
actual connecting tree. It loses the tree's actual ratios and their joint
valuation data. The minimum is the correct direction for that information
loss. A next precise object is a tree carrying compatible prime-valued
potentials; a favorable independent edge assignment is not its realization.

## 6. Reproduction and validation

Run the standalone Python source
`04-computation/continuing8_20260906_lrc_minimum_tree.py` normally and with
`-O`. It reconstructs the ratio atlas, full profile domains, every compatible
edge minimum, and all positional tree tests. NumPy performs exact bounded
integer operations, with no floating-point weights. The certificate records
every pair table, complete-domain word hashes, all minimum owners and the
entire hostile survivor list. The output and JSON are compared as raw LF
bytes between modes; full gate counts are recorded in the final manifest.

The independent referee reconstructs the atlas by a separate multiplicative
method, its interval geometry and endpoint sweep independently, the complete
full-word bank by unpruned native enumeration, and trees by Kruskal instead
of vectorized Prim. Verification of the whole finite universe is separate
from the analytic tree proof. This result does not claim any of the remaining
clocks is itself a counterexample.
