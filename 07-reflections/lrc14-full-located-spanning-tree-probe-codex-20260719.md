# LRC(14) full located spanning-tree probe

**Verdict: SOUND, with two essential qualifiers.**  The six-comb cover must
be by the actual **strict** danger teeth, and the tree must be extracted from
the consecutive handoffs of an irredundant tooth chain.  It is not valid to
attach `gcd(d_i,d_j)/(14d_i d_j)` to an arbitrary component of
`(D_{d_i} intersect D_{d_j}) intersect G`: a slow-gap endpoint may clip that
component, exactly as the counterexample in THM-1199 shows.  The private-mass
input from THM-1198 is what permits a special five-edge tree for which no
edge is clipped.

The resulting theorem is an existence statement:

> If the six strict danger combs of
> `c<d_1<...<d_6` cover a complete closed `c`-safe gap `G`, then there is a
> spanning tree `T` on the six fast labels such that
>
> ```text
> H:=sum_i 1/d_i
>   >=1/c+(49/6) sum_({u,v} in T) gcd(u,v)/(14uv)
>   =1/c+(7/12) sum_({u,v} in T) 1/lcm(u,v).          (1)
> ```

This is the label-gcd sharpening of THM-1178 for one specially located tree.
It also pays all five Hunter edges locally, whereas THM-1244 pays two edges
on the smaller slowest-spoke component and leaves two periodic errors.

## 1. Conventions and the private-label input

Put

```text
D_d={t in R: ||dt||<1/14},
G=[A,B]=[(14k+1)/(14c),(14k+13)/(14c)],
L=|G|=6/(7c).                                         (2)
```

The strict convention is topologically material here.  Measure bounds are
unchanged if teeth are closed, but a finite closed-interval cover can hand
off at one point and have zero overlap.  In contrast, a finite strict-tooth
cover of a compact interval has positive overlaps at every handoff.

THM-1198 is proved against the stronger closed teeth.  After deleting label
`i`, it supplies a survivor of the other five closed combs of measure at
least

```text
1/(49c).                                               (3)
```

Because the original six **strict** combs cover `G`, this survivor lies in
the strict comb `D_{d_i}`.  Thus one may take a measurable private-provider
set

```text
Q_i subset G intersect D_{d_i},
Q_i intersect D_{d_j}=empty for j!=i,
|Q_i|>=1/(49c).                                       (4)
```

Only nonemptiness is needed for the base theorem, but the quantitative
length will give a multiplicity strengthening in Section 7.

## 2. The irredundant interval-chain lemma

The periodic strict teeth meeting `G` form an open cover of the compact set
`G`.  Choose a finite subcover and then delete teeth until it is irredundant.
Write its distinct open intervals as

```text
J_a=(alpha_a,beta_a),     a=1,...,N,                  (5)
```

ordered by increasing left endpoint.

There are no ties among the `alpha_a`: with equal left endpoints the shorter
interval is contained in the longer and is redundant on `G`.  More
generally, irredundancy gives

```text
alpha_1<...<alpha_N,       beta_1<...<beta_N.          (6)
```

Indeed, `alpha_a<alpha_b` and `beta_a>=beta_b` would imply
`J_b intersect G subset J_a intersect G`, so `J_b` could be deleted.

Exactly one selected interval contains `A`.  If two did, the one extending
less far to the right would be contained in the other after restriction to
`G`.  Dually, exactly one contains `B`.  Consequently

```text
alpha_(a+1)>=A,       beta_a<=B       (1<=a<N).        (7)
```

The inequalities may be equalities, but the intervals themselves are open.
Coverage forces

```text
alpha_(a+1)<beta_a.                                  (8)
```

If equality held, that common point would be omitted by every selected open
interval; if the left endpoint were larger, the intervening open segment
would be omitted.  Monotonicity (6) rules out rescue by a nonconsecutive
interval.

Equations (7)--(8) show that the **entire raw tooth overlap**

```text
W_a=J_a intersect J_(a+1)
   =(alpha_(a+1),beta_a)                              (9)
```

is a positive open interval contained in `int(G)`.  This is stronger than
saying that the two combs intersect somewhere in `G`; it is the fact that
removes the slow-wall endpoint owner.

Every private set `Q_i` must be covered by selected teeth carrying label
`d_i`, because no other label covers any of its points.  Hence all six labels
occur in the word

```text
lambda_1 lambda_2 ... lambda_N,
lambda_a in {d_1,...,d_6}.                            (10)
```

Two consecutive labels are distinct: different teeth of one fixed speed
are disjoint.  The graph whose edges are the unordered consecutive pairs
`{lambda_a,lambda_(a+1)}` is connected, since (10) is a walk visiting all
six labels.  Choose a five-edge spanning tree `T` of this graph, and for each
tree edge retain one consecutive occurrence witnessing it.

This proves the full-rank topological claim.  Private mass is used only to
force every label into the word; connectedness of the word then gives rank
five.

## 3. Exact gcd quantum on every chosen edge

Suppose a retained handoff goes from the tooth of speed `u` and address `n`
to the tooth of speed `v` and address `m`.  Because right endpoints and left
endpoints cross in chronological order, its full length is

```text
omega
 = (14n+1)/(14u) - (14m-1)/(14v)
 = [v(14n+1)-u(14m-1)]/(14uv)>0.                     (11)
```

(If the chronological labels are reversed, rename `u,n` and `v,m` in
(11).)  The numerator is a positive integer divisible by `gcd(u,v)`, hence

```text
omega>=gcd(u,v)/(14uv)=1/[14 lcm(u,v)].               (12)
```

By (9), the whole interval of length `omega` lies in `G`, not merely its
intersection with `G`.  Therefore

```text
mu(G intersect D_u intersect D_v)
  >=gcd(u,v)/(14uv)       for every {u,v} in T.        (13)
```

This is exactly the positive-divisor lemma already isolated in the Lean core
of THM-1244.  Addresses should be represented as integers in a full formal
proof, since a real lift of `G` may meet negative-address teeth.

## 4. Hunter with actual located credits

Let

```text
A_i=G intersect D_{d_i},
Q_T=sum_({u,v} in T) mu(G intersect D_u intersect D_v).
                                                               (14)
```

For any active label set `S`, the subgraph induced by `S` in a tree has at
most `|S|-1` edges.  Thus the pointwise forest-Hunter inequality is

```text
1_(union_i A_i)
 <=sum_i 1_(A_i)-sum_({u,v} in T)1_(A_u intersect A_v). (15)
```

The six combs cover `G`, so integration gives

```text
L+Q_T<=sum_i mu(G intersect D_{d_i}).                  (16)
```

THM-1166's sharp singleton interval discrepancy, also used in THM-1237 and
THM-1244, is

```text
mu(G intersect D_d)<=L/7+6/(49d).                     (17)
```

Putting `H=sum_i 1/d_i`, equations (16)--(17) yield

```text
L+Q_T<=6L/7+6H/49,
L/7+Q_T<=6H/49.                                       (18)
```

Since `L=6/(7c)`,

```text
H>=1/c+(49/6)Q_T.                                     (19)
```

Finally sum (13) over `T` and insert it in (19), proving (1).

No pair-period discrepancy `rho(1-rho)/gcd` occurs.  Those errors enter when
one transports a **global average** pair mass to a specified interval.  Here
each of the five credits is already an actual interval lying inside `G`.

## 5. Endpoint and quantifier audit

The following distinctions are necessary.

1. **Strict versus closed cover.**  The theorem uses a strict open cover.
   Closing the teeth is harmless for (17), but not for the handoff lemma: two
   closed teeth may cover through a zero-length meeting.
2. **A special tree, not every nerve tree.**  THM-1199 gives the exact row
   `(c,d_i,d_j,k)=(7,10,34,2)`, where a clipped intersection component has
   endpoint owners `(7,34)` and length `1/3332`, smaller than
   `gcd(10,34)/(14*10*34)=1/2380`.  Thus an arbitrary label-intersection edge
   cannot receive the label gcd.  Consecutive interior handoffs avoid this
   counterexample.
3. **Full raw overlaps, not positive intersections after clipping.**  Merely
   knowing `mu(G intersect D_u intersect D_v)>0` does not give (12).  Its
   endpoints can include a wall of speed `c`.  Equation (9) is the missing
   topological statement.
4. **Existence, not simultaneous switching.**  One tree is selected once
   from the chronological quotient graph and then used pointwise in (15).
   There is no phase-by-phase tree choice.
5. **Repeated labels are allowed.**  The chronological word need not be a
   Hamiltonian path.  Projecting its connected walk to six labels and then
   choosing a tree is valid; it does not assume that each label occurs once.

With these qualifications, there is no open-endpoint or inequality gap.

## 6. Lean-ready algebraic core

It is useful to clear all divisions before formalization.  Put

```text
S_T=sum_({u,v} in T) gcd(u,v)/(uv).                   (20)
```

The analytic/topological providers give

```text
c>0,
7cL=6,
49L<=42L+6H-49Q_T,          -- Hunter + six singleton bounds
S_T<=14Q_T.                 -- five located tooth quanta
                                                               (21)
```

The desired conclusion is the polynomial inequality

```text
12+7c S_T<=12c H.                                    (22)
```

Indeed the middle relation in (21) gives
`7L+49Q_T<=6H`.  Multiply by `2c`, use `14cL=12`, and use
`7cS_T<=98cQ_T`.  Division by `12c>0` gives (1).

A robust Lean skeleton is therefore:

```lean
theorem full_located_tree_debt
    {c L H Q S : ℝ} (hc : 0 < c)
    (hgap : 7 * c * L = 6)
    (hhunter : 49 * L ≤ 42 * L + 6 * H - 49 * Q)
    (hquant : S ≤ 14 * Q) :
    12 + 7 * c * S ≤ 12 * c * H := by
  have hc0 : 0 ≤ c := le_of_lt hc
  have htree : 7 * L + 49 * Q ≤ 6 * H := by
    linarith
  have htreeScaled :
      2 * c * (7 * L + 49 * Q) ≤ 2 * c * (6 * H) :=
    mul_le_mul_of_nonneg_left htree (by positivity)
  have hquantScaled : c * S ≤ c * (14 * Q) :=
    mul_le_mul_of_nonneg_left hquant hc0
  nlinarith [htreeScaled, hquantScaled]
```

The discrete provider can reuse THM-1244's
`positive_overlap_numerator_quantum`; the new paper-facing provider is the
finite irredundant interval-chain lemma in Section 2.

## 7. Stronger consequences already latent in the proof

### 7.1 Common-gcd and ordered-star debt

Let

```text
g_0=gcd(d_1,...,d_6).
```

Every tree-edge gcd is at least `g_0`.  Root `T` at `d_6`.  Each nonroot
vertex `d_i` contributes an edge to a parent no larger than `d_6`, so

```text
sum_({u,v} in T) gcd(u,v)/(uv)
 >=(g_0/d_6) sum_(i=1)^5 1/d_i
 =(g_0/d_6)(H-1/d_6).                                 (23)
```

Hence

```text
cH-1
 >=(7c g_0)/(12d_6) (H-1/d_6).                       (24)
```

This strictly strengthens THM-1178's ordered tree-free seam debt whenever
the six fast speeds have a nontrivial common gcd.  It remains scale-sensitive
when `g_0` is small; no uniform LRC(14) contradiction follows from (24).

### 7.2 Private mass forces a long owner word

Let `n_i` be the number of selected teeth carrying label `d_i`.  They cover
the private set `Q_i`, and one full `d_i`-tooth has length `1/(7d_i)`.
Therefore

```text
n_i/(7d_i)>=|Q_i|>=1/(49c),
n_i>=ceil(d_i/(7c)).                                  (25)
```

Thus

```text
N>=sum_i ceil(d_i/(7c)).                              (26)
```

This is more structural than mere label occurrence: high-speed owners must
recur many times in the chronological tooth word.  It is a direct possible
interface to the finite relative-address/owner-reuse program.

### 7.3 Repeated handoffs can be credited on one Hunter edge

For an unordered label pair `{u,v}`, let `m_uv` be its number of consecutive
occurrences in (10).  The corresponding handoff intervals are pairwise
disjoint.  If two such intervals met, their intersection would lie in two
distinct teeth of at least one of the speeds `u,v` (two handoffs can share
one middle tooth, but then their other two teeth have the same label), which
is impossible.  Hence

```text
mu(G intersect D_u intersect D_v)
 >=m_uv gcd(u,v)/(14uv).                              (27)
```

Give the complete graph on six labels edge weight

```text
w_uv=m_uv gcd(u,v)/(uv),                              (28)
```

with weight zero when the pair never occurs.  Its positive support is
connected.  A fixed edge belongs to exactly one third of the `6^4` labelled
spanning trees of `K_6`; averaging their weights therefore gives a tree
`T_*` with

```text
sum_(uv in T_*) w_uv >=(1/3)sum_(u<v)w_uv.            (29)
```

A maximum-weight tree has no zero edge because the positive support is
connected.  Applying Hunter with (27) strengthens (1) to

```text
H>=1/c+(7/12)sum_(uv in T_*)m_uv/lcm(u,v)
 >=1/c+(7/36)sum_(u<v)m_uv/lcm(u,v).                  (30)
```

Since `sum m_uv=N-1` and every pair has
`1/lcm(u,v)>=g_0/d_6^2`, equations (25)--(30) imply the completely scalar
but weaker corollary

```text
H>=1/c
 +[7g_0/(36d_6^2)]
    [sum_i ceil(d_i/(7c))-1].                         (31)
```

The multiplicity upgrade (27)--(31) is not needed for the base theorem and
should receive its own short referee/Lean audit before canonization.  Its
proof uses only disjoint same-speed teeth, Cayley's spanning-tree count, and
THM-1198's six private-length floors.

## 8. Structural and tournament reading

The ordinary speed-order tournament is transitive and loses the theorem.
The faithful first-layer vertices are individual tooth occurrences, with
chronological path edges.  Quotienting occurrences by speed label preserves
connected graphic rank and handoff multiplicity but destroys tooth addresses
and overlap positions.  The minimum sufficient carrier is

```text
(irredundant labelled tooth word;
 interior handoff intervals;
 quotient edge multiplicities m_uv;
 overlap numerators and their gcd divisors).           (32)
```

Orienting a quotient edge by chronological first occurrence gives a
relation, but repeated handoffs may occur in both orientations, so it is not
a tournament.  Its useful fingerprints are connected rank five, the
weighted maximum spanning tree in (29), repeated-edge counts, and the owner
word.  This challenges the runner-vertex assumption directly: runner labels
carry the final Hunter tree, while tooth occurrences carry the reason its
five edges are genuinely located.

## 9. Exact scope

The proposed lemma is valid and closes a real gap between THM-1178/1199 and
THM-1244: there exists a full five-edge label-gcd tree of actual overlaps on
the original slow gap.  It does **not** prove six-comb noncoverage or
LRC(14).  The correction can still shrink under large relatively coprime
speeds.  The highest-leverage next use is not another scalar relaxation, but
to combine the long private-owner word (25) with the centered blocker's
relative-address/positive-holonomy cycle: repeated owner occurrences must
either reuse an oriented germ or spend additional located pair mass through
(27).
