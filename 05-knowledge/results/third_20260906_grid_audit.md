# Independent audit of translated-grid excess and the balanced scale cutoff

**Status: PROVED relative to the named inherited and CITED suppliers;
FINITE-EXACT independent audit passed.** The translated-grid mechanism and
the proposed coarse scale bound are correct with strict bad arcs and weak
safety. An independent optimization of the complete actual inert pair atlas
improves the sufficient integer cutoff to **`t>=97,097`**. This does not
prove LRC(14), classify a failing row, or prove that the sufficient cutoff
is sharp for actual decoder entries. No theorem ID or external priority is
claimed. This note accepts the full proof in the
[producer theorem](third_20260906_grid.md), including its alternative origin
credit, optional forest credit, unit-case refinement, and physical phase
transfer. Its standalone source imports no producer code.

## 1. Inheritance and the first hostile

The entry geometry and the actual 5,855-pair atlas come from **THM-3818**,
[scaled inert cube-class support-two packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md),
Section1. The hypothetical-failure gcd ceilings come from the independently
audited [joint-shadow theorem](lrc14_joint_shadow_empty_core_next_sep06.md).
Its lower-runner existence input is **CITED**, as recorded in that supplier;
this audit does not re-prove its external computer-assisted theorem.

The pair-overlap supplier is the main identity and envelope in **THM-739**,
[pairwise-coprime bad-overlap exact Bernoulli form](../../01-canon/theorems/THM-739-pairwise-coprime-bad-overlap-exact-bernoulli-closed-form.md).
This slug is essential because the repository also has another THM-739.
Its microscopic addendum, as inherited at this audit's start, omitted
containment clipping: at `(p,q)=(1,2)` it gives `3/28`, whereas the true
overlap is `1/14`. That addendum is not a dependency. Section4 below gives
the repaired elementary formula, independently matching the main Bernoulli
identity on every actual atlas row. The root was notified immediately so
the demonstrated correction can be carried into the maintained canon.

The inherited windowed addendum has the same missing operation: each
resonant arc must be intersected with the window, not merely retained if
wholly inside it. Its unqualified far-ratio bulk claim also fails: for
`p=1,q=14n`, the overlap on `[0,1/14)` is exactly `1/98`, whereas one
fourteenth of the full-circle overlap is `1/686`. The script retains n=1;
the general example follows by integrating n complete q-periods inside
the window. Neither windowed claim is used here.

The live concepts are labelled grid multiplicity, strict endpoints, nested
origin excess, hereditary subset gcds, exact pair-arc components, and actual
decoder ratios. The decisive tests were the containment hostile, exact grid
phase walls, and independent reconstruction of the entire inert atlas.
The map preserves one six-body safe phase through all of its lifts; it
forgets the initial phase when bounding the remaining seven masks. Its
sidecars are coprime scales, a common labelled grid, and one actual pair.

## 2. The universal translated-grid count

Let `U=(u_1,...,u_b)` be nondecreasing positive integers, let `t>=1`, and
consider every translated grid

```text
G(theta,t)={theta+j/t mod1 : j=0,...,t-1}.
```

Define bad sets by `B_i={x: ||u_i x||<1/14}`. Weak-safe points have every
clearance at least `1/14`. Put

```text
d_i=gcd(t,u_i),
B=sum_i d_i ceil(t/(7d_i)),
R0=sum_(i=2)^b (ceil(t/(7u_i))-1).
```

Then every translate has at least

```text
t-B+R0                                             (1)
```

weak-safe points. The bound may be negative, in which case it asserts no
existence.

**Proof.** Multiplication by `u_i` maps the t-grid onto a translated
`t/d_i`-grid with multiplicity `d_i`. An open interval of length `1/7`
contains at most `ceil(t/(7d_i))` points of that quotient grid. Thus the
total number of bad incidences is at most B.

For every `i>=2`, the origin interval

```text
C_i={x: ||x||<1/(14u_i)}
```

is contained in every `B_j` for `j<=i`. Since the `C_i` are nested,
pointwise bad multiplicity minus its union indicator is at least
`sum_(i=2)^b 1_(C_i)`. An open interval of length ell contains at least
`ceil(t*ell)-1` t-grid points. Summing these origin layers gives excess
at least R0, proving (1). The open-interval correction minus one is needed
when both endpoints are grid points. Repeated u_i are allowed in this
general counting lemma.

For a selected pair, one may instead credit its entire overlap. These are
alternative lower bounds on the same multiplicity excess: **R0 and the
full pair credit must not simply be added**. A further combined credit
would need a pointwise forest or another non-double-counting certificate.

## 3. Exact excess ceilings in the seven-mask critical case

Under a hypothetical primitive thirteen-speed failure, write a partition as

```text
tV union gU,  |V|=6, |U|=7,
gcd(V)=gcd(U)=gcd(t,g)=1.
```

For any k-subset I of U, the physical subset consisting of all six tV
labels and the k selected gU labels has gcd

```text
gcd(t,(u_i)_(i in I)) = gcd((d_i)_(i in I)).
```

The equality uses both primitive V and coprime t,g. The inherited ceilings
therefore imply, for `k=1,...,7`,

```text
gcd((d_i)_(i in I)) <= (90,30,9,4,2,1,1)_(k).       (2)
```

The final entry is primitiveness of the full row. These are necessary
conditions under hypothetical failure, not restrictions on all safe rows.

Let `E=B-t`. If `7` does not divide t and `tau=t mod7`, then every d_i is a
unit modulo7 and direct division with remainder gives

```text
7E=sum_i d_i*((-tau*inverse(d_i,7)) mod7).            (3)
```

A value d can occur at most as many times as there are bounds in (2) at
least d. This relaxes all joint subset restrictions to per-value
multiplicity caps. Selecting the seven largest retained weights in (3)
gives exact residue-wise ceilings

```text
tau       1    2    3    4    5    6
E max   429  426  438  435  447  445.               (4)
```

The independent output records the seven selected divisors for each row.
Each row also passes all127 genuine subset checks in (2). Moreover take
`T=lcm(d_i)` and multiply T by a suitable element of `1,...,6` to obtain
the desired nonzero residue tau modulo7. With `u_i=d_i`, each is exactly
`gcd(T,u_i)`, realizing its value of E. Thus (4) is sharp even for this
abstract gcd-cap problem. This realization does **not** establish an
actual connected decoder component or a failing physical thirteen-row.

If `a=v_7(t)>=1`, only divisors having `v_7(d_i)=a` contribute to E.
For a=1 there are at most three such labels by the four-subset bound in
(2). For each possible residue of `t/7`, the relaxed top-three ceilings
are respectively

```text
151,134,152,135,146,122.
```

For a=2 there is at most one such label, and the singleton bound forces
`d_i=49`, giving E at most42. For a>=3 there are none. Consequently the
uniform bound **`E<=447`** is valid in every case.

## 4. Exact pair measure, components, and the grid error

Let `1<=p<q` be coprime, write `S=p+q`, and let

```text
R=ceil(S/14)-1,             J=2R+1.
```

Then `B_p intersect B_q` has exactly J open circular components and measure

```text
mu(p,q)=[p+sum_(r=1)^R min(2p,S-14r)]/(7pq).        (5)
```

**Proof.** Bad arcs centered at `j/p` and `k/q` meet precisely when their
circular center difference has an integer lift r satisfying `|r|<S/14`.
Coprimality gives one center pair for each r modulo pq; the retained range
has no repetitions. These intersections are disjoint because each family's
bad arcs are disjoint. Their individual lengths are

```text
min(1/(7q), (S-14|r|)/(14pq)).                      (6)
```

At r=0 the narrower arc is contained in the wider one, giving `1/(7q)`.
The positive and negative r contributions yield (5). If `S/14` is integral,
the two extreme touching pairs have no open intersection and are excluded;
this proves the displayed ceiling formula for J, including its boundary.

For actual speeds `u=hp,v=hq` on a translated t-grid, multiplication by h
reduces to a `t/e`-grid repeated e times, where `e=gcd(t,h)`. Apply the open
interval lower bound separately to the J components. Since
`ceil(z)-1>=z-1`, the literal overlap count is at least

```text
t*mu(p,q)-e*J.                                    (7)
```

Combining this one pair credit with total bad incidence at most B gives
at least `t*mu-E-eJ` weak-safe grid points. Under (2), e is the selected
two-subset gcd and is at most30. Thus the strict inequality

```text
t > (447+30J)/mu(p,q)                             (8)
```

forces at least one weak-safe point. Equality in this real-valued bound
does not by itself force a point.

The original coarse bound also passes audit: the main THM-739 envelope
gives `mu>=1/49-1/(4pq)`, which exceeds `1/91` for `pq>=27`; the36 coprime
pairs with `pq<=26` give `mu>=1/91` by (5). Every pair with S<=356 has
J<=51, hence `t>179907` suffices without optimizing the ratio.

## 5. Complete actual-atlas optimization and physical consequence

Reconstruct the actual atlas independently by trial division: retain every
coprime `1<=p<q`, `p+q<=356`, whose sum has only prime factors congruent
to2 modulo3, each with exponent at most two. This gives exactly **5,855**
pairs and **94** sums, agreeing with the inherited declared atlas universe.

Formula (5) equals the main Bernoulli expression exactly on every row.
The exact maximum is uniquely

```text
max_atlas (447+30J)/mu = 6019965/62,
(p,q)=(5,348),       J=51,       mu=62/3045.
```

Since `97096 < 6019965/62 < 97097`, condition **`t>=97097`** is a uniform
sufficient integer bound for (8). This is the exact optimum of the declared
independent E<=447, e<=30, one-pair component-error relaxation. It is not
a lower bound on any genuine counterexample and need not be the best
consequence after retaining correlations among these quantities.

To obtain an actual common time, the cited lower-runner result supplies
a V-safe phase y. Every physical time `x=(y+j)/t` preserves that six-body
safety. The gU phases are a translated t-grid because g is a unit modulo t.
If t>=97097 and U has even one pair in the actual atlas, (2)--(8) contradict
a hypothetical failure by producing a U-safe lift. This proves weak safety
for the full thirteen-speed row.

In particular the proposed balanced actual decoder equality entries close
above this scale: their seven-component has an actual edge. In fact the
proof after that edge is supplied uses neither the full decoder equality,
the finite height box, nor full U connectivity. The native theorem applies
to any stated primitive six/seven partition with at least one such U-pair.
No strict full-row safety is exported from the weak-safety count.

## 6. Reproduction and exact scope of verification

The independent [source](../../04-computation/third_20260906_grid_audit.py)
and [retained output](third_20260906_grid_audit.out) enumerate every phase
wall and open phase cell for1,818 grid rows, giving67,884 phase cases. The
universe is all distinct positive subsets of1,...,7 of sizes1,...,4, plus
the full seven-label row and repeated-label controls, for t=1,...,18.
It separately checks82,364 pair-grid phase cases for coprime p<q with
p<=7,q<=12, common factors1,...,3, and clocks1,...,12. Literal arc sweeps
check270 pairs, including the exact optimizer and large-component controls.

Every actual atlas pair is checked by both clipped resonance and Bernoulli
arithmetic. All subset gcds of every cap-maximizing row are retained.
The full script has **303,923 active requirements** and uses no assertions
that disappear under optimization. Finite controls challenge the analytic
lemmas; they do not replace their universal quantifiers.

```text
python3 -B 04-computation/third_20260906_grid_audit.py
python3 -B -O 04-computation/third_20260906_grid_audit.py
```

The retained semantic SHA-256 for the exact ordered atlas is
`892c0c48d83510e123ebda349719b93706996d47c23da1df96c006903d29c5c6`.
Normal and optimized outputs agree byte for byte. Raw LF-byte SHA-256:

```text
source 41e6e0f0f89ffe0554e2d977fee4eb5f76ba24709306acbbf0d39f087a105b25
output 2c6e58e08398aec7f61f383f3900b17806406d750714549f01a0dda6534a4471
```
