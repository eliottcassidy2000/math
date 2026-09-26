# Auditing the finite-poset barrier and quantifying the escaping boundary

**Status:** **SOURCE CLAIM / INTERNAL AUDIT PASSED, WITH EXTERNAL DEPENDENCIES RETAINED:** the attached `KSBFT_v7.pdf`, *Breaking the infinite barrier in the 1/3–2/3 conjecture*, by Max Aires, Swee Hong Chan, Igor Pak, and Greta Panova, dated September 25, 2026. Its claim is an absolute positive improvement over the BFT constant, **not** the full 1/3–2/3 conjecture. No fatal error was found in the internal derivations audited below. The cited BFT, log-concavity, window, Shepp, Haqi, and bounded-width/large-range inputs were not all independently reproved in this lane. The paper explicitly omits the derivation of its displayed global epsilon numeral; that numeral is **NOT independently certified here**.

**PROVED, elementary:** the finite Fibonacci-poset matching formulas, boundary localization, and norm/limit separation in section 3. **FINITE-EXACT:** independent controls for 5230 naturally labelled posets on 2..6 elements and the stated auxiliary universes. **OPEN:** a descent implication for a prescribed positive Collatz integer. This note supplies an exact obstruction to a careless transfer, not a Collatz proof.

Source pin: user-supplied `C:/Users/Eliott/Downloads/KSBFT_v7.pdf`, 31 pages, SHA256 `48bfecea2d4b47698861bd3eeb544f350617e5828a018e1511c67efdf1ddb759`. Text extraction and rendered pages were read; material formulas on PDF pages 2,3,8,11,13–14,16–19,21–25,29–31 were visually checked. This local source pin makes no assertion about publication or independent external acceptance.

Igor Pak's [September 25, 2026 author announcement](https://igorpak.wordpress.com/2026/09/25/how-to-make-small-improvements-and-break-barriers/) independently corroborates the title, authors, and claimed positive improvement. The attached file remains the proof source audited here; byte identity with the linked public PDF was not established.

Reproduction: [script](../../04-computation/experiments/crossroads_poset_20260926_barrier.py), [output](crossroads_poset_20260926_barrier.out). The script uses only exact integer/rational decisions and no disabled assertions. Decimal columns are display only.

## 1. Inheritance and the precise statement

For a finite poset `P`, choose a uniformly random linear extension `f`. Its balance is

```text
delta(P;x,y)=min(Pr[f(x)<f(y)],Pr[f(y)<f(x)]),
delta(P)=max_(x!=y) delta(P;x,y).
```

Write `C=(5-sqrt(5))/10`, approximately 0.2763932. The attached Theorem 1.2 claims that some fixed `epsilon>0` gives `delta(P)>=C+epsilon` for **every finite nonchain**. It does not claim `delta(P)>=1/3` in general.

The theorem is assembled from three separately quantified results:

1. **Bounded range:** if each element is incomparable with at most `D>=2` others, then `delta(P)>C+eta_D`, where the visually verified formula is `eta_D=(D+1)^(-4(D+1))/4096`. The number 4096 is a denominator, not an exponent.
2. **Large width:** sufficiently large width gives `delta(P)>e^(-1)-10^(-100)>1/3`. Its proof is in sections 5–6 of the attachment.
3. **Bounded width, unbounded range:** the imported Aires–Kahn theorem [AK25b, Theorem 1.6] says that for every fixed width bound `K` and positive error, sufficiently large range gives balance arbitrarily close to 1/2.

The combination is logically uniform: choose a fixed large-width threshold, apply the third result with a slightly larger strict width bound, and let its range threshold determine one fixed `D` for the first result. This yields some positive absolute gap. Making that gap numerically explicit additionally needs effective constants through all inputs; the attachment does not display the arithmetic producing its enormous epsilon numeral.

The closest repo mechanism is the distinction between a finite residue family and a prescribed integer in [the orbit-family note](crossroads_family_20260926_orbits.md), sections 3,5–8. A genuinely successful but different finite-to-infinite construction is [the pairing Bellman note](crossroads_family_20260926_flow.md), section 8: its integer-address selector and uniform boundary/reset estimates establish a natural density in the **modified pairing model**. Its conclusion does not transfer to the original Collatz map.

The canonical hostile here is the infinite Fibonacci poset, defined by `x_i<x_j` iff `j-i>=2`, for `i,j` in the integers. Its probabilities are defined in the paper by limits of finite restrictions; we do not assume an undefined uniform measure on all infinite linear extensions. The corrected near miss is interchanging a maximum over moving labels with a limit at fixed labels. The least-used useful coordinate is the displacement from the ordered list of expected ranks.

| Live object | Question / operation | Essential extra information |
|---|---|---|
| Finite expected ranks | Select an extremal displacement | An attained maximum and the finite endpoints |
| Near-BFT triples | Quantitative stability | Positive local-event probability, not just formal possibility |
| Fibonacci restrictions | Take a limit / discard boundaries | Moving labels and the norm being measured |
| Poset modification | Add very likely relations | Conditional law, variance, and original versus modified heights |
| Collatz parity families | Regard words as linear extensions | The actual source integer, carry, height cutoff, and conditional measure |

## 2. Internal proof audit

### 2.1 Local probability cannot be arbitrarily small at bounded range

For a conjunction `E` of comparisons involving a set `S` of `m` elements, the attachment's Lemma 3.1 gives

```text
Pr(E)>0  =>  Pr(E)>=(D+1)^(-m(D+1))=:q_m.             (1)
```

Its proof survives audit. Add the comparisons to obtain `Q`; let `T` consist of `S` and all elements incomparable with something in `S`, so `|T|<=m(D+1)`, and set `U=X\T`. Added relations do not change the induced order on `U`: any newly proposed route through `s in S` would already be oriented the same way by the old comparisons with `s`. Thus `Q[U]=P[U]`.

Every extension of an induced subposet extends to the full finite poset, giving `e(Q)>=e(P[U])`. Once the order on `U` is fixed, each `t in T` has at most `D+1` possible ranks, because its predecessors and incomparable elements bound its allowable rank interval. Hence `e(P)<=e(P[U])(D+1)^|T|`, proving (1).

This is a uniform-in-size **local-event floor**, not a bound on arbitrary events or arbitrary changes of measure. Its dependence on `D` is load-bearing.

### 2.2 Near equality first forces an exact structure

A BFT triple has ordered mean heights `h(x)<=h(y)<=h(z)<=h(x)+2`, and does not form a chain. Cases other than the `2+1` configuration are excluded under the near-BFT hypothesis by the cited BFT bounds and the attachment's Appendix A. In the remaining case `x<z`, with `y` incomparable to both, the proof retains a slack variable `R>=0` that earlier BFT estimates could discard:

```text
max(delta(P;x,y),delta(P;y,z)) >= C+R/(5+3sqrt(5)).    (2)
```

Here `R` dominates the probability that either reversed adjacent pair is separated by at least one intervening element. The normalized optimization in Claim 2.5 and its convex-chord estimate were checked algebraically. Appendix A's two packed-height cases also pass the internal algebra audit, subject to the packed-height and monotonicity inputs quoted from BFT.

Under `delta(P)<=C+eta_D`, (2) gives `R<q_3`. Therefore every three-element event implying a separated reversal is impossible by (1). The resulting exact structure is:

- `z` covers `x`;
- `y` is incomparable only with `x,z`;
- writing `L={w:w<y}` and `U={w:y<w}`, the whole poset partitions into `L`, the triple, and `U`, with `L<z`, `x<U`, and `L<U`;
- the two incomparability edges at `y` are bridges.

The cover argument does not assume a false converse: if something lay strictly between `x,z`, then every extension would have gap at least 2, and the positive event placing `y` between them would make its expectation strictly larger than 2. The comparison-conjunction consistency steps likewise preserve their directions.

### 2.3 The finite boundary selects the right triple

Order the elements `v_1,...,v_n` by increasing expected rank, and define

```text
d_i=h(v_i)-i,    M=max_i |d_i|.
```

For a connected incomparability graph on at least three vertices, Lemma 4.1 selects a consecutive BFT triple with `max(|d_k|,|d_(k+2)|)=M`. The finite endpoint inequalities `d_1>=0`, `d_n<=0`, `d_1+d_2>=0`, and `d_(n-1)+d_n<=0` force a maximal positive displacement far enough from the right edge, or a maximal negative displacement far enough from the left edge. The resulting mean span is at most 2.

If the triple were a chain, that span would force all three ranks to be consecutive in every extension. Any incomparability linking it to another element supplies an extension with an extra interior element, a contradiction. Thus connectedness excludes the chain case. This is an actual extremal selection argument; choosing an arbitrary central triple would lose the boundary signal.

### 2.4 Stability at the rigid triple

Set `L'=L union {x}`, `U'=U union {z}`, and `m=|L'|`. With independent uniform extensions of these blocks, let `r` be the probability that `x` is last in the left block, and `s` the probability that `z` is first in the right block. Let `Z=1+r+s`.

There are exactly three allowed placements of `y`: the ordinary middle placement, an extra placement before `x` when `x` is last, and an extra placement after `z` when `z` is first. If both conditions occur there are three, not four, placements. Thus

```text
e(P)=e(L')e(U')Z,
p=Pr(y before x)=r/Z,   p'=Pr(z before y)=s/Z.        (3)
```

Writing `t=E[m-f_(L')(x)]`, `u=E[f_(U')(z)-1]`, the two remaining mean-rank contributions are `(1+s)t/Z` and `(1+r)u/Z`. These weighted expectations, their normalization, and all three rank identities were independently checked on 289 glued finite block models.

Since `t>=1-r`, `u>=1-s`, and the mean span is at most 2, one obtains

```text
r+s+2rs>=2.                                        (4)
```

The extremal algebra has `r=s=rho=(sqrt(5)-1)/2`. If the local balance is `a=C+epsilon<=0.2764`, the displayed estimates give

```text
|r-rho|,|s-rho| <=5sqrt(epsilon),
M<=21sqrt(epsilon),
a>=C+M^2/441.                                      (5)
```

All inequalities used in this stability estimate were checked, including the factors 11,50,20,21. The remaining excess rank contributions are nonnegative and sum to at most `11epsilon`. The conclusion concerns the **selected** triple carrying the global maximum, not every BFT triple indiscriminately.

Finally, connectedness supplies an element incomparable with `v_1`. Its positive inversion event and (1) give

```text
M>=d_1=sum_(j>1)Pr(v_j before v_1)>=q_2.
```

Together with (5), this yields a gap at least `q_2^2/441`, strictly larger than `eta_D=q_2^2/4096`, contradicting the near-BFT hypothesis. The endpoint probability and the attained maximum are exactly the finite ingredients.

Rationality alone would not do this: finite extension probabilities are rational whereas `C` is irrational, but strict inequality at each finite size does not imply a uniform positive gap across sizes.

### 2.5 Large width: internal conditioning and constants

The large-width proof uses the continuous uniform order-polytope variables `Z_x=(n+1)F_x`. It retains two different quantities: the mean gap `Delta(x,y)=|h(x)-h(y)|` and the standard deviation `d(x,y)` of `Z_x-Z_y`.

Under `delta(P)<=e^(-1)-epsilon`, the quoted log-concavity and window estimates give `d<=Delta/epsilon`, separated expected ranks, and a bound on width in terms of the largest expected-rank gap `G`. Shepp's quoted continuous correlation inequality yields the **squared** triangle inequality

```text
d(x,z)^2<=d(x,y)^2+d(y,z)^2,
```

and hence the improved global estimate `d(x,y)^2<=G Delta(x,y)/epsilon^2`. This distinction is essential to the following tail estimate; an ordinary triangle inequality would not supply the same bound.

The proof adds relations at expected-rank distance at least `ell=(288/epsilon^2)G(log G)^2` from a maximal gap. They are oriented by the original heights, so the resulting `Q` is consistent. The log-concave tail sum and rank separation show the total exceptional probability is at most `G^(-4)`. Conditioning on the added relations is exactly uniform measure on `O(Q)`, with the same normalization `n+1`. Cauchy–Schwarz bounds the perturbation of each cross-cut mean gap using the variance of the **difference**, not the potentially much larger variance of either coordinate.

The resulting canopy antichains are large by the imported Haqi–Kahn ideal inequality but lie in an original-height window of size `ell`. The window estimates force

```text
G <= (18432sqrt(3)/epsilon^3)(log G)^2,
```

which bounds `G` and then the width. The floor in the tail-sum lower index, all directions of conditioning, the use of original versus modified heights, and the final elementary inversion of this inequality pass the internal audit. Inverting `G<=A(log G)^2` uses the large explicit `A`; it is not a universally valid bound without that size condition.

**Dependency boundary:** this audit does not certify the cited general Shepp, window, Haqi, or bounded-width/large-range theorems merely because their statements work algebraically here. The finite controls below are checks, not substitutes for those source proofs. The unspecified effective constants behind the global epsilon numeral remain outside this audit.

## 3. Exact boundary escape in the Fibonacci family

The following elementary reconstruction makes the paper's mechanism quantitative without relying on the new theorem.

Let `F_N` be the poset on `1,...,N` with `i<j` whenever `j-i>=2`. Let Fibonacci numbers satisfy `F_0=0,F_1=1`.

### 3.1 Linear extensions are path matchings

An extension can reverse only adjacent labels. Two adjacent reversed pairs cannot share a vertex, since that would reverse a pair at distance 2. Conversely any collection of disjoint adjacent swaps gives a valid extension. Thus extensions correspond bijectively to matchings of the path on `N` vertices, and

```text
e(F_N)=F_(N+1).
```

If the edge `(j,j+1)` is selected, the remaining matchings live independently on the `j-1` vertices to its left and the `N-j-1` vertices to its right. Therefore its reversal probability is

```text
p_j=F_j F_(N-j)/F_(N+1),       1<=j<N,
p_0=p_N=0.                                         (6)
```

The expected displacement of label `j` is consequently

```text
d_j=p_j-p_(j-1).                                    (7)
```

For the symmetric restriction `N=2m+1`, relabel by `i=-m,...,m`. The Fibonacci determinant identity gives

```text
|d_i|=F_(2|i|)/F_(2m+2),
M_m=F_(2m)/F_(2m+2).                               (8)
```

The maximum occurs at the two endpoints. Also `F_j F_(N-j)<=F_(N-1)`, by the Fibonacci addition formula, so (6) shows that the endpoint adjacent pairs maximize the balance too. Thus for these odd-sized restrictions,

```text
delta(F_(2m+1))=M_m --> (3-sqrt(5))/2.
```

For every fixed labelled adjacent pair in the two-sided exhaustion, (6) instead tends to `C=(5-sqrt(5))/10`; all nonadjacent pairs are comparable. Hence

```text
lim_m max_(pairs in F_(2m+1)) delta(pair) = (3-sqrt(5))/2,
sup_(fixed pairs) lim_m delta(pair) = (5-sqrt(5))/10.  (9)
```

The difference is `1-2/sqrt(5)>0`. This is an exact failed interchange of operations.

### 3.2 Large displacement occupies a bounded boundary strip

If a vertex has distance `r=m-|i|` from the nearer endpoint, then

```text
|d_i|=F_(2m-2r)/F_(2m+2)<=2^(-r-1).                (10)
```

Indeed `F_(a+2)>=2F_a` for `a>=1`, and the central numerator `F_0` is zero. Iterating this inequality proves (10). Therefore any fixed threshold `tau>0` can be exceeded only within a number of endpoint layers depending on `tau`, independently of `m`.

For a fixed boundary distance `r`, Binet's formula gives the sharper limiting profile

```text
|d| --> phi^(-2r-2),    phi=(1+sqrt(5))/2.
```

By summing the even-index Fibonacci numbers, one also gets an exact total-mass identity:

```text
sum_(i=-m)^m |d_i|
 =2(F_(2m+1)-1)/F_(2m+2) --> sqrt(5)-1.             (11)
```

In particular the mean absolute displacement tends to zero, although its maximum tends to about 0.382. More generally, (10) bounds every unnormalized positive power sum independently of `m`, so every normalized finite positive moment tends to zero. The nonzero displacement has escaped to moving endpoints rather than disappeared uniformly.

This proves a concrete warning stronger than a visual analogy: average or fixed-coordinate control does not control an extremal witness unless one also controls where that witness can move. It does **not** assert that all finite-to-infinite passages fail.

## 4. What can and cannot transfer to Collatz

The root lane supplies an exact separate map from fixed-total slope-positive parity words to linear extensions of a width-two poset. This audit does not duplicate that construction. Its width stays at most two, so the large-width theorem above cannot be invoked on that object. Its incomparability range may grow with the horizon, so the small-range gap is not automatically uniform either.

The most promising transferable ingredient is the local-event floor (1). To use it after imposing an ordinary-integer height cutoff, one would need a faithful statement about the **conditioned** extension law. It is not enough that a comparison event has positive probability under the unconditioned uniform word law. Conditioning can erase it, and it can make the surviving order a chain. General arithmetic height cuts need not themselves be conjunctions of poset comparisons.

The root's cheapest explicit hostile is the two growing words with three odd letters and one even letter, `1110` and `1101`, whose least positive residue representatives modulo 16 are 7 and 11. The comparison has probability 1/2 under the two-word law, but an initial-height cutoff at 8 retains only the first word. This is a measure change, not a counterexample to a poset theorem.

| Source → target | Candidate map | What is preserved | What remains missing |
|---|---|---|---|
| Finite Fibonacci poset → infinite local limit | Keep fixed labels while expanding the interval | Each fixed pair probability | The maximizing boundary pair, lost in (9) |
| Bounded-range poset → local comparison events | Restrict to a finite set of labels and count extensions | Positive event floor (1) for the uniform full law | A similar floor after arithmetic conditioning |
| Finite parity-word family → width-two extension poset | Root lane's chain/cross-relation encoding | The finite uniform word law | Absolute source height and which word belongs to the prescribed source |
| Extremal poset displacement → Collatz obstruction | No descent-preserving map established | An instructive extremal-selection mechanism | Carry, source identity, measure, and a descent predicate |

A witness moving to later **time** is not itself a defect in a Collatz argument: the desired descent time for one fixed initial integer is allowed to depend on that integer. The invalid step would be changing the starting integer or the parity branch while passing to longer horizons. The endpoint analogy must keep this distinction.

The operational lesson is to seek an observable whose maximum can be selected on actual prefixes of one source, together with a lower floor that survives the relevant arithmetic conditioning. A density estimate or a two-adic limit alone supplies neither condition. This is a research requirement, not a claimed reduction of Collatz to the new poset theorem.

## 5. Exact controls, adversarial coverage, and stopping point

The standalone script enumerates every transitively closed naturally labelled poset on 2..6 vertices: respectively 2,7,40,357,4824, for 5230 total. This is an explicitly restricted labelled universe; it contains a representative of every isomorphism type at those sizes but is not a claim to enumerate all labelled posets. It enumerates 252265 linear extensions across that universe.

It then checks:

- 4376 maximal-displacement BFT-triple selections for connected incomparability graphs;
- 1089374 continuous squared-variance triangle inequalities, using exact order-statistic moment formulas;
- 91196 ideal-cut inequalities and 96426 antichain window inequalities;
- 289 independently enumerated rigid-triple block gluings, including all three height identities;
- 1230 rational feasible triples in the Case D optimization and 101 Appendix A boundary values;
- direct Fibonacci extension counts and moments through 13 vertices, and 63000 exact boundary-localization instances through `m=250`.

The near-BFT global hypothesis is not populated by small examples, so the script does not pretend to validate its difficult conditional theorem by a vacuous finite census. Instead it checks the unconditioned selection lemma, exact structural identities, algebraic feasible region, and the hostile limiting family independently.

The continuous difference variance used in the script is reconstructed directly from discrete extensions. If `g=|f(x)-f(y)|` and `mu=E[f(x)-f(y)]`, then

```text
Var(Z_x-Z_y)=((n+1)/(n+2)) E[g^2+g]-mu^2.
```

This follows by conditioning on the uniform order-statistic gaps inside each extension simplex. Thus the variance check is not a simulation of the cited correlation inequality.

Run normally and with `-O`:

```text
python 04-computation/experiments/crossroads_poset_20260926_barrier.py
python -O 04-computation/experiments/crossroads_poset_20260926_barrier.py
```

**Audit conclusion:** the internal finite-poset mechanism and its constants passed the stated scoped checks. The displayed global epsilon and the complete upstream proof graph remain outside this lane's certification. The durable new object is the exact boundary-escape profile (8)–(11), with a typed explanation of why it does not preserve a fixed Collatz source. No canon promotion, shared-document edit, or Collatz closure is made by this lane.
