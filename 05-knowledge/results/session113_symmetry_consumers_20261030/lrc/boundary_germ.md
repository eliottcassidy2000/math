# LRC: an exact boundary-sensitive mass germ and finite relaxation gate

**Session 113, 2026-10-30. Status: elementary proofs plus finite exact verification. LRC(14) remains open.** This is a source-faithful reduction, not a new universal lower bound or a claim of literature priority.

## Result first

Let V be a nonempty finite set of positive integer speeds, M=max V, and δ=p/d in lowest terms with 0<δ<1/2. Write

S(r)={t in R/Z : ||vt|| >= r for every v in V}.

Then there is an explicit endpoint-owner coefficient B(V,δ)>=0 such that

**μS(δ−ε)=μS(δ)+ε B(V,δ) for every 0<ε<1/(2dM).**  (Restrict also ε<δ; the applications below satisfy this automatically.)

Consequently the following are equivalent:

1. S(δ) is nonempty.
2. μS(δ−1/(4dM)) is positive.
3. μS(δ−1/(4dM)) >= 1/(2dM²).

For thirteen speeds and δ=1/14 this gives the exact finite target

**μS(1/14−1/(56M)) >= 1/(28M²).**

A proof of this inequality for every source row would prove LRC(14); we do not supply it. Unlike an unguarded fixed relaxation or positive mass at the original radius, this target neither invents near-miss witnesses nor discards isolated closed witnesses. It replaces an existing all-ε compactness bridge by one explicitly guarded finite-radius inequality. It does not imply a uniform small clock, a universal original-radius mass floor, or entry into an inherited dyadic/decoder gate.

## Provenance and exact source contract

Read against the current LRC frontier and:

- `05-knowledge/results/lonely_runner_conditional_survival_20261011e.md`, §§5–6: tight AP zero mass; complete endpoint clock d lcm(V); all-ε limiting bridge.
- THM-3731: owner/phase/word data cannot be replaced by valuation symmetry alone.
- THM-4453: conditional component/mass entry is not universal entry.

The underused sidecar is the **signed active endpoint-owner set**, including every tie. The operation is physical radius relaxation on the same named row and same time circle. No independent phase, synthetic density, reordered owner, or auxiliary clock is substituted. The two mandatory large-clock hostile rows are checked below. No historical 'closed LRC' claim is inherited.

## 1. The owner coefficient

As δ>0, S(δ) has a unique presentation as disjoint maximal closed intervals [a_j,b_j] inside (0,1), allowing a_j=b_j. Define

P(a)={v in V : va mod 1=δ},
N(b)={v in V : vb mod 1=1−δ}.

At every left endpoint P(a) is nonempty; at every right endpoint N(b) is nonempty. Otherwise all active constraints would permit continuation to that side. Put

B(V,δ)=Σ_j [1/max P(a_j) + 1/max N(b_j)].

The MAX is essential. Under relaxation by ε a positive-slope constraint with speed v moves its left wall from a to a−ε/v. Among tied walls the most restrictive is the one with largest v. Similarly the right wall moves from b to b+ε/max N(b). For a singleton both contributions are present. This is not a choice of one arbitrary owner; ties remain in the certificate, and the extrema are derived from them.

For sufficiently small ε the old intervals therefore expand exactly by these amounts. The next section gives a uniform source-height lower bound on 'sufficiently small'. If S(δ) is empty the coefficient is zero.

A small tie hostile makes this load-bearing: at δ=1/3 both rows {1,2} and {1,2,4} have exactly the same two singleton witnesses. Their coefficients are respectively 3 and 3/2. At t=1/3 the latter positive owner set is {1,4}; retaining only owner 1 gives the wrong expansion.

## 2. Event-free radius: proof

The safe intervals for one speed v have walls

(k+r)/v and (k+1−r)/v, 0<=k<v.

Their intersection is governed by the ordering of finitely many affine wall functions. An ordering can change only where two distinct functions agree. An isolated agreement has radius

r=(v l−w k)/(w s−v t), with s,t in {−1,+1},

and denominator, after reduction, at most v+w<=2M. Identical functions cause permanent ties, not new events. Fixed circle-cut walls 0 and 1 introduce only integer-radius events. Birth, death, merger, endpoint-owner change and singleton creation are all included in these ordering events.

If r0<δ=p/d is any such event, rational separation gives

δ−r0 >= 1/(d denominator(r0)) >= 1/(2dM).

Thus the entire open radius interval (δ−1/(2dM),δ) contains no event (above radius zero). Its combinatorial type is constant. Continuity from δ identifies the expanding endpoint walls exactly as in §1; empty sources cannot acquire a component before an event. Summing interval lengths proves the displayed affine identity.

### Scale-invariant refinement

The stronger event height

H(V)=max({2} union {(v+w)/gcd(v,w): v,w in V, v!=w})

can replace 2M in the event-free bound. In the wall-intersection formula the numerator is divisible by gcd(v,w), so the reduced denominator divides (v+w)/gcd(v,w) or |v−w|/gcd(v,w); the same-speed event is r=1/2. Thus the exact affine identity holds throughout 0<ε<min(δ,1/(d H(V))). Unlike M, H is unchanged by common dilation. The checker uses this stronger interval. The simpler M-based gate below remains a convenient universal corollary. More sharply, with M0=max(V)/gcd(V), nonemptiness is equivalent to the scale-invariant finite gate

μS(δ−1/(2d H(V))) >= 1/(d H(V) M0).

This follows from the same proof after dividing all speeds by their gcd. No independence of phases is used.

## 3. The single finite inequality: proof

Choose ε=1/(4dM). If S(δ) is nonempty, choose a witness t. Each function ||v·|| is v-Lipschitz, so an arc of radius ε/M about t survives at δ−ε. Its length is 2ε/M=1/(2dM²). Conversely, if S(δ) is empty then both its mass and B are zero, and the affine identity gives zero relaxed mass. This proves all three equivalences.

An unguarded relaxation is false: for V={1,2}, δ=17/50>1/3 the source has no witness. Relaxing by 1/100 crosses its maximum-clearance event at 1/3 and creates positive mass. The guarded ε=1/400 does not.

## 4. What symmetry does and does not preserve

Reflection t -> 1−t swaps P and N and preserves each component's weight. If V=gW, physical pullback under t -> gt creates g copies of every component, while every owner speed is multiplied by g. Therefore

μS_V(δ)=μS_W(δ), and B(V,δ)=B(W,δ).

The unweighted number of isolated witnesses multiplies by g; the boundary-weighted carrier is invariant. For the finite inequality it is preferable to divide by gcd(V) first and use the smaller primitive M, retaining the physical lift of any resulting witness.

**Closest hostile: tight AP.** For V={1,...,n−1}, δ=1/n, the source survivor is exactly {a/n:gcd(a,n)=1}; its original mass is zero. At a/n, the positive owner is the inverse r of a mod n in {1,...,n−1}, and the negative owner is n−r. Hence

B=Σ_(r unit mod n) [1/r+1/(n−r)] = 2 Σ_(r unit mod n) 1/r.

At n=14, the six witnesses form one orbit under multiplication by units modulo 14. Nevertheless their physical boundary weights are

| numerator a | positive, negative owners | weight |
|---|---|---|
| 1,13 | (1,13), (13,1) | 14/13 |
| 3,11 | (5,9), (9,5) | 14/45 |
| 5,9 | (3,11), (11,3) | 14/33 |

Thus B=23324/6435. Multiplying the first representative's weight by orbit size gives the **wrong** value 84/13. The unit action preserves the closed grid predicate but not the off-grid physical derivative. Modular speed congruences are not equality of speeds under perturbation. This pinpoints exactly why a transitive symmetry orbit still needs an owner/phase sidecar when imported into mass or entry estimates.

## 5. Verification and remaining target

`check_boundary_germ.py` reconstructs closed survivors by rational interval intersection, derives owner sets from actual phases, and independently recomputes relaxed survivors at three radii. Checks remain active with `python -O`.

- 1,275 rows/radii: all speed subsets of {1,...,9} of sizes 1–4, at five rational radii.
- Tight APs n=3,...,20, including zero mass and exact owner-weighted coefficient.
- AP14 dilations by 2,3,7 and reflection ownership.
- V=26·{1,...,12} union {339}: mass 755507/25913160, B=45067831/6786780, 368 closed components.
- V={1,...,12,5460}: mass 301/10296, B=43103/25740, 192 closed components.
- Explicit unguarded-relaxation false positive.

Normal and optimized execution both pass **18,542 checks** with byte-identical output. Reproduce with `python check_boundary_germ.py` and `python -O check_boundary_germ.py`. These are audits of the proved formula, not evidence for an unproved universal bound.

**Next structural obligation:** bound the guarded relaxed survivor mass from below on the actual source row, using the marked owner carrier. A modular quotient is useful only if it transports the reciprocal-owner weights (or a valid lower bound for their total), not merely witness-orbit cardinality. This preserves tight-boundary information absent from Haar averaging but does not itself resolve the physical entry/owner-transport bottleneck.
