# Signed fixed-half-path extremal sequence: exact extension through seven pairs

## Status and scope

Let F_m be the maximum number of Hamiltonian paths fixed by a specified standard anti-involution, over all signed controller tournaments on m vertex-pairs, with the fold-compatible sign normalized to plus. This is not an identification with the older undefined repo symbol M_k, not the ordinary Hamiltonian-path maximum, and not automatically the maximum number of blue tiles after dividing by centralizers.

The candidate was tested rather than assumed. The complete finite computation gives

m: 1, 2, 3, 4, 5, 6, 7
F_m: 1, 3, 9, 39, 195, 1185, 8295.

The first five values come from the preceding exact session29 controller work; this directory extends the calculation at six and seven pairs. The two new maxima agree with the candidate sequence C_m, but no formula for F_m at arbitrary m, no universal bound F_m <= C_m, and no asymptotic for F_m/m! has been proved.

## Search protocol, including negative initial evidence

The initial deterministic sample used std::mt19937 seed 290063. It included q=0 (transitive), q=2 (initial directed triangle followed by an ordered tail), disjoint consecutive directed triangles with all cross-arcs forward, a cyclic-style controller and its converse, then distinct uniformly sampled labelled controller masks until reaching the requested total.

- 200 sampled m=6 controllers: best 1177, below C_6=1185.
- 30 sampled m=7 controllers: best 8253, below C_7=8295.

These samples did not justify rejection. Ordinary isomorphism reduction is cheap at these orders: scan labelled bit masks, and for each unseen mask mark its orbit under all m! permutations. This leaves exactly 56 controller representatives at m=6 and 456 at m=7. All representatives were then optimized exactly. This is not enumeration of all full tournaments on 12 or 14 vertices.

For each controller, fix the first outgoing label at each nonsink to plus. Enumerate the other labels. This visits every outgoing bipartition modulo a simultaneous flip of each star. For each permutation of pair vertices, multiply the factors 1-b_vu*b_vw at every internal vertex controlling both neighbors. Summing these weights gives A. Sum weight*b_vu into B_v when final vertex v controls its predecessor u. Independent star flips then give the exact optimum A+sum_v |B_v|. The emitted final negative-label mask incorporates all maximizing star flips; the printed B vector is before those flips.

Total bipartitions examined in the complete reduced runs:

- m=6: 34,816 across 56 representatives.
- m=7: 8,388,608 across 456 representatives.

Each order has exactly one maximizing controller isomorphism class. Representative files contain every class's exact local optimum, including nonmaximizers. This establishes the finite upper bounds computationally under the already proved phase formula; the new witness values are independently checked without that formula by signed endpoint/pair-subset DP.

## Literal witnesses and controller structure

Edges are indexed lexicographically by unordered pairs i<j, from zero. Bit q_e=0 means i controls j, bit 1 means j controls i. Negative-label bit e=1 means b=-1; b is the required sign when the controller occurs after the controlled vertex in a half-path ordering.

m=6: q=140, negative-label mask=26896, F=1185.
A=960, B=(45,45,45,45,45,0).
Controller is C5{1,2} -> one terminal vertex.

m=7: q=268, negative-label mask=39375, F=8295.
A=7800, pre-phase B=(-78,-78,-78,-78,-78,105,0).
Controller is C5{1,2} -> a transitive two-vertex tail.

The same witnesses lift to odd order by adjoining a center, with i+ -> center -> i-. In even order i+ -> i-. They therefore give fixed-path values at orders12/13 and14/15, respectively. Do not infer blue tile values without the centralizer calculation.

`witnesses.json` supplies literal even-order tournament outmasks and sigma, with sigma(i)=i XOR 1. `verify.py` builds the directed graph independently and counts paths with a DP indexed only by visited pair subset and signed endpoint, finally summing plus endpoints. It verifies 1185 and8295.

## Existing sequence, exact generating function, and the missing theorem

The matching sequence is OEIS A130905:
https://oeis.org/A130905

Its established EGF is exp(x^2/2)/(1-x). The entry also gives the interpretation as permutations with each transposition allowed either of two colors. These properties belong to C_m independently of any tournament conjecture.

A directly useful alternate interpretation is a partial matching plus a linear order of all unmatched vertices. If there are j matched pairs, choose their2j vertices, match them, and linearly order the remainder. Thus

C_m = sum_j binom(m,2j) (2j-1)!! (m-2j)! = m! sum_j 1/(2^j j!),

where 0<=j<=floor(m/2). This proves the candidate EGF without borrowing a tournament claim. Equivalently C_0=C_1=1 and C_m=m C_(m-1) for odd m, while C_m=m C_(m-1)+(m-1)!! for even m.

The numerical asymptotic C_m/m! -> sqrt(e) is rigorous for this candidate alone. For r=floor(m/2)+1, its deficit satisfies

1/(2^r r!) <= sqrt(e)-C_m/m! <= [1/(2^r r!)]/[1-1/(2(r+1))].

This follows by bounding the positive exponential-series tail by its first term and a geometric series.

The missing structural statement is an injection from signed-controller fixed half-path configurations into partial matchings with ordered unmatched vertices, or another proof of the same upper bound, together with an all-order attaining family. The local source/wedge expansion supplies natural pair interactions, but a source's two incident edges are not already a disjoint matching. Replacing those wedges by arbitrary independent matched pairs without proving disjointness and injectivity would be an invalid shortcut.

The optimizer sequence already changes controller type: transitive at m4, C3->T2 at m5, C5->T1 at m6, C5->T2 at m7. This suggests testing larger odd regular cores with one or two terminal vertices, but it does not specify all future optimizers. The transitive controller alone is strictly suboptimal from m5 onward in the computed range.

A repository text search for A130905 and the adjacent1185/8295 values found no prior Markdown hit in the available checkout. This is only a scoped text-search result, not a claim of first discovery. The generating function itself is established mathematics; the putative extremal tournament identification is the new unproved question.

## Reproduction

Compile sample.cpp with g++ -O3 -std=c++17.
Run `./sample 6 200` or `./sample 7 30` for the initial seeded samples.
Run `./sample 6` and `./sample 7` for all ordinary-controller classes.
Run `python verify.py` after both class files exist.

No Git or repository mutation was performed by this worker.

## A bounded analytic bridge, not a Gaussian tournament model

Let X have the exponential distribution of mean1 and let Z be an independent standard real Gaussian. The moment generating function of X+Z is exp(z^2/2)/(1-z), for z<1. Consequently C_m=E[(X+Z)^m]. Expanding this moment gives the same partial-matching formula: Gaussian moments pair selected positions; the exponential moment orders the remaining positions. This makes a Wick-style expansion a precise candidate proof target. It does not assert that a tournament is a Gaussian field or identify this variable with any actual rank88 or Lorentzian Hessian operator.

At m3 the transitive controller has exact polynomial
f=6+x+y+z-2xy,
where x,y are the two outgoing labels of its source and z is the remaining label. Since x,y,z are each +/-1, z<=1 and x+y-2xy<=2, proving f<=9 directly. Equality uses opposite x,y and z=1. The unsigned relaxation11 is therefore genuinely too large. A prospective general matching bound must retain these label incompatibilities instead of dropping the negative interaction terms. A deletion/insertion recurrence with a perfect-matching residual is another concrete target, but no such tournament recurrence is proved here.
