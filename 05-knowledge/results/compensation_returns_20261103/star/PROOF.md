# A complete fixed-star compensation census

2026-11-03. **PROVED one-apex formula and guards; FINITE-EXACT complete
named-star census.** Neither the 1379 debt state nor the 1361 debt state
can regain its requested count1409 or1427 by changing only edges incident
with inner mark0 while retaining strongness and outerroot99 end defect0.
The complete search has4096 words,1024 strong words, and512 words satisfying
both guards. It produces323 distinct guarded counts. This packet does not
claim an arbitrary-order or arbitrary-deleted-graph impossibility.

The input is the first two actual triangle repairs in the inherited
[entry-aware receipt](../../c4_port_returns_20261103/entry_selector/receipt.json).
The closest proved mechanism is the retained-mark reader and its terminal
entry law. The hostile is the C2 child's consumption of q=1 into q=3. The
corrected near miss is to infer a neutral repair from a restoration of q
alone: the actual strong triangle reversals cost30 and66 Hamilton paths.
The least-used sidecar is the distinction between e=E_inner0(1), which
remains1, and q=H(T-inner0), which changes. The concept board is exact named
star / deleted order / root guard / compensation count / next entry / e-versus-q.

## 1. The actual common carrier and exact word formula

Both debt graphs have the same physical names and the same non-inner0
edges. Deleting inner0 leaves the transitive order

    99,10,11,2,1,3,4,12,200,401,400,201.

Let a be inner0, let d0,...,d_(m-1) be this order, and put x_i=1 iff a→d_i.
For any m>=1, every Hamilton path is an increasing subsequence of D,
then a, then the increasing complementary subsequence. Therefore

    H(x)=x_0+1-x_(m-1)
         +sum_(i<j) x_i(1-x_j) 2^(j-i-1)
         +sum_(i=0)^(m-2) (1-x_i)x_(i+1).                 (1)

The first two terms have the apex first or last. For a nonempty path on
each side, if the after-apex minimum i is below the before-apex maximum j,
the intermediate j-i-1 vertices independently choose a side. If the
before-apex maximum is below the after-apex minimum, they must be adjacent.
These cases are disjoint and exhaustive; all coefficients count actual
paths. A word determines every actual row after the deleted order and
names have been fixed. There is no isomorphism substitution.

For m>=2, strongness holds iff x_0=1 and x_(m-1)=0. Necessity follows from
the otherwise global source d0 or sink d_(m-1); sufficiency follows from
the directed Hamilton cycle a,d0,...,d_(m-1),a. In every word q=H(D)=1.

If strongness holds, outerroot r=d0 has sole in-neighbor a. Concatenating
an ordered two-path cover whose first path ends at r gives a Hamilton
path, except when the second path starts at a. In that exceptional case
the first path must be the singleton r. Hence

    E_r(2)-H(T)=S_a(T-r,1)=x_1.                          (2)

The last equality uses the unique order on D-r. Thus the exact typed
guard is x_0=1,x_1=0,x_(m-1)=0, leaving2^(m-3)=512 words at m=12.
Here E_r(2) uses an ordered first slot ending at r; empty second slots
are allowed. The start-entry defect is S_v(2)-H in the dual convention.

## 2. Complete census and nearest alternatives

`census.json` lists all4096 masks, their words, exact H, strongness, guarded
root defect, and Hamming distance from each actual debt star. Every pair
outside the star is fixed. `typed_profiles.json` independently recomputes
all512 guarded graphs by literal subset path DP, retaining their actual
rows, names, E1,S1,E2,S2, full ordered endpoint-joint matrix J, and all
start/end defects. A separate reachability computation checks strongness
for every one of the4096 graphs.

| Actual debt | Desired H | Closest below | Closest above | Multiplicities below/above |
|---:|---:|---:|---:|---:|
| 1379 | 1409 | 1407 | 1417 | 2/2 |
| 1361 | 1427 | 1425 | 1433 | 4/2 |

There is no exact compensation witness. The two nearest-below1407 words
each require five star flips from the actual1379 state. For1361, two of
the nearest1425 states require just one flip, respectively {0,12} and
{0,1}; the other two require four. Every nearest lower/upper state lacks
start defect1 or2 at every vertex. These alternatives do not solve the
requested compensation and are not replacements for an equal-count source.

`nearest.json` saves all ten lower/upper witnesses, their literal star
cuts and flip sets, full joint profiles, and fresh named golden/F49/F343
observers. These observers preserve their original definitions; no graph
automorphism or native transport is inferred from a Frobenius observation.

The full guarded family does contain entries:256 states have defect1 at
physical201. One of these, mask5 with H2559, also has defect2 at inner0.
Exactly one guarded state has an end defect1/2 anywhere: mask1,H2049,
at physical10. Thus the no-compensation conclusion is a count obstruction,
not the assertion that all nearby q1 graphs have lost every positive entry.

## 3. The source resource that survives

For the inner mark0, write the full seed carrier as

    (h,e,S2,F,p_deleted,q,P2,S3).

Here e is the number of Hamilton paths ending at0; F=E_0(2);
p_deleted is the ordered two-path cover count after deleting0;
q is its Hamilton count; S3 is the ordered three-slot start count at0.
The fresh profiles are

| Graph | Full eight-field carrier |
|---|---|
| Original child1409 | (1409,1,3633,2373,6144,3,279796,354169) |
| Repaired debt1379 | (1379,1,2423,2327,4096,1,281044,290659) |
| Original child1427 | (1427,1,3537,2469,6144,3,281716,325657) |
| Repaired debt1361 | (1361,1,2359,2359,4096,1,274900,267331) |

Thus e=1 survives the C2 step and the triangle repair, independently of
whether q is3 or1. The inherited lone-reader-zero terminal law can use e=1
for a new construction, but that construction changes the graph and its H;
it is not a neutral restoration established here. `controls.json` retains
all five actual stages for both families, all their ports and full joints,
and fourteen large weak/guard-destroying controls.

## 4. Scope, controls, and reproduction

The complete universe is all star masks on the two fixed named thirteen-
vertex debt states. They have the same non-star graph, so a single mask
census covers both, with distinct actual-source distances. The exact
formula also makes the scalar count set independent of the names assigned
to this fixed-length transitive order. Relabelling alone does not preserve
the full named J matrix or observer roles.

The independent literal profiler checks all512 typed graphs,254 small
stars (every mask through total order8, including weak cases), fourteen
large hostile masks, and ten retained actual stages:790 full profiles.
The fresh eight-field helper independently overlaps all existing H,e,S2,F
and P2 coordinates, and computes the deletion and three-slot coordinates.
No Python assert is used. Counts fit signed64-bit at the declared order13;
even ordered three-slot covers are bounded by13!*choose(15,2).

Reproduce from the repository root:

    python 05-knowledge/results/compensation_returns_20261103/star/reproduce.py

The normal and `-O` runs must produce byte-identical artifacts and output.
`manifest.json` pins the proof, scripts, artifacts, and inherited input
dependencies. No old receipt is modified. No all-order two-branch return,
native sequence, or {7,21}-only theorem is claimed by this packet.
