# Why fixed-root doubling and the moving-root tower do not rearm these digits

This note concerns the nonzero synchronized outputs of the [first ternary splice](../ternary_deposit/README.md). It
does not preclude changing the source/profile; in particular it does
not conflict with the literal second digit with a different required
continuation. Write h=H, E,O for the two weighted root parities,
g=H(T−root), B=E+O, and e for the full Hamilton paths ending at the
root. The signed root-count identity is O−E=e. A root is called
normalized if e=g. The ready equation of this particular source family
is Δ=0 where
```
Δ=I+g=h−2O+3g,    I=h−2O+2g.
```
At normalized roots this is equivalently Δ=h−B+2g=h−2E+g.
The nonzero first digits have h=252R−15u+d u, E=126R,
O=126R+15u, e=g=15u, B=252R+15u, and Δ=d u for d=2,4.
Here u=H(M)>0 and R=P_M(2)≥2u; empty M can instead use the
literal values, which satisfy the same h/E/O/g formulas with R=u=1.

## Fixed-root doubling and opaque append obstruction

In a uniform near-sink block C, all old nonroot vertices beat C
and C beats the root. If h_C=H(C), q_C=P_C(2), the exact rooted
transfer is I'=h_C I, g'=h_C g, and therefore Δ'=h_C Δ.
For an ordinary sink block, all old vertices beat C and
h'=h_C h, g'=h_C g, E'=O'=q_C(E+O)/2. Consequently
```
Δ'−h_C Δ = −q_C E−(q_C−2h_C)O ≤ −2h_C E.
```
Both types of nonempty block preserve or increase E by retaining
old rooted paths with a complement of C. At the first nonzero
output, 0<Δ≤4u<2E. Therefore no near-sink block word can zero
it, and the first ordinary sink block sends it strictly negative,
from which either block keeps it strictly negative. Ordinal
prefixes scale all the coordinates, not restoring zero. The
fixed-root argument needs no finite bank or ordered-size cap.

These transfer identities are the literal odd/even block split of
the later nonsparse-reader packet. The antecedent scalar operations
are in Session195's `full-binary-return/PROOF.md` and companion
`binary-return-barrier/PROOF.md`. That theorem studies its own binary
condition e=g, B=h+g±1; it must not be quoted as a classification
for Δ=0. The preceding defect calculation is the adaptation.

## Turning the root after any binary sink word still fails

For a strong tournament, a root w has e_w=g_w iff the terminal SCC
of T−w entirely beats w. Indeed every Hamilton path meets the
ordinal SCCs in order, and each vertex of the terminal SCC can end
one. In particular a normalized root of a nontrivial strong
tournament is a strong articulation vertex.

On each literal first nonzero output, deleting any vertex except root
7 leaves a strong tournament; deleting the root leaves a terminal
module which beats it. This is checked on the finite quotients
(including the empty case) in the replay. Inflating the marked slot
by any nonempty tournament preserves strongness whenever the quotient
is strong and of order ≥2; if deleting one vertex empties the slot,
use the empty quotient. Thus there is no immediate change of root
to a different normalized port.

Now append an arbitrary word of ordinary sink and fixed-root
near-sink singletons to the first nonzero output. If it ends in
an ordinary sink, the only normalized root is the last sink,
whose defect is 2h>0. If it ends in a near-sink y, deleting any
candidate other than original root or y leaves y reachable from
all other vertices: all nonroots beat y and the root can use one
of its at least two remaining old outneighbors. Every such candidate
beats y, so cannot normalize. At y, the one-exit formulas give
e_y=h_prev−e_prev, g_y=h_prev, and B_y=B_prev+h_prev, where e_prev,
B_prev are read at the original root. If e_prev>0, the new root
does not normalize. If e_prev=0, the one-exit formula gives
H_y=h_prev+B_prev=B_y, so its normalized defect is 2h_prev>0.
The original root cannot zero its defect by the fixed-root proof.

## The THM-4865 moving-newest-root tower

In this tower the first newcomer x1 beats the old root x0; each
subsequent x_j beats only x_(j−1) among earlier vertices. The
inherited theorem gives
```
h_(i+1)=h_i−e_i+B_i,  e_(i+1)=h_i−e_i,  B_(i+1)=B_i+h_i,
```
where h_i is the count at height i and e_i,B_i read the newest root.
It gives the scalar tribonacci recurrence after height two.

At height k≥2, the only normalized roots are x_(k−1) and x_(k−2).
The newest root has e_k<h_(k−1)=g_k. On deleting any old nonroot,
all surviving old vertices reach the tail (the original root reaches
x2), and every old nonroot beats the tail. Deleting x_i for 0≤i<k
leaves the prefix before the tail in the SCC order. If the tail
has at least three vertices, it is strong: adjacent backward
steps and forward jumps of size two generate reachability. In that
case not every tail vertex beats x_i. With one or two tail vertices,
the terminal SCC is exactly the singleton x_(i+1), which beats x_i.
At height one, only the original root normalizes, by the binary
word argument above.

For i≥1, the normalized defect of x_i as penultimate root is
2h_(i−1). Indeed B_i=h_i+e_(i−1) and e_i=h_(i−1)−e_(i−1)
follow directly from the one-exit updates. After one stationary
near-sink step the defect is
h_i−B_i−e_i+3h_(i−1)=2h_(i−1).

For any root of a base with profile (h,g,B), append two moving
births while ultimately reading that same initial root. The resulting
count is h+2B from THM-4865; deletion count is still g. The new
unsigned rooted path sum is 4B−g: an old rooted path can end with
neither new vertex, with the first alone (except the old singleton),
with the second alone, or with the suffix second→first.
Its defect is therefore h−2B+3g. At x_i for i≥1 this is
`2h_(i−1)−e_(i−1)−B_(i−1)=D_(i+1)`; at original x0 it is −252R+d u, nonzero.

Finally put D_i=B_i−3e_i. Its first terms are
```
D0=B−3e, D1=B−2h+3e, D2=2h−e−B, D3=B−e, D4=B+e.
```
The linear updates above have tribonacci characteristic polynomial,
so D_i follows that recurrence as well (or subtract the recurrences
for B_i,e_i). For the displayed nonzero first digits, D0>0, D1<0,
D2>0, D3>0 and D4>0 using R≥u. Thus all subsequent D_i remain
positive, and none is zero. Every tower height and every root
fails at least one of normalization and Δ=0.

## Archive implications, in order of relevance

- `01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci-h-k-equals-h-k-1-plus-h-k-2-plus-h-k-3-for-every-rooted-base-the-lineage-c2-law-at-any-vertex-and-the-three-root-response-numerations-attain-every-odd-number-in-69-702189.md` and
  `07-reflections/the-root-response-one-object-behind-the-c2-grammar-the-near-join-and-the-tribonacci-tower-klein-20261104.md`:
  keep the full B readout and e, not only H. Their height-two value
  h+2B is already a 3h+2σ mechanism, but the ready equation here
  asks B−h=2g on a normalized root. That is not identical to
  the old σ=1,2 readiness. A proposed handoff should state which
  normalization and digit equation it uses.
- `session195/full-binary-return/PROOF.md`: the earlier all-word
  no-go even for its named one-arc gadget switches targets
  e=g,B=h+g±1; adapting it to Δ=0 would need new inequalities.
- `05-knowledge/results/session189_response_frontier_20261104/returns/PROOF.md`:
  the four-vertex joint-entry adapter carries a named J=1 port
  through a physical three-birth dial and breaks mod-3 confinement,
  but proves exact dyadic spacing prevents +2 among wrappers of a
  fixed component. This is a distinct witness for changing source.
- `01-canon/theorems/THM-4882-signed-one-apex-responses-slope-and-parity-classification-session198-clock-merged-positive-count-barrier.md`:
  the signed aU−bH clock can span all integers of a parity in
  specified ratios, but an actual H insertion cannot run its
  negative shift indefinitely. Pursue it only as a new readout
  with a physical continuation/cancellation, not a scalar substitute.
- `session197/minimum-positive-clock/PROOF.md`: transitive root
  cloning already realizes 72·3^m+1 from a different profile;
  its run-matching condition gives a genuinely different growth
  axis, but not consecutive fresh digits or a rearm of this source.
- `05-knowledge/results/ternary_ports_golden_defect_20261018j.md`:
  neighboring Petersen/Tait material shows how full labelled tensor
  support can block a gluing even when both scalar sides are present.
  This motivates preserving the named marked port and adding a
  continuation sidecar, not a claim about tournament counts.
