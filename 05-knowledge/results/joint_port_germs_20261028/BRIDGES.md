# What the connections preserve

**CITED inherited mechanisms / PROVED elementary comparisons / OPEN applications.**
This note compares actual maps and operations. Equal numerical values alone
are not evidence of a common classification.

## 1. The Borwein ladder measures the dimension of first contact

Write a_j=1/(2j+1), O_m=1+sum_(j=1..m) a_j. The corrected
[THM-4837 Borwein account](../../../01-canon/theorems/THM-4837-the-owners-borwein-ladder-first-nonvanishing-of-the-odd-sinc-cosine-transform-at-2-3-and-2pi-closed-form-first-failure.md)
identifies the first crossings of 2,3,2pi at

| Threshold xi | First m | Last denominator 2m+1 |
|---|---:|---:|
| 2 | 7 | 15 |
| 3 | 56 | 113 |
| 2pi | 40249 | 80499 |

These indices and the first-bite formula are published results, not new
discoveries here. See [Basel--Baillie](https://arxiv.org/html/1510.03200);
the repo's [directed-interval receipt](../minimal_blue_selectors_20261028/bridge-receipt.json)
rechecks the crossings. MISTAKE-692 in the THM-4837 lineage corrects credit
and normalization.

There is a useful geometric comparison. Let U_j be independent uniforms on
[-a_j,a_j], S=sum U_j, xi>=1, m>=1, and epsilon=O_m-xi. In the first-contact
range 0<=epsilon<=2 min a_j, the event S>=xi-1 is a corner simplex:

    x_j=a_j-U_j >= 0,       sum x_j <= epsilon.

The other bound S<=xi+1 is automatic in this range. The simplex volume is
epsilon^m/m!, so the known cosine transform is

    F_m(xi)=pi Pr(xi-1<=S<=xi+1)
           =pi epsilon^m / (m! product_(j=1..m) 2a_j).

The orders of first nonzero contact are therefore 7,56,40249. The
extraordinarily small deficits come from both the slack and this exponent.

By contrast, [L1-L2](LRC.md) keep a one-dimensional time circle. Every
nonempty isolated LRC witness opens linearly when the common threshold
decreases. Its first contact order is one, irrespective of the number of
runners. The coefficient records intersecting endpoint constraints.

**Elementary contact-order rule.** If positive quantities have expansions
A(epsilon)=a epsilon^r+o(epsilon^r) and B(epsilon)=b epsilon^s+o(epsilon^s),
with a,b>0, their product has order r+s and their sum has order min(r,s).
This explains independent product corners and disjoint positive
contributions. Overlapping unions need inclusion/exclusion; signed sums
can cancel the leading term. Reflection traces fall into this last warning.

The useful common object is a positive measure under an explicit
perturbation, with its order and coefficient. The dimension being integrated
must remain part of the interpretation. Borwein's convolution does not
identify 7 independent box variables with seven LRC constraints on one time.

## 2. Signed even graphs, trijections, and operational completeness

Here "even graph" means the sign-even notion in
[Royle et al.](https://arxiv.org/abs/2204.01947), not an Eulerian graph with
all vertex degrees even. The Walsh transform sends tournament orientation
functions to graph supports with signs from their stabilizers. Taking
invariants gives equinumerosity. It supplies a linear correspondence and
counts, not automatically a literal equivariant bijection of the objects.

The inherited
[signed operation spectrum, S1-S2](../odd_sign_spectra_trijections_20261028/THEORY.md)
retains more: degree, free-complement trace and the flip operator spectrum.
The present [J5](PORTS.md) shows its next boundary. For the rigid seven-pair,
every marked r-vertex cross spectrum is (1+z)^(7r), whereas an actual
selected one-vertex extension separates them. A trace has forgotten which
old arc each allowed operation interrogates.

For an oriented carrier with |E|=|V|-1, the complement has the right number
of pairs for a tournament of one smaller order. That dimension equality is
only the starting point. Audited
[THM-4835, Euler-characteristic criterion](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md)
proves that a pair-coordinate trijection exists exactly when every
subgroup-fixed graph has Euler characteristic one, using the FULL carrier
automorphism group. The underlying graph can be disconnected. Its
cycle rank then equals its number of components minus one.

The explicit end-edge map for trees and the root-swap map for functional
orientations satisfy this stronger criterion. Bigons exhibit a weaker
class correspondence without the same coordinate spectrum. An odd-order
Euler defect is not by itself a proof of failure of weak trijectivity:
the theorem's corrected seven-vertex examples must remain in view.

Thus the levels are different: dimension/count; signed character; actual
coordinate transport; and closure under the selected operations. J2-J3
prove that the last level can demand every arc bit even when the earlier
levels coincide.

## 3. Three kinds of oddness, with different mechanisms

1. A tournament automorphism group has odd order, because an involution
   would reverse the arc of one of its two-cycles. This holds at even as
   well as odd vertex orders.
2. The aligned-pair model of odd-order self-converse tournaments has the
   free Klein action in
   [P2-P3 of the inherited signed packet](../odd_sign_spectra_trijections_20261028/THEORY.md).
   It explains four equal parity sectors and divisibility by four from
   order five. Eplett's class correspondence does not retain a chosen
   anti-involution's centralizer and its port transports for free.
3. Reflection on the LRC time circle fixes only 0 and 1/2. Since 0 is
   forbidden at a positive threshold, its fixed witness channel sees the
   all-odd speed sector. Mixed-parity witness sets may be nonempty while
   reflection traces vanish. Positive growth in L1 detects those sets.

There is a new concrete reminder of the distinction. J1's ordered register
can have a terminal C3 kernel. It recovers the pointed core class, but no
equivariant choice of one vertex of a freely rotating triangle exists.
If a later port names that original sink, the triangle rotation is killed.
The exact exception is

    Aut(register)=Aut(pointed core) x C3

when the first port is the universal sink and no later port repeats it;
otherwise the kernel is trivial. Counting classes, reconstructing a core,
and preserving its automorphisms are three separate requirements.

## 4. A common test algebra, and the limits of the analogy

For a finite state space, fix both legal future operations and readouts.
Identify two states only when every legal future readout agrees. The
resulting behavioral quotient is the least sufficient deterministic state
space. This is the standard observation-algebra principle already used in
the repository, not a newly invented universal classification theorem.

| Source and operation | Target/map | Preserved test | Lost information and needed sidecar | Decisive check |
|---|---|---|---|---|
| Hamilton presentations; movable dents/C3 | Recognizable ordered ports | Every selected old arc | None at the discrete minimum; module/triangle kernels remain relevant for groups | J1-J3; seven-pair and81/27 |
| LRC sets; common threshold decrease | Positive germ (L,D) | Weak witness existence | Absolute phase is lost | L5's speed7 test |
| Tight LRC sets; adding a speed | Phase/owner records, reflection-folded | Future witnesses and uniform germ | Full directional rates need all owner sets | L6,9240 additions |
| Borwein box sum; moving support | Corner mass | First nonvanishing and exponent | Individual box coordinates and later dependent constraints | Known simplex formula |
| Marked graph completions; pair transport | Signed-coordinate cube | All coordinate flips only when transport is supplied | Weak orbit counts do not supply coordinates | Bigon and Euler-fixed-subgraph tests |
| Arithmetic square-class selections | Binary parity-check kernel | Square-total predicate | Valuation size and integer lifts are lost | Fano/regularity hostiles in the inherited arithmetic packets |

The arithmetic row is inherited from
[relation gluing and integer lifts](../relation_gluing_7_15_20261028/README.md).
The [signed-incidence comparison](../consecutive_flags_order10_20261028/README.md)
of graceful and square-sum edge equations likewise preserves an equation
with its incidence signs, not either conjecture's global existence claim.

Concurrent work provides a useful further test:
[Session91's golden observer and guarded flip operators](../session91_packing_golden_carriers_20261029/README.md)
shows that an arithmetic ring value can lose a named-toggle state, and that
a guard leaves a positive boundary operator. Its two-register compatibility
condition is another actual joint-data requirement. Our port register is a
combinatorial realization of the same need for an operation-complete state;
our LRC phase hostile shows that even an entire local response function can
still lose the next operation.

## 5. What changed on the concept board

- Spectral action: still useful for symmetry and parity, insufficient for
  selected response on rigid cores.
- Persistent selectors: arbitrary mobility forces full information; a
  restricted selector menu is now the meaningful place to seek smaller K.
- Strong modules: provide intrinsic port names using only C3, with known
  internal group factors.
- Positive germs: repair parity cancellation and strict/weak measure loss;
  phase must be restored for continuation.
- Contact order: explains tiny Borwein tails and linear LRC boundaries,
  while blocking an unsupported identification of their dimensions.

Next precise targets are minimal carriers for fixed successor selectors,
bounded substitution depth, and phase quotients for a restricted future-speed
alphabet. Proving LRC(14) still requires a uniform positive inequality for
all thirteen-speed rows. The height-adapted Gibbs equivalence in L3 supplies
an exact target, not that missing inequality.
