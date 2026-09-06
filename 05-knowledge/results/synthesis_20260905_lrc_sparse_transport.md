# Sparse interval transport removes the all-height flow obstruction

**Status: PROVED ANALYTICALLY; FINITE-EXACT controls; independent proof and
raw-address/clipping audit PASS. No theorem ID has been reserved.** This extends the finite-only observation
in [THM-4409](../../01-canon/theorems/THM-4409-lrc14-third-sheet-component-network-certificate.md)
that its degree-zero max-flow equals the sum of edgewise minima. The equality
holds at every height. In fact every graph involved is a disjoint union of
stars. A second exact identity identifies the network relaxation with a
one-pair roof on the raw carrier lattice. Neither result proves the universal
`6/77` upper bound, physical entry, synchronization, or LRC(14), all **OPEN**.

## Inheritance and portfolio

The independent raw-address review is recorded in
[the cross-lane audit](synthesis_20260905_moments_lrc_audit.md); it reconstructs
the three sheet preimages and the clipped-tooth minimum independently.

The closest proved mechanism is THM-4409's component-incidence capacity,
with THM-4386's raw component address as its arithmetic carrier. The canonical
strict-loss hostile is `(1,19,79)`, where the best network is `8/553` but the
physical mass is `108/10507`. The equality control is `(1,5,11)`.
The corrected near miss is to replace a local maximum by an average, or to
move optimization through an aggregation without retaining its selector.
The relevant current warnings are ACTIVE-GUARDRAILS 7, 12, 26 and the
quantifier-order firewall; no matching correction to THM-4409 was found.

The least-used relevant sidecar is **parent separation**:
[THM-4100](../../01-canon/theorems/THM-4100-residual-component-three-outlier-lrc-compiler.md)
and [THM-4112](../../01-canon/theorems/THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families.md)
use gaps to control components of unions. Here the dual operation,
intersection, preserves an even simpler gap envelope. That operation makes
the apparent transportation optimization degenerate.

The anchor was the all-height THM-4409 carrier. The niche was the forgotten
ancestry/separation mechanism. The wildcard was to exchange max-flow for
raw-lattice roofs, then identify exactly which optimization cannot be moved
inside the sum. The concept board and its final update are:

| Concept | Original obstruction | What the new result changes |
|---|---|---|
| Pair/third network | Possible high-height Hall competition | None at degree zero: explicit optimal star flow |
| Component ancestry | A component might split repeatedly | The shortest-tooth parent cannot split under intersection |
| Raw affine-defect carrier | Graph and lattice seemed distinct | They compute the same one-pair roof relaxation |
| Equality / crossing geometry | Graph loses exact mass | Loss is exactly the sum of crossing protrusions |
| LRC consumer / entry | Local bounds lack the other speeds | Unchanged; no chart-entry implication is supplied |

The META-PATTERNS cards used were “Search the statement before the method,”
“Correct the object before sharpening the technique,” and “Type every
analogy and every implication.” This report does not promote a new method
card from a single thread.

## 1. A sharp sparse-interval transport lemma

Work on a fixed real interval; components touching only at endpoints do not
form contact edges. A finite disjoint interval family `A` is
`(delta_A,gamma)` separated when every member has length at most `delta_A`
and every gap between consecutive members has length at least
`gamma delta_A`. The two shores may have different `delta` values.
Empty families have zero capacity and may be omitted.

**Proposition 1 (PROVED).** Suppose `A,B` are respectively
`(delta_A,gamma)` and `(delta_B,gamma)` separated, with `gamma>=2`. Form the
bipartite graph of positive-length contacts. With vertex capacities equal to
interval lengths, its maximum flow is

```text
Cap(A,B;1) = sum_(I in A, J in B, |I intersect J|>0) min(|I|,|J|).       (1)
```

Every connected component of the contact graph is a star or an isolated
vertex. An optimal flow is given explicitly by the summands in (1).

**Proof.** First, if `I` contacts at least two members of `B`, then it spans a
gap of length at least `gamma delta_B`, so `|I|>gamma delta_B>=|J|` for
every adjacent `J`. Thus adjacent vertices cannot both have degree at least
two: that would force both `|I|>|J|` and `|J|>|I|`. This proves the star
claim, even for `gamma>=1`.

Now suppose `I` contacts `d>=2` intervals of `B`. They must occur
consecutively, and `I` spans all `d-1` gaps between them. Hence

```text
|I| > (d-1) gamma delta_B >= 2(d-1) delta_B >= d delta_B
     >= sum_(J adjacent to I) |J|.                                     (2)
```

Assign each contact `x_(I,J)=min(|I|,|J|)`. A degree-one vertex respects its
capacity directly. At a higher-degree vertex, (2) proves the capacity
constraint. The same argument applies on both shores, so all edge minima
are simultaneously feasible. Conversely every feasible edge load is at
most each endpoint's capacity, hence at most its minimum. Summing proves
optimality and (1). No max-flow/min-cut theorem is needed. QED.

**Sharp boundary.** The coefficient `2` cannot be lowered for arbitrary
interval families. For any `0<=c<2`, take unit leaves
`J_1=(0,1), J_2=(1+c,2+c)` and
`I=(1-epsilon,1+c+epsilon)`, with sufficiently small positive `epsilon`.
Then `0<|I|=c+2epsilon<2`, while two edge minima sum to
`2 min(1,|I|)>|I|`. The exact saved control uses `c=3/2`, `epsilon=1/10`:
the edge sum is `2`, and the max-flow is `17/10`.

This sharpness concerns the **general interval lemma**, not a claim that the
specific integer-comb application first fails at radius `1/6`.

## 2. Intersections preserve the parent separation

**Proposition 2 (PROVED).** Fix `gamma>=2`. A finite intersection of
interval families, each `(delta_i,gamma)` separated, has connected
components that are `(min_i delta_i,gamma)` separated.

**Proof.** Choose a constituent family with the smallest `delta_*`.
One of its intervals has length at most `delta_*`; it cannot meet two
different intervals of any other family, since their intervening gap has
length at least `gamma delta_i>=2delta_*`. Thus each chosen parent interval
contributes at most one interval to the intersection. Distinct output
components inherit the gaps between distinct parent intervals. QED.

For the ordinary shifted comb

```text
D_(a,theta) = {x mod 1 : ||a(x+theta)|| < lambda},
0 < lambda <= 1/6,
```

tooth lengths and gaps are

```text
delta_a=2 lambda/a,
gap_a=(1-2lambda)/a = gamma delta_a,
gamma=(1-2lambda)/(2lambda)>=2.                         (3)
```

Cut the circle at a fixed point and retain the clipped first/last pieces.
Clipping can only shorten teeth and does not shrink any interior gap.
Propositions 1 and 2 therefore apply to **every partition of any finite
collection of positive integer speeds into two nonempty subcores**, with
arbitrary phases and no height bound. Neither oddness, ternary-unit status,
nor primitivity is needed for this structural fact.

At THM-4409's `lambda=1/14`, `gamma=6`. For each of its six permutations and
each chosen pair, the pair intersection and third comb consequently form a
star forest. Its degree-zero certificate is exactly

```text
U_ij = sum_pi sum_contacts min(pair-component length, third-tooth length), (4)
```

at **every height**, extending THM-4409 Section 4's previously finite-exact
observation. Searching for a high-height transportation bottleneck cannot
improve this degree-zero certificate.

The same inheritance applies to the literal antipodal danger combs of
THM-4100/4112. Their tooth length is `1/(7v)`, and the gap/tooth ratio is
`6` for even `v` and `5/2` for odd `v`. Hence finite intersections and every
two-subcore transport split satisfy Proposition 1, even with mixed parity.
This is a valid transfer of their parent-gap sidecar; it is not a transfer
of their sufficient LRC supplier conclusion.

## 3. Weighted boundary and exactness

**Proposition 3 (PROVED).** Under Proposition 1, let a common weight satisfy
`0<m<=phi(x)<=M` on all intervals, with `M/m<=gamma/2`. Replacing interval
lengths by their `phi` integrals still gives

```text
Cap(A,B;phi) = sum_contacts min(integral_I phi, integral_J phi).          (5)
```

Indeed, for a degree-`d>=2` star center, (2) and
`M d <= m gamma(d-1)` imply that the sum of all leaf integrals is at most
the center integral. Matching components are automatic. At radius `1/14`,
the sufficient distortion ratio is `3`. This includes the constant case,
and for example the remainder weight `1-g` whenever `0<=g<=2/3`.

**REFUTED extension: arbitrary nonnegative weights.** With
`I=(9/10,71/10)`, `J_1=(0,1)`, `J_2=(7,8)`, the leaf gap/tooth ratio is
six. Set `phi=1/100` on `I` and `phi=1` outside `I`. The center mass is
`31/500`, each leaf mass is `901/1000`, and the two edge minima sum to
`31/250`. This exceeds the center capacity. The precise failed implication
is geometric length domination to weighted mass domination. Thus the
unweighted theorem does not remove THM-4409's nonzero-Fejer flow problem.

For constant weight, (1) gives the exact information-loss identity

```text
Cap(A,B;1)-|union A intersect union B|
 = sum_contacts [min(|I|,|J|)-|I intersect J|].                          (6)
```

Every summand is nonnegative. It vanishes exactly when the two contacting
intervals are nested. Consequently the capacity is the true intersection
measure **iff every positive-length contact is nested**. Matching is not
required. For a crossing edge the summand is the smaller interval's
protruding length. The missing coordinate is crossing endpoint position,
not a hidden Hall constraint.

## 4. An exact raw-lattice formula for the network relaxation

Use the primitive distinct odd ternary-unit assumptions of
[THM-4386](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md)
and [THM-4392](../../01-canon/theorems/THM-4392-lrc14-raw-carrier-box-spline-fourier-poisson-duality.md).
Write `r=3/14`, `Lambda_w={C in Z^3:C.w=0}`, and, for `{i,j,k}={1,2,3}`,

```text
q_i = 2r/w_i,
d_ij(C) = r/w_i+r/w_j-|C_k|/(w_i w_j),
L(C) = max(0,min(q_1,q_2,q_3,d_12(C),d_13(C),d_23(C))),
Omega_w = {C in Lambda_w : every C_i nonzero mod 3 and L(C)>0},
K_ij(C) = min(q_1,q_2,q_3,d_ij(C)),  C in Omega_w.                       (7)
```

**Proposition 4 (PROVED relative to the raw-address theorem).** The exact
degree-zero network and physical mass are

```text
U_ij = sum_(C in Omega_w) K_ij(C),
mu(F_w) = sum_(C in Omega_w) L(C),
L(C) = min_(ij) K_ij(C).                                                 (8)
```

**Proof.** The raw-address theorem identifies a component with the lift
class `n mod Zw` and `C=w cross n`. In its real lift, the three intervals
have lengths `q_i`, and the distance between the `i,j` centers is
`|C_k|/(w_i w_j)`. Their pair-intersection length is
`min(q_i,q_j,d_ij(C))`; it is positive on `Omega_w`. The edgewise minimum
with the third interval is exactly `K_ij(C)`.

The change of variable `y=3x` gives three `x` preimages for each raw
component, each with every interval length scaled by `1/3`. The distinct
owner residues choose exactly the three cyclic sheet assignments across
these preimages. Reflection accounts for the other orientation through the
corresponding raw carriers; no extra factor of two is introduced. Thus
summing all physical edges across all six sheets gives the first equation
of (8). Proposition 1 identifies that edge sum with the network capacity.
The second equation is THM-4386, and the third is immediate from (7).

There is a cut-at-zero detail. No pair component for distinct sheets crosses
zero, since a nonzero sheet excludes zero. If the third tooth does cross
zero, it has center zero. Any nonzero-sheet tooth at speed `a` lies at
distance at least `(1/3-lambda)/a` from zero, while its length is
`2lambda/a`. At `lambda=1/14` the distance/length ratio is `11/6>1`.
Therefore any pair component contacting a clipped third half is shorter
than that half. Clipping cannot change the relevant minimum. This verifies
that (8) uses exactly THM-4409's clipped interval convention. QED.

The support `Omega_w` in (7) is load-bearing. It retains **all three**
positive-overlap conditions. Replacing it by the weaker support
`K_ij(C)>0` is a different and potentially much larger lattice sum.

The previously concealed quantifier is now explicit:

```text
mu(F_w) = sum_C min_ij K_ij(C)
       <= min_ij sum_C K_ij(C) = min_ij U_ij.                              (9)
```

Equality holds iff at least one pair attains the pointwise minimum for
every raw carrier. That pair is a global selector; a carrier-specific
selector always exists but recovers the exact physical computation.
No mean of the three pair capacities proves the desired maximum bound.

For the canonical hostile `(1,19,79)` the three values are

```text
U_12=12/553, U_13=8/553, U_23=184/10507,
mu=108/10507, min(U)-mu=44/10507>0.                                      (10)
```

Thus no global pair realizes every carrier's cheapest roof, despite the
absence of any flow competition. Formula (9) connects the network to the
existing affine-defect roof classifiers with the exact coordinate loss
exposed.

Only two antipodal carrier classes are needed to expose that failure:

| Representative (also include its negative) | `L` | `K_12` | `K_13` | `K_23` |
|---|---:|---:|---:|---:|
| `(-16,5,-1)` | `5/1501` | `3/553` | `3/553` | `5/1501` |
| `(-7,17,-4)` | `1/553` | `3/553` | `1/553` | `3/553` |

The first class needs pair `23`; the second needs pair `13`. This is a
concrete selector conflict, not merely a numerical gap between bounds.

## 5. Connection contracts

| Source -> target | Map and preserved predicate | Destroyed information / needed sidecar | Decisive test |
|---|---|---|---|
| THM-4112 ancestry -> intersection carrier | Choose the shortest-tooth parent; preserve its component gap | Union ancestry and its LRC supplier are not preserved; retain the actual operation | Check each parent tooth meets at most one tooth in every other family |
| THM-4409 network -> star edge list | Keep every contact and both lengths; preserve exact degree-zero capacity | Overlap endpoints remain absent; crossing protrusions restore mass | (1,19,79) has strict loss even though the flow is explicit |
| Raw carrier -> network roof | `C=w cross n`, then `K_ij`; preserve every pair certificate | Primitive direction alone loses raw scale; retain full `C` and `Omega_w` | Independent integer intervals versus raw roof on every primitive triple through 43 |
| Three pair roofs -> actual mass | Take pair minimum separately on each carrier | A single global selector may not exist | Strict gap `44/10507` in (10) |
| Constant -> weighted transport | Integrate common density | Length domination need not imply mass domination; retain distortion or true flow | The exact ratio-six hostile in Section 3 |

## 6. Reproducible checks and remaining frontier

The source is
[synthesis_20260905_lrc_sparse_transport.py](../../04-computation/synthesis_20260905_lrc_sparse_transport.py),
with matching output
[synthesis_20260905_lrc_sparse_transport.out](synthesis_20260905_lrc_sparse_transport.out).
It computes integer-scaled interval endpoints, selects certificates before
using physical overlap measures, checks every vertex load and star edge,
and independently replays the augmenting max-flow implementation on four
controls. The raw-lattice routine imports neither interval nor graph code.

```powershell
python -B 04-computation/synthesis_20260905_lrc_sparse_transport.py
python -B -O 04-computation/synthesis_20260905_lrc_sparse_transport.py
```

The interval universe independently reproduces all `2,910` THM-4409 triples
through height `79`, including `1,747` exact best capacities and the sole
target equality `(1,5,11)`. A separately declared deterministic stress set
has `430` triples through height `3419`; it is **not** a complete height
census and gives no additional universal claim. It finds no sharp-target
hostile. The direct raw-roof replay covers every primitive distinct odd
ternary-unit triple through height `43`, plus three additional controls.
The final script executes `217,719` explicit checks. Its `457` raw-roof
controls and four full augmenting-flow controls pass, and normal/optimized
streams agree byte for byte. The source and output SHA-256 values, using
raw LF bytes, are respectively:

```text
8faa67e5cb1e713951c5fda82d327068d75dcb9381a35adf0595e8bb9805e2da
e23754d48ffb2ad64e0c622b61e33cfb6e1502dc01d0b7c816c0ee4c9551ee1c
```

**Direct progress:** the all-height flow/edge-min identity and star-forest
structure are proved. **Niche progress:** sparse interval families form a
class closed under intersection, including the parity-dependent antipodal
families. **New bridge:** the component network is an exact one-pair raw
roof relaxation. **Refuted extension:** unrestricted weight transport.

The next nonredundant target is the arithmetic inequality

```text
min_ij sum_(C in Omega_w) K_ij(C) <= 6/77                                  (11)
```

for every primitive distinct positive odd ternary-unit triple. A proof
must control the raw carrier population and the global pair selector,
or quantify crossing debt sharply enough to work directly with the mass.
Additional max-flow machinery at degree zero cannot address that task.
Even proving (11) would leave chart entry and synchronization with the
other speeds, so **LRC(14) remains OPEN**.
