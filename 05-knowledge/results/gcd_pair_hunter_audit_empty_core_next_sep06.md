# Independent CRT-overlap audit and six hereditary gcd ceilings

**Status: PROVED ELEMENTARILY / RELATIVE TO CITED LRC THROUGH THIRTEEN
RUNNERS + FINITE-EXACT complete profile certificate; independent audit PASS.**
For a primitive thirteen-speed strict counterexample `M(V)<1/14`, every
`(13-d)`-subset has gcd at most

| Complement size `d` | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| Maximum possible gcd under the proved cuts | 1 | 2 | 4 | 9 | 32 | 96 |

These are necessary conditions, not unsafe realizations. In particular,
a primitive row containing eight speeds of gcd greater than 32, or seven
speeds of gcd greater than 96, is `1/14`-safe. Common dilation replaces
the bounds by 32 or 96 times `gcd(V)` without primitivity. The final two
ceilings use the finite exact classification below after an analytic
reduction to explicit bounded clocks; there is no bounded speed-height
extrapolation.

## 1. Inheritance and the exact added input

[THM-1217, mixed-period beat-mask Hunter tree](../../01-canon/theorems/THM-1217-mixed-period-beat-mask-hunter-tree.md)
already proves lifting to one common clock and maximum-spanning-tree
intersection credit. [THM-860, primitive Hamming-six reduction](../../01-canon/theorems/THM-860-primitive-hamming-six-finite-ramification-reduction.md)
contains its fibre version and pointwise forest proof. Both were read;
neither is being asserted as a new abstract theorem.

[THM-4447, composite-clock capacity](../../01-canon/theorems/THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction.md)
gives the current exact open-arc count and divisor absorption. The
[nine-subset companion](gcd_nine_audit_empty_core_next_sep06.md) verifies
the preceding ceiling nine and its equality-phase hostile. The added
application here is a uniform lower bound on intersections of padded
capacity blocks, derived from balanced residue counts and CRT, and its
use in the complete divisor-profile classifier. Targeted prior recovery
found the Hunter mechanisms above; no external priority claim is made.

The live concepts are common clocks, maximal bad blocks, balanced residue
counts, tree overlap, and divisor inheritance. The corrected near miss is
subtracting overlap forced only after padding from unpadded singleton
sizes. The smallest hostile appears in Section 3.

## 2. Sharp pair-intersection lower bound for periodic blocks

Let `q,r` divide a common positive clock `c`. A unit-step cyclic block of
length `k<=q` in `Z/qZ` and one of length `ell<=r` in `Z/rZ` are pulled
back to subsets of `Z/cZ`. Put

```text
a=gcd(q,r),        k=Aa+u,        ell=Ba+v,
0<=u,v<a.
```

Reducing the first block modulo `a`, every residue occurs `A` or `A+1`
times, with exactly `u` high residues. The step remains a unit modulo
`a`, so this follows by dividing the block into full periods of length
`a` and a final distinct-residue segment. The second block has the
corresponding counts `B` or `B+1`, with `v` high residues.

CRT says a pair of residues modulo `q,r` has a common lift exactly when
they agree modulo `a`; in that case it has `c/lcm(q,r)` labelled lifts.
If `U,V` are the two sets of high residues, the intersection cardinality is

```text
(c/lcm(q,r)) * [aAB+A v+B u+|U intersect V|].
```

Hence it is at least

```text
L_c(q,r;k,ell)=
(c/lcm(q,r)) * [aAB+A v+B u+max(0,u+v-a)].             (1)
```

This lower bound is sharp over independent choices of unit steps and
block starts. Take both steps one, let the first extra segment be
`{0,...,u-1}`, and start the second extra segment at `u`. Their overlap
has size `max(0,u+v-a)`. Sharpness here concerns arbitrary periodic block
pairs. It does not assert that all pairwise minima can occur at one
physical body phase for a fixed tuple of tails.

For the LRC application use `k=ceil(q/7)` and `ell=ceil(r/7)` and write
`L_c(q,r)` for `(1)`. It includes `q=1` or `r=1`; no support or parity
filter is used. For coprime orders, the formula is the exact CRT product
`c*k*ell/(qr)`.

## 3. Padding and the lawful tree bound

At a fixed body phase, a tail of effective order `q` has a strict bad
set which is a pullback of a unit-step cyclic block of length at most
`k=ceil(q/7)`. Enlarge this block to length exactly `k`, retaining its
step. This is possible even for an empty block. Call its pullback `B_i`;
the actual bad set is contained in `B_i`, and

```text
|B_i|=(c/q_i) ceil(q_i/7)=:C_i.
```

For any spanning tree `T` on the tail indices,

```text
|union_i B_i| <= sum_i |B_i| - sum_(ij in T)|B_i intersect B_j|.
```

Indeed a point occurring in `s>=1` sets belongs to at most `s-1` edges
of the induced forest. Combining this with `(1)` gives the uniform upper
bound on the **actual** bad union

```text
H_c(q_1,...,q_d)=sum_i C_i
                   - max_(spanning trees T) sum_(ij in T)L_c(q_i,q_j). (2)
```

The maximum uses only proved pair lower bounds. They need not be jointly
attainable: each one holds for every padded pair in the actual chosen
family. Thus the tree optimization remains valid. If `(2)<c`, there is a
safe common label at every body-safe phase.

Both restrictions are load-bearing. Three identical singleton sets show
that subtracting all three pair overlaps, a cycle, can give zero for a
union of size one. At `c=6`, `y=1/2`, tails `(2,3)` have **empty actual
bad sets**, while their padded order-three/order-two blocks must intersect
once. Subtracting that padded credit from the actual singleton sum gives
`-1`, an invalid union bound. Formula `(2)` consistently uses padded
capacities and has neither error.

Named arithmetic controls missed by the scalar sum are:

```text
c=36, g=(1,4,4,6,9): scalar sum37; the q9/q4 pair alone forces overlap2;
c=27, g=(1,1,1,3,9): scalar sum27; the q27/q3 pair forces overlap1;
c=288,g=(1,1,4,4,18,32): scalar sum290; q16/q9 forces overlap12.
```

In every case the full maximum-tree bound is strictly below `c`; the
independent sidecar records its exact value. These are exclusions of
arithmetic profiles, not asserted actual non-lonely rows.

## 4. All-divisor absorption and child inheritance

For an exact body gcd `c` and tails with `g_i=gcd(c,w_i)`, choose any
`a|c`, `a>=2`. Absorb all tails with `a|g_i` into the pack. The remaining
tails have gcd coordinates `gcd(a,g_i)` and effective orders
`a/gcd(a,g_i)`. Full primitivity prevents every tail being absorbed;
the enlarged pack therefore has at most twelve speeds. Cited
lower-runner LRC supplies a phase with clearance strictly above `1/14`.
Apply `(2)` at clock `a`. A hypothetical failure must pass every such
cut.

Adjoining tail `i` to the core produces the child clock `g_i` and other
coordinates `gcd(g_i,g_j)`. For any divisor `a|g_i`, this child's absorbed
residual set equals the parent's residual set at divisor `a`: tail `i`
is already absorbed, and
`gcd(a,gcd(g_i,g_j))=gcd(a,g_j)`. Thus all parent divisor cuts include
all child divisor cuts. The child is primitive because
`gcd(g_1,...,g_d)=1`.

This gives an inductive exact classifier. Start with the zero-tail
primitive clock set `S_0={1}`. For `d<7`, any surviving profile has some
order `q_i` with `beta(q_i)>=1/d`. To see this, an order-one tail already
suffices; otherwise the clock-`c` cut `(2)>=c` implies the scalar sum at
least one. Since `beta(q)<=(q+6)/(7q)`,

```text
q_i<=Q_d:=floor(6d/(7-d)),
c=q_i g_i<=Q_d max(S_(d-1)).                            (3)
```

The child argument gives every `g_i in S_(d-1)`. Consequently enumerate
every clock through the explicit bound in `(3)`, every sorted divisor
word in `S_(d-1)`, require gcd one, and test all divisor cuts. No large
clock is omitted. Conversely, every accepted flat word has all child
cuts; induction identifies its child profiles in the previous complete
bank. The source checks this membership literally as an additional
independent control.

## 5. Independent complete classifier and evidence

The [source](../../04-computation/gcd_pair_hunter_audit_empty_core_next_sep06.py)
imports no primary producer. It loops over **every integer clock** through
`1,2,8,32,135,1152`, rather than generating clocks as products. It computes
the maximum spanning tree by Prim vertex growth rather than sorting edges.
The complete results, including clock one, are:

| `d` | Clock search bound | Primitive candidate words | Accepted profiles | Distinct gcds | Maximum gcd |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 2 | 2 | 2 | 2 |
| 3 | 8 | 16 | 5 | 4 | 4 |
| 4 | 32 | 208 | 19 | 7 | 9 |
| 5 | 135 | 3,592 | 110 | 16 | 32 |
| 6 | 1,152 | 474,844 | 1,217 | 43 | 96 |

The [full JSON bank](gcd_pair_hunter_audit_empty_core_next_sep06.json)
contains every accepted profile and gcd set. All six layers agree with
root's independent producer. Beyond classifier agreement, all unit/start
pairs of capacity blocks for `1<=q<=r<=18` are constructed literally;
their exact minimum intersection equals `(1)`. Common-clock multiplicity
is checked separately. Every unordered triple of subsets of a four-point
space checks the tree union inequality, and both hostiles in Section 3
are retained.

The [output](gcd_pair_hunter_audit_empty_core_next_sep06.out) reports
1,579,350 explicit gates. Normal and optimized outputs are byte-identical.

```sh
python3 -B 04-computation/gcd_pair_hunter_audit_empty_core_next_sep06.py
python3 -B -O 04-computation/gcd_pair_hunter_audit_empty_core_next_sep06.py
```

Semantic digest:
`b48461a5ad4384872372d48cdf1a31d28736367f8cc87a835a6f78df1688e84a`.
Raw SHA-256 source/output/JSON:

```text
dc5b2af3ae2df980081379a10a1b5517e272839b6332bb05a6c441592319865a
5cf8f668064b7321dcd7341b725f6cd53efc90373fa9e6381b049df7aa2c1c1c
dde323ce07917e530c562c2e9994386e005b5e1e42fa3fb7a9a93f93d0544a49
```

The source-to-target map retains common labels, divisor gcds, block
lengths and forced pair overlap. It discards the full unit and phase
arrangement. Thus surviving profiles remain abstract necessary data.
The old clock-nine fully-spoiled phase shows that this remaining
coordinate can matter even at small clocks. No LRC(14) closure or
extension of the finite bound to seven or more tails is claimed.

## Final primary-proof acceptance

The complete [primary proof](lrc14_recursive_gcd_empty_core_next_sep06.md)
and its standalone source have now passed final mathematical and type
audit. In particular, its sharper pivot-generated candidate universe is
complete, every child divisor cut appears in the parent absorption test,
and its named maximal signatures have the stated exact safe realizations.
The entire sorted signature arrays and gcd sets in all six primary JSON
layers agree exactly with the independent rectangular-universe JSON;
this comparison checks every retained tuple, not just counts or maxima.

The primary proof's seven-tail stopping result is also correct. For
`c=p` prime and seven gcd coordinates equal to one, the only nontrivial
divisor is `p`, no tail is absorbed, and total padded capacity is
`7ceil(p/7)>=p`. Every pair credit is zero because
`2ceil(p/7)<=p`: check `p=2`, and for `p>=3` use
`ceil(p/7)<=(p+6)/7<=p/2`. Thus the declared divisor/tree relaxation
retains arbitrarily large clocks. Formula `(7)` in the primary proof,
with six body speeds and all tail gcds one, realizes every such profile
by a primitive row safe at `1/(14p)`. This proves a boundary of the
arithmetic relaxation and supplies no unsafe-row or full-LRC claim.

All primary raw SHA-256 values were independently recomputed and matched:

```text
primary source
30d9e767cc502ec3bb818edf07d8ce8adf1546a4224abe731aaa0ee37d287063
primary output
86ad4e7f2a61c4efaa5e538646132ae9305c4fab849e295fa20e9809fae958c5
primary complete profile JSON
7af5425db35516037e03bb2d114c34bd1f2c8971cdee91866dc77b754b13f557
```

No mathematical correction was needed. A minor wording clarification
makes the forest proof refer explicitly to labels with positive
membership; labels in no padded set contribute zero to both sides.
