# A Hilbert–Gray atlas with a permanent seal and an explicit seed-frame reset

**PROVED:** a source-retaining sealed reader, its four-coordinate transfer,
a marked Hilbert–Gray atlas, an exact local edge response, and a physical
seed-promotion/frame-retarget template. **FINITE-EXACT:** small named graph,
phase and native-guard controls in `check.json`. No all-integer coverage,
prime/7p induction, recursive219-key closure or minimal-dimension claim.

## 1. Recover the right predecessor and its corrections

The closest proved mechanism is **H1–H4, `return701_algebra_20261103/SPACE_FILLING.md`**.
It labels Hilbert cells by certified binary arithmetic sources, keeps the
source clock and finite orientation controller, and identifies ancestry
with cell inclusion. It explicitly does not identify scan adjacency with
a native Collatz edge. Its H4 separates a finite controller from an unbounded
integer register. We retain that separation here.

The canonical hostile is loss of a named joint attachment after marginal
compression. The corrected near miss is **MISTAKE-065** and **THM-452,
`01-canon/theorems/THM-452-doubling-tiling-coset-law.md`**: converse, reversal
of the Hamilton path, and tile complement are different operations; reversing
both graph and path cancels one putative bit negation. The least-used relevant
sidecars are the fixed boundary, two mirror halves and their defect in
**THM-830, `01-canon/theorems/THM-830-b3-deletion-deck-mirror-current-calculus.md`**,
and the marked carry in `gray_selector_recursion_20261018e.md`. Its chart union
is not achieved by gluing all alternatives onto the same ports.

The historical S136/S160 Gray scans are useful provenance, but several prose
claims are not imported. Flipping a→b changes degrees by `(d_a-1,d_b+1)`;
the sorted score multiset is unchanged exactly when `d_a=d_b+1`, not when
the two degrees are equal. A directed triangle is a smallest hostile to
that latter claim. Also a transitive4 with its long arc reversed has H5
instead of H1, contradicting S160's asserted `|ΔH|<=Hmax(3)=3`. The proposed
contraction is not necessarily a tournament. The exact ordered-cover current
is the correct replacement. A finite Gray traversal is a discrete address
schedule; it is not, by that fact alone, a continuous space-filling theorem
or an intermediate-value theorem for H.

The concept board is **actual seed; permanent seal; finite frame; arithmetic
register; named mixed attachment; acceptance condition**. The new map below
connects the first five. The sixth remains a separate induction obligation.
No external literature result is needed: we use the previously proved finite
Hilbert recursion and prove the new graph transport directly.

## 2. The permanent zero seal makes reader ancestry literal

Let U be any nonempty tournament with a marked vertex w. Add a seal z that
beats exactly w and loses to every other U vertex. Add a transitive prefix
A=TT_k, k>=1, dominating U and z, and an outer root r that loses exactly to w.
Call the result `C_k(U,w,z)`. All these are physical names.

**Strong-collar lemma.** The full C is strong for every U and every w,
including a singleton U or a mark that does not reach U.

Indeed, `w→r→z→w` is a directed triangle. Every remaining vertex v receives
an arc from r and sends one to z, so it reaches and is reached from that
triangle. The argument even permits an empty A, but k>=1 is retained for
the zero-ready/dyadic-clock interface used here. The inner component U+z
need not be strong. A sufficient condition for that stronger claim is
`|U|>=2` and w reaching all U vertices: then w reaches a second vertex,
which reaches z, and every non-w vertex reaches z→w. Do not confuse these
two strongness statements.

A free reader letter is inserted immediately **before z** in logical reader
order. Give its newborn x exactly the following old out-neighborhood:

```
bit1: {z},
bit0: {z,w}.
```

Thus x loses to the prefix, r and every old unsealed vertex except possibly
w. Its restriction to U is exactly the old D/O reader: a global sink for
bit1, or a vertex beating only w for bit0. It points to the retained z.
Consequently the same old seal stays the final zero-letter in logical
reader order, and every old name and old arc survives every new letter.
The actual parent is an induced subtournament, not another representative
with the same H. If the original seed is nontrivial and its mark reaches it,
every inner sealed component remains strong as well.

## 3. Four exact scalar registers and the converse command audit

Use the unsealed U fields from the previous reader:

```
h=H(U), e=E_w,U(1), q=H(U-w), t=S_w,U(2),
F=E_w,U(2), p=P_(U-w)(2), M=2^k.
```

Define the four-coordinate carrier

```
(a,b,c,d)=(h+MF,e,q,t+Mp).
```

The sealed full graph has exactly

```
H(C)=a+d-b.
```

To prove this, a final zero reader letter changes `(h,F)` into
`(h+t-e,F+p)`; wrapping in TT_k gives the formula. Applying the established
D/O reader laws to the unsealed U gives

```
D / bit1: (a,b,c,d) -> (a,0,c,2d),
O / bit0: (a,b,c,d) -> (a+d-b,c,c,2d-c).
```

This is a four-coordinate exact carrier, not a claim of minimality or a
finite set of numerical states. Its registers are unbounded integers.

There is also an exact older observer interpretation. In `C^op`, z is a
one-inport child of w, and the old parent `(C-z)^op` has rooted profile
`(h,s,q,F)=(a,b,c,d)`. Its fixed ordered live ports are `(w,z)`. A free bit1
becomes a newborn beaten only by z; a free bit0 becomes a newborn beaten by
both w and z. Both old live ports are retained. In the universal convention
`3*incoming_mask+discard`, the commands are therefore exactly

```
bit1 -> mask2, discard2 -> command8,
bit0 -> mask3, discard2 -> command11.
```

No path reversal or tile complement is silently composed with this graph
converse. H is unchanged because reversing every Hamilton path is a
bijection. `check.py` verifies literal adjacency and the exact coefficient
identities `L_bit F = F A_command` on the inherited four-to17 root-child
factor in the universal static query basis. This connects the small reader
to the actual operation language, not merely two equal dimension counts.

An independently retained spare u is a physical vertex throughout, but is
not a third live slot of this two-port language. Its later promotion needs
its joint data or the full actual graph. The earlier ten-coordinate phase
carrier and mixed three-bin counts remain necessary for their stated
larger observation tasks; this four-coordinate reduction only targets H.

## 4. Hilbert cells can carry both induced ancestry and one-arc scan moves

For depth d, let `g(i)=i XOR (i>>1)`, padded to2d bits and read most-significant
bit first. Those bits are the free reader word, followed logically by the
same retained zero seal. Put its actual named graph at the finite Hilbert
cell `point_d(i)` using H2's recursion:

```
point_0(0)=(0,0),
q,j = divmod(i,4^(d-1)),  s=2^(d-1),  (x,y)=point_(d-1)(j),
point_d(i) = [(y,x),(x,y+s),(x+s,y+s),(2s-1-y,s-1-x)]_q.
```

The inherited induction gives a bijection onto the2^d square grid, unit
Manhattan scan joins, and

```
floor(point_d(i)/2)=point_(d-1)(floor(i/4)).
```

The needed new address compatibility is exact:

```
g(i)>>2 = g(i>>2).
```

Therefore deleting the final two free reader vertices gives the actual
named parent at index `floor(i/4)`. The permanent seal is crucial: appending
a new zero-seal each time would shift or delete an old boundary name.

Consecutive Gray words differ in exactly one bit. At fixed length, all
reader vertices form the same transitive tail, and every incidence is
unchanged except that bit's mark–tail arc. Thus consecutive Hilbert scan
cells differ by exactly one old arc. By the proved one-helper theorem in
`adaptive_metagraph_reformulation_20261103/METAGRAPH.md`, that single arc
change has a literal realization using one helper birth, two guarded native
complements and helper deletion. All old physical names survive. Strongness
holds at the cell endpoints; it is not asserted at every native intermediate.

This is a genuine chart of actual sources. It does not make every neighboring
cell a3H+2/3H+4 return, and it does not make every nearby scalar value reachable.
For the C3 seed, k=2, the initial sealed value is33. The four first-level
Hilbert cells have values `141,107,91,127`; neither101 nor103 occurs.

The orientation matrices are generated by coordinate swap S and -S, hence
belong to `{I,S,-I,-S}`. Converting an address digit `(x,y)` to Gray letters
requires the previous binary bit c: output `(x XOR c,y XOR x)` and retain y.
Thus eight reachable frame/carry states suffice. Coordinates, depth, the
four arithmetic registers, the seed identity and physical names are separate
unbounded data. This is compatible with H4 and with the new arithmetic
agent's nonregular synchronized address/value graph; no finite Boolean
selector or image-language consequence is inferred from that result.

## 5. The exact local scalar response prevents a false intermediate-value step

Take state `(a,b,c,d)` just before a free letter. Compare bit0 with bit1,
then apply a common suffix v of length l and finally the retained seal.
Immediately after the differing letter, the state difference is

```
(d-b,c,0,-c).
```

If l=0, the final H difference is `d-b-2c`. If l>=1, it is

```
ΔH = d-b-c [2^l + sum_(i:v_i=0) 2^i + 1_(v_0=0)].
```

Proof: the difference in c is zero; the difference in d doubles at every
suffix step; the first suffix step erases the difference in b. Each zero
suffix letter adds the current d difference to a, with the extra first-step
b term if that first letter is zero. Applying output `a+d-b` proves the
formula. This is a specialization of the already retained arc-cover
current, not an independently invented general arc-derivative theorem.
For a reversed old arc `w→x`, the same difference is `C_xw-C_wx`, with
C counting first-path end / second-path start.

A value step of +2 or +4 is therefore a specific equation in these retained
registers and suffix bits. Unit geometric distance does not enforce it.
For C3, k=2, take2d free bits, all1, and change just the last free bit to0.
They are consecutive in the Gray scan and adjacent Hilbert cells, but

```
ΔH = 19*2^(2d-1)-2.
```

This grows without bound. Local geometry and one-arc locality cannot support
a discrete intermediate-value inference for the Hamilton spectrum.

## 6. A seed switch has an explicit physical frame-reset template

Let the full current unsealed graph be U, including **all old free reader
vertices**. Choose another vertex w' of U. Change exactly the four pairs

```
{z,w}, {z,w'}, {r,w}, {r,w'}.
```

The result is precisely `C_k(U,w',z)`. All U internal arcs, all prefix arcs
and all physical names are unchanged. The full endpoint is strong by the
collar lemma, independently of w' reachability inside U. Four applications
of the inherited arc compiler give eight guarded native complements,
four helper births and deletions, and no deletion of an old vertex. If
w'=w, no reset is needed. The old active mark becomes a retained passive
vertex; the new mark may be the old spare. Their roles are explicitly changed.

A one-vertex seed promotion can precede this reset. Birth y with any desired
out-neighborhood inside U, with y beating z and losing to r and A. This
preserves the entire old graph and realizes `C_k(U+y,w,z)`. Then retarget
w if required. This pays the physical realization of a proposed new seed;
it does not produce the target arithmetic equation or a bounded seed bank.

**Absorption is essential.** The new seed is the entire old U, optionally
plus y. It is generally false that the original small seed with the same
old word can now be read at w'. For the C3 seed and one free zero-letter,
retargeting0→1 while reinterpreting that letter would additionally change
the two internal pairs `{0,free}` and `{1,free}`. The four-edge reset leaves
both untouched. Its honest state space can grow; it is not closure of the
original219-key dictionary.

The receipt gives a literal one-vertex promotion followed by four arc
compilers, retaining the old outer root6, seal3, prefix4→5 and old spare1,
which becomes the new mark. It records actual source/target golden and
F49/F343 data. These finite-field records are current observer codecs, not
automorphisms or a claim of phase conservation through native intermediates.

## 7. Connection contract and remaining acceptance obligation

| Source -> map -> target | Preserved | Lost / retained sidecar | Cheapest decisive test |
|---|---|---|---|
| Actual named seed + free word -> sealed births -> full wrapper | Old induced source, names, exact H register and strong full endpoint | Independent spare responses require full graph/mixed counts | Literal birth cuts and counterpart command8/11 |
| Binary address -> Gray word + Hilbert cell | Prefix ancestry and unit scan adjacency | Arithmetic registers, frame, seed identity, physical names | `g(i)>>2=g(i>>2)` plus one marked-arc difference |
| Rooted parent profile -> natural-pair17 | All future H under the specified language | Arbitrary seed rewrite, other live names, phase | Exact matrix intertwining on all four coordinates |
| Current unsealed U -> optional one birth + four-star reset | Every old physical name and U internal arc | Original word is absorbed; new-root joint fields must be computed | Actual native guards and full endpoint matrix |
| Integer reader address -> arithmetic inverse | Exact represented target equation | A covering/reselection theorem is still needed | The independently found fixed-bank miss13835447 |

The fixed219-key miss is consistent with this construction. The chart
organizes existing words; the reset allows a proposed new seed to be
realized with its names. A constructive induction still needs a rule that
selects such a seed/mark/clock for every required prime or7p target, proves
termination or decreasing complexity, and retains the required acceptance
and attachment data. A Hilbert frame cannot replace that arithmetic rule.

## 8. Reproduction and bounded scope

Run `python3 .../hilbert/check.py` and the same command with `-O` from the
repository root. The script uses standard-library Python and the pinned
existing exact reader, endpoint-gap and arc-compiler libraries. It checks
all85 C3-atlas cells through depth3, every actual small seed/mark through
order3 for strong-collar controls, exact four-to17 identities, local edge
responses, the15 Gray scan arc compilations at depth2, and one explicit
promotion/four-edge frame reset. It performs no large seed or target census.
All checks use explicit exceptions. The script and source dependencies are
pinned in `check.json`; normal and optimized receipts agree byte-for-byte.
