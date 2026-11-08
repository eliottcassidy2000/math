# Adaptive scalar switches and their exact propagation boundary

**Status: PROVED, independently audited.** The identities and
all-order sensitivity statement below follow from the inherited physical graft
law. The discovery census is separately **FINITE-EXACT**. No radix coverage,
cofiniteness, improved unrestricted scalar cutoff, or Chan–Pak proof is claimed.

## 1. Inheritance and the constructive test

Use the rooted tournament state

    x(A,a)=(D,E,U,d),

where D=H(A), E counts Hamilton paths ending at a, d=H(A−a), and U counts
ordered spanning two-path covers whose designated nonempty path starts at a.
The other path may be empty, and H(empty)=1.

The physical maps and the general graft matrix are inherited from
`session225_graft_recovery_20261108/PROOF.md`, Sections 1–3 and 5. Here O
adds a global sink retaining the root. Z adds a vertex which beats the old
root and loses to every old nonroot, retaining the old root. J makes the same
addition and marks the new vertex. L grafts a directed edge, marks its source,
and directs both new vertices into the old root and all old nonroots into
both new vertices. C is the analogous rooted directed-triangle graft.
The exact L matrix and its rooted-edge interpretation are also inherited from
`session227_chan_pak_constructions_20261108/boundary_escape/PROOF.md`,
Sections 5–6; they are not new claims in this packet.
Words are chronological, acting left to right. The initial rooted edge has
state t=(1,0,2,1).

The first exploratory route sought small physical radix returns and did not
find a coverage mechanism. The different constructive test here asks for two
literal words which produce different D but exactly the same (E,U,d).
Such a pair could potentially store a digit while leaving the quantities
used by future O/Z moves unchanged. The missing issue is whether that digit
can subsequently survive a root-changing construction without entering the
other coordinates. Sections 3–4 give the exact answer for a common context.

The closest inherited mechanism is the four-state closure. The canonical
hostile is replacing a physical predecessor by a formal cone point. The
corrected near miss is treating a two-point scalar fiber as a recursively
filled interval. The useful sidecar is the full physical D-column, rather
than the scalar slope alone. These are local construction results, not a
classification of all adaptive representations.

## 2. Two exact physical switches

Suppose a physical input has E=0. Every input can be put in this form by an
initial O, although that initial O changes U and is not silently omitted.
Direct multiplication of the inherited matrices gives:

    JOLO(D,0,U,d) = (5D+7U, 0, 18D+20U, 3D+3U),
    JJLO(D,0,U,d) = (5D+5U, 0, 18D+20U, 3D+3U).           (1)

Thus the two literal constructions have the same complete three-coordinate
sidecar (E,U,d) and a scalar gap 2U. There is a second pair, using only the
older one-vertex operations:

    JOZJO(D,0,U,d) = (5D+7U, 0, 12D+14U, 3D+3U),
    JJJJO(D,0,U,d) = (5D+6U, 0, 12D+14U, 3D+3U).         (2)

Its scalar gap is U. Equation (2) is included as a control and companion
identity, not as a claim of priority over the existing O/Z/J literature.
Every map displayed is a graph construction, so the formulas prove physical
realizability for every physical E=0 input. Neither identity prescribes the
value of U freely. In particular, its parity is even on such physical inputs.

For t, the first pair gives

    t.JJLO=(15,0,58,9),  t.JOLO=(19,0,58,9).

Appending Z yields (73,9,107,9) and (77,9,107,9). These are explicit physical
examples, not newly discovered unrestricted spectrum values. The maps before
this last Z have a global sink, so strongness is not asserted for them.

### Exact two-point fiber bound

For the first pair, if the common output sidecar is (0,V,r), then

    U=(V−6r)/2,  D=(20r−3V)/6.                           (3)

For the second pair the inverse is

    U=(V−4r)/2,  D=(14r−3V)/6.                           (4)

Indeed the two maps from (D,U) to (V,r) have determinant −6. The old d is
forgotten, but it does not occur in either output scalar. Thus, across an
arbitrary collection of physical E=0 inputs, applying the two words in (1)
produces at most two D values for each fixed output (E,U,d). The same statement
holds separately for (2).

The independent proof check sharpens this to arbitrary mixtures of all four
words. Their physical sidecar images are disjoint, since

    V−6r = 2U > 0                  for pair (1),
    V−6r = −6D−4U < 0              for pair (2).

Here D,U>0 on every nonempty physical input. Thus the final sidecar identifies
which pair was used, and that pair's inverse recovers the preceding (D,U).
Across arbitrary physical E=0 seeds and arbitrary positive-depth iterations
mixing all four displayed words, each final sidecar therefore has at most
TWO scalar D values. Every intermediate output has E=0, so the final-step
argument applies regardless of the preceding choices. Arbitrary other moves,
and the union with arbitrary zero-step seeds, are outside this bound.

Consequently neither these pairs separately nor their four-word union stacks
a growing interval of scalar values in one fixed sidecar fiber. This is an exact
failure of the proposed induction, not a conclusion from a finite coverage
experiment.

## 3. The D-column of every graft-containing context is physical

Let e_D=(1,0,0,0)^T. For every nonempty rooted gadget (B,r), including a
singleton, the inherited graft matrix K_B satisfies

    K_B e_D = (H(B), e_r(B), U_r(B), H(B−r)) = x(B,r).   (5)

This is exactly its first column. For a singleton, J e_D=(1,1,1,1), the
singleton's physical state. Also

    O e_D = Z e_D = e_D.                                (6)

**Physical sensitivity theorem.** Let W be any fixed finite context made
from O, Z, and arbitrary nonempty rooted grafts. The gadgets and root choices
may vary from step to step. If W contains no graft, its D-column is e_D.
If it contains a graft, its D-column is the state of an actual rooted
tournament T_W.

**Proof.** All O/Z moves before the first graft fix e_D by (6). The first
graft sends it to x(B,r) by (5). Every subsequent move has the proved physical
state law and can be performed on this actual rooted gadget. Their resulting
rooted tournament is T_W. Linearity proves W e_D=x(T_W). This works even if
the gadget root is a sink or B is a singleton. ∎

The theorem is a Session227 consequence of the inherited graft law; the
matrix formula itself is not new. It concerns one prescribed common context,
not arbitrary rerooting depending on the input graph or a general graph
operation outside this language.

## 4. Arithmetic consequence and the stopping reason

Let two physical input states differ only by a nonzero scalar amount:

    x_+ − x_- = Δ e_D,  Δ≠0.

Applying the **same** context W gives

    W x_+ − W x_- = Δ e_D                 if W has no graft,
    W x_+ − W x_- = Δ x(T_W)              otherwise.       (7)

Every finite tournament has a Hamilton path; this includes the empty deletion
under the convention H(empty)=1. Its U-coordinate is also positive, since
one can take the root singleton as the reserved path. Therefore a common
context containing any graft produces nonzero U and d differences of the
same sign as Δ. Their magnitudes are at least |Δ|. In particular no later
O/Z/graft suffix can erase those two differences again.

Rédei parity, inherited as `01-canon/theorems/THM-001-redei.md`, makes
H(T_W) odd. Hence the scalar gap in (7) is always a
positive odd multiple of Δ, preserving its sign and exact 2-adic valuation.
Before any graft the multiplier is 1. This is about the difference of two
Hamilton-path counts, not the parity of either count alone.

This explains the limitation of (1)–(2) without restricting the length of the
common continuation or its gadget alphabet. A pure-D switch is retained by
O/Z tails, but it cannot pass through a common graft and remain a pure-D
switch. A binary expansion in which an already stored gap must be doubled by
a common context is also impossible: its multiplier is odd.

The result does **not** exclude branch-dependent continuations, newly created
switches after a graft, intervals with varying sidecars, a state-dependent
Euclidean algorithm, or general tournament coverage. None of those missing
constructive mechanisms was supplied in this lane. An explicit next design
question is whether a genuinely many-to-one sidecar return can merge preceding
scalar choices, unlike the invertible sidecar maps (3)–(4). This is an OPEN
requirement for that proposed induction, not a theorem excluding other methods.
This is the precise
stopping reason, rather than a claim that adaptive arithmetic is impossible.

## 5. Finite discovery universe and independent controls

The canonical checker `verify.py` uses only Python's standard library and
explicit exception gates, which remain active under `python -O`.

- It checks the exact switch matrices, the inverse identities (3)–(4),
  the disjoint sidecar signs, and all four-word mixtures through depth six.
- Every labelled tournament of orders 1 through 4 is used with root 0:
  75 source presentations. All five literal O/Z/J/L/C moves are checked by
  independent Hamilton subset DP (375 checks). Every E=0 source also tests
  all four displayed words literally (124 checks).
- For all 75 labelled rooted-at-zero gadgets through order 4, four literal
  graft evaluations on a physical spanning basis independently recover the
  first column and compare it with the gadget's own profile. The signed
  combination used in this verification is linear algebra, not a claim that
  signed states are physical.
- All 3,906 O/Z/J/L/C words through length 5 check the sensitivity identity.
  Positivity in arbitrary contexts comes from the proof, not this finite list.

The additional finite discovery search starts with the E=0 input subspace,
examines **every word of at most seven letters** over O/Z/J/L/C, and adds a
final O. It deduplicates the resulting linear transfer matrices, not graph
isomorphism types. There are 76,543 distinct returned matrices, 353 sidecar
buckets with different D rows, and 355 unordered pairs of distinct scalar
rows in those buckets. None has D-row difference depending only on the old
d, and each branch of these colliding pairs has rank at most two on the
three-coordinate E=0 subspace. These last observations are FINITE-EXACT only.
The all-order sensitivity theorem is not inferred from the search.

Earlier scratch searches of affine eigenvalues and finite state fibers were
exploratory and do not enter the theorem or its receipt. In particular no
absence of a short radix loop is promoted to absence of a longer one.

Reproduce the canonical packet with:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

Chan–Pak completeness remains **OPEN**.
