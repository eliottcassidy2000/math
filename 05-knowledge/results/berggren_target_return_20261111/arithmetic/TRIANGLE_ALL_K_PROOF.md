# All births: half-contact budgets for triangle-chain extension content

**PROVED IN NOTE; independently reviewed within this session.** This note
classifies the complete arbitrary-gap, arbitrary-k extension-content
filtration of the stated triangle-chain family. It began as a proof candidate;
the root reviewer checked the full grouping and budget argument, and the
Berggren lane independently passed the analytic proof plus exhaustive and
seeded abstract-budget controls in
[ALL_K_AUDIT.md](../berggren/ALL_K_AUDIT.md). The earlier independently reviewed
C2 theorem remains a separate specialization. Status was promoted only after
both reviews passed; no external referee or literature-priority claim is made.

Let

    T=TT_(g0) -> C3 -> TT_(g1) -> ... -> C3 -> TT_(gr),

with r>=1 and all g_i>=0. Let nu be the maximum matching size in the path
on the r triangles whose edge i exists iff g_i=0. Equivalently, if the
maximal zero-gap runs have lengths L_1,...,L_a, then

    nu=sum_j floor(L_j/2).

For the gcd C_k(T) of all k-vertex extension Hamilton counts, the claim is

    C_k(T)=3^max(r-k-min(k,nu),0),
    L(T)=r-nu.                                         (1)

The constructor upper bound and the weak divisibility
`3^max(r-2k,0) | C_k(T)` are already proved in
[TRIANGLE_FILTRATION.md](../fibres/TRIANGLE_FILTRATION.md). It suffices to
establish

    3^max(r-k-nu,0) | C_k(T).                           (2)

## 1. Coefficient groups are products of local cover counts

Use a Boolean variable for every old-new and new-new edge orientation.
The Hamilton count is a multilinear polynomial. Its coefficient gcd is
its Boolean-value gcd by Möbius inversion, as in the reviewed
[C2 proof](TRIANGLE_C2_PROOF.md). For any coefficient, partition contributing
Hamilton paths by the order of the k new vertices, which old segments
between them are empty, the old block containing each old-new boundary
contact, the selected/unselected flag of every orientation factor, and the
selected endpoint labels.

An old block is either a whole C3 or a singleton. A contact is one used
old-new path edge: it names an end or a start of an old path segment.
Contacts at the same old vertex or in the same triangle are retained
separately. New-new edge factors and the signs from expanding1-X are fixed
within a group. **Unselected old endpoint labels are summed within their
specified blocks.** A group therefore has the form

    sign * product_(old blocks B) c_B,                  (3)

where c_B is a nonnegative integer local ordered-path-cover count. The
ordering of the old SCCs makes the local choices independent: each old
path visits those components in their ordinal order. Its required start
and end blocks determine where it is available. A group with c_B=0 may
be discarded. All valuations below refer to positive local counts.

Let m_B be the number of contacts in block B and let m=sum m_B. If e_NN
is the number of used new-new path edges and e_end is the number of full
path endpoints occupied by new vertices, incidence counting gives

    m+2e_NN=2k-e_end,
    m/2+e_end/2=k-e_NN<=k.                             (4)

The argument will allocate half a unit to every contact, plus half a unit
to each new full-path endpoint. These are numerical proof budgets, not
assumed graph symmetries.

## 2. Triangle losses and the only unpaid local case

For a triangle B write

    loss(B)=1-v_3(c_B).

If no selected coefficient variable names one of its vertices, cyclic
rotation of its labels acts freely on its local covers, including any
unselected endpoint labels. Hence 3 divides c_B and loss(B)<=0.

If m_B>=2, regardless of selected labels, integrality gives

    loss(B)<=1<=m_B/2.                                 (5)

Thus repeated contacts need no detailed enumeration: their two half-units
already pay for losing the triangle's one baseline factor3.

The only potentially unpaid case is a marked triangle with exactly one
contact. Its endpoint is fixed by the selected variable. Suppose q old
paths are available in this triangle, counting the path with that boundary
event. Its exact local count is

    c_B=q^2-q+1.                                       (6)

Indeed, if the distinguished path contains1,2,3 triangle vertices, the
numbers of choices for the remaining vertices are respectively
`(q-1)^2`, `q-1`, and1. The first number includes assigning both remaining
vertices to one of the q-1 other paths, or assigning them to distinct paths.
The cyclic orientation makes every indicated internal order unique.

The number in (6) is divisible by3 exactly when q=2 mod3, and then its
3-adic valuation is exactly1. In that case loss(B)=0 and no extra payment
is needed. Otherwise loss(B)=1: its one contact pays1/2 and the triangle
requests one additional half-unit.

Call these latter triangles requesters. Each requester has exactly one
boundary event. If q_- and q_+ are the available-path counts just before
and just after it, then |q_+-q_-|=1 and q=max(q_-,q_+). When q is not2
mod3, exactly one of q_-,q_+ is divisible by3. Attach its request to the
corresponding region immediately before or after the block.

If R is the total number of requesters and m_T is the number of contacts
in triangle blocks, the preceding bounds give

    sum_triangle loss(B) <= m_T/2+R/2.                 (7)

## 3. Regions and three ways to pay requests

A region lies between successive blocks that contain contacts, or between
the first/last such block and the boundary of the old source. Intermediate
blocks have no contact, so the number q of available paths is constant.
The requester rule attaches a request only to a region with q=0 mod3.
There are at most two requests at a region, one from each endpoint block.
A region with q=0 containing any old vertex makes the local product zero;
only empty q=0 regions survive in a nonzero group.

**A free singleton pays.** If a region contains any old singleton, that
singleton has no contact and may be assigned to any one of the q available
paths. Its local factor is q. If q=0 the whole group vanishes; otherwise
this factor contributes at least one power of3, paying both possible
half-unit requests. Regions are disjoint, so the same singleton factor
is never spent twice.

**A singleton contact uses its own budget.** If an unpaid region ends in
a singleton block that contains a contact, charge the neighboring request
to that singleton's contact budget. A singleton with one contact changes
availability by one, so at most one of its adjacent regions has availability
divisible by3. It can therefore receive at most one half-unit charge.
A singleton with two contacts can receive at most two such charges, covered
by its budget1. There cannot be more than two contacts at one old vertex.
The opposite endpoint of any charged request is the requesting triangle,
so no ambiguity about splitting a request between singleton endpoints arises.
Writing m_S for all singleton contacts, this payment uses at most m_S/2.

**A terminal request uses endpoint slack.** At the left boundary the number
of available old paths is1 if the full Hamilton path starts old, and0 if
it starts new. Thus a terminal region with availability divisible by3 has
q=0 and requires a new full-path endpoint. The same holds on the right.
Each terminal region receives at most one half-unit request, so the total
terminal payment is at most e_end/2. This is precisely the slack retained
in (4); terminal regions must not be silently treated as positive gaps.

Classify requests in the displayed order, so these payments are disjoint.

## 4. The remaining requests are bounded by the zero-gap matching

Every remaining request is attached to a region whose two endpoint blocks
are triangles, with no singleton in between. There are no contact blocks
inside a region. Any intermediate triangle therefore has no contacts;
all singleton gaps along this interval are zero. Its endpoints and all
intermediate triangles belong to one maximal zero-gap run.

Fix a run of length L and let B be the number of its triangles whose
requests remain unpaid and attached internally to this run. Every triangle
requests at most one half-unit, so B<=L. If L is even, this immediately
gives B/2<=floor(L/2).

If L is odd, the only possible failure would be B=L. Then every triangle
in the run has exactly one contact and requests internally. There is no
contact-free triangle left in the run, so each unpaid region joins two
consecutive triangles. Its availability is a multiple of3. Each endpoint
has one event, so this is the unique multiple-of3 side of both endpoint
events: **both triangles request into the same region.** The requesters
are therefore paired by disjoint adjacent edges. An odd number L cannot
all be paired. Thus B<=L-1, again yielding B/2<=floor(L/2).

Summing over all maximal zero-gap runs shows that the remaining requests
cost at most nu. This step handles repeated events in other triangles:
such a triangle is not a requester, so it makes B<L rather than creating
an unmatched half-unit. No local start/end order inside a cyclic triangle
is imposed.

## 5. Sum the budgets

The free singleton valuations pay their assigned requests directly.
Writing V_S for the sum of the free-singleton local valuations, the explicit
request inequality is

    R/2 <= V_S+m_S/2+e_end/2+nu.

Equations (7) and this inequality imply, for every nonzero
coefficient group (3),

    v_3(group) >= r-m_T/2-m_S/2-e_end/2-nu
               >= r-k-nu.                             (8)

The second inequality is (4). Therefore every coefficient is divisible by
the integer power in (2): summing signed coefficient groups cannot destroy
a common divisibility factor. Combining (2) with the weak rotation bound
gives exponent

    max(r-k-nu,r-2k,0)=max(r-k-min(k,nu),0).

The existing disjoint H41/H5 constructors attain the reverse divisibility,
so this completes the proof of (1), including the lifetime r-nu.

## 6. Hostile checks and scope boundaries

The naive alternative of multiplying a separate one-event factor3 for
every contact fails. In C3 on0->1->2->0, prescribe one old path to start0
and end1, with one other path available. The sole local cover is
`(0,1)` and `(2)`, of count1. Its two events cannot be treated as two
independent one-event triangles. Inequality (5) handles it correctly.

A triangle need not be a consecutive block in the full Hamilton path even
when none of its vertices is adjacent to a new vertex. The rotation action
is on its labels and the local covers, not on an assumed consecutive block.
The concrete hostile and repaired proof are recorded in the sibling lane.

All sixteen complete C2 polynomial certificates in TRIANGLE_receipt.json
agree with (1). The k3 samples remain finite controls, not the proof of
the all-k statement. This theorem concerns the gcd of every extension
count, not existence of a specified count such as3H+2 or3H+4. It does not
close the target-preserving Chan–Pak return gate.

For a concrete target-specific boundary, the complete one-birth spectrum
of the fixed source C3->C3 is `{9,15,31,41,45,61,75,99}`. Its content is1,
but the +2 target29 is absent; the +4 target31 occurs, for instance when
the new vertex beats old vertices0,3,4 (mask25). TRIANGLE_local.py checks
all64 cuts by independent permutation enumeration and also checks (6)
for q=1,...,12 and the repeated-contact hostile of count1. This statement
does not exclude29 from other sources or from larger extensions.
