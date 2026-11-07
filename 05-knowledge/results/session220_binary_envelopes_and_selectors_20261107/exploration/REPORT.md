# Single cross-edge switch: exact recovery, not a completeness advance

Status: elementary proved identities below; 841 FINITE-EXACT literal controls. **No novel gate or Chan–Pak proof is claimed.** No independent proof audit yet. The source snapshot is acbad6c8636cf8cbb675b4c947fa19e33ac0df23. All output is outside that read-only checkout.

## Inheritance and concept board

Closest mechanism: arbitrary ordered root amalgamation, including the seven-state formula in Session216 `paired-reset/FINITE_SEARCHES.md`. Closest hostile example: the Session217 adaptive uniform-module obstruction. Corrected near miss: Q/Z are commuting counters, not independent binary digits. Underused sidecar examined here: one selected cross-edge and its two actual endpoints.

The live concepts were: (1) marked cross-edge surgery, (2) two-path cover products, (3) root-retaining transfer, (4) root-amalgamation provenance, (5) span-two/three ribbon collective occupancy, and (6) odd-cycle packing. The fourth comparison shows that the proposed surgery is already subsumed by known machinery. The ribbon check reached THM-4880/HYP-9473 rather than a new invariant; its finite coverage does not supply all-target occupancy. No new OCF construction was established.

## Physical operation and exact formula

Let A,B be disjoint nonempty tournaments, with a in A and b in B. Start with A=>B and reverse only a->b. Call the resulting tournament T and retain root a.

For a rooted tournament U,v define h=H(U), e=end_v(U), d=H(U-v), and R the number of ordered spanning two-path covers whose designated path starts at v, allowing the other path to be empty. Use H(empty)=1. For B,b write u=H(B), s=start_b(B), E=end_b(B), l for the analogous root-ending two-cover count, r for its root-starting two-cover count, and c=B_B(2), counting two labeled paths with empty bins allowed.

Then

    h(T) = u h + l R - s e,
    R(T,a) = c R - r d,
    e(T,a) = E d,
    d(T,a) = u d.

The formulas include singleton inputs without exceptions. The four-state column (h,R,e,d) therefore has an exact linear transfer for every actual B,b. This is a stateful construction, not a scalar gate depending only on H(A).

### Proof

Old Hamilton paths destroyed by reversing a->b necessarily consist of an A-Hamilton path ending at a, followed by a B-Hamilton path starting at b. There are es of these.

New Hamilton paths use b->a exactly once. Deleting that seam uniquely gives the pattern A1 B1 b a A2 B2, where A1 and B2 may be empty. Equivalently, choose a two-cover of A with its designated path starting at a, and a two-cover of B with its designated path ending at b; concatenate the nondesignated A path, designated B path, designated A path, and nondesignated B path. The cross joins are all legal. This bijection gives Rl new paths and proves the first formula.

In a two-cover of T with designated path starting at a, b->a cannot be used: a is already the designated start. Before the flip, such covers of A=>B are counted by Rc, by independently splitting A and B between the two labeled bins and joining corresponding bins. The destroyed covers use a->b at the start of the designated path. Its A portion must then be the singleton a; the other A portion has d choices, and the B two-cover has r choices. This proves the second formula. Paths ending at a must end with b->a; the preceding path on A-a=>B contributes dE. Root deletion leaves A-a=>B, contributing du.

## Exact equality boundary

The switch never decreases H. It preserves H if and only if a is a sink of A and b is a source of B; otherwise its positive increment is even and at least 2.

Indeed R>=d>=e and l>=H(B-b)>=s. The increment is

    Rl-es = (R-e)l + e(l-s).

Here l>0. Equality forces R=e, hence R=d=e. If a had an outneighbor x, a designated path (a,x) and any Hamilton path of the remainder would be an additional cover beyond the d singleton-root covers. Thus R=d implies a is a sink. Conversely a sink has R=e=d>0. Equality then requires l=s. The reversed argument proves that this holds exactly when b is a source. Tournament Hamilton-path parity gives the even increment. The sink/source case merely swaps adjacent singleton boundary components, so its equality also has a direct graph interpretation.

## Why this is inherited rather than a new construction law

Let W be B plus a new root a, with b->a and a->B-b. Form the ordered root amalgamation of (A,a) and (W,a), orienting all A-a vertices toward B. This is exactly T. Thus the inherited arbitrary root-amalgamation transfer already supplies the operation. The separate universal seam transfer in `05-knowledge/results/seven_return_native_recursions_20261029/SEAM.md` and THM-4875 also cover the underlying one-arc join, with converse orientation conventions.

When B is a singleton, the displayed four-state transfer becomes precisely inherited Z:

    (h,R,e,d) -> (h+R-e, 2R-d, d, d).

When B is transitive with b at position j (zero based), l=2^j, r=2^(|B|-1-j), c=2^|B| and u=1. This is exactly j global-sink additions, one Z, then |B|-1-j global-sink additions, retaining a throughout. Nontransitive B changes the deletion count multiplicatively, d'=u d; that coordinate cannot be dropped to assert binary induction.

## Verification and remaining boundary

Run `python3 verify_switch.py`. It independently computes endpoint subset dynamic programs and convolutions on literal graphs, then compares all four transferred coordinates for every pair of rooted labeled input tournaments of orders 1,2,3. There are 29 rooted inputs, 841 ordered input pairs, and 81 zero-increment cases; all pass. The zero classification and even nonnegative increments are also checked. These finite checks support the written identities; they do not establish any asymptotic coverage.

The precise missing positive result is a target-dependent choice of physical inputs and marks, or a closed reachable-set invariant, that solves u h+l R-s e=N while maintaining the needed continuation coordinates. No such rule was found. Neither the switch identity nor the ribbon finite census proves that every odd N>=23 occurs. This recovery packet should not be promoted as a new Chan–Pak milestone.
