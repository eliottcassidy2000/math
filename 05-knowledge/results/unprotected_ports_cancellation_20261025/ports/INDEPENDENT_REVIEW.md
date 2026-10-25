# Independent review of the unprotected-port and two-triangle repairs

2026-10-25. Reviewer: the counting/recurrence agent, independent of the
port producer. **PASS:** the structural conclusions in [README.md](README.md)
and [ENDCAP_NEXT.md](ENDCAP_NEXT.md) follow as stated. All four maintained
normal/optimized checker runs passed. No demonstrated proof gap or required
producer correction was found.

## Scope accepted

For n>=12 the endpoint-four family I={h,n-3,n-2},5<=h<=n-4 has three
directed black moves to the exact labelled K_P(T), with no middle filter.
The general {n-3,n-2}-port theorem has precisely its stated sufficient
degree condition: k odd, or n outside {2k+1,2k+3}; the consecutive tail
at n=2k+3 is separately covered by the inherited theorem. The h=4,
1->3 branch is also proved. The addendum proves h=3,3->1 in the exact
endpoint-four family. These are local marked-path replacements. They
do not prove unrestricted black-core connectivity or remove guards from
the inherited definition of the endpoint-four frontier.

## Forward chord ledger

For 1->3, the A words have all the stated joins. Moving vertex2 keeps
the insertion predecessor h-1, unlike the earlier removed-vertex failure.
The first reversal sends the ordered pair (a,z) to the reverse-image pair
(b,x), giving the closing-edge contradiction. For A2, a is at position
h-2 and z at position n-2, so r_A2(z)=3. Reflecting the position of a
gives exactly w=n-h+2 when2h<=n+2 and w=n-h+3 when2h>=n+4.
The P-reflections of w are respectively h-3 and h-4, below the first
port. For h>=5 neither tested zw pair is the retained cz or zx edge;
the two values in equation (2) are therefore correct.

The omitted equality2h=n+3 fixes a under A2. Its degree in K_A1(T)
is k+1, so anti-reflection forces n=2k+3, then h=k+3. Cardinality
forces the complete consecutive tail. The last A word has the precise
insertion form required by the previously audited rotation lemma and its
fixed-a refinement. The same exact tail is its only remaining exception.
The inherited all-k repair starts from the original labelled T and applies
to that entire source, so no partially executed path is being spliced to
an unrelated certificate.

At h=4 the main A2 pair witness is correctly replaced: w=c would have
been retained. The literal reversal instead compares S(3,4)=1 with
S(y,z)=0. The final product has two disjoint four-cycles at n>=12.
The added h=4 statement is therefore supported by its own witness,
not by extending the old witness beyond its range.

## Backward chord, degree cocycle and even-order antipodes

For 3->1 the D ledger's new joins use a->4, z->b, a->x and3->b as
claimed. Its first witness T(a,x)=1 versus T(4,3)=0 is literal.
For its final witness at B, original P-blue reflection gives
T(n-4,c)=T(b,3)=0; this is a nonpath pair and is complemented to1.
The h>=5 hypothesis also gives B(a,3)=0. The correct final target is
established independently by the unordered path-edge XOR.

Comparing K_D1(T) with B changes only F={01,34,04,1z}. The endpoint
degree formula for a path-complement move gives delta_3=-1,delta_z=+1
and zero elsewhere. If D2 is blue, its reversal r and the P reversal
give the degree equation (6) exactly. Direct position calculation gives
the useful complete formula

    g(v)=v+4 for0<=v<=n-5,
    g(n-4)=2, g(n-3)=3, g(n-2)=0, g(n-1)=1,
    where g=r r_P.

It produces the displayed increments at0,n-5,n-4,n-3 and the two odd
cycle tables. The sums of degree offsets are -n for n=1 mod4 and0
for n=3 mod4. Applying the average degree identity in B forces respectively
n=2k+1 and n=2k+3, with k even. This step does not incorrectly assume
that g is already an automorphism of B: it uses the proved degree cocycle.

The even-order argument is also sound as an arc argument. When n=4t,
g has the two parity cycles of length2t. When n=4t+2, its cycles are the
four residue classes modulo4 of lengths t+1,t+1,t,t. In each case the
only within-cycle pairs from F are04 and1z, and both are adjacent pairs.
The remaining two F pairs cross cycles. At n>=12 the selected even cycle
has length at least4, so its antipodal pairs cannot be adjacent.

The identity r_P g r_P=g^(-1) follows because r and r_P are involutions.
It sends an entire g-cycle to a g-cycle of the same length, reversing
its cyclic direction; therefore it also sends antipodal pairs to
antipodal pairs. Both flip terms in equation (8) vanish throughout the
chosen pair orbit. After half a cycle, g exchanges the endpoints of an
unchanged arc, a contradiction. This avoids the invalid shortcut of
treating g as an automorphism on pairs that meet the defect set F.

The retained n13 and n15 examples show that the two odd critical degree
possibilities are feasible blue intermediate words. They are not hidden
empty cells and are not advertised as disconnected metagraphs.

## The h=3 two-triangle addendum

For I={3,y,c},3->1, the backward ledger remains legal: a->4 and a->x
are still true in this exact row. Its first two words are black at n>=12
by the same k=4 argument, whose odd critical orders9 and11 are below
scope. The former final-word witness can fail, as correctly retained.

For D3=(x,a,b,3,...,z), direct reversal composition is
(a,x,b)(z,c,y), fixing every middle vertex. Its automorphism condition
in B forces each triangle's incidence to every fixed middle vertex to
be uniform. The complete a-row gives the stated left-to-middle pattern;
reflection gives the right-to-middle pattern. Since a points to every
right-triangle vertex, the simultaneous rotations force all L-to-R arcs
forward. The path edges close each induced triangle cyclically. Conversely
these incidence rules make the two rotations an automorphism. Thus the
carrier is exact, with the actual blue middle retained; it is not an
ordered sum assumption about outer strongly connected components.

The carrier forces T(b,4)=1. On the opposite 4->b branch, the differing
B(a,4) and B(b,4) immediately make D3 black. On b->4 the E words are
legal for every source in the branch, even without conditioning on old
D3-blue. Their new nonpath pairs b4,zx,(n-4)a,ay occur twice, and all
other path-edge multiplicities have exactly the P parity. Their actual
reverse-image pairs are

    E1: (a,x) -> (b,3),
    E2: (b,x) -> (4,3),
    E3: (a,z) -> (b,4).

At the stated sources their values are respectively1 versus0,0 versus1,
and1 versus0. This verifies all three blackness arguments without a
new assumption on any middle arc. Both chord directions are exhausted.

## Reproduction and independent controls

Both modes of both maintained scripts passed during review:

    python check.py
    python -O check.py
    python endcap_check.py
    python -O endcap_check.py

The main checker certified1,044 endpoint-four chord cells for n12..40,
1,662 declared general cells for n12..17,1,008 h4 forward cells, and
degree/permutation/antipodal formulas through n300. Four source routes
use the inherited tail certificate. It retains the explicit n13/n15
blue-step sources and the h4 illegal-first-join and old h3-blue-final-word
hostiles.

The addendum checker certified104 complete affine chord cells at n12..60
and80,100,140, plus52 conditioned carrier cells. It verifies that the
opposite b4 chord is infeasible in the carrier and replays two literal
source models per chord cell. Its affine inconsistency certificates
quantify all reflected middle completions; the literal models supplement
them rather than defining the universe.

A separate no-import review script reconstructed reversals directly from
the words, independently verified the explicit formula for g, checked
the internal-defect adjacency and both antipodal flip terms for all251
even orders12..512, and checked the E-word XOR and the three reverse-image
pairs at all501 orders12..512. It passed normally and with `-O` using
explicit exceptions. The all-order scope rests on the structural proofs
above, not extrapolation from these finite intervals.

That independent control is preserved as [check_independent.py](check_independent.py).
Reproduce it with `python check_independent.py` and
`python -O check_independent.py`; it imports no producer implementation.

The resulting residual within the inherited endpoint-four remainder is
h=3 with1->3, or h=4 with3->1. No surviving cell is asserted to be
disconnected, and the global black-core problem remains open.
