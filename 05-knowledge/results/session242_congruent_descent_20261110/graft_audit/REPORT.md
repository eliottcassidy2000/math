# Independent hostile graft audit

Status: VERIFIED mathematical argument and FINITE-EXACT literal witnesses.
The final producer `exception/PROOF.md` has now received a line audit and passes. This is
not a proof of Chan–Pak or a classification of arbitrary normalized states.

## Claim, quantifiers, and proof

A proper graft here has a nonempty rooted source A with a nonsink root a,
and any nonempty tournament B. All original source vertices remain. Each
vertex of B beats a, and every vertex of A other than a beats every vertex
of B. Let H and Q count the gadget Hamilton paths and ordered two-bin path
covers (empty bins permitted); g=gcd(H,Q). Then its count N satisfies N/g !=7.

Use inherited source state (D,E,U,d), with D=H(A), E the number of paths ending
at a, U the number of ordered covers whose first path starts at a, and
d=H(A-a). For f=U-D-d+E, the scalar identity is

    N=Hd+Q(D-E)+(Q-H)f.

Physical inequalities give d>=1, D-E>=1, f>=0, Q>=2H. Thus N>=H+Q.
Equality holds exactly when d=1, D-E=1, f=0, since all three coefficients
are positive. For a nonsink, E<=(D-1)/2 then forces D=1,E=0,U=2.
This equality is attained by the transitive two-vertex source rooted at its
initial vertex, for every gadget B.

For any tournament, the Hamilton-path count is odd (the classical parity
identity, inherited throughout this repository). For |B|>=2, pair each subset
S with its complement in

    Q=sum_S H(B[S]) H(B[V-S]).

Every summand is odd, and there are 2^(|B|-1) complementary pairs. Each pair
contributes 2 modulo4, and their number is even. Therefore 4 divides Q.
Because H and g are odd, h=H/g is odd and q=Q/g is divisible by4.
Also q>=2h. If N/g=7, positivity gives h+q<=7. Together these bounds force
(h,q)=(1,4). The inequality q>=2h is indispensable: h+q<=7 alone also
permits (3,4).

The pair (1,4) is physically realized by TT2. Replacing B by TT2 gives the
same rooted source and an actual tournament with seven Hamilton paths, which
is impossible by the inherited forbidden-seven theorem. No division of a
rooted state and no assumption that an arbitrary normalized pair is physical
occurs here. This one small forced pair is physically realized explicitly.

There is also a direct obstruction for this non-singleton case. Parity
U=E modulo2 makes f even. In

    7=d+4(D-E)+3f,

f>=2 would give at least11. Thus f=0. Positivity forces D-E=1 and d=3;
but D>=d and E<=(D-1)/2 give D-E>=2, a contradiction. The singleton gadget
B has (H,Q,g)=(1,2,1), so this case instead follows directly from the inherited
unrestricted absence of seven. The 4-divisibility proof must not be applied
to singleton gadgets.

For N=7p with p prime, g divides N and g<=H<=N/3. Hence g cannot be7p.
The case g=p would yield N/g=7 and is excluded. Therefore g is1 or7. The
statement remains correct for p=7 (duplicated divisor cases collapse).
This is necessary, not sufficient, and does not say all gadgets have gcd1.

## Sharp boundary and physical controls

The nonsink hypothesis is essential, not a technical convenience. With a
singleton rooted source and gadget rows [22,60,8,81,4,29,55], the output has
49 paths, H=49,Q=812,g=7 and N/g=7. Thus the theorem is false for unrestricted
root positions. In general sink roots give N=Hd.

For the same gadget, rooted cyclic-triangle source gives a strong ten-vertex
output with1673=7*239 paths. Therefore g=7 really occurs for prime cofactor
targets. The pair (H/g,Q/g)=(7,116) need not itself be a realizable gadget:
H=7 is forbidden. This explicitly blocks a general physical-normalization
argument even though the quotient-seven obstruction above is valid.

The second strong source has rows [126,180,120,2,168,72,154,45], root3, and
(D,E,U,d)=(259,70,266,77), of common content7. Appending a sink and then a
singleton graft at the same original port gives a strong ten-vertex graph
with791=7*113 paths. State content is not an automatically removable physical
factor. In particular E/7=10 and d/7=11 coexist with D/7=37,U/7=38;
the audit does not assert that this divided state is physical or nonphysical.

## Independence, scope, and reproducibility

`verify.py` imports no producer code. It uses memoized first-vertex recursion
for all induced counts, independently compares full-vertex permutation
counts on the two seed graphs and both ten-vertex outputs, verifies strong
connectivity, and exhausts all labelled source tournaments through order4,
all roots, and all labelled gadgets through order3. Literal graph outputs
check the scalar row, proper lower bound, normalized quotient exclusion and
sink identity. Checks raise exceptions and remain active under optimization.

Run from this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -O -B verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

`--fail` is an intentional failing control and must exit nonzero in both modes.
These finite tests are not the all-order proof above.

## Priority and dependency audit

The scalar row, positivity bound H+Q<=N, strict scalar descent and prime-target
gcd1 condition are inherited, explicitly present in Session229
`inverse/PROOF.md`, with earlier Sessions225/227 attributed there. They are
not new Session242 discoveries. The present new contribution is the
quotient-seven obstruction for every proper graft and its sharp7p consequence,
plus explicit content-seven hostile witnesses. A targeted exact-phrase and
formula search found no earlier copy of this obstruction; this is a repository
provenance check, not a claim of literature novelty.

Session229's proof correctly distinguishes a finite sound inverse from
existence of a successful branch, and retains full rooted state rather than
scalar D alone. Its scalar inequalities and gcd consequence check out.
However its reference to
`05-knowledge/results/session225_graft_recovery_20261108/PROOF.md` is absent
from the audited checkout. Its copied `inverse/inherited.py` is present.
This audit does not re-certify every absent dependency of the finite physical
membership oracle, nor promote its universal branch existence to a theorem.

## Final producer line audit

`exception/PROOF.md` was read after completion. Its normalization principle
explicitly requires an actual replacement gadget, handles singleton gadgets
separately, retains q>=2h, proves the local parity obstruction correctly, and
separates root content from gadget content. Its cut-injection argument is
correct, and its finite-bank content filter is necessary only, as stated.
Both the all-order quotient-seven theorem and seven-prime corollary PASS.
The exact witness statements PASS independently. No global inverse coverage
or congruent-number-to-tournament map is asserted.
