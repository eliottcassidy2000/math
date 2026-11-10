# Physical content normalization and the seven-prime descent gate

**PROVED elementary consequences of inherited graft identities; FINITE-EXACT witnesses.**
Chan–Pak completeness remains OPEN. No new unrestricted coverage range or
literature-priority claim is made. The two-prime-lane reduction itself is
inherited, not a result of this note.

## Inheritance and scope

Closest mechanism: the unrestricted proper-graft inverse of Session229,
`05-knowledge/results/session229_combined_constructions_20261109/inverse/PROOF.md`,
Sections1–3, and its independent audit. Canonical hostile: a numerical gcd
quotient need not be an actual tournament interface; 7 and21 are not full
Hamilton counts. Corrected near miss: scalar strict descent does not construct
a physical predecessor. Least-used sidecar: the normalized **two-coordinate
gadget pair**, not just its normalized Hamilton count.

The board is: multiplicative atoms; full physical source state; gadget content;
fixed-port content; congruent-number square class. The direct target is the
seven-times-prime lane. The niche is primitive gadget normalization. The
wildcard is the precise contrast with the square-AP witness31,41,49.

Dependencies, each with the relevant concrete source:

* THM-4094, `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`:
  ordered joins, exceptional atoms, and the p/7p reduction. The exact filename
  in this checkout is listed in SOURCE_FILES.txt if a slug differs.
* THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`:
  the two full-count exclusions.
* Session229 inverse proof above: proper graft, physical inequalities, parity,
  two-input strict scalar descent and the finite physical inverse.
* THM-5092, `01-canon/theorems/THM-5092-four-test-tail-content-and-invisible-reroot-escape.md`:
  fixed-port S/B transfer and exact four-test content.

## 1. Definitions, and the normalization rule that is actually legal

Let B be a nonempty tournament. Write H=H(B) and

    Q = sum_{X subset V(B)} H(B[X]) H(B[V(B) minus X]),

where H(empty)=1. Thus Q counts ordered pairs of disjoint directed paths
covering all vertices, allowing an empty path. This Q is not Session237's
binary pair-energy statistic.

Let (A,a) be a physical rooted source with a a nonsink. Use its inherited
state (D,E,U,d), in which D=H(A), E counts Hamilton paths ending at a,
d=H(A-a), and U counts ordered covers with the first path starting at a.
The literal graft adds B with B→a and A-a→B. Set

    f=U-D-d+E.

Inherited physical constraints and the exact output are

    D,d positive odd; 1<=d<=D; 0<=E<=min(d,(D-1)/2);
    f>=0; U=E mod2; Q>=2H;
    N=Hd+Q(D-E)+(Q-H)f;  H+Q<=N; H<=N/3.            (1)

In particular f is even: D+d is even and U+E is even.

**Physical normalization principle.** Suppose (H,Q)=g(h0,q0), with g a
positive integer. If there exists an actual nonempty tournament B0 with
(H(B0),Q(B0))=(h0,q0), grafting B0 onto the SAME source (A,a) gives exactly
N/g. Consequently N/g cannot equal7 or21.

Proof: (1) is homogeneous and linear in the two gadget coordinates. The
source is unchanged, and the replacement gadget is explicitly physical.
This supplies the witness, not merely divisibility. The conclusion does
not hold by this argument when only the normalized pair of integers exists.

## 2. A general seven-cofactor obstruction

**Theorem.** For every proper graft and g=gcd(H,Q),

    N/g is odd, and N/g != 7.                            (2)

If B has at least two vertices, also N/g>=5. The singleton gadget has g=1
and may give N=3, so this last qualification is necessary.

First, H is odd. For |B|=n>=2, complement pairs partition the subsets into
2^(n-1) pairs. Every H(B[X])H(B[X^c]) is odd by Rédei parity, including empty
parts. Hence Q is twice a sum of an even number of odd integers, and

    4 divides Q.                                        (3)

As g divides odd H, g is odd, so (3) survives division by g. Put h=H/g,
q=Q/g. Then h is positive odd, q>=2h, and 4 divides q. Formula(1) gives
N/g>=h+q>=5. If N/g=7, these inequalities force exactly

    (h,q)=(1,4).                                        (4)

This is the actual pair of the two-vertex transitive tournament TT2. The
normalization principle would produce an actual tournament with count7,
contradicting the inherited exclusion. If B is a singleton, g=1, and the
same inherited exclusion directly applies. This proves(2).

There is also a proof at the local physical-inequality level, not using the
global exclusion for the forced pair(4). Its normalized equation is

    7=d+4(D-E)+3f.

If D>=3, then D-E>=(D+1)/2>=2, so the right side is at least9. If D=1,
then d=1 and E=0, so the right side is5+3f, with f even; it cannot be7.
This isolates the mechanism: complement parity forces the normalized TT2
interface, and physical endpoint parity excludes the prospective source.

**Corollary (seven-prime gate).** If N=7p for an odd prime p!=7, every proper
graft representation satisfies

    gcd(H,Q) in {1,7}.                                  (5)

Indeed g|7p, while H<=N/3 rules out g=7p, and(2) rules out g=p. For p=7,
the same final set follows simply from g|49 and H<=49/3. The case p=3
has no physical output at all, by the inherited exclusion of21.

This improves the generic statement g|N on precisely the composite atomic
lane. Ordinary prime targets force g=1 by the old H<=N/3 bound. The new
point is that the OTHER prime p cannot be hidden in gadget content while
leaving the forbidden cofactor7. It does not prove that an allowed-content
representation exists.

## 3. The surviving content-seven branch is real, not a formal remainder

Use seven-vertex gadget B with out-neighbor bit rows

    [22,60,8,81,4,29,55].

Exact counting gives

    (H,Q)=(49,812), gcd(H,Q)=7, normalized pair=(7,116).

The gadget is strongly connected. Its normalized pair is NOT physically
realizable: its first coordinate is the forbidden7. Thus generic division
by gcd cannot be substituted for the normalization principle.

Take A=C3 rooted at any vertex. It has

    (D,E,U,d)=(3,1,3,1), f=0.

The literal proper graft is strong and has

    N=49+2*812=1673=7*239,

with239 prime. Therefore content7 in(5) occurs in an actual representation, and banning
all nonprimitive gadgets would discard actual seven-prime constructions.
This does not prove that their scalar targets lack alternative primitive-gadget
representations.
No novelty of the attained integer1673 itself is claimed.

A different, independently counted eight-vertex rooted seed C is

    rows=[126,180,120,2,168,72,154,45], root=3;
    (h,R,e,d)=(259,266,70,77), gcd=7.

It is strong, and h=7*37. Here R=U, and e=E, but the coordinate order differs
from(1). All fixed-port S/B descendants have content divisible7, by the
inherited four-test theorem. Adding a sink and then a bridge gives

    H(C_SB)=h+2R=791=7*113,

again strong, with113 prime. The final bridge gadget is a singleton with
gcd(H,Q)=1, so this also exhibits the g=1 branch of(5). **Rooted source
content7 and gadget content7 are different invariants**: C itself has
(H(C),Q(C))=(259,4456), whose gcd is1.

For any finite fixed-port seed bank, let c_i=gcd(h_i,R_i,e_i,d_i). If a tail
has H=7p and p>max_i c_i is prime, then its seed must have c_i=1 or7,
since c_i|7p. For prime outputs p above the bank contents, only c_i=1 is
possible. The C example shows why the seven-prime lane has a genuinely
additional permitted root-content stratum. This finite-bank filter is a
necessary condition, not a target coverage theorem.

## 4. Auxiliary size control and the infinite-descent comparison

Every Hamilton path of an n-vertex gadget can be cut at n+1 positions into
an ordered covering pair. The path and cut can be uniquely recovered by
concatenating the ordered pair and reading its first length. This injection
proves the elementary bound

    Q>=(n+1)H, and hence N>=(n+2)H.                      (6)

This sharpens a prospective inverse's order-dependent scalar budget. It
is a cut injection, not a classification of equality, nor a new general
path-cover theorem. We make no priority claim for this elementary bound.

The exact congruent-number arithmetic is

    41²-31²=49²-41²=720=5*12².

Dividing those squares by12² preserves an equal-gap square AP, with common
difference5. The associated rational right triangle has legs3/2 and20/3,
hypotenuse41/6, area5. Here the map and its preserved predicate are explicit.
This witness proves existence for area5; it is not itself a decreasing map
on a set of counterexamples.

The useful transfer to Chan–Pak is the **normalization obligation**: division
must keep the object in the target category. Our theorem identifies one
case where arithmetic normalization is forced to be physical (TT2), and
one sharp surviving case where it is not (normalized gadget pair(7,116)).
Unlike scaling a rational triangle, dividing a tournament's two integers
need not yield a tournament. Nor does dividing a realized square count
preserve realization:49 is realized but7 is not. Ordinary square-class
normalization of values is therefore not a valid tournament descent rule.

A least missing allowed target is inherited to lie on p or7p after the
finite carry atoms are supplied. On p use only primitive proper-graft
gadgets; on7p use only gadget contents1,7. What remains OPEN is existence
of a successful full-state physical predecessor in those restricted lanes.
The filtered, terminating inverse is still a decision procedure for its
image, not a proof that the image contains every target. This is where
an honest minimal-counterexample proof must supply new mathematics.

## 5. Reproduction and finite scope

Run from this directory:

    python3 -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -B -O verify.py | cmp - CHECK.json

The executable uses explicit failure gates, not Python assert. Its
`--fail` option must exit nonzero in both modes. CHECK.json records the
consequence objects: actual strong graft values, primes, full state and
the empty source-solution set for normalized quotient7.

The finite universe for the parity/cut regression is **every labelled
tournament on1,...,5 vertices**, totaling1099 objects. The two seven/eight-
vertex seed witnesses are checked both by subset DP and by independent
literal permutation enumeration over all induced subsets. The two final
10-vertex graphs are counted by literal subset DP and checked strong;
the results are compared with the independently evaluated transfer formula.
No randomized scouting output is used as a universal negative certificate.
The all-size statements above have proofs independent of the finite census.
