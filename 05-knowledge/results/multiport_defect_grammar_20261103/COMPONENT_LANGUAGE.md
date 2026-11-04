# Exact acceptance in the singleton/triangle component language

2026-11-03. **PROVED:** a complete finite inverse and word-count formula for
every fixed positive marginal f in the declared atom language.
**FINITE-EXACT:** all first three levels of the common-clock +2/+4 targets
from H_t=-5+7*2^t. The language accepts both first letters, but only 22 and
42 at level two and only 422 at level three. **NOT CLAIMED:** a native
program, a physical transformation from any preceding source, unrestricted
target nonexistence, or an arbitrary-depth global classification.

The inherited mechanism is the seven-coordinate
[ordered-component transfer](../session182_physical_returns_20261104/induction/PROOF.md),
read at `b486271d24`. The closest positive control is the
[two-digit entry](transport/TWO_DIGIT_ENTRY.md). The hostile examples here
are the first missing words 24 and 44. The least-used sidecar is the joint
correction j: the complete inverse fixes f but still has to test h-j.
This is an acceptance computation in a realized component language, not an
inference from equal Hamilton counts to interchangeable physical sources.

## 1. Actual atoms, orientation, ports, and signatures

Use ordered sums of these nonempty atoms, every vertex of an earlier atom
beating every vertex of a later atom:

- T_b is a singleton with port bit b in {0,1}.
- C_k is a directed triangle with k ports, for k in {0,1,2,3}.

The final atom must have at least one port. This is exactly the condition
that the port set meet the suffix's terminal strongly connected component.
Consequently attaching any nonempty transitive prefix and the usual root
produces an actual strong tournament. Each triangle's vertex orientation
is cyclic; its three rotations are transitive on the subsets of each
cardinality k. Thus all actual k-port choices have the same signature.
The words classify those choices up to internal rotations; arbitrary
physical names and the required transports of those names remain extra
data. The checker uses consecutive names and chooses the first k names
of each triangle as a canonical actual representative.

In coordinates sigma=(h,p,f,g,j,e,s) of the cited transfer, the atoms are

    T_b: (1,2,b,b,0,b,b),
    C_k: (3,12,3k,3k,k(k-1),k,k).                       (1)

The directed triangle has one Hamilton path ending at each vertex and one
starting there. For a pair of distinct ports, the two-slot joint count is
two in the orientation of their triangle edge and zero in the reverse
orientation; hence the total joint correction is k(k-1). These facts give
all seven entries in (1), also checked independently by literal subset DP.

For consecutive components with signatures (h,p,f,g,j,e,s) and
(H,P,F,G,J,E,S), ordered assembly gives

    (hH, pP, fH+pF, gP+hG,
     jH+hJ+eS+gF, hE, sH).                             (2)

The cross term gF is retained. All outputs in the receipt are realized by
this literal operation and then independently recounted from adjacency.

## 2. Complete inverse for a fixed marginal f

For a nonempty component word write c for its number of triangle atoms
and a for its number of singleton atoms. Then

    h = 3^c,       p/h = 2^a*4^c.

Set q=f/h and rho=p/h. Appending an atom with radix r=2 or 4 and port
digit k updates

    q_new = q + rho*k,       rho_new = rho*r.            (3)

Therefore q is an integer in mixed radices: a singleton supplies a base-2
digit 0 or 1, and a triangle supplies a base-4 digit 0,1,2, or 3. The first
component supplies the lowest digit. The final digit is nonzero by the
terminal-port condition.

**Exact inverse theorem.** Fix a positive integer f. For every integer
c >= 0 such that 3^c divides f, put q=f/3^c. Enumerate all mixed-radix
expansions of q with exactly c radix-4 positions and a nonzero highest
digit. Replace radix-2 digits by T_b and radix-4 digits by C_k. This lists
every admissible atom word with marginal f, and no other words.

Proof. Necessity follows from (3) and h=3^c. Conversely each such expansion
has the prescribed Hamilton count factor h and marginal hq=f, and its
nonzero final digit satisfies the exact strongness condition. The literal
ordered sum realizes it. For a chosen first radix r, its low digit must be
q mod r and the remaining word must represent floor(q/r); the triangle
budget decreases exactly when r=4. Thus recursion over the two possible
first radices is complete and introduces no duplicates. If the quotient
is zero, stop exactly when the low digit is positive and the remaining
triangle budget is zero. In particular q=1 with one triangle remaining
has the valid one-atom expansion C_1: a final radix-4 digit may use a zero
padding bit above the top binary bit. Excluding that case would already
lose the initial uniform +2 witness. QED.

The inverse is finite with no arbitrary bound on tournament order. There
are only v_3(f)+1 possible values of c, and each division reduces a positive
remaining q by a factor of at least two until it vanishes. No additional
high zero atoms are allowed: those would violate the terminal-port guard.

There is also an exact all-order size formula. Put
L_c=bitlength(f/3^c). For fixed c, a radix-4 position joins two neighboring
binary positions; its start may be the final bit, using the single padding
bit just described. Choose c nonconsecutive starts from positions
0,...,L_c-1. The number of choices is binomial(L_c-c+1,c), interpreted as
zero when its top entry is smaller than c. Hence the number of admissible
words at fixed f is exactly

    N(f) = sum_(0 <= c <= v_3(f)) binomial(L_c-c+1,c).     (4)

The subset-of-domino-starts construction is also an independent inverse
implementation in the checker. Equality with the quotient recursion and
with (4) is checked for every f used below. Formula (4) counts actual
component words modulo each triangle's internal rotations; it is not a
count of native programs or of rooted physical transformations.

## 3. Exact acceptance through three levels

For any realized suffix signature, the full common-clock attachment has

    H_t = h-j+2^t f,       t >= 1.                      (5)

Start with T_1 T_1 T_1, whose (f,h-j) is (7,-5). For a digit word
d_1...d_l in {2,4}, a uniform return target at the same prefix clock must
have

    f_l = 7*3^l,
    a_0=-5,        a_i=3a_(i-1)+d_i,
    h_l-j_l = a_l.                                     (6)

The exact inverse finds every possible atom word at the prescribed f_l,
then applies the full signature law to test the intercept a_l. There are
4, 7, and 13 candidate words at levels one, two, and three, respectively.
All are included in `component_language.json`, with actual rows, port
sets, signatures, and the target-to-witness index lists.

| Level | Target f | Accepted digit words | Missing digit words |
|---|---:|---|---|
| 1 | 21 | 2, 4 | none |
| 2 | 63 | 22, 42 | 24, 44 |
| 3 | 189 | 422 | 222, 224, 242, 244, 424, 442, 444 |

The accepted witnesses, using the atom notation of Section 1, are:

| Digits | Actual suffix word | sigma=(h,p,f,g,j,e,s) |
|---|---|---|
| 2 | T_1 T_1 C_1 | (3,48,21,39,16,1,3) |
| 4 | T_1 T_0 T_1 T_0 T_1 | (1,32,21,21,12,1,1) |
| 22 | T_1 T_0 T_1 C_2 | (3,96,63,66,40,2,3) |
| 42 | T_1 T_0 C_1 T_1 | (3,96,63,57,34,3,3) |
| 422 | T_1 T_0 C_1 C_1 | (9,576,189,333,100,3,9) |

The complete intercept sets at the successive levels are

    f=21:  {-15,-13,-11},
    f=63:  {-61,-43,-39,-37,-33,-31,-19},
    f=189: {-183,-181,-139,-137,-129,-91,-69,-57,-49}.

Thus the first missing target intercepts -35 for 24 and -29 for 44 are
excluded by a complete fixed-f classification in this atom language,
including every internal triangle port subset and every word length.
This is stronger than failing to find them within an order cutoff.

Acceptance still means that a fresh actual word has the required uniform
count; it does not construct a map from the preceding physical tournament.
The cases 22 and 42 cannot be treated as reusable operations merely because
their target signatures are present. Conversely, a missing word is not a
missing Hamilton count in arbitrary tournaments. Allowing other strong
components, another prefix clock or clock-dependent output, or a different
attachment model is outside this fixed-f test. In particular the present
negative conclusions concern the common prefix clock in (6); they are not
claimed invariant under output-clock shifts.

## 4. Independent checks and the next structural obligation

Run:

    python 05-knowledge/results/multiport_defect_grammar_20261103/check_component_language.py
    python -O 05-knowledge/results/multiport_defect_grammar_20261103/check_component_language.py

Both modes use explicit RuntimeError guards. They write the same
deterministic `component_language.json` receipt. The inverse recursion is
compared with the separate binary-domino enumeration and exact word count.
For all 25 candidates including the initial source, a literal directed
subset-path DP reconstructs all seven coordinates without using (2).
For each, actual attachments at t=1 and t=2 independently pass strongness
and Hamilton-count checks, giving 50 full-attachment controls. The final
one-bit C_1 boundary is an explicit positive control. Normal and optimized
stdout and JSON receipts are byte-identical.

The smallest new obstruction is already a joint-intercept obstruction at
the second +4 return, not an inability to reach the required marginal f.
The seven-coordinate carrier evaluates these words exactly, but its
realizable atom set does not supply the needed intercepts. A larger
accepted-source grammar must introduce another operation, atom, or clock
scheme and then separately prove physical source compatibility. This
calculation identifies the first precise equations such an extension has
to repair; it does not assert that any particular extension will do so.

A concrete next atom bank is the familiar strong four-vertex tournament
with adjacency rows [8,1,3,6]. Independent literal counting gives h=5,
p=32, and one-port f values [7,9,6,5]. In particular p/h=32/5 and three
of those one-port ratios f/h are nonintegral. Such an atom leaves the
integral mixed-radix restriction proved above. The next target is its full
port bank inside ordered assembly, with complete acceptance and physical
compatibility kept separate. This is a specified extension to test, not a
claim that it supplies the missing second +4 returns. No search in that
extended language is performed in this receipt.
