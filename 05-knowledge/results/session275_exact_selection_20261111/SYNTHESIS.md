# Beyond a universal grammar: physical arithmetic controllers and count corridors

**Session275, 11 November 2026. Chan–Pak remains OPEN.**

This session went beyond Session274's geometric generation theorem by asking
which choices can actually select counts, how they interact, and which extra
state makes a failed construction usable. It produced new all-length physical
controllers, an exact classification of a normalized hole, and a concrete first
step in a different connected-corridor route. It did not produce an eventual
threshold or a proof that all allowed odd counts occur.

## 1. Independent physical digits, not just equal counts or a gcd

Two explicit six-vertex pointed tournaments, one with5 paths and one with9,
have the same root-deletion and receiving coordinates but different sending
coordinates. Their4×4 fusion matrices differ only in the first row.

This gives the exact algebraic law

    D_left C D_right = 0

for two such same-direction trade differences and EVERY intervening physical
fusion context C. Consequently arbitrary independent binary choices have NO
mixed interactions, even with arbitrary fixed spacers. They form real
tournaments at every word length, not formal signed solutions.

For a length-m pure word, the final count is exactly a baseline plus a subset
sum of explicit weights. The weights are superincreasing, so there is a complete
greedy inverse: it either constructs the requested member of this family or
certifies that the target is outside this family. The generating polynomial
factors as a product of binary digit polynomials. This is an all-size exact
selector, but for a sparse language, not all odd integers.

The investigation did not stop at its first obstruction:

- The first H-preserving pair had a growing factor3^m. Allowing the inputs'
  counts to differ removed that inherited common-prime debt.
- The resulting controller still had a fixed difference gcd4. Adding a second
  physical one-way bank and a terminal triangle reduced the gcd to exactly2
  at EVERY length, while keeping independent choices.
- The allowed digits remain too widely spaced. Removing congruence obstacles
  is not the same as obtaining overlapping integer intervals.

There is also a proved obstruction to the cheapest repair: a nontrivial
one-way trade cannot have root-deletion count1, because the remaining receiving
coordinate encodes the complete transitive root-incidence word.

[Proof, literal seeds, exact decoder and gap theorem](switches/PROOF.md)

## 2. Turning the controls around gives a paid multiplication gate

Taking converses reverses the direction of the trade. Same-direction controls
still do not interact; a forward control followed by a reverse control does.
The interaction is explicit and positive. Every third and higher mixed
interaction vanishes, with arbitrary baseline spacers retained.

For example an actual edge-spaced construction implements

    H(x,y)=125+80x+100y+72xy,  x,y in{0,1}.

More generally two independently encoded root-response integers u and v,
separated by a nontrivial physical spacer S, have mixed response Q_S uv.
The baseline and linear costs are retained, not subtracted as though they
were tournaments. This supplies a genuine constructive multiplication
ingredient in the same reusable interface algebra.

The next target is now precise: an interval-overlap theorem for these paid
quadratic forms, or a compatible offset bank that controls their baseline and
linear terms without losing the rooted state. Neither the multiplication gate
nor its primitive coefficient gcd by itself supplies that theorem.

[Quadratic controller, parity repair and multiplication](switches/EXPANSION.md)

## 3. Two-sided seams have an exact divisor decoder and a renewable +2 pair

Repeated source and sink seams at a fixed root have the all-length form

    H(m,n)=q·2^(m+n)+beta·2^m+gamma·2^n+delta,  m,n≥1.

Four physical observations recover the packet. Its normalized determinant
kappa=(q·delta−beta·gamma)/q is unchanged by continuing either fan.
Prescribing N becomes the exact factor equation

    (q·2^m+gamma)(q·2^n+beta)=qN−(q·delta−beta·gamma).

The factors have explicitly checkable dyadic forms. This is a terminating
exact inverse problem, not just a residue test. The zero-fan initial cases
are separate, and a finite fixed bank is provably sparse.

Changing the marked root changes the packet. A single explicit seven-vertex
source has two balanced-fan outputs that differ by exactly2 at every length:

    168·4^m−241·2^m+109,
    168·4^m−241·2^m+107,                 m≥1.

Thus the construction yields297/299,1831/1833, and an unbounded family of
physical neighboring odd counts. A cut switch plus a bounded number of old
edge repairs transports the root. This is NOT yet a low-count-jump path
between successive scales; its intermediate repairs need not have small
count changes.

The next extension preserves whole families, not just two endpoints. Apex-cut
families with the same two boundary-flow charges acquire the SAME additive
count translation under balanced fans, preserving all their old-edge count
differences. Every vertex-transitive deletion core has all ports balanced, so
its entire apex-cut cube satisfies this condition. The actual61/63 corridor edge
lifts to367/369,1747/1749,7579/7581 and onward, still differing in the same
single old free edge. The translated corridor does not grow in count width;
inter-scale overlap and attachment to the initialized component remain open.

A hostile example also shows that identical ENTIRE one-sided scalar fan
sequences can conceal different two-sided menus. The root packet, not just
its past outputs, must be retained.

[Two-sided selector, boundary cases and root transport](seams/PROOF.md)

## 4. Normalized21 really occurs, and its physical lift can escape

The prior proper-graft theorem excluded normalized cofactor7 but left21 open.
This session resolves that particular question exactly:

    N/gcd(H(B),Q_B(2)) = 21

holds precisely when Q_B(2)=20H(B) and the source is transitive with its root
penultimate. Every such gadget has at most19 vertices. An explicit seven-vertex
gadget has(H,Q)=(29,580), giving an actual strong nine-vertex tournament with
H=609=21·29. The normalized packet(1,20) is NOT physical, so division by29
would be an invalid constructive descent.

The next experiment retained the graph and moved its root instead. An internal
root gives an actual H1567 output, escaping factor29. Its complete renewable
fan packet is

    (D_k,E_k,U_k,d_k)=(825·2^k−83,241,825·2^k+241,241), k≥1.

Every exact29-adic valuation occurs infinitely often in this real family,
with a unique explicitly lifted bad exponent class at each modulus. The
valuation controller does not cover arbitrary absolute counts, but it proves
that the content trap can be escaped physically, without division and without
losing the state required for continuation.

[Exact normalized21 classification](descent/PROOF.md) ·
[Physical escape and full-state renewal](descent/ESCAPE.md)

## 5. A different route: make one count corridor unbounded

Instead of generating every tournament, connect physical tournaments by moves
whose counts change by0 or2. A connected component then visits every odd
integer between any two of its count values. Therefore an unbounded-count
component containing the explicit23 seed would prove all odd counts≥23.
Unbounded vertex count is not enough; neutral padding already gives that.

The experiment uses fixed-path tile flips, cut gauges, complement and neutral
reframing. Through order7, the combined moves connect every whole maximal
consecutive-odd spectral run. Reframing is decisive: at order7 there are111
fixed-chart components but only10 after neutral frame changes. Single-tile
moves alone leave159 reframed components. This is a structural comparison of
operations, not a larger count-coverage census.

There is a constructive first growth step. The new63 birth structure is joined
by count steps63→61→59 to a neutrally padded seven-vertex graph, and by63→65
to another. Every move has a common Hamilton path, so it is a legal free-tile
move in a suitable frame. It therefore stitches the previous23..61 and65..105
corridors into one variable-order component through63. This repairs the first
birth obstruction by reconfiguration before neutral deletion.

The OPEN continuation problem is to make these corridor extensions recur at
all heights. The renewable high-count+2 siblings in Section3 are possible
bridge ingredients, but no attachment of them to this component is proved.

[Corridor criterion, exact finite structure and63 portal](landscape/PROOF.md)

## 6. What unifies the advances, and what still does not follow

Three concrete principles emerged:

- **One-way nilpotent control:** arrange physical differences so their mixed
  products vanish; reversing one information flow creates a controlled product
  rather than an uncontrolled interaction cascade.
- **Lift-sensitive descent:** normalize arithmetic only while keeping the
  integer and rooted data that make the normalized state physically realizable.
- **Count corridors instead of universal generation:** geometric completeness
  may be much stronger than needed; unbounded low-jump continuation of one
  correctly initialized component would suffice.

These are names for the mechanisms above, not new universal meta-theorems.
Random and deliberately distant repository sources supplied useful changes of
view: old Morse landscapes, tensor parity, sharp context divisors, interval
contexts, Berggren causal horizons, Rule30 moving observers, and Zeckendorf
carry/decoder work. Each connection's mathematical use and limitations are
recorded in its lane. Historical claims were not re-audited wholesale.

Concurrent S93 work arrived during the session with a DRAFT3/4 counting-exponent
argument. It suggests a complementary additive-energy lens for the quadratic
selector. Its reported status is not promoted here, and many distinct outputs
would still not give exact prescribed-target coverage.

See [research board](RESEARCH_BOARD.md), [remaining obligations](OPEN_QUESTIONS.md),
[reproduction](REPRODUCE.md), and [focused review](review/REVIEW.md).
