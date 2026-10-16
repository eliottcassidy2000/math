# Golden elliptic curves as marked extension counters

2026-10-16. **PROVED scoped:** the golden model and twist identities, the
signed three-anchor extension formula, and the stated bounds and losses.
**CITED:** the classical Legendre invariant, Hasse bound, and the
Paley/elliptic character-sum method. **FINITE-EXACT:** the independent
enumerations below. **OPEN:** the planar Jacobian conjecture, universal
Collatz convergence, and the plane six/seven alternative. No priority is
claimed for the classical arithmetic ingredients.

## 1. Inheritance and the live connection

The closest proved bridge is
[THM-4602, tournaments as chirotopes and Legendre frieze strips](../../01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md).
It already retains the quadratic-character orientation gauge, and separates
real total positivity from finite-field quadratic residues. The newer
[THM-4653, the Paley K6-minus-P6 Legendre reduction](../../01-canon/theorems/THM-4653-paley-k6-minus-p6-terminal-reduction-to-a-legendre-bilinear-form.md)
likewise makes an elliptic sum an actual pattern-counting object. Its
stronger modular identification remains a separate hypothesis.

The corrected near miss is to identify two configurations by a shared
modular invariant and then infer identical extension support. The hostile
below has the same elliptic j-invariant and the same ordered three-vertex
tournament, but one configuration admits a common out-neighbour and the
other does not. The least-used necessary sidecar is the quadratic twist
together with the actual ordered roots.

Board: golden roots; marked elliptic model; quadratic-character relation;
full extension support; unit-distance realization; differential and native
clock. The source is an ordered field triple. The target is its eight
labelled extension counts. The map is a character expansion. It preserves
every three-port existence predicate and count, but does not preserve
metric geometry or the history of a Collatz source. Those require their
own exact maps, examined in the companion packets.

## 2. One concrete golden curve

Put phi=(1+sqrt(5))/2. For lambda^2=lambda+1, the Legendre curve

\[
 E_\lambda:v^2=u(u-1)(u-\lambda)
\]

has

\[
 j(E_\lambda)=256\frac{(1-\lambda+\lambda^2)^3}
 {\lambda^2(1-\lambda)^2}=2048.
 \tag{1}
\]

Indeed the numerator's inner factor is 2 and
lambda(lambda-1)=1. A rational model in this geometric isomorphism class is

\[
 E_0:Y^2=(X+2)(X^2-X-1)=X^3+X^2-3X-2.
 \tag{2}
\]

Its invariants are c4=160 and discriminant=2000=2^4*5^3, giving j=2048.
It has the golden conjugates and -2 as its three finite branch points.
All Legendre parameters with this j are displayed by the exact factorization

\[
 (1-\lambda+\lambda^2)^3-8\lambda^2(1-\lambda)^2
 =(\lambda^2-\lambda-1)(\lambda^2-3\lambda+1)
   (\lambda^2+\lambda-1).
 \tag{3}
\]

Thus the six cross-ratio parameters lie in Q(sqrt(5)); this is a concrete
arithmetic link, not an inference from the number 2048 alone.

There is an especially relevant independent occurrence. Francis and
Whitcher's [five-colour hypercube Adinkra paper, §3.1](https://arxiv.org/html/2410.11137)
gives five elliptic factors with parameters phi+1, phi and -phi, all of
j-invariant 2048. Their graph-height moves become signed additions of
specified divisors. This is a useful model of a combinatorial operation
with an authenticated arithmetic image. Its Jacobian is a Jacobian
*variety*, not the Jacobian conjecture. No reduction between those two
meanings is asserted here.

### Twists must survive the translation

Let s^2=5 and lambda=(1+s)/2. The change

\[
 X=\lambda-su,\qquad Y=-sv
\]

turns (2) into

\[
 v^2=-s\,u(u-1)(u-\lambda).                 \tag{4}
\]

It is therefore a quadratic twist of the displayed Legendre equation.
For an odd prime p different from 5 with 5 a square, write
a_p(E)=p+1-#E(F_p), and let chi be the quadratic character, chi(0)=0.
Then

\[
 a_p(E_0)=\chi(-s)a_p(E_\lambda),\qquad
 a_p(E_{1-\lambda})=\chi(-1)a_p(E_\lambda). \tag{5}
\]

The second identity follows from u -> 1-u, which negates the cubic.
It can reverse the trace even though j is unchanged. At p=5 the rational
model is singular and s=0; (4) cannot be used there.

There is also a type distinction. The relation a -> b iff chi(b-a)=1
is a tournament precisely when chi(-1)=-1. For prime fields containing
the golden roots this requires p=11 or 19 modulo 20. If 5 is nonsquare,
the golden roots instead lie in F_(p^2), where -1 is a square: the same
quadratic-character rule gives an undirected graph. In particular the
F_9 golden construction is not a Paley tournament.

For such inert odd primes the base-change and twist identities give

\[
 a_{p^2}(E_\lambda)
 =\chi_p(-5)\big(a_p(E_0)^2-2p\big).        \tag{6}
\]

Here Norm(-s)=-5, so the character of the twist over F_(p^2) is chi_p(-5).
This is a genuine squaring recursion with its correction and twist retained.

## 3. Elliptic trace measures all eight extension patterns

Take distinct ordered anchors a_1,a_2,a_3 in F_p and signs
epsilon_i in {-1,1}. Let

\[
 N_\epsilon=\#\{t\notin\{a_1,a_2,a_3\}:
                     \chi(t-a_i)=\epsilon_i\ (i=1,2,3)\}.
\]

Set A=a_p(y^2=prod_i(x-a_i)) and

\[
 B_\epsilon=\sum_{x\in\{a_1,a_2,a_3\}}
                 \prod_{i=1}^3(1+\epsilon_i\chi(x-a_i)).
\]

Then the complete extension signature is

\[
 \boxed{8N_\epsilon=p-\sum_{i<j}\epsilon_i\epsilon_j
             -(\epsilon_1\epsilon_2\epsilon_3)A-B_\epsilon.}
 \tag{7}
\]

*Proof.* Expand the product and sum over F_p. Each singleton character
sum is zero, each distinct-root pair sum is -1, and the triple sum is -A
by the elliptic point-count formula. Subtract the three anchor contributions
before dividing by 8. Outside the anchors the product is precisely 8 on
the requested cell and zero elsewhere. This also proves that the eight
counts sum to p-3. The graph version of this method and its use of Hasse
are classical; see [Brouwer and Martin, Triple intersection numbers for
the Paley graphs](https://aeb.win.tue.nl/preprints/p3g.pdf).

Each B_epsilon belongs to {0,4,8,12}. Hence Hasse's |A|<=2sqrt(p) gives

\[
 8N_\epsilon\ge p-15-2\sqrt p>0\quad(p>25). \tag{8}
\]

Thus every signed three-port pattern occurs for every odd prime p>25.
This includes both the graph and tournament cases, with their distinct
interpretations. It is a local extension theorem, not a Hamiltonian-path
balance theorem or a Euclidean unit-distance realization theorem.

For ordered transitive tournament anchors, B_(+++)=B_(---)=4. In that
case the common out- and in-neighbour counts simplify to

\[
 N_{+++}=(p-7-A)/8,\qquad N_{---}=(p-7+A)/8. \tag{9}
\]

The golden triple is always transitive when p=3 modulo 4: the identity
lambda(lambda-1)=1 gives chi(lambda)=chi(lambda-1), so the third point is
either before both 0,1 or after both. Its abstract tournament therefore
has very little information; the trace in (9) supplies the missing count.

### Minimal retained-data hostile

At p=11 compare the ordered triples (0,1,4) and (8,0,1). Both are
transitive in the displayed order, both have j=2048, and their traces are
-4 and +4. Their common out-neighbour counts are respectively **1 and 0**.
The existence predicate changes. A model retaining only j or only the
three-vertex tournament cannot be a lossless interface for this task.
Keeping the ordered roots and twist repairs the failure exactly by (7).

## 4. What this supplies to the three proof frontiers

The companion packets test actual consumers rather than infer a reduction
from a shared elliptic vocabulary.

| Frontier | Positive transfer | Required data and surviving obligation |
|---|---|---|
| [Planar Jacobian](jacobian_elliptic_bridge_20261016c.md) | Every finite j has a smooth period-five Lyness fibre after a parameter extension of degree at most 12, with an explicit Weierstrass map and differential scale. | At j=2048 the Lyness parameter has exact degree 12, although the Legendre parameter has degree 2. Matching j forgets the torsion marking, twist, punctures and differential. The Lyness logarithmic two-form cannot become an exact rational coordinate area form, even after finite algebraic pullback. |
| [Hadwiger--Nelson](hn_elliptic_bridge_20261016c.md) | The inherited positive-rank elliptic contact family gives exact unit geometry; a cyclotomic eleven-cycle attachment escapes its old field with a fully computed boundary relation. | It multiplies every old six-colour extension fibre by 1,627,604 while deleting none. The new theorem preserves every old six-colouring under at most 17 exterior vertices in the specified degree-five extension. A forcing connector there needs at least 18 coupled exterior vertices. |
| [Collatz](collatz_curve_bridge_20261016c.md) | A source-marked elliptic coordinate lift faithfully stores every supplied native word and commutes with its exact signed-terminal projection. | It is a representation, not new descent. All fixed-source models are isomorphic, and even infinitely many grounded sources carry nontorsion points. A new observable must prove an actual source payment or signed-energy bound while retaining the carry, native cylinder and clock. |

For the Jacobian lane, the degree-2 versus degree-12 contrast is an exact
example of representation cost: choosing labelled two-torsion roots and
choosing the particular period-five Lyness marking are different problems.
Neither costs nor the preserved differentials can be read from j alone.
The same packet explains the five-torsion marks through an exact Frobenius
action at three: orbit lengths are 6 for cyclic subgroups, 12 for points
modulo sign, and 24 for nonzero points. These counts retain different data.

The geometric bound has an equally explicit small structure. Its excluded
sixteen-vertex core splits into layers of sizes 1, 5 and binomial(5,2)=10.
The forced planar rhombi yield an impossible rank-five Gram matrix. The
same abstract graph has an exact rational unit-distance realization in
five dimensions. Thus the geometric dimension is a necessary part of
this extension test; the incidence pattern alone does not carry the proof.

Incoming [elliptic cross-edge work, G7](elliptic_crossedges_20261016c.md)
then supplies a new mixed unit direction after a real cubic extension.
The geometric companion proves that the seventeen-vertex bound survives
every finite odd-degree extension of a conjugation-stable old field
without sqrt(-3). In particular it survives this new degree-fifteen tower.
The resulting concrete search target is a coupled configuration of at
least eighteen exterior points using the new direction, with full unit
incidence and boundary support authenticated. The direction itself does
not assert such a connector exists.

For the Collatz lane, Chang's golden Möbius recurrence does admit a genuine
cyclic complex-torus quotient. Replacing that cyclic scaling by both 2 and
3 makes the positive scaling group dense and the quotient non-Hausdorff.
The companion gives the proof. Discarding more trajectory coordinates does
not repair this quotient. One possible alternative is a discrete marked
cover, which would still have to price the original integer trajectory.

The reusable working principle is therefore concrete: transport the
*operation and its observable*, then audit what the quotient erased.
Elliptic trace transports three-anchor counts; geometric field equations
transport complete colouring support; marked differentials transport
area-form information; exact native coordinates transport Collatz payment.
Only a further map preserving the required predicate would turn one of
these bridges into a proof of the corresponding universal conjecture.

## 5. Exact controls and reproduction

[Script](../../04-computation/experiments/golden_legendre_ports_20261016c.py)
and [receipt](golden_legendre_ports_20261016c.json).

The universe is all 1,109 normalized triples (0,1,lambda), lambda!=0,1,
over odd primes through 101. For each, direct character-cell enumeration
is compared with (7), including all eight patterns. Affine changes include
the character orientation gauge. Literal x,y point enumeration supplies
an independent trace check for p<=19. The golden specialization checks
22 parameters, both twist signs, exact coordinates, and the moduli
factorization (3). Inert cases p=3,7,13,17,23 are checked by explicit
quadratic-field arithmetic. Altered traces, root order, boundary terms,
boolean counts and invalid source types are rejected by the receipt audit.

```text
python -B -X utf8 04-computation/experiments/golden_legendre_ports_20261016c.py
python -B -O -X utf8 04-computation/experiments/golden_legendre_ports_20261016c.py
```

Normal and optimized runs give 7,892 exact checks and identical
LF-normalized JSON, SHA256
`4e2001dfa2eadb4d7e954f23e474a024662fa76888c5ad07264d02afd4f30ea7`.
The family-algebra and Collatz-memory agents independently audited the
proof, source/type guards, and optimized replay against the saved output.
The latter also enumerated 815 ordered transitive triples at p=3,7,11,19
as an independent check of (9). Proof identities and Hasse imply the
unbounded scoped statements; the finite census alone does not.
