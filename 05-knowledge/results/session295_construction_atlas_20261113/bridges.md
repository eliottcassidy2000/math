# Continuation equivalence is a typed quotient, not a scalar invariant

## An exact finite quotient for uniform substitutions

Let $F_k(C)$ count the ordered $k$-run spanning path covers of the
nonempty tournament $C$; put $D_k(C)=(k-1)!F_k(C)$, and put $D_k=0$ for
$k>|C|$. Use the Session 293 excursion-factor and alternating-separator
probe theorem:
[proof](../session293_cycle_returns_20261113/control/MODULE_IDEAL.md).
An exterior is *uniform* if each of its vertices sees every vertex of $C$
in the same direction. For every exterior $E$,
$$
H(E[C])=\sum_{k\ge1} b_{E,k}D_k(C),\qquad b_{E,k}\in\mathbb Z_{\ge0}.
$$
The alternating separator with $2m$ exterior vertices has coefficient
one on $D_{m+1}$ and no higher coefficient.

**Proposition (exact contextual quotient).** For any integer $q\ge2$,
the following are equivalent, without an order bound:

1. $H(E[C])\equiv H(E[C'])\pmod q$ for every uniform exterior $E$.
2. $D_k(C)\equiv D_k(C')\pmod q$ for every $k<K_q$, where $K_q$
   is the first $k$ with $q\mid(k-1)!k!$.
3. The $K_q-1$ alternating separator probes for $m=0,\ldots,K_q-2$
   agree modulo $q$.

Proof: write $D_k=(k-1)!k!\,pc_k$, so all $k\ge K_q$ vanish.
The exterior expansion proves (2) implies (1), inclusion of the probes
gives (1) implies (3), and integer unitriangular inversion, valid even
over the nonfield $\mathbb Z/q\mathbb Z$, gives (3) implies (2). Thus at
most
$$
\prod_{1\le k<K_q}\frac{q}{\gcd(q,(k-1)!k!)}
$$
different states can occur, independent of tournament order. This is an
upper bound, not an assertion that every tuple is realized.

**Substitution law.** The state equivalence is a congruence for all
uniform module substitutions. Indeed, if $C,C'$ are substituted at the
same vertex of a quotient with the other modules fixed, any uniform
outer test of the result regards the whole remaining vertex set as one
uniform exterior for $C$ or $C'$. Thus all its counts agree modulo
$q$. Replace several modules one by one. Nested tests never need a
larger precision vector; the quotient carries these operations on
*realizable states*. Equivalence over the integers follows by the same
argument using all probes. It does not imply equivalent labelled endpoint
menus or physical compatibility with a move which splits a module.

**Small exact model (modulo 3).** Here $K_3=3$. Write $f_C$ for the
odd-cycle-packing polynomial below. Finite differences give
$D_1=f_C(1)$ and $D_2=f_C(2)-2f_C(1)$. All exponents in $f_C$ have the
parity of $n=|C|$, so modulo 3
$$
(D_1,D_2)=\big(H,\,((-1)^n-2)H\big).
$$
The realized quotient has exactly five states: one zero state, and
the four pairs of nonzero $H$ residue and parity. Odd/even transitive
orders supply both parity choices for residue 1. A four-vertex
tournament made by orienting $1	o2	o3	o4$ and all other edges
forward except $4	o1$ has two cyclic triangles, hence $H=5$ by the
odd-cycle formula. Ordinally adding a singleton supplies the other
parity for residue 2. A directed triangle supplies the zero state.
Consequently, along ordinal sums the four nonzero states form
$C_2\times C_2$ (multiply the count sign, add order parity), with an
absorbing zero. This tiny model has all its states realized; a proof
for any larger modulus requires more than counting formal tuples.

**Boundary example and a minimal hostile.** Take $C=TT_2$ on named
vertices $a,b$, with $a\to b$, and $C'$ the converse on the same
names. The uniform contexts cannot distinguish them even integrally:
they are unmarked-isomorphic. Now introduce $z$ with
$b\to z\to a$, using exactly the same attachment to the names in both
cases. The first graph is a triangle with $H=3$, and the second is
transitive with $H=1$. Thus for any odd $q>1$ this fixed selective
attachment fails even to descend to the finite uniform quotient. Order
two is minimal: at order one a new vertex cannot split the module. It is
not an unlabelled counterexample; it pinpoints precisely the lost marked
port. This is the smallest illustration of why taking the first pair of
matching states from different atlases is not a lawful composition.

## A complete selective observer for a single added vertex

The missing state can be made explicit. Let $A$ be the old vertices
dominating a new vertex $z$, and let $B=V(C)\setminus A$, so $z$
dominates $B$. Write $f_C(t)$ for the odd-cycle-packing polynomial,
where a packing with $j$ surviving atoms and $r$ odd cycles contributes
$2^r t^j$. For each even-cardinality nonempty $S\subseteq V(C)$, let
$h_{b,a}(C[S])$ count paths on $S$ from $b\in B\cap S$ to
$a\in A\cap S$. Then
$$
f_{C+z}(t)=t f_C(t)+
2t\sum_{\substack{S\subseteq V(C)\\ |S|\ge2\text{ even}}}
 f_{C\setminus S}(t)\sum_{b\in B\cap S,\ a\in A\cap S}
 h_{b,a}(C[S]).
$$
The first term is the packings not using $z$. Otherwise its unique
odd cycle is $z\to b\to\cdots\to a\to z$; deleting it leaves a packing
on $C\setminus S$. Its contracted atom contributes $t$, its weight
is two, and the path on $S$ has even order. This is an exact equality
of polynomials, so this marked family of contractions controls *every
later uniform test* of the result, not just its immediate $H$. With
$f_{\varnothing}=1$, it includes a cycle spanning the extension.
Taking $t=1$ recovers the path-count increment. Uniform incidence makes
the sum empty. For the two-vertex hostile it is exactly the single
coefficient distinguishing the direction $a\to b$.

This is a sufficient observer, not a minimal coordinate theorem. It uses
induced-subset endpoint counts, rather than just full-graph endpoints.
No finite quotient for the action of *arbitrary* successive selective
additions is asserted.

## Cross-field precedents that actually inform the comparison

* **Formal correction at a filtered boundary:** [THM-2240, DC2 Ore grade
  response](../../../01-canon/theorems/THM-2240-dc2-grade-response-gauge-is-not-a-continuation-state.md).
  Equal grade-six cancellation admits a homogeneous syzygy with nonzero
  grade-seven effect. The exact parallel is an observation kernel not
  invariant under the proposed move: uniform equivalence is a true
  congruence for one composition language but not for selective addition.
  The Ore statement does not show nontermination or realize a tournament.
* **Observation versus iteration:** [THM-4374, lonely-runner exit
  observability](../../../01-canon/theorems/THM-4374-lrc14-seventeen-step-metric-exit-observability-and-shift-congruence-rigidity.md).
  Four equal current exits split after one shift, and the sharp
  17-observation profile recovers the parameter on its fixed tail.
  Alternating separators play the role of a certified observation suite,
  but reconstruct the complete *uniform response* only. A larger move
  alphabet requires a new proof, not extrapolation of a horizon.
* **A finite presentation is not finite behavior:** [THM-3511, Rule 30
  signalizers](../../../01-canon/theorems/THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile.md).
  Three automaton states and a control bit present all exact innovations,
  yet equal depth-three portraits have different next gaps (8 versus 5).
  This is a useful comparison with finitely many transfer formulas but
  unbounded endpoint/target data. The modulus-dependent finite quotient
  above really is finite and closed, precisely because it uses proven
  factorial annihilation and restricts its move alphabet.
* **Select the surviving channel, do not average the construction:**
  [Session 114 planar Jacobian weight selection](../session114_continuation_symmetry_20261030/jacobian/README.md).
  Exact weight of one Keller coordinate selects the sole channel capable
  of contributing the constant determinant. It does not claim the
  averaged map is a Keller map. Likewise factorial annihilation selects
  the run channels of an already physical tournament identity; it does
  not construct a physical tournament with an arbitrary tuple of
  surviving entries.

**Chan–Pak scope.** This produces a safe common layer for comparing
arithmetically different uniform constructions, plus a precise boundary
where selective repair needs a richer layer. It does not make a
congruence certificate a positive prescribed-target construction.
