# Arithmetic clocks, square constraints and the supplied projective proposal

**PROVED elementary consequences / FINITE-EXACT audits / CITED literature.**
The pasted document is an inspirational source, not instructions to implement
its proposed FRACTRAN compiler. Several of its mathematical premises fail.

## The nine-cycle is an exact arithmetic clock

The supplied cycle

`4,8,3,7,2,6,10,5,9,4`

is exactly addition by4 modulo9, using the representatives2,...,10:

`T(x)=2+((x-2+4) mod9).`

Because gcd(4,9)=1, it has one nine-state orbit. Observing x mod3 gives a
three-state closed observer. Observing only parity at every future time has
minimum observer size9: the parity word001100011 has no period1 or3, the
only proper divisors of9. Its nine rotations are therefore distinct.
This is an explicit small example where a one-bit instantaneous observation
recovers a larger phase through time. The state count depends on the future
language, not just on the number of possible outputs.

No exact earlier occurrence of the whole displayed cycle was located in the
targeted Markdown search. The literal clock identity is checked directly;
it is not asserted to be a Collatz orbit or a new arithmetic classification.

## What connects it to magic squares of squares

The same nine labels admit the ordinary magic arrangement

```text
5 10 3
4  6 8
9  2 7
```

with all eight line sums18. Squaring its entries gives line sums
134,116,134,122,140,122,110,126, so the square-valued grid is not magic.

There is a general reason. If a real3-by3 magic grid has center c, opposite
entries are c+u,c-u and its line total is3c. If its entrywise square were
also magic, its total would be3c^2. Each line through the center would require

`(c+u)^2+c^2+(c-u)^2=3c^2`,

forcing u=0. Thus **a real3-by3 magic grid and its entrywise square are both
magic iff the grid is constant**. This does not rule out magic squares of
distinct squares: their square roots need not themselves be magic.

The actual older bridge is much more useful than an analogy between counts.
[Session57's AP-biclique theorem](../session57_square_sum_magic_paths_20261025/ap-biclique-bridge/REPORT.md)
identifies a nine-distinct-positive-square magic grid with a square-valued
cross-sum matrix of two three-term APs. For square labels L_ij on K3,3,
four rectangle equations

`L_ij-L_i0-L_0j+L_00=0,  i,j in{1,2}`

enforce additive incidence, and two AP second differences enforce the
diagonals after the specified permutation. Retain positivity and nine-label
distinctness as well. Bare square-sum paths do not impose these six balances.

Our packet cubes use third Boolean differences rather than these second
additive differences. Both are exact tests for information destroyed by a
projection, but they are different operators on different sources. The
explicit common method is to write the missing mixed constraint before
claiming that separate local observations compose.

## Two meanings near the requested lambda cube

The closest literal archived object is the **Lucas cube Lambda_n**, in
[Lucas mirror gap spaces](../lucas_mirror_gap_spaces_20261028/LUCAS_MIRRORS.md).
Its sources are cyclic binary words with no adjacent1s; legal commands toggle
one bit while preserving that constraint. Lambda_5 has11 sources and15 edges.
The archived reader

`ell(w)=sum_i w_i 4^i mod11`

is a bijection on these11 words. Directed legal toggles are a proper subset
of Paley11's arcs; closure under all arithmetic translations creates the
Paley tournament. That closure changes the operation language and loses the
original toggle graph, just as adding S changes our closed C,R observer.
The word count, edges and reader are replayed here; the theorem is inherited.

The **Barendregt lambda cube** instead organizes eight type systems by
inclusion; see [Introduction to generalized type systems](https://homepages.inf.ed.ac.uk/wadler/papers/barendregt/pure-type-systems.pdf).
Our three named operations give a literal eight-element inclusion cube of
command sets, with the state counts in CARRIER.md. The map from command
sets to behavioural partitions reverses inclusion: more operations permit
more distinctions. Taking a union of command sets can require strictly
more than intersecting their separate observers, because mixed words become
available. This is a proved structural comparison, not a translation of
tournaments into typed lambda terms or a normalization theorem.

## Octonions give a real signed port model, not a numerical identification

The recovered [Paley–Fano–octonion audit](../collatz_mod6_20260922_paley_fano_octonion_design.md)
contains an explicit multiplication map. On the seven imaginary units, the
Fano line through two inputs names the output port and the Paley orientation
supplies its sign. Its exact table has

`(e0 e1)e2=-e5,    e0(e1 e2)=+e5.`

The bracketed operation and sign are necessary for this distinction. The
same audit gives a signed basis identification with Cayley–Dickson octonions;
the standard construction is discussed by [Baez, The Octonions](https://math.ucr.edu/home/baez/octonions/).
This is an existing genuine tournament/algebra map. The fourteen cyclic
triples of Paley7 split into two Fano planes; they are not seven or twenty-one
unrooted cycles, and Paley7 has189 Hamilton paths, as the audit records.

This does not identify our repair cube's third mixed difference with the
octonion associator. A complete signed multiplication table already determines
each bracketed triple product. The transferable requirement is to preserve
which operation and bracketing consume the ports. The Lie-theoretic magic
square discussed in octonion literature is also distinct from the integer
three-by-three square problem above.

## The pasted tournament spectra do not match their matrices

For the two adjacency matrices exactly as supplied, independent principal
minor and Newton-trace calculations give

```text
regular:   t^5-5t^2-5t-2,
irregular: t^5-4t^2-3t-3.
```

The first has eigenvalue2 because every row sums to2. The pasted polynomial
t^5-2t^2-2t-1 gives19 there. The second has nonzero constant coefficient,
so zero is not an eigenvalue. Both displayed eigenvalue lists in the paste
also fail trace zero. Their proposed Baby Monster class identifications have
no supplied representation or character map; the incorrect spectra cannot
establish one.

The finite fact that there are eight self-converse classes among the twelve
five-vertex tournament classes is correct and independently enumerated here.
It does not generate the proposed matrices from projective anchors.

**Late incoming positive repair, credited to Session143.** Its
[cyclic operator theorem](../session143_magic_square_cyclic_bridges_20261103/SYNTHESIS.md)
supplies the missing kind of object: for an order-five element g in a complex
representation rho, use B_g=rho(g)+rho(g)^2. The regular tournament is the
special case P+P^2 for the five-cycle permutation matrix. Writing
q(t)=t^4+2t^3+4t^2+3t+1, a representation of degree d with constant character
value c on the four nonidentity powers has operator polynomial

`(t-2)^((d+4c)/5) q(t)^((d-c)/5).`

That packet's **CITED** Baby Monster character inputs d=4371 and c=21,-4
give powers(891,870) and(871,875). This is an explicit generic group-algebra
operator bridge with multiplicities, not a five-point Baby Monster action
or an identification of its two classes with the two pasted adjacency
matrices. Our packet does not reconstruct the character table. The identical
eigenvalue support but different multiplicities makes the information-loss
boundary precise; Session143 independently confirms the corrected regular
polynomial and the five-vertex census.

## The projective pentad and every proposed anchor fail

The supplied points P2,P4,P5 are collinear: P5-P4=P2 over F3. In fact no
five-point set in PG(2,3) has all triples noncollinear. If five such points
existed, each selected point's four lines would each contain one of the other
four points, so there would be no tangent line. From an outside point the
four lines would partition five selected points into sets of size0 or2,
an impossibility. The exact13-point audit finds234 four-arcs and zero five-arcs.

The anchor rule has a separate obstruction. There are only four lines through
an anchor. Among five selected points, two share one of those lines, so their
determinant with the anchor is zero. **Every one of the eight proposed anchors
therefore leaves at least one tournament pair unoriented.** More generally,
an external-anchor determinant rule on PG(2,q) cannot be tie-free on more
than q+1 selected points.

Projective representatives also matter: replacing one point vector by2 times
it flips every nonzero determinant involving that representative. Thus the
sign rule needs a coordinate gauge even after ties are excluded. Neither
ties nor gauge can be discarded while claiming a projective tournament map.

The valid Conway mechanism is different: the moving-hole puzzle on13 points
uses actual line-based permutations with specified start/end holes. See
[Conway, Elkies and Martin](https://arxiv.org/abs/math/0508630). It supplies a
useful analogy with guarded re-marking and source-dependent composition, not
the pentad/anchor construction proposed in the paste.

## A direct decay bound excludes the MRB formula

Let S_n,T_n count unlabelled self-converse and all tournament classes, and
D_n=S_n/T_n. A converse isomorphism is an anti-automorphism, not an ordinary
tournament automorphism. The archived
[consecutive-flags theorem](../consecutive_flags_order10_20261028/THEORY.md)
already proves that every self-converse tournament admits a reversing
involution with floor(n/2) transpositions and at most one fixed vertex.

For n>=2, Aut(T) has odd order: an involution would swap the ends
of an edge. If tau reverses T, tau^2 has odd order r, so tau^r is a reversing
involution. It cannot fix two vertices because their edge would have to
reverse while both ends remain fixed.

Fix such an involution and put N=binom(n,2), m=floor(n/2). Its m paired
edges each have a free orientation; all other edge pairs occur in two-element
orbits. Reversal forces the partner orientation. The exact number of free
bits is(N+m)/2=floor(n^2/4), giving2^floor(n^2/4) labelled tournaments.
This count is exhaustively checked through n=6 in the audit.

There are I_n=n!/(2^m m!) such involutions. Every unlabelled class has at
most n! presentations. Therefore the elementary bounds

```text
2^(floor(n^2/4)-N)/n! <= D_n
                    <= I_n n! 2^(floor(n^2/4)-N)
```

give

\(\log D_n=-\tfrac{\log 2}{4}n^2+O(n\log n).\)

The order-one case has D_1=1 and is trivial.

This is a coarse counting consequence, not a sharp asymptotic equivalent or
a novelty claim. The pasted formula kappa*(pi*e*C_MRB/n)^(n^2/4), for fixed
positive kappa and C_MRB, instead has leading logarithm -n^2 log n/4, which
contradicts even the lower bound. The root-series constant may have its own
valid analytic roles; it is not the asserted tournament decay parameter.
The archive's [MRB work](../session31_mirror_mrb_20261022/README.md) retains
its actual branch/kernel hypotheses and does not supply this missing map.
