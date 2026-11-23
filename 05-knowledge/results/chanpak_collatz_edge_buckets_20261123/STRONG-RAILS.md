# Reusable binary width, and literal \(3h+2\) rails

**PROVED via finite templates.** Expand a distinguished vertex into a
nonempty uniform module \(U\). Each outside vertex either beats all of
\(U\) or loses to all of \(U\). Fix an order of the outside vertices and
select the occupied gaps for visits to \(U\). Feasibility is exactly
that the resulting word is a directed walk in the small template,
with one named placeholder for each occupied gap. Each such pattern
lifts independently to an ordered cover of \(U\) by that many nonempty
paths. Every Hamilton path has one such decomposition.

In the templates below no pattern occupies more than two gaps.
Write \((c_1,c_2)\) for the two pattern counts. The resulting read is
\(c_1 H(U)+c_2 F_2(U)\), where \(F_2(U)\) counts ordered covers into
two nonempty paths. For \(U\) transitive of width \(w\ge1\),
\(H(U)=1,\ F_2(U)=2^w-2\): choose any nonempty proper subset for
the first path. Expansion of a strong template remains strong:
uniform contacts let a path enter or leave through any chosen
module member, and a cycle through the old port connects any two
distinct module members via the exterior.

## Exact strong fork in counts

Port 0, out-neighborhood masks and occupied-gap counts:

| Row | Masks | \((c_1,c_2)\) |
|---|---|---:|
| parent | (2,8,3,5) | (5,2) |
| lower | (8,17,3,6,13) | (13,6) |
| center | (4,41,2,5,15,29) | (15,6) |
| upper | (8,41,3,4,15,29) | (17,6) |

For common \(U\), these outputs are \(h\) and
\(3h-2H(U),3h,3h+2H(U)\). For a transitive module this gives
\(h=2^{w+1}+1\) and three strong children. These are separate
templates, not one single-arc scaffold or a renewable port.

## Exact one-arc \(3h\leftrightarrow3h+2\) on unbounded heights

Again with port 0:

| Row | Masks | \((c_1,c_2)\) |
|---|---|---:|
| parent | (2,12,9,1) | (5,3) |
| center | (2,36,57,51,35,1) | (15,9) |
| upper | (2,36,41,51,39,1) | (17,9) |

For each common \(U\), the last two rows differ by the literal
reversal of the outside pair \(\{2,4\}\), changing \(3h\) to
\(3h+2H(U)\). For every transitive width \(w\ge1\) this is a
single-arc seam \(3h\leftrightarrow3h+2\) at
\(h=3\cdot2^w-1\), attained in the first row.
The flip changes the one-visit coefficient and preserves the
two-visit coefficient. This gives unbounded heights using
**different** strong hosts, not an unbounded chain of unit moves
from one graph.

A different strong rail has center masks (50,188,25,17,32,76,31,125),
port 7 and flipped pair \(\{3,5\}\). Its reads are (171,84) and
(173,84), giving \(84\cdot2^w+3\) and that number plus 2.
The factor \(h=28\cdot2^w+1\) is attained. For \(w=1,2\),
witnesses are (80,9,3,5,14,31,62) and (48,33,67,7,14,28,59).
For \(w\ge3\), take \(D=(8,1,3,102,15,23,55)\) and a transitive
block of width \(w-2\) oriented into \(D\). Add apex \(b\) oriented
out to the block, beaten by vertices 4 and 5 of \(D\), and
beating the other vertices of \(D\). For this family the exterior
template at width one has \((c_1,c_2)=(225,112)\), also checked
below, giving \(112\cdot2^{w-2}+1=28\cdot2^w+1\).
A second slope 132 occurs for masks (242,76,9,17,38,78,28,126),
port 7, flip \(\{3,6\}\); this assertion does not need an
independent realization of the quotient.

Run **python3 replay.py** for the bounded gap counts and connectivity
checks. The module argument then holds at arbitrary widths.
