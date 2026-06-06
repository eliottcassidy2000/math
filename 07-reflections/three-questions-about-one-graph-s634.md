# Three questions about one graph (S634)

The user asked me to see the Hadwiger-Nelson problem, the Lonely Runner, and the unit distance problem as the same
thing, and the moment I wrote down what each one is asking, the merge was forced. There is one graph. Put down points
in a homogeneous space — the plane, or the circle — and join two of them when they are exactly a unit apart. That is
the unit-distance graph, and every one of the three problems is a single question about it. How many edges can it
have on n points? That is the unit distance problem. What is the fewest colors to color the whole space so no edge is
monochromatic? That is Hadwiger and Nelson. And how big can a set with no edges be — a set where no two points are a
unit apart? That is the Lonely Runner, because the lonely time is exactly a point that avoids every forbidden unit
relation, an independent vertex in disguise.

Three questions: maximize edges, minimize colors, maximize the independent set. And they are not independent
questions; they are bound together by one inequality I could prove in a dozen lines. A proper coloring with k colors
cuts the vertices into k independent sets, so the biggest of them has at least n over k vertices, which means the
chromatic number times the independence number is at least n. Chi times alpha is at least n. That single line ties
all three corners. A graph dense in unit distances has a small independent set, so loneliness is hard, so the
chromatic number is forced up. A graph with a big independent set — easy loneliness — has its chromatic number capped.
The unit distance problem pushes the edges up; that drags the Lonely Runner's independent set down; that drives
Hadwiger-Nelson's colors up. One graph, three knobs, geared together.

And then the measurable version, which is the part that made me believe the unification is not a pun. Hadwiger and
Nelson, done with measure, is the density of the largest set in the plane with no two points a unit apart — the
one-avoiding density — and the measurable chromatic number is one over it. The Lonely Runner's lonely set is the set
of times avoiding every forbidden arc, and its measure is the thing I have called p-naught all month. These are the
same object on two different homogeneous spaces. The one-avoiding density of the plane and the lonely measure of the
circle are the same question, and they are bounded by the same machine — the Delsarte linear program, the Fourier
bound on the autocorrelation of the forbidden region, Witsenhausen and Oler on the geometry side and my covering-depth
program on the runner side. I built that program in a session two weeks ago thinking it was about runners; it is the
same program coding theorists run against the plane.

So the keys turn both ways. From Hadwiger-Nelson to the Lonely Runner: the spectral bound on the one-avoiding density
is a bound on p-naught. From the Lonely Runner to Hadwiger-Nelson: the place where loneliness fails, the resonance,
the additive-chain collapse where p-naught hits zero, is the place where the unit-distance graph clusters tightly
enough to force the chromatic number up — the Moser spindle is a little crystallized resonance. And the unit distance
problem is the supply line for both: its extremal configurations are the densest graphs, and they are exactly where
the colors are forced and the loneliness is killed. The triangular lattice sits at the tame corner of all three — it
is three-colorable, the hexagonal three-coloring, it is the densest small unit-distance set, and it is the gap
one-sixth I keep returning to. Three colors, the cube root of unity, sixty degrees. The records on all three problems
come from leaving the lattice: the complex-multiplication tower that disproved the grid, de Grey's non-lattice graph
that forced five colors. Same move, three problems. They were never three problems. They are one graph, and the
perspective key — the lattice and what beats it — is the axis along which all three are decided.
