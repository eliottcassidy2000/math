# Erdős 61: the chromatic-4 targets are a tower from C₅ (S726)

Erdős's sixty-first problem is infinite graph theory, the kind the cluster has no finite machinery to
crack: two graphs each with uncountable chromatic number `aleph_1`, must they share a common subgraph of
chromatic number four, or even countably infinite chromatic number? Erdős's own guess pointed the way —
probably every `aleph_1`-chromatic graph already contains every chromatic-four graph of large girth, in
which case any two such graphs share all of them and the problem dissolves into a question of
unavoidability: what substructure is forced inside a graph whose chromatic number is uncountable? The
paper the prompt sent me to, Huang–Ju–Zhou on Erdős–Hajnal beyond the five-vertex path, is the finite
cousin of exactly this instinct — H-free graphs are forced to contain a large clique or independent set —
and the shared theme is the one the cluster keeps meeting from the other side: complexity forces
substructure. I cannot solve the infinite problem, but I can say precisely where its targets come from,
and they come from a place the cluster already lives.

The targets are a tower, and the tower starts at `C₅`. The canonical chromatic-four graph with no
triangles is the Grötzsch graph, and the Grötzsch graph is the Mycielskian of the five-cycle. Mycielski's
operation takes a graph, adds a shadow of every vertex and one apex, and returns a graph with chromatic
number one higher and no new triangles. Iterated from `C₅` it climbs chromatic number three, four, five and
upward, and I verified the first rungs: `C₅` at chromatic three and girth five, Grötzsch at chromatic four
on eleven vertices, the next Mycielskian at chromatic five on twenty-three, the independence numbers two,
five, eleven, the independent-set counts eleven, a hundred three, seven thousand four hundred seven. So the
chromatic-four obstruction that Erdős wants forced inside every uncountably-chromatic graph is the first
nontrivial rung of a Mycielski tower, and the foot of that tower is `C₅` — which is the Paley graph on five
vertices, the exact object S725 found at the bottom of Erdős's ordinal Ramsey problem 592. The two
infinite problems, one about colouring ordinals and one about colouring uncountable graphs, rest on the
same five-vertex stone.

The Mycielski tower is a transfer, and that is the cluster's home ground. Each step is a fixed graph
operation that advances chromatic number by exactly one — an additive ladder, the cold end of the
temperature spectrum S720 mapped out, with chromatic number playing the role of the running coordinate and
triangle-freeness the frozen invariant. Independence tracks up the rungs the way the cluster's master
object `H = I(Omega,2)` tracks everything else, and the common-subgraph question becomes, in this language,
whether two towers share a rung. The finite skeleton of problem sixty-one is therefore a transfer-spectrum
problem, and the cluster has spent a month learning to read transfer spectra.

What is genuinely hard is the same thing that was hard in the last three sessions: girth. Mycielski caps
the girth at four — it manufactures four-cycles at every step — so the Mycielski tower delivers
chromatic-four triangle-free graphs but not the large-girth ones Erdős' "probably" really wants. Large
girth means locally tree-like, the laminar ultrametric structure S725 built for the ordinals and S722–S723
met as the no-four-concyclic and convex-hull constraints. Chromatic four with large girth is a
high-temperature target — sparse yet stubbornly chromatic, the probabilistic and Ramanujan constructions,
not the Mycielski one. And the host is constrained too: an uncountably-chromatic graph necessarily contains
a four-cycle, so it has small girth itself even while it must, conjecturally, contain large-girth subgraphs.
The tension between the host's forced short cycles and the target's required long ones is the heart of
sixty-one, and it is the same chromatic-versus-girth tension that runs through the whole distance corner of
the repository.

The creative move, then, is to read `aleph_1` as a limit and run the S725 renormalization. The finite
levels of the chromatic tower are individually avoidable — you can build a graph that dodges any one fixed
chromatic-four graph — but the uncountable chromatic number is a fixed point that should be unable to dodge
the whole tower at once, exactly as a limit ordinal restores the partition relation that every finite
exponent below it breaks. The program that falls out is concrete even if the execution is not mine to
finish: prove that an `aleph_1`-chromatic graph contains the Mycielski targets by an elementary-submodel or
compactness tower argument, treating the uncountable colouring as the limit that forces the finite
obstruction; and attack the large-girth case through the laminar structure, where the cluster's ultrametric
and shell-tower tools live. The sharpest single question to settle first is the smallest one: is there one
fixed chromatic-four graph — Grötzsch, say — contained in every uncountably-chromatic graph? If even one
fixed chromatic-four graph is unavoidable, the chromatic-four form of sixty-one is immediate, because both
`G₁` and `G₂` would contain it and so share it. That is a finite target inside an infinite problem, and it
sits at the first rung of a tower whose foot we have already named.
