# Exact closure under complement and a selected involutive marking

**PROVED all-size observer theorem; FINITE-EXACT counts at order seven.**
The operation language is named and specified below. This is not closure
under all re-markings, arbitrary attachments, or an unmarked metagraph.

## The source and legal operations

Let X_n be the Boolean chart of tournaments on0,...,n-1 containing the
directed path0,1,...,n-1. Coordinates are the lexicographically ordered pairs
(i,j) with j>=i+2; bit1 means i->j. Output H counts directed Hamilton paths.
C reverses every free arc and retains the selected path. Thus C^2=1.

Fix a nonidentity vertex permutation sigma with sigma^2=1. Let D be the
sources for which sigma(0),...,sigma(n-1) is a directed Hamilton path. On D,
R re-marks along this path: new vertex i carries old vertex sigma(i). It
preserves the actual ordinary tournament up to the displayed relabelling.
Commands outside their domain fail observably; legality is part of the output
language. The permutation is retained, not silently replaced by an arbitrary
representative of the same unmarked class.

**Guard separation.** R(D)=D and R^2=1 on D, while D and C(D) are disjoint.
The first assertions follow from sigma^2=1 and the old standard path. If D
is nonempty, the nonidentity selected path must use a free pair: a path using
only the edges of the standard chain is either that chain or its illegal
reverse. On that free pair C reverses the required orientation. Consequently
the selected path cannot be legal both at x and Cx. If D is empty the result
below reduces to the complement-only observer.

## The minimum deterministic observer

Define a phase epsilon(x):1 on D,2 on C(D),0 elsewhere. Put a=H(x),b=H(Cx).
The following state is the **coarsest exact deterministic observer** for all
finite words over{C,R}, with output H and command legality:

```text
epsilon=0:  (0,a,b),
epsilon=1:  (1,a,b,d),  d=H(C R x),
epsilon=2:  (2,a,b,d),  d=H(C R C x).
```

Here displayed operator compositions act rightmost first. In command-word
notation, the last two queries are respectively RC and CRC.

The complete transitions are

```text
C: (0,a,b)   -> (0,b,a),
C: (1,a,b,d) -> (2,b,a,d),
C: (2,a,b,d) -> (1,b,a,d),
R: (1,a,b,d) -> (1,a,d,b).
```

R is illegal in phases0 and2. The current H is always the second coordinate.

**Sufficiency.** These formulas follow from C^2=R^2=1, H(Rx)=H(x), and guard
separation. They give every next class and output, hence every future word.
The actual transition components have at most four sources: an R-pair in D
with its two C-partners, or just a C-pair outside D union C(D). No actual
source is fixed by the nonidentity involution R: that would be an order-two
tournament automorphism, reversing an exchanged pair's edge. A quotient
state can nevertheless be fixed by R when b=d. The inverse transitions
retain their legal source.

**Necessity.** Legality of R separates phase1; legality of command word CR
separates phase2 from phase0. The empty word and C separate different a,b.
The word RC in phase1 or CRC in phase2 observes d. Thus any two different
displayed states have a distinguishing word of length at most3. Every exact
observer must refine this one. No abstract reset or arithmetic averaging
enters the argument. QED.

The extra information is one integer response and a three-way legality phase,
not one bit in general. This phase records path guards, not Frobenius transport.

## The actual seventeen, seventy-one and thirty-three marking

For n=7 choose sigma=(0,5,6,3,4,1,2), in one-line notation. It swaps the first
and third of the three ordered two-vertex blocks, fixing vertex0. D is cut
out by exactly three free-edge requirements, so it has4096 sources.

The inherited hostile replays literally:

```text
x=24511, R x=31,
H(x)=H(Rx)=17, H(Cx)=71, H(CRx)=33.
```

The new states are(1,17,71,33) and(1,17,33,71), exchanged by R. The ordinary
tournament class remains the same; the selected marking changes its actual
complement target. The source is inherited from
[Boolean realization](../boolean_realization_71_17_20261102/THEORY.md).

The complete32768-tile calculation gives

| Language | Minimum states |
|---|---:|
| H only |77|
| C |3008|
| C,R |10455|

The last total is2955 phase-zero states and3750 in each active phase. The
source masses are24576,4096,4096 respectively. These are observer states,
not tournament isomorphism classes.

An independent endpoint-subset DP checks H at every tile. A second observer
calculation lists every command word through length3; distinct response
vectors certify necessity, and direct transition congruence certifies
sufficiency for arbitrarily long words. It agrees with the closed formula
and with the producer's partition refinement.

## A second selector and the cube of operation languages

Add S with permutation(0,3,4,1,2,5,6), swapping the first two two-vertex
blocks. It also has4096 legal sources and is an involution there. The
permutations R,S generate the simultaneous S3 action on the three blocks,
but their commands still have source-dependent path guards.

The complete minimum-state counts for all eight operation sets are

| Commands | States |
|---|---:|
| none |77|
| C |3008|
| R |146|
| S |151|
| C,R |10455|
| C,S |9545|
| R,S |431|
| C,R,S |15874|

For the last row the refinement counts are77,6439,12771,15849,15874,15874.
Independent full word observation through length4 plus a congruence check
certifies the15874 result. It is a complete finite theorem, not an all-size
formula for two selectors.

The joint refinement of the separate C-only and R-only observers has4680
classes, strictly fewer than10455. Combining two operation families requires
closure under alternating words, not merely collecting their separate
summaries. This is a precise memory-growth mechanism.

For example, tiles74 and76 share the first closed state(0,91,37), but command
word CS separates them once S is permitted. An operation-cube analogy with
type-system inclusion preserves the order of allowed languages; it does not
identify this cube with a lambda-calculus semantics or an octonion algebra.

## A Boolean inverse with an explicit reachability boundary

Multi-source BFS from the68 actual H17 tiles supplies shortest legal words,
with every intermediate Boolean orientation retained. Under{C,R}, precisely
148 tiles reach H17, at distances0,1,2,3 with counts68,68,6,6. Under{C,R,S},
152 reach it, with distance counts68,68,6,8,2 through distance4.

The executable [constructor.py](constructor.py) returns the shortest word,
literal masks, full endpoint adjacency, and DP-checked counts. For sources
outside the basin it exhausts their actual operation component. Rejection
means no inverse in this language, not no tournament realization of the target.

Thus more accurate memory does not itself supply broader constructive
coverage. A uniform Boolean inverse still requires an actual operation atlas
whose reachability covers the desired targets. Neither this observer theorem
nor the inherited81/27 hostile is a full metagraph classification.
