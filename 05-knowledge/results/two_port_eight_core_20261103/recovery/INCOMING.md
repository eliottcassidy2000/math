# Incoming mechanisms and their exact scope

Inheritance pass at baseline `1b50a1d763`, incorporating the mathematical
signal of incoming checkpoints `0a463a024f`, `3f2fa00280`, `1b50a1d763`.
The maintained proof files below, rather than a chronological log, determine
the claims. This note adds no universal two-port theorem, recursive spectrum
closure, or eight-private-point extension theorem.

## The source graph is part of the arithmetic state

[recursive_seed_switching_20261103/README.md](../../recursive_seed_switching_20261103/README.md)
and its [transport/FOLLOWUP.md](../../recursive_seed_switching_20261103/transport/FOLLOWUP.md)
give the closest exact repair and the hostile to scalar-only closure.

**REFUTED:** the fixed 219-key seed bank is closed under all required digit
steps. The actual chain

```
170807 -> 512423 -> 1537271 -> 4611815 -> 13835447
```

ends at a count absent from that bank at every reader length and prefix
clock. The complete inverse is based on binary sections of the reader,
not on an invalid monotonicity assumption. This is a bank obstruction;
the count itself is not forbidden.

**FINITE-EXACT repair:** promote two old tail names into an enlarged
eight-name seed and change the selected physical mark, preserving all 23
old names. The resulting reader and prefix give

\[
(h,F)=(594839,827538),\qquad h+16F=13835447.
\]

The physical repair changes 18 arcs using 36 native complements. Its next
two actual children, 41506343 and 41506345, return to the old bank; the
enlarged seed alone accepts neither. Their recorded native programs use
40 and 46 complements. Four specified ordered port pairs retain full
17-query observers and independent endpoint/phase checks. The all-order
selection rule for such seed switches remains **OPEN**.

**Exact connection to two ports:** a transition must map the named source,
not just its integer $H$, to the named target. Equal Hamilton counts can
have different cut or endpoint responses: the inherited 81/27 pair has
the same $H=462915$ but cut-256 responses 1061775 and 1395279. The
needed sidecar is the paired endpoint state and physical port identity.
Native legality, strongness, scalar acceptance, and equality of complete
boundary responses are separate predicates.

## Finite boundary, unbounded path-component depth

[multiport_defect_grammar_20261103/transport/PROOF.md](../../multiport_defect_grammar_20261103/transport/PROOF.md)
proves an exact interface with paired component endpoints at active names.
Forgotten endpoints become `*`; a polynomial variable counts components
whose two endpoints are forgotten. The formal number of coordinates at
$p$ active ports obeys

\[
g_0=1,\quad g_1=4,\quad
g_p=4g_{p-1}+2(p-1)g_{p-2}.
\]

This bounds the number of polynomial coordinates, not their degree, their
linear rank, or the number of reachable states. Exact joins act on
distinct components. A flip acts as $I+J_{reverse}-J_{forward}$; products
of flips require mixed terms.

**Preserved:** Hamilton words and path-cover counts under the specified
join, flip, forget, and birth operations. **Lost:** the identities of
forgotten names. A later birth must treat every forgotten vertex uniformly;
arbitrary later attachments require retaining more names or the full
squarefree endpoint state. Thus this grammar does not by itself prove
universal two-port reconstruction.

The cheapest hostile to truncating future transport by the current
observer is already exact: a three-flip mixed Hamilton/end response of
396 requires the five-component layer. The full two-birth response of a
fixed core uses at most three old paths, but repeated operations can draw
on deeper layers. These statements have different quantifiers.

## Two uniform digits exist; a common recursive chart does not follow

[transport/TWO_DIGIT_ENTRY.md](../../multiport_defect_grammar_20261103/transport/TWO_DIGIT_ENTRY.md)
proves that the same infinite source family has both uniform entries:

\[
H_0(t)=-5+7\,2^t,\quad
H_2(t)=-13+21\,2^t=3H_0(t)+2,\quad
H_4(t)=-11+21\,2^t=3H_0(t)+4.
\]

The +2 output uses a transitive pair followed by a partially marked
directed triangle; the +4 output uses a transitive five-vertex suffix.
The +2 signature's marked endpoint entry is **1, not 3**. Its zero-start
physical new port is distinct from the outer root, whose end defect is
16. Iterating the existing two-birth channel there triples the count;
it does not add another 2. The +4 construction includes a specified root
edge reversal, with strongness proved; it has no claimed native-guard
implementation in this entry theorem.

[TRANSITIVE_CHART.md](../../multiport_defect_grammar_20261103/TRANSITIVE_CHART.md)
gives a **PROVED obstruction** for a fixed transitive suffix chart. Write

\[
H_t(f)=1-j(f)+2^t f,
\quad
j(f)=\sum_{i<k\text{ selected}}
       \bigl(2^{k-i-1}+[k=i+1]\bigr).
\]

A fixed output chart and fixed clock shift satisfying
$H_{t+k}(f')=3H_t(f)+d$ for all legal clocks must obey

\[
2^k f'=3f,\qquad d=3j(f)-j(3f)-2.
\]

Thus one fixed transitive chart cannot implement both uniform digits.
For $f=7$, the digit is 4; at its $f=21$ output the next uniform
digit is -28. Clock-dependent output charts and nontransitive suffixes
are outside this obstruction. The scan of odd $f<2^{16}$ finding no
+2 chart is **FINITE-EXACT**, not an all-order exclusion.

## Exact component language: acceptance versus a physical arrow

[COMPONENT_LANGUAGE.md](../../multiport_defect_grammar_20261103/COMPONENT_LANGUAGE.md)
keeps the seven-coordinate signature $(h,p,f,g,j,e,s)$. Its ordered
composition law is

\[
(hH,\ pP,\ fH+pF,\ gP+hG,\ jH+hJ+eS+gF,\ hE,\ sH).       \tag{1}
\]

The term $gF$ is the retained cross-component joint, not a dispensable
scalar correction. A singleton atom with marking $b\in\{0,1\}$ has
signature $(1,2,b,b,0,b,b)$; a triangle with $k\in\{0,1,2,3\}$ marked
vertices has $(3,12,3k,3k,k(k-1),k,k)$.

For a word with $c$ triangles, $h=3^c$. The ratio $q=f/h$ is an
integer mixed-radix expansion with singleton radix 2 and triangle radix 4.
The last digit must be nonzero for strongness. At fixed $f$, $c\le
v_3(f)$, so the inverse word search is finite and complete, including
the possible zero binary padding bit in a leading radix-four digit.

The source chart with $f=7$ and intercept -5 admits both one-letter
entries. At common clock, the accepted digit words through depth three
are exactly

```
depth 1: 2, 4
depth 2: 22, 42
depth 3: 422
```

There are respectively 4, 7, 13 candidate atom words at $f=21,63,189$.
These are **FINITE-EXACT complete acceptance results for this grammar and
clock**. They do not supply a physical arrow between successive word
representations, exclude a target Hamilton count globally, or address
other clock choices. [incoming_check.py](incoming_check.py) independently
enumerates the bounded word universe from the atom signatures and verifies
these values; it does not reuse the producer's inverse search.

The related [session182 induction proof](../../session182_physical_returns_20261104/induction/PROOF.md)
confirms the same ordered law and the terminal-component strongness test.
The seven scalars close under ordered assembly; the full paired forest
state is needed for arbitrary finite operation words.

## Consequence for the eight-private-point lane

The common principle has an exact, limited form: retain the data needed
by the **next operation**, then test realizability separately. In tournament
transport the forgotten pairing or component depth can invalidate an
integer recurrence. In the geometric lane an abstract bad list can fail
to arise from actual unit contacts with the fixed seven-point carrier.

The map used here is explicit: private unit graph plus a fixed boundary
colouring maps to the allowed list at each private vertex. The inherited
contact bound gives lists of size at least three. A nonzero graph-polynomial
coefficient from [LIST_CERTIFICATE.md](LIST_CERTIFICATE.md) proves that every
such list assignment colours that graph, eliminating a candidate regardless
of geometry. Zero capped support preserves the candidate for separate
list and metric analysis. It is not a six-colour obstruction.

The live board is therefore: paired endpoint state, anonymous component
depth, physical source identity, graph-polynomial list certificates, and
reciprocal geometric incidence. The next decisive tests are complete
coefficient verification on the 38 supplied private-core types, exact
geometry of the eight residuals, and full endpoint/three-path tests on any
modular tournament survivors. Neither a numerical grammar acceptance nor
a necessary metric filter may be promoted to its target predicate.
