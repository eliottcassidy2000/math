# Target languages, neutral reservoirs, and exact all-height two-apex fibres

**Status: PROVED from the inherited positive transfer, SCC multiplicativity,
induced monotonicity, and the cited Busch floor. Finite core polynomials and
automaton sizes below are FINITE-EXACT. Chan–Pak remains OPEN.**

The standard finite-semiring and finite-automaton facts used here are not
claimed as new. The contribution is their target-sensitive physical use and
the complete three-reservoir fibre normal form, rather than a new count census.

## 1. The question changed: retain the whole future language of an absolute target

The four-letter word model has named apexes `u -> v`, a transitive spine
`s_1 -> ... -> s_m`, and letters

```
l_i = a_i + 2 b_i,
a_i = [u -> s_i], b_i = [v -> s_i].
```

Write `T(w)` for this actual tournament and `H(w)` for its Hamiltonian-path
count. The empty word is the two-vertex tournament, so `H(empty)=1`.

For a prescribed integer `N>=1`, do not identify two prefixes merely because
both have the same H. Identify them only when they have exactly the same
successful suffixes:

```
p ~_N q  iff  for every word z,
                 H(pz)=N  iff  H(qz)=N.                 (1)
```

This is a target-specific right congruence. Its exact finite realization below
is a usable replacement for scalar deduplication. It does NOT say arbitrary
prefix replacement is a physical local operation on an existing graph: the
replacement reconstructs a graph from its word and preserves this specified
future-target predicate, not all root coordinates or all other target counts.

### Theorem 1. Exact capped continuation automaton

The entire language `L_N={w:H(w)=N}`, allowing every word length, is recognized
by an effectively constructible finite deterministic automaton. A simultaneous
Moore automaton records every H at most N and a single overflow state.

Use the inherited 36-state positive run-assignment transfer in
`04-computation/experiments/s93_chanpak_wave_20261111/seed/two_apex_probe.py`
(the positive transfer also appears in Session270). Its initial vector has
ones in coordinates0 and18 and zeros elsewhere; each of four letter matrices
has nonnegative integer entries and the final functional is a sum of selected
coordinates. Replace every arithmetic value x by `min(x,N+1)` after each
addition and multiplication. There are at most `(N+2)^36` capped vectors.
Whenever the final functional exceeds N, collapse that state to one dead state.

**Proof.** For nonnegative integers, truncation at N+1 commutes with truncated
addition and multiplication. Thus every capped coordinate equals the
truncation of its actual value after every word. The final value is exact
whenever it is at most N; otherwise overflow is correctly detected. Inserting
one new spine vertex gives an induced extension, so H never decreases.
Consequently no overflow word can subsequently return to the target. The
finite reachable-state closure therefore recognizes the entire infinite
language, not just a bounded prefix. Standard deterministic minimization gives
(1); minimizing with every count as output gives the stronger simultaneous
continuation-count equivalence. QED.

The nonnegative 36-state carrier is essential to this elementary saturation
proof. A smaller rational minimal linear representation can have negative
coordinates; blindly capping that representation is not justified. Likewise,
capping H alone discards hidden future contributions.

The threshold language `{w:H(w)<=N}` is downward closed under taking
subsequences: deleting spine names is induced deletion, hence cannot increase
H. This gives a second view through finite forbidden-subsequence bases
(Higman's lemma), but no use of that general lemma is needed for the explicit
positive automaton or the next, sharper structural normal form.

## 2. Every infinite fibre has only three neutral reservoirs

Define finite core sets as follows.

* `A_h` consists of binary cut words for STRONG one-apex tournaments of count h.
  A nonempty word starts1 and ends0. Put `A_1={empty}` separately, representing
  the singleton apex; for h>1 there are no empty words.
* `C_h` consists of four-letter words whose whole two-apex tournament is
  STRONG and has count h. In particular `C_1` is empty.
* Define the alphabet embedding `phi(0)=1`, `phi(1)=3`.

### Theorem 2. Unique SCC fibre normal form

Every word of count N belongs to exactly one of the following alternatives:

```
(I)   0^p c 3^r,             c in C_N, p,r>=0;
(II)  0^p a 1^q phi(b) 3^r,  a in A_d, b in A_e,
                              de=N, p,q,r>=0.           (2)
```

In each alternative the core words, their counts, and every clock are unique.
Thus (2) is an unambiguous constructive grammar for the entire all-order fibre.
It includes empty cores and N=1, rather than applying a strong bound to neutral
padding.

**Proof.** Every nontrivial strong component contains an apex because a
subtournament of the spine is transitive. A strong component contains a
contiguous interval of spine names: if it contains `s_i,s_j`, with `i<k<j`,
then `s_i -> s_k -> s_j` places `s_k` in that same component.

If u and v are in the same strong component, all other components are spine
singletons. Earlier singletons dominate both apexes, hence have letter0;
later ones lose to both, hence letter3. This gives (I), with the core strong.

If they lie in different components, their condensation order is u's component
before v's, since `u -> v`. All intervening non-apex components are singleton
spine names. Before u's component they have letter0; between the components
they have letter1; after v's component they have letter3. Inside u's component
every spine name beats v, so its letters are0/1 and form the one-apex word a.
Inside v's component every spine name loses to u, so its letters are1/3 and
form phi(b). A singleton apex gives the empty core. Multiplicativity across
SCCs gives de=N. SCC decomposition is unique, so are all data in (2).
Conversely every displayed word has exactly the specified SCC sequence and
count; all arcs are actual tournament arcs. QED.

The finite-core/SCC mechanism also underlies THM-5129
`one-vertex-sources-finite-supply-complete-lists-permanent-exceptions-line-counts`,
whose retained observable is a source slope with a dyadic tail. Here the
retained object is the entire fixed-apex word language and its exact ordered
reservoir locations, not a new proof of generic SCC finiteness.

This generalizes the special H13 core normalization in Session191
`two_apex/APPEND_OBSTRUCTION.md`; it does not claim that the SCC method itself
is new. That earlier special case is a useful hostile control: scalar H13
admits six different strong core presentations, not one universal root menu.

### Effective finite-core bound

The cited Busch lower bound, in the weaker convenient form
`5^(n-1)<=H(T)^3` for every nontrivial strong tournament, gives, with

```
B(N)=floor(log_5(N^3)),
```

at most B(N) spine names in a one-apex strong core of count<=N, and at most
B(N)-1 in a two-apex strong core of count<=N. No floating logarithm is required:
compare integer powers. The singleton A1 is handled separately.

Moreover the TOTAL number of nonneutral spine names in alternative II is at
most B(N): if both cores are nontrivial their two inequalities multiply; if
one is empty the surviving inequality suffices. In alternative I the bound is
B(N)-1. This is a global all-order normal form with only logarithmically many
nonneutral names, not a bound on the possibly arbitrarily long word itself.

## 3. Exact generating functions and the Dirichlet convolution of the two cores

Let the same symbols denote length polynomials

```
A_h(z)=sum_(a in A_h) z^|a|,      C_h(z)=sum_(c in C_h) z^|c|.
```

Put `D_N(z)=sum_(de=N) A_d(z) A_e(z)`, summing ORDERED factor pairs.
Then the ordinary generating function for the number of named word
presentations of count N is EXACTLY

```
F_N(z) = C_N(z)/(1-z)^2 + D_N(z)/(1-z)^3.             (3)
```

This counts words, not isomorphism classes and not distinct count values.
It follows directly from the unambiguous physical grammar (2). In particular,
for every `m>=B(N)` the fibre size is the explicit quadratic polynomial

```
[w length m, H(w)=N]
 = sum_j c_(N,j)(m-j+1)
   + sum_j d_(N,j) binom(m-j+2,2).                    (4)
```

Its quadratic coefficient is half the Dirichlet convolution
`sum_(de=N) A_d(1) A_e(1)`. If this vanishes and C_N is nonzero, the fibre is
eventually linear. If both vanish, the fibre is empty at EVERY order.
No oscillating quasipolynomial or new eventual core can appear: every
unbounded direction is one of the displayed neutral clocks.

The degree is a physical coexistence test: quadratic growth occurs exactly
when N factors as the product of two attainable one-apex core counts (allowing
1); linear growth means all realizations force both named apexes into the
same strong core. For instance13,29,41 and43 have linear, not quadratic,
word fibres in the complete inventory. This is a statement about named
presentations and shared strong components, not a symmetry classification
of ordinary tournament classes.

Examples proved by the finite core inventory plus the all-order bound:

```
F_1 = 1/(1-z)^3;
F_3 = z/(1-z)^2 + 2 z^2/(1-z)^3;
F_5 = 6 z^2/(1-z)^2 + 4 z^3/(1-z)^3;
F_13 = 6 z^3/(1-z)^2;
F_7 = F_21 = 0.
```

Thus the word multiplicities for H1,H3,H5 are respectively
`binom(m+2,2)`, `m^2`, and `2(m^2-1)` (the last for m>=1), whereas the H13
multiplicity is `6(m-2)` for m>=3. These multiplicities do not promote the
absence of 7 or21 in this grammar to a new global tournament omission proof;
those global omissions are inherited separately.

The two-apex version of Chan–Pak has now been reduced EXACTLY to the positivity
question `C_N(1)+D_N(1)>0` for every allowed odd N. Formula(3) classifies all
later repetitions once that core positivity is known. It does not establish
that positivity uniformly. In particular more neutral word presentations
cannot repair a missing core.

## 4. Recurrent automata are padding, not new strong structure

For fixed N, the grammar (2) has a bounded number of nonneutral letters and
only the three single-letter reservoirs0*,1*,3*. Any infinite family within
one fibre must therefore have unbounded neutral padding along a subsequence.
A putative fixed nonconstant word period cannot remain inside the fibre at
all repetition depths: its repeated changes of letters cannot be absorbed
into three constant reservoirs with bounded cores. A period containing2
cannot remain either. This is the physical explanation of the all-height
collapse, not merely a fitted recurrence of numerical counts.

For the threshold language, the finite union over counts<=N gives the same
conclusion. Its minimal DFA, and the exact-fibre minimal DFA, have only
singleton strongly connected components. To see the last statement carefully,
a nontrivial recurrent component would contain two different states and a
closed walk. If a closed walk has two different letters or the letter2, its
arbitrarily repeated accepted context contradicts (2). Otherwise every closed
walk uses one constant letter in{0,1,3}. A cycle of period greater than1 would
force acceptance to depend periodically on the length of that constant run.
But for a fixed surrounding word, induced insertion makes H nondecreasing,
so equality to N and the threshold predicate each occupy an interval of run
lengths and cannot alternate periodically. Distinct states on such a unary
cycle are consequently right-equivalent and are merged in the minimal DFA.
Only the dead state may carry all four loops without an accepting continuation.

This supplies a practical structural diagnostic for a computed target
certificate: after minimization its productive SCCs cannot hide an additional
counter or recurrent mixing mechanism. It does NOT mean target automata have
bounded size uniformly in N.

## 5. What the different old subjects contributed

* The old additive-basis/Zeckendorf discussion
  `07-reflections/goldbach-additive-basis-normal-forms-s494.md` distinguishes
  abundance from a canonical legal representation. Here the canonical object
  is the SCC/core/clock decomposition, and abundance is exactly (4), not a
  promise of additional target values. The new Dirichlet convolution is
  multiplicative because SCC path counts multiply; replacing it by additive
  convolution would lose the actual graph operation.
* THM-4057 `stern-brocot-depth-pullback-and-rational-edge-tournament-gauge`
  retains scale when identifying rational objects with discrete paths. Its
  useful methodological transfer here is to retain the neutral clock rather
  than quotient it out before asking about continuation. No Stern/Berggren
  theorem is imported into tournament coverage.
* THM-3183 `factorial-hecke-lattice-square-and-oriented-wedge-continuant`
  distinguishes a rationally compressed state from its positive/integral
  realization. Here this distinction selects the positive36-state carrier,
  not the smaller signed linear quotient, for exact saturation. The map is
  state compression; the lost predicate is nonnegative truncation validity;
  the sidecar is an actual positive run-assignment state.
* Session153 `rewriting/NOTE.md` proves that local affine residue coverage can
  lose supply on return. The target automaton keeps the entire successful
  continuation language and detects that loss exactly for each N. Its finite
  construction does not itself supply a uniform-in-N return theorem.

These are mechanism-level analogies with explicit retained and lost data,
not claims that these mathematical systems are isomorphic.

## 6. Artifacts and the remaining constructive obligation

`capped.py` constructs and minimizes complete positive closures. `fibres.py`
constructs every permitted strong core under the exact integer Busch cutoff,
forms (3), and independently compares all words through length7 and full
subset-DP through length4. `compare.py` compares exact/threshold continuation
classes across target changes rather than only tabulating attained values.

At N21 the closure has1469 live capped states,76 simultaneous count-Moore
states,42 threshold-DFA states, and a one-state rejecting exact-H21 DFA. This
is a finite proof about ALL word lengths. The smaller threshold and exact
machines must not be confused with the76-output machine.

The missing positive theorem is now sharply located: generate the finite
strong-core numerators `C_N,A_N` recursively in the TARGET N with a verified
nonempty branch for every allowed odd N, while preserving the actual apex
interfaces. A per-N finite procedure, rational fibre formula, state count,
or large multiplicity does not supply this target-uniform positivity.
