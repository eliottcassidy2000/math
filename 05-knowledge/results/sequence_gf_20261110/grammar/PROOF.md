# Recurrence compilers: which modules grow, and which generating function remembers support?

**2026-11-10. PROVED IN THIS NOTE, with exact new-formula examples.**
No independent canon audit or mathematical-priority claim. Chan–Pak remains
OPEN. The companion is `compiler.py`; its output is `results.json`.

## 0. Inheritance and the question being changed

The inherited mechanisms are already strong:

- **THM-3121**, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`,
  gives the complete quotient-walk / ordered-path-cover substitution law.
- **THM-3166**, `01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md`,
  makes the path-colour polynomial multiplicative under order join and
  proves finite exponential sums at each fixed path-cover depth.
- **THM-3181**, `01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`,
  gives half-grid recovery and fixed extreme-depth recurrences.
- **THM-3202**, `01-canon/theorems/THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction.md`,
  proves that making ALL THREE cyclic factors grow destroys C-finiteness.
- **THM-3248**, `01-canon/theorems/THM-3248-q4-paired-owner-stirling-compiler.md`,
  retains the complete walk numerator, not merely its denominator, and
  explicitly separates observer symmetries from tournament automorphisms.
- **THM-4878** and
  `05-knowledge/results/spine_grammars_cutwidth_apexes_hankel_seams_20261104.md`
  already give nonnegative finite-state representations for bounded
  spine grammars, with the one-apex four-coordinate transfer inherited
  from earlier sessions. Their finite coverage is not a completeness theorem.
- Session24,
  `05-knowledge/results/session24_research_synthesis_20261021.md`, already
  gives rational fixed-H **unlabelled-class** multiplicities using ordered
  SCCs. We do not rebrand that theorem as new.
- Session248,
  `05-knowledge/results/session248_recursive_symmetry_20261110/SYNTHESIS.md`,
  distinguishes symmetry of a physical tournament, a marked hierarchy,
  and an observer. Its free action on ordered path covers supplies the
  symmetry charge in §4 below.

Closest mechanism: fixed-depth Newton inversion plus the walk kernel.
Canonical hostile: `C3[P_r,P_r,P_r]`, whose singleton-cover contribution
is factorial-cubic. Corrected near miss: a rational quotient-walk kernel
does NOT make every sequence obtained by growing its factors C-finite.
Least-used relevant sidecar: the subset of quotient vertices whose factors
actually grow; the fixed outside vertices bound how often they can be revisited.

The active concept board has five entries:

1. Growing quotient / transitivity / bounded runs versus factorial returns.
2. Path-colour evaluations / Newton inversion / finite exponential spectrum.
3. Marked grammar / nonnegative matrices / noncommuting operation order.
4. Value support / positive catalytic series / exact fixed-target automata.
5. Module symmetry / free ordered-cover action / arithmetic divisibility.

This note applies the repository's “Keep observer, recurrence class, and
finite head”, “Controlled forgetting and unlabeled quotients require a
sidecar”, and “Separate local support from bounded-height coverage” patterns.
The new connections are explicitly mapped in §9.

## 1. Definitions

Let `H(T)` count literal directed Hamiltonian paths of a nonempty tournament.
Let `pc_T(c)` count unordered spanning covers by `c` nonempty directed paths,
and set

    F_T(c) = c! pc_T(c),
    Q_T(t) = sum_c F_T(c) binom(t,c).

Thus `Q_T(1)=H(T)`. All `F_T(c)` are nonnegative integers.
For a nonempty seed `A` and `r>=1`, let `A^{▷r}` be the ordered join of
`r` disjoint copies. The inherited identities are

    Q_(A^{▷r})(t) = Q_A(t)^r,
    F_(A^{▷r})(c) = sum_(j=1)^c (-1)^(c-j) binom(c,j) Q_A(j)^r.   (1)

Formula (1) remains valid when `c>|A|r`: it then equals zero.
The omitted `j=0` term is zero because `r>=1` and `Q_A(0)=0`.

Fix a quotient tournament `Q`, a nonempty subset `G` of its vertices, and
one nonempty seed `A_i` for each `i in G`. At the other quotient vertices
place fixed nonempty tournaments `B_v`. Define

    T_r = Q[ A_i^{▷r} for i in G; B_v for v outside G ],  r>=1.    (2)

Let `g=|G|` and `b=sum_(v outside G) |B_v|`. The growing subset `G`, not
the entire quotient, is the diagnostic.

## 2. The exact growing-quotient dichotomy

**Theorem 1 (PROVED IN NOTE).** For (2):

**A.** If `Q[G]` is transitive, there are explicitly computable integers
`beta_j`, indexed by `j=(j_i) in {1,...,b+1}^g`, such that

    H(T_r) = sum_j beta_j ( product_i Q_(A_i)(j_i) )^r, r>=1.      (3)

Consequently the sequence has a rational ordinary generating function and
a constant-coefficient recurrence of order at most `(b+1)^g`.

**B.** If `Q[G]` is not transitive, then

    H(T_r) >= 3 (r!)^3,  r>=1.                                 (4)

It has no eventual constant-coefficient recurrence and no rational
ordinary generating function.

Thus for this precise class, transitivity of the growing quotient is an
if-and-only-if criterion for C-finiteness. The theorem does NOT assert an
if-and-only-if criterion for D-finiteness / polynomial-coefficient recurrences.

### Proof of A: outside vertices are separators

Expand every fixed `B_v` into its actual vertices, retaining their arcs.
Keep each growing block as one quotient label. For a Hamiltonian path of
`T_r`, contract each maximal consecutive run lying in one growing block
to its block label. Keep each of the `b` fixed vertices individually.

Two consecutive occurrences of the same growing label cannot occur:
they would be one maximal run. Delete the fixed vertices from the word
by cutting at them, not by reconnecting across them. The resulting at
most `b+1` intervals consist only of growing labels. Within each interval,
every step is an arc of the transitive tournament `Q[G]`, hence strictly
increasing in its unique order. Therefore each growing label occurs at
most once per interval, and its TOTAL occurrence count satisfies

    1 <= c_i <= b+1.                                           (5)

Let `w_c` be the number of these quotient words using each fixed actual
vertex once and growing label `i` exactly `c_i` times. It is a finite
nonnegative integer, computed solely from the fixed context.
For each word choose an ordered `c_i`-path cover of growing block `i` and
assign its paths to that label's occurrences in encounter order.
Concatenation is legal, and maximal-run splitting is its inverse. Thus

    H(T_r) = sum_(1<=c_i<=b+1) w_c product_i F_(A_i^{▷r})(c_i).   (6)

This is a finite specialization of THM-3121, with the boundedness now
explained by growing-quotient transitivity. Inserting (1) gives (3), with

    beta_j = sum_(c_i>=j_i) w_c
                    product_i [(-1)^(c_i-j_i) binom(c_i,j_i)].   (7)

The compiler enumerates quotient words by `(used fixed-vertex mask,
growing occurrence counts, last label)`. It stops by word length
`b+g(b+1)`. It does not enumerate Hamiltonian paths of the large `T_r`.

### Proof of B: the growing cycle can avoid all fixed context

A nontransitive tournament contains a directed triangle. Take its three
vertices in `G`. In each corresponding `A_i^{▷r}`, select one vertex
from each copy. The resulting induced subblock is `P_r`, the transitive
`r`-vertex tournament. Hence `T_r` contains induced `C3[P_r,P_r,P_r]`.

In that induced tournament, alternate cyclically among the three blocks,
using each block vertex as a singleton run. There are three initial block
choices and `(r!)^3` assignments of actual vertices to the occurrences.
This gives `3(r!)^3` distinct Hamiltonian paths.

Hamiltonian-path count is monotone under adding vertices: insert each
new vertex into each old path by a fixed deterministic legal insertion
rule. Distinct old paths remain distinct because deleting the inserted
vertices recovers the old order. Thus the bound passes to `T_r`.
Every eventual C-finite sequence grows at most polynomial times
exponential, contradicting (4). This also excludes a rational OGF. QED.

### Boundary and exact recurrence

- `b=0`: the growing quotient must itself be transitive, and (3) reduces
  to the familiar `H(T_r)=(product_i H(A_i))^r`.
- `g=1`: every fixed context gives C-finiteness, even if the whole
  quotient is highly cyclic.
- `g=2`: the growing quotient is always transitive. A THIRD growing
  module is the first possible source of the factorial obstruction.
- Equal exponential bases in (3) must be combined and zero coefficients
  removed. Write the result `sum_(lambda in L) B_lambda lambda^r`.
  Its minimal characteristic polynomial is exactly
  `product_(lambda in L) (X-lambda)`, by Vandermonde independence.
- With `G_T(z)=sum_(r>=1) H(T_r)z^r`,

      G_T(z)=sum_(lambda in L) B_lambda lambda*z/(1-lambda*z).   (8)

- Nothing in (3) automatically extends to `r=0`. Empty growing blocks
  change the quotient. The omitted `0^0` Newton terms then matter.
  One may separately prepend the actual `r=0` value to the OGF.

## 3. Three exact examples and a multivariate extension

For `C3[A^{▷r},1,1]`, the only kernel weights are `w_1=3,w_2=1`.
Therefore

    H = Q_A(1)^r + Q_A(2)^r.                                (9)

Taking `A=C3`, where `Q_A(t)=t^3+2t`, gives

    H_r=3^r+12^r,
    H_(r+2)=15 H_(r+1)-36 H_r,
    H_1,H_2,H_3 = 15,153,1755.

For `C3[A^{▷r},B^{▷r},1]`, the kernel weights on `(c_A,c_B)` are
`3,1,1,1` at `(1,1),(1,2),(2,1),(2,2)`. If
`a=Q_A(1),a'=Q_A(2),b=Q_B(1),b'=Q_B(2)`, then

    H_r=(a'b')^r-(a'b)^r-(ab')^r+3(ab)^r.                 (10)

For `A=B=K1` this is `4^r-2*2^r+3`, with sequence
`3,11,51,227,963,...` and characteristic polynomial
`X^3-7X^2+14X-8`.

For `A=C3,B=K1` it is `24^r-12^r-6^r+3*3^r`, beginning
`15,423,11961,309987,...`, with characteristic polynomial
`X^4-45X^3+630X^2-3240X+5184`.

The third growing singleton module changes the behavior completely:
`H(C3[P_r,P_r,P_r])=3,45,2721,421425,...`, with the factorial floor (4).
This is the inherited THM-3202 hostile, used to locate the exact boundary.

The same proof allows independent positive repetition parameters `r_i`:

    H(T_(r_i)) = sum_j beta_j product_i Q_(A_i)(j_i)^(r_i),

with multivariate OGF

    sum_j beta_j product_i [Q_(A_i)(j_i) z_i/(1-Q_(A_i)(j_i)z_i)]. (11)

Thus one fixed-context compiler supplies whole multisequences without
guessing recurrences from sampled terms.

## 4. Symmetry charge survives the sequence compiler

**Corollary 2 (PROVED, inherited free-action mechanism).** Let `K_i` be
any actual automorphism subgroup of `A_i` and `D=product_i |K_i|`.
Then every `Q_(A_i)(j)` is divisible by `|K_i|`, and

    D^r divides H(T_r).                                    (12)

Proof: an automorphism fixing an ordered path cover fixes each ordered
vertex word, hence every vertex. The action on ordered covers is free,
so `|K_i|` divides every `F_(A_i)(c)`. The binomial evaluation of `Q`
preserves that divisibility. Formula (3) gives (12) in the transitive
case; in general the product group acting independently on the copies
extends to `T_r` and acts freely on its Hamiltonian paths.

This is a direct use of Session248's arithmetic symmetry charge, not a
claim of a newly found group action. In the mixed example, all bases are
divisible by3 and `H_r=3^r(8^r-4^r-2^r+3)`. A family can have rapid,
beautiful recurrence growth while being arithmetically poor for covering
odd targets. Observer coordinate symmetries, such as THM-3248's paired
owners, cannot be substituted for actual `K_i` here.

## 5. Why a finite list of these recurrences cannot cover the spectrum

**Corollary 3 (PROVED).** A finite union of the one-parameter transitive-
growing families in Theorem 1 contributes only `O(log X)` distinct
Hamiltonian-path values at most `X`.

After combining bases, either the largest base is1, giving a constant
sequence, or it is an integer `lambda>1`. Its coefficient is positive,
because `H(T_r)>0` for every `r`. Therefore
`H(T_r)~B_lambda lambda^r`, and only `O(log X)` indices can contribute.
Finite unions preserve the bound.

This is not a prohibition on recurrence methods. It says where to change
the carrier: allow a growing context or a BRANCHING grammar, rather than
another finite catalog of periodic rays. Within a fixed ordinal module,
different seed letters act diagonally via `Q_A(j)` and commute. The
ordered seed word is invisible to this particular outside context. A
marked boundary whose matrices do not commute can retain the word order.

The exact next object should therefore be a positive series that remembers
the output VALUE, not just a rational series summing values by word length.

## 6. A general positive support-series compiler

Suppose an inherited finite grammar has deterministic control states `q`, legal letters `a`,
and a nonnegative integer column state `v in N^d`. A letter acts by
`v -> M_a v`; affine maps are included by appending a constant coordinate.
The final observable is `h=ell^T v` for a nonnegative integer `ell`.
The grammar control may be carried separately or incorporated as blocks.

For each control state define the catalytic series

    F_q(z;X)=sum_(legal words ending at q) z^|w| X^(v(w)),
    X^v=product_i X_i^(v_i).

For a transition `p --a--> q`, define the monomial substitution

    Phi_a(X)_j = product_i X_i^(M_a[i,j]).

Then the EXACT formal functional equations are

    F_q(z;X)=seed_q(X)+z sum_(p --a--> q) F_p(z;Phi_a(X)).       (13)

Proof: `X^(M_a v)=Phi_a(X)^v`; decompose every nonempty word by its
last letter. There is a unique solution in the `z`-adic formal series
ring, obtained recursively by word length. All coefficients are
nonnegative, so no cancellation can hide a witness.

Putting `X_i=u^(ell_i)` gives

    P(z,u)=sum_w z^|w| u^(h(w)).                              (14)

Now the target `h` is attained if and only if some coefficient
`[z^n u^h]P` is positive. This is a literal support-preserving map. It is
not the sum-of-values series `sum_w z^|w| h(w)`, which retains only one
moment of the value distribution. When the substitutions include powers
such as `y->y^2`, (13) is a multivariate Mahler-type system; no claim
that all such systems have the classical univariate Mahler form is needed.

### Fixed-value fibers are always regular in a nonnegative representation

**Theorem 4 (PROVED).** Fix an integer target `h>=0`. Replace every
coordinate of every state vector by `min(value,h+1)` after each letter.
The resulting finite deterministic automaton accepts EXACTLY the legal
words with output `h`. Hence their length OGF is rational, with an
effective constant-coefficient recurrence and an effective witness search.

Proof: clipping at `h+1` commutes with addition and multiplication by
nonnegative integers after clipping the result. Thus it commutes with
every matrix update and with the final nonnegative linear readout.
The clipped readout equals `h` iff the true one does. There are at most
`|control states|*(h+2)^d` vector/control states. Reachability supplies
the finite transition matrix `R_h`, and the fiber OGF is

    e_start^T (I-z R_h)^(-1) f_h.                            (15)

The actual reachable set may be much smaller. Breadth-first search gives
a shortest legal word, which the inherited grammar decodes into an actual
tournament. This constructs witnesses without enumerating all tournaments.
If an accepting word exists, one has length at most the number of reachable
states minus1: a shortest accepting path cannot repeat a state. For a fixed
grammar this is a pseudo-polynomial decision-and-construction algorithm
in the numerical target h, not a polynomial-time claim in its bit length.

**Necessary boundary:** use a NONNEGATIVE integer realization, not an
arbitrarily rank-minimized rational/signed representation. Negative entries
can cancel large hidden coordinates, making this clipping invalid.

For a proposed Chan–Pak grammar, (13) is a positive all-target support
problem and (15) is an exact finite-target decision-and-witness compiler.
Uniform coverage still requires proving nonemptiness for every allowed
odd `h`; a succession of finite automata is not that proof.

## 7. Explicit one-apex Mahler system and rational fibers

Use the inherited one-apex family: start with a transitive chain with
vertices in increasing order, then add an apex `x`; bit `b_i=1` means
`x->v_i`, and `b_i=0` means `v_i->x`. A word of length `n` represents an
`n+1`-vertex tournament. This counts marked WORDS, not isomorphism classes.

The standard transfer can be written in three affine coordinates:

    initial (H,S,p)=(1,0,0),
    bit0: (H,S,p) -> (H+S+p, 2S, 0),
    bit1: (H,S,p) -> (H,     2S+1,1).                     (16)

Here `p` is the previous bit and
`S=sum_(i<j) b_i 2^(j-i-1)` immediately before reading position `j`.
The matrices after appending a constant coordinate are nonnegative.
The two matrices do not commute: words `10` and `01` give H3 and H1.

One derivation of (16), using inherited odd-cycle collections, is that
every cycle must pass through the apex, so at most one odd cycle can be
chosen. A new bit0 closes cycles beginning at previous bit1 positions.
Their total contribution to H is2 for an adjacent `10`, and `2^(j-i-1)`
for separation `j-i>=2`; this is exactly `S+p`. Equivalently, this is a
compressed form of the previously proved one-apex transfer.

Let `F_p(z;x,y)` count words ending in bit `p`, with monomial `x^H y^S`.
The empty word belongs to state0. Formula (13) becomes

    F_0(z;x,y)=x+z F_0(z;x,xy^2)+z x F_1(z;x,xy^2),
    F_1(z;x,y)=z y F_0(z;x,y^2)+z y F_1(z;x,y^2).         (17)

The positive support enumerator is `F_0(z;x,1)+F_1(z;x,1)`.
This retains `H` exactly, unlike a length-only mean recurrence.

### Every one-apex value fiber has an especially small denominator

For `h>1`, remove all initial0s and all terminal1s from a word with
`H=h`. These are respectively source and sink padding, and preserve H.
The remaining unique core starts with1 and ends with0.
For a core of length `n`, its last zero contributes at least `2^(n-2)`,
so

    h >= 1+2^(n-2),
    n <= 2+floor(log_2(h-1)).                              (18)

Let `C_h(z)` be the finite polynomial counting those cores by length.
Then the exact all-length marked-word fiber is

    sum_(n>=0) #{words length n with H=h} z^n
      = C_h(z)/(1-z)^2,          h>1.                       (19)

For `h=1`, the words are exactly `0^a1^b`, so the answer is
`1/(1-z)^2`. Uniqueness of trimming prevents overcounting in (19).
This is a marked-family application of the general support compiler,
not a rediscovery of Session24's unlabelled SCC multiplicity theorem.

The exact companion finds

    C_3=z^2, C_5=2z^3, C_9=3z^4, C_11=z^4,
    C_13=0, C_15=2z^5, C_21=0, C_25=z^6.

Thus missing13 is certified for this family at ALL orders, not inferred
from a finite prefix of vertex sizes. The general spectrum does contain13;
the small family is insufficient. Direct one-apex support through511 has
145 values, without adding order-join products. This is a reproduction/
example, not a new census record.

Finally, along the periodic word `(10)^r`, (16) gives

    H_r=(4^(r+1)+6r+5)/9, r>=0,
    sum_(r>=0) H_r z^r=(1-3z)/((1-4z)(1-z)^2).             (20)

The sequence starts `1,3,9,31,117,459,1825,...`. The repeated pole at1
comes from a Jordan block in the marked affine grammar. The pure-
exponential conclusion of Theorem1 does not apply to this different,
nonuniform attachment operation.

## 8. What this does and does not give toward Chan–Pak

**Constructive progress:**

1. A fixed substitution description now tells us BEFORE fitting terms
   whether rational GF / constant-coefficient recurrence search is lawful.
2. In the lawful case, formulas (6)-(8) generate exact coefficients,
   the minimal denominator, and the recurrence from small seed profiles.
3. For bounded-interface grammars already available in the repository,
   (13) converts their nonnegative matrices into a positive all-values
   generating object. Formula (15) converts each target into a terminating
   exact decision and actual tournament-witness search.
4. Theorem4 can be applied to the inherited cutwidth2 or two-apex
   nonnegative representations. A successful all-target proof should give
   a uniform reduction or invariant for these positive equations, not
   conclude support coverage from their moments, spectral radii, or many
   periodic examples.
5. The symmetry charge (12) lets a compiler reject arithmetically
   unsuitable seed palettes for a proposed target before expanding words.

**Still open:** a uniform proof that a particular realizable branching
grammar has support every odd integer except7 and21. The automata for
different h do not come with a uniform accepting-state induction here.
No claim of density or all-value coverage is made for a new family.

**A000568 is a different sequence.** It counts unlabelled tournaments by
vertex order. In (3) one parameter indexes growth of one explicit family;
in (14) words are counted with their grammar presentations; in (19) the
value is fixed and presentation length varies. Quotienting presentations
to unlabelled tournaments requires unique decomposition and orbit counts.
Repeated values or repeated encodings are harmless for a support witness,
but cannot be treated as A000568 coefficients.

## 9. Typed connections and next cheap tests

- **Fixed context -> recurrence.** Source: module substitution. Map:
  maximal-run splitting then Newton inversion. Preserves literal H.
  Loses internal isomorphism information. Sidecar: which blocks grow and
  their `Q_A(j)` values. Test: the three finite-context examples versus
  literal vertex subset DP, with the third-growing-module hostile.
- **Nonnegative grammar -> support series.** Source: exact path-state
  matrices. Map: `v -> X^v` and monomial substitutions. Preserves every
  output value and legal word. Loses nothing about value support, but
  does not remove presentation multiplicity. Test: (17)'s coefficients
  versus direct one-apex counts and capped fibers.
- **Symmetry -> arithmetic filter.** Source: actual module automorphisms.
  Map: free action on ordered covers. Preserves divisibility, destroys
  subgroup identity after more general symmetry-breaking operations.
  Sidecar: actual group, not observer symmetry. Test: C3 seed gives
  `3^r` divisibility in both compiled examples.
- **Fixed-target finite state -> all-target proof obligation.** Source:
  exact capped automata. Map: accept output h, then decode an accepted
  word. Preserves realizability. Lost uniformity: state space depends
  on h. Cheapest useful next test: generate the same support-series
  equations from the already-proved positive cutwidth2 or two-apex
  matrices and look for a finite residue/inequality classification of
  accepting-state witnesses. Merely minimizing their signed Hankel
  representation is not a safe substitute.

## 10. Reproduction and scope of computation

Run from the repository root:

    python3 05-knowledge/results/sequence_gf_20261110/grammar/compiler.py

The dependency-free program computes the context kernels and exponential
coefficients, derives characteristic polynomials, and checks only the new
formulas against an independent vertex subset DP. It includes nine context
instances up to13 vertices, the all-growing cyclic control through12
vertices, and all127 one-apex incidence words of length at most6. It then
compares nine exact capped-fiber automata to the independently enumerated
core polynomials, and evaluates the periodic example.

These are cheap positive and hostile controls for newly derived claims,
not an audit of previously proved repository results. No new canon theorem number is assigned in this note.
