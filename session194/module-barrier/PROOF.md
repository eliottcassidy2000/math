# A sharp intact-module threshold for ternary return

Session 194, 2026-11-05. **PROVED elementary theorem; finite controls below.**
The unrestricted tournament Hamilton-path spectrum outside `{7,21}` is still
**OPEN**. This is an architecture obstruction, not another absent count.

## 1. Inheritance and the changed question

The Session 193 complete marked-clock operator
(the symplectic packet linked in the published [Session 193 guide](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/00-navigation/S193.md), Theorem 1) can be iterated on full path-cover
polynomials, but its transported fan has increments `2H(S),4H(S)` rather than
literal `2,4`. Its Proposition 2 further rules out literal ternary children
when the intact clock has profile content greater than one. That divisor
test leaves primitive clocks, notably the H13 clock, untouched.

This packet closes that loophole for **every context retaining the whole
input as a module**, including primitive inputs and arbitrarily complicated
contexts. It uses positivity and a small induced cyclic context, not another
content calculation. It does not constrain constructions that split the
input module, modify its internal arcs, use its vertices as separate
substitution factors, or reconstruct a different witness.

The closest inherited identity is the C3 substitution kernel, THM-3121,
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
We prove the needed specialization directly. The hostile is the primitive
five-vertex H13 tournament. The corrected near miss is that content one
would suffice to permit small affine returns. The useful sidecar is the
positive two-path-cover coordinate, rather than its gcd with H.

## 2. Definitions and exact cyclic minimum

Let H(T) count directed Hamiltonian vertex words, with H(empty)=1. A nonempty
vertex set M is a **module** if each vertex outside M either beats all of M
or loses to all of M. Write S=T[M], n=|S|, h=H(S). Contract M to one vertex
v, retaining each other vertex separately; this gives a quotient tournament Q.

Let F2(S) count ordered pairs (P,R) of nonempty directed paths partitioning
V(S). Explicit row lists below are outgoing bitmasks on named vertices
0,...,n-1: bit j of row i is one exactly when i->j. With P_S(d) the labelled-bin profile used in Session 193,

    F2(S)=P_S(2)-2h
         =sum_(empty != A proper-subset V(S)) H(S[A]) H(S[V(S)-A]).

Every induced tournament has a Hamiltonian path (the usual insertion proof),
so F2(S)>=2^n-2. Splitting each Hamilton word at each internal position gives
another injection, proving F2(S)>=(n-1)h. Neither bound is substituted for
the exact F2 in the theorem below.

**Lemma 1 (exact two-singleton cyclic extension).** Let C(S) have three
blocks S,{a},{b}, with S->a->b->S. Then

    H(C(S))=3h+F2(S)=h+P_S(2).                         (1)

Furthermore C(S) is strong for every nonempty S.

**Proof.** Hamilton paths using S in one maximal run have exactly the three
block orders S,a,b; a,b,S; and b,S,a, contributing 3h. A path using S in two
runs must be P,a,b,R, and is uniquely an ordered nonempty two-path cover of S.
There is no third run because only two exterior vertices exist and b->a is
absent. This proves (1). Every S vertex reaches a then b then every S vertex,
which proves strongness. QED.

**Theorem 2 (module-context dichotomy and sharp threshold).** For every
tournament T containing S as a module, exactly one of the following
structural cases occurs.

1. The quotient vertex v is a singleton strong component of Q. Then S is
   traversed contiguously in every Hamilton word of T, and h divides H(T).
2. The quotient vertex v lies in a nontrivial strong component of Q. Then

       H(T)>=3h+F2(S)=h+P_S(2).                       (2)

The threshold in the second case is exact for each individual S: C(S)
attains it. The numerical alternatives need not be disjoint; the structural
cases are disjoint.

**Proof.** Strong components in a tournament are totally ordered. In case 1,
all components preceding v beat S and all those following v lose to S.
Every Hamilton word has a Hamilton word of S as one contiguous block, with
independent preceding and following paths. Hence H(T)=H(X)hH(Y).

In case 2, v belongs to a directed triangle in its strong component. Here is
an elementary proof of that fact: v has nonempty out-neighbor set B and
in-neighbor set A inside the component. If there is no arc B->A, no path from
B can return to A or v, contradicting strongness. Hence some b in B and
a in A give v->b->a->v. The induced subtournament on S union {a,b} is C(S).

Hamilton paths in any induced subtournament inject into Hamilton paths of
the whole tournament: insert each missing named vertex at the first legal
position, in a fixed order. A legal position always exists, and deleting
the inserted vertices recovers the input path. Therefore H(T)>=H(C(S)).
Apply Lemma 1. Its explicit strong extension proves sharpness. QED.

This proof uses neither a strong-order cutoff nor Camion's theorem, a finite
tournament census, or an assumption that S itself is strong.

## 3. No literal ternary return retains a nontransitive input module

**Corollary 3.** If h=H(S)>1, there is no tournament T containing S as a
module with H(T)=3h+2 or H(T)=3h+4.

**Proof.** Hamilton-path parity for tournaments gives h odd, so h>1 implies
h>=3 and n>=3. In case 1, h divides H(T), but an odd h>1 cannot divide 2 or
4. In case 2, F2(S)>=2^n-2>=6, so H(T)>=3h+6. QED.

Only the familiar oddness theorem is used in the divisibility branch. If
one prefers not to invoke parity, the conclusion applies verbatim whenever
h is coprime to 2; this includes every tournament input.

More generally, every context value below 3h+F2(S) is a multiple of h. Thus
there is an exact input-dependent gap before a cyclic context can change
the multiplicative arithmetic. This gap persists after arbitrarily many
substitution steps whenever the original S remains an intact final module.
It also persists if the context contains a new isomorphic copy of S rather
than the original labels.

**Necessary exception h=1.** A tournament has H=1 precisely when transitive.
For any n, `TT_n -> D5`, with D5=(2,12,9,1) and H(D5)=5, retains TT_n as a module and has
H=5=3H(TT_n)+2. The divisibility alternative therefore cannot exclude h=1.
If the quotient component containing TT_n is nontrivial, however, (2) still
gives the sharp threshold 1+2^n. For n>=3 this exceeds both 5 and 7; for
n=2, C(TT_2) attains 5; for n=1, C(K1) attains 3. We claim no h=1 universal
exclusion and do not infer anything about the globally forbidden 7 here.

**Primitive hostile.** Take

    D=(2,12,9,17,7),       H(D)=13, F2(D)=78.

Its profile content is one: H(D)=13 and P_D(2)=104 alone do not show this
(their gcd is 13), but P_D(3)=471 is coprime to 13. The verifier computes the
actual bin values rather than assuming content one. Any context in which D
lies cyclically has at least 117 Hamilton paths. Every smaller context count
is a multiple of 13. Hence 41 and 43 are impossible in **all** intact-module
contexts, despite primitive profile. These scalar targets are attainable by
other constructions; they are not new spectrum gaps.

## 4. Consequence for the reusable grammar

S193's profile-content obstruction was necessary but not sufficient: even
g(S)=1 never permits a literal ternary step on a nontransitive intact input.
The complete polynomial state and rank-four substitution compiler remain
correct, but no choice of later cores, ports, or finite/unbounded nesting
repairs this problem while retaining S as one module.

The workable interfaces must therefore expose **individual old vertices**,
change internal old arcs, delete/replace an input, or use independent
witnesses for different branches. The S191 physical fan does expose vertices
of its changing core; it retains only a transitive clock as a module and
does not contradict Corollary 3. Likewise S193's H41 old-arc/birth repair
and its branch-specific product remain available.

## 5. Verification

`verify.py` is a standalone standard-library certificate. It enumerates all
labelled input tournaments on one through four vertices, plus the explicit
H13 input and transitive five-vertex exception. For each it checks every
labelled quotient of order one through four, substituting the input at the
fixed quotient vertex zero. No isomorphism or value filtering is used.
The theorem, including arbitrarily large inputs and contexts, is proved
above; this finite universe is a hostile/positive control only.

The verifier separately uses literal permutation enumeration for the cyclic
extension identity on every input; subset endpoint DP for context counts;
direct subset splits for F2; and independent graph reachability for quotient
SCC membership. It checks both quantitative F2 bounds, the exact minimum,
the primitive H13 profile, the transitive exception, and the ternary
exclusion in the stated finite universe. Checks use exceptions under both
normal Python and `python3 -O`.

Reproduce from this directory with `python3 verify.py` and
`python3 -O verify.py`. Compare stdout and `receipt.json` bytes.
