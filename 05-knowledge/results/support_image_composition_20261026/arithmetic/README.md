# Residue clocks: exact finite path lifting and the missing integer guard

**Status: PROVED elementary clock and guard results; independently audited
[PASS](../audit_arithmetic/PROOF_AUDIT.md). Computation FINITE-EXACT.** The parity-cylinder mechanism is classical; the
contribution here is its explicit arbitrary-modulus arrow interface, the
odd-prime versus dyadic refinement criterion, and the bounded-source guard.
No claim of Collatz convergence, twin infinitude, or a new complexity
reduction follows.

## Inheritance and the live board

The closest mechanism is the actual-image quotient criterion in
[Session25, §§1–2](../../session25_spectrum_wiring_20261021/session25-paper/PAPER_CONNECTIONS.md),
together with the classical finite parity-vector bijection. The canonical
hostile is a quotient that joins two arrows using different representatives
of the middle state. The corrected near miss is treating an existential
residue relation as an induced deterministic map. The least-used sidecar is
the exact positive source inside a finite residue cylinder.

Our board is: integer orbit; labelled clock; finite path cylinder; prime
wheel; native time; positive-source guard. Anchor: exact composition.
Niche: distinguish spatial prime refinement from temporal dyadic precision.
Wildcard: the same signed boundary −1 in two different arithmetic carriers.
The research cards used are “Separate local support from bounded-height
coverage” and “Controlled forgetting and unlabeled quotients require a
sidecar”; no new maintained card is proposed.

Recent ancestry is explicit:

- [THM-4808, branch character and every-modulus obstruction](../../../../01-canon/theorems/THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md)
  is an incoming **pre-audit** result. We use only the independently proved
  elementary finite-residue family in §6, not its
  spectral claims. [HYP-9418](../../../hypotheses/HYP-9418-legal-path-sub-observers-berggren-sieve-levels-and-the-branch-character.md)
  remains OPEN: no density model or spectral extrapolation is a dependency.
- [Growing sieve W1–W3](../../growing_sieve_path_20261026/sieve/README.md)
  supplies the prime-equality exceptions and height-completeness gate.
- [The F8/abc observer](../../growing_sieve_path_20261026/field/README.md)
  preserves literal Paley tournaments after retaining its determinant bit,
  but loses the abc-hit predicate. [Joint Frobenius phase](../../strong_atom_continuation_vitali_20261026/FROBENIUS-ABC.md)
  is a different sidecar; neither bit is identified with Collatz parity.
- [Root27 machines](../../root27_machine_bridges_20261026/anchor/README.md),
  C3–C5, already has exact carry cylinders, ordered composition, and an
  OPEN fixed-source certificate condition. Its ceil(3z/2) map and macro
  credit are different from the shortcut Collatz map below.

## 1. The precise clock is a relation

Let C be shortcut Collatz on integers:

    C(n)=(3^e n+e)/2,  e=n mod2 in {0,1}.

On positive integers this is a total map. For M≥1 define a labelled clock
graph on Z/M by

    r --e--> s  iff  2s = 3^e r+e (mod M).              (C1)

Vertices are residues, not time indices. Loops and distinct arrows with
the same endpoints are allowed. This graph is neither a tournament nor in
general a functional graph. Its relation is exactly the existential image
of one integer transition, even when the source is required positive.

**Proof.** Necessity is immediate. For sufficiency choose integer
representatives r,s. Write D=(3^e r+e−2s)/M. An integer n=r+Mk gives
the required arrow precisely when D+3^e k is even. This chooses exactly
one k modulo2, forces the required input parity, and gives C(n)=s moduloM.
Adding 2M preserves everything and supplies positive sources. QED.

Every vertex has exactly two labelled outgoing arrows. If M is odd each
label has one successor. If M is even, the congruence itself forces
e=r mod2, and there are two successors differing by M/2. Thus no residue
map C-bar on Z/M can commute with C for M>1: the outgoing relation has
more than one successor somewhere. For odd M>1 the two successors coincide
only at r=(M−1)/2; every other vertex has two. M=1 is the trivial exception.

At M=7, r=3 has two labelled arrows to s=5: source 10 uses e=0 and
source 3 uses e=1. Forgetting the label loses native elapsed time: ordinary
Collatz takes respectively one and two steps to 5. When M is even, the
vertex already determines e, so that particular loss is repaired.

## 2. Exact finite composition, including a source certificate

**C2 (all M≥1, L≥0).** Every length-L labelled path π starting at r has
as its set of integer sources exactly one residue class

    n=a_π (mod Q),  Q=M·2^L,  0≤a_π<Q.                 (C2)

Consequently there are exactly 2^L paths from each r, in bijection with
the 2^L lifts of r modulo M·2^L. Every finite clock path has positive
integer realizations, and the relational power R_M^L is exactly the
existential projection of genuine L-step integer orbits. The chosen
representative may depend on the entire path.

*Proof by lifting.* For L=0 the class is r moduloM. Assume a prefix of
length j has source class a modulo M2^j and contains o odd steps. Iterating
the branch formulas shows that changing the source by M2^j k changes the
j-th state by 3^o M k. To realize a prescribed next arrow, the one-step
calculation in §1 specifies one parity of this k, since 3^o is odd.
Exactly one lift modulo M2^(j+1) survives. Conversely every source in
that class has the whole specified trace. Induction proves all assertions.

This exactness is special, not a generic property of quotient relations.
For example f(a)=b, f(b)=b, f(b')=c, f(c)=c and q(b)=q(b')=B,
q(a)=A,q(c)=C produce A→B→C in the quotient relation. No f-orbit starting
over A realizes that two-step path.

For a parity word w of length L with o ones, write

    C_w(n)=(3^o n+b_w)/2^L.
    b_empty=0; b_(we)=3^e b_w+e·2^L.
    b_(uv)=3^(#ones(v)) b_u+2^|u| b_v.                 (C3)

These are ordinary ordered affine composition identities. They provide
the endpoint evaluator only on the word's parity cylinder; arbitrary
substitution must retain that guard. Ordinary, unshortened Collatz duration
along a realized word is exactly L+o. This is not a native Turing-machine
transition clock.

## 3. Odd refinement is a covering; dyadic refinement is not

Let M|N and project vertices and arrows moduloM.

**C4.** Every path in R_M has a unique path lift from every prescribed
starting residue above it in R_N if and only if N/M is odd, equivalently
v2(N)=v2(M). Here “covering” means this outgoing-arrow/path lifting
property; no assertion about incoming local degrees is intended.

*Proof when N/M is odd.* For a path π of length L, its source cylinder is
a modulo M2^L. A prescribed fine start is b moduloN, with b=a moduloM.
The gcd of N and M2^L is M. CRT therefore gives exactly one class modulo
N2^L in their intersection. C2 at level N gives the unique lifted path.

*Failure when N/M is even.* If M is odd, parity at the fine start fixes
one label and excludes the other coarse label. If M is even, both fine
successors have the same label and differ by N/2, which is a multiple
of M. They project to the same coarse arrow; the other coarse successor
has no lift from this fine start. QED.

In particular, adding an odd prime to a primorial wheel is a covering of
these clock graphs, whereas adding a dyadic precision bit is not. This
explains two different roles of a “larger modulus.” In the deterministic
finite-horizon diagram the time map instead runs between different levels:

    Z/(2^L M) --C--> Z/(2^(L−1) M),  L≥1.              (C4a)

This map is well-defined because adding 2^L M preserves input parity and
changes the output by an integer multiple of 2^(L−1) M. Odd-prime reduction
commutes with this time map. There is no time map at the same fixed level
unless one keeps the existential arrows.

## 4. Bounded realization and the exact sieve product

For a specified path π and positive height H the exact number of realizations is

    N_π(H)=floor((H−a_π)/Q)−floor(−a_π/Q).             (C5)

Thus a formal path need not have any source in a given finite interval.
If Q≥H there is at most one such source, which can be read off and tested.
Example: the e=0 loop at residue0 modulo3 has positive sources divisible
by6, so it has no realization at height H=5.

For a prime cutoff B≥2 let M_B be its primorial. On this clock, the strict
twin-wheel predicate gcd(r(r+2),M_B)=1 is a well-defined vertex predicate.
The corrected predicate E_B(n), allowing a divisible entry only when it
equals the tested prime, depends on the actual integer n. It does not
descend merely to residues: at B=3, sources3 and9 agree modulo6 but E_B(3)
holds and E_B(9) fails. The equality sidecar is necessary.

Let q be the least prime>B. From inherited W2–W3, if H+2<q² then
E_B(n) is equivalent to “n,n+2 are prime” on 3≤n≤H. Therefore the exact
twin-start path count is the count of n in the class C2, in [3,H], satisfying
E_B(n). When Q≥H it is decided by at most one candidate and those explicit
tests. This composes a finite path certificate with a complete finite-height
sieve; it does not say C preserves twins, nor infer an infinite twin set.

## 5. Every finite prefix lifts; the infinite positive source can fail

For every M and L≥1, the e=0 loop at residue0 repeated L times has source
class 0 moduloM2^L. The e=1 loop at residue−1 repeated L times has source
class −1 moduloM2^L. Their least positive sources are respectively M2^L
and M2^L−1. Thus every finite prefix is positively realized.

No positive integer realizes either infinite labelled loop. For all-zero
parity it would be divisible by every power of2, forcing n=0. For all-one
parity the identity

    C^L(n)+1 = 3^L(n+1)/2^L

forces n+1 divisible by every power of2, hence n=−1. These are actual
nonpositive integer fixed points; the obstruction is not that every
compatible path lacks an integer representative.

The growing strict twin wheels also admit the coherent class −1 at every
level. There it corresponds under the *different* Pythagorean map T(x) to
the degenerate signed triple (−1,0,1). This common boundary illustrates the
same quantifier failure, “for each finite depth some positive source”
versus “one positive source for every depth.” It is not a conjugacy of
Collatz with the sieve or the Berggren tree.

## 6. Connection contract and limits

| Source → target and map | Preserved | Lost | Guard / decisive test |
|---|---|---|---|
| Integer orbit → R_M, reduce states and keep parity | Every finite transition and labelled composition | Exact source, height | C2 cylinder and C5 interval count; all-one loop |
| R_N → R_M, reduce modulo divisor | Arrows; unique lifts exactly for odd N/M | Extra prime data; dyadic step choice | C4; compare 3→15 and 3→6 |
| Clock cylinder + E_B → twin-start traces | Exact finite-height primality when H+2<q² | Infinite realization or frequency | Height and prime-equality gates |
| Berggren word → vector modM | All branch transitions | Determinant sign if bit omitted | Root versus A D^(M−1) U D, as THM-4808(c) |
| Root27 machine carry word → affine endpoint | Guarded arithmetic macro endpoint | Native tape time, positivity/credit | Keep fixed-source cylinder and pending-slot guard |

The THM-4808 family has pair (16M+2,4M+1), congruent to root(2,1),
and opposite determinant. Indeed the root-first word successively gives
(5,2), (4M+1,2), (8M,4M+1), (16M+2,4M+1), with one A. This proves
the all-M claim directly; no pre-audit dependency is needed and no novelty
is claimed. The lost observation is the Euclid vector (or its triple),
**not the entire transport matrix**: the latter retains determinant ±1
modulo every odd M>1 and therefore distinguishes the signs. Adding
the bit yields the inherited Paley tournament observer; its arrows act
on tournaments and are not themselves their tournament arcs. Similarly
the F8 joint-degree sidecar and abc radical/height sidecars remain typed.

We do not derive a halting equivalence for any Turing machine from this
clock graph. The actual guarded machine configuration and positive native
duration are required by the inherited acceleration theorem. No claim
that classical Collatz itself is Turing-complete is used. The mathematical
finite path lifting theorem is a certificate interface, not a solution to
an infinite fixed-source problem.

## Reproduction

Run `python3 -I -B check.py` and `python3 -I -B -O check.py` here. Both must
match the stored normal/optimized receipts. The checker compares the graph
paths with independently generated integer traces, tests all refinement
ratios in a finite grid, checks affine/native clocks and both infinite-path
hostiles at finite depths, and composes a height-complete sieve with trace
classes. Exact universes and source hashes are in the accompanying files.
