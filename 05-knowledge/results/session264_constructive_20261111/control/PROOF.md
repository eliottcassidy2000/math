# Physical modular control and a hidden affine clock

**Status: PROVED from the inherited THM-5130 transfer law, with separately
labelled FINITE-EXACT extensions. Chan–Pak and integer target coverage remain
OPEN.** These are new controller deductions, not a re-audit of that law.

## 0. Inheritance and the question changed

The accepted input is THM-5130,
`01-canon/theorems/THM-5130-sign-flip-towers-gadget-transfer-law-and-seeds-of-both-signs.md`,
parts A–C, at file commit `7c14d6741eab74550c9337f6143751399d36053e`.
Session263's `transfer/PROOF.md` supplies a free F3/F5 state orbit, and its
`ONE_QUERY_EXTENSION.md` observes the complete state in one bounded-depth
count. Freeness and observability do not imply target controllability. That
is the corrected near miss here; the hostile is a free thin integer orbit.
More sharply, `session199/profinite-source-escape/PROOF.md` already proves
scalar odd-compatible modular density for EVERY modulus in another physical
architecture, and proves that its witnesses can escape every finite target.
No novelty is claimed for scalar modular density by itself. The new subject
is the initialized moving-mark grammar, exact unit offset, and independently
controllable full-state residues.

The least-used sidecar is E, the selected ending-path count. It will first
be discarded by an actual operation that no longer depends on it, then
retained to expose a new affine translation. The diverse old sources and
typed maps are in `SOURCES.md`.

Live concepts: physical moving marks / Euclidean shears / affine Jordan
clock / finite positive inverse words / integer-versus-modular target.
Anchor: source supply for Chan–Pak. Niche: complete modular controller.
Wildcard: a Pell-style hidden clock inside a rapidly expanding positive
counting system.

For a selected tournament (B,O), put r=H(B)-W_B(O), M=U_B(O),
h=H(B), q=h+r, and E=End_B(O). A unit source means r is exactly +1 or -1.
All vertices below are actual tournament vertices. Congruences never turn
an arbitrary formal tuple into an assumed literal initial source.

## 1. A two-coordinate controller valid at every modulus coprime to 30

There are two legal operations, valid on every literal selected tournament:

- D adds one sink, keeps the old selected set, and sends
  (r,q,M,E) to (r,q,2M,0).
- A takes the closing C(B,O;TT3), selects the SECOND and THIRD vertices of
  TT3, and forgets the old selection. It adds four vertices and sends
  r to -r, q to 8M, and M to 3q+20M. Its ending coordinate is
  E'=6M+2r-q, but no later A or D needs that coordinate.

These follow from THM-5130(A), Remark 2 and (B), with i=2, delta=1, m=0.
Thus forgetting E is an exact closed quotient for THIS grammar, not a claim
that the unrestricted source algebra has only two coordinates. On x=(q,M),

    D = [1 0; 0 2],          A = [0 8; 3 20].                 (1)

### Theorem 1 (every slope residue, with an exact offset sign)

For every positive integer n coprime to 30, every a modulo n, and each
r in {+1,-1}, there is a literal source of EXACT offset r and M=a modulo n.
One can additionally require q=1 modulo n. There is an explicit finite
constructor producing at most 8n²-2 vertices for n>1.

*Proof.* Work over R=Z/nZ. Both matrices in (1) are invertible. Enlarge
the finite matrix group by the exact sign: A flips r and D preserves it.
Every inverse in that finite group is a nonnegative power of its original
positive physical generator. For A use its order in C2 x GL2(R), so a
modular inverse is never incorrectly interpreted as restoring the old
vertices, decreasing an integer count, or having the wrong offset parity.

Set L=ADA^(-1). Direct multiplication gives

    L = [2 0; 5/2 1],
    D L D^(-1) L^(-1) = [1 0; 5/4 1].                       (2)

Both words preserve r exactly. Since 5/4 is a unit, powers supply every
lower shear E21(t). The matrix

    W=E21(-5/2) A = [0 8; 3 0]

satisfies W E21(t) W^(-1)=E12(8t/3), giving every upper shear, again with
offset preserved. These are positive physical words after finite-group
inverse replacement.

For r=-1 start with TT2 and both vertices selected: (q,M)=(0,3).
Apply E12(1/3), then E21(a-3), obtaining (1,a).
For r=+1 start with one selected vertex: (q,M)=(2,1).
Apply E12(-1), then E21(a-1), again obtaining (1,a).
This proves the assertion without assuming an arbitrary primitive vector
has already been realized.

For an effective size bound, breadth-first search the graph of residues
(r,q,M), using precisely the two positive edges A,D. There are at most
2n² vertices, so a reachable target has a word of length at most 2n²-1.
Each letter adds at most four vertices; the larger initial object has two.
This gives 2+4(2n²-1)=8n²-2. Rebuilding that word gives the actual tournament.
For n=1 either initial source suffices and no bound from the nontrivial
search is needed. ∎

**Corollary 2 (actual scalar modular coverage).** For every n coprime to30,
every odd residue b modulo 2n occurs as H(T), with T on at most 8n² vertices
when n>1. Indeed choose a=(b-r)/2 modulo n and close the resulting source
with one transitive tail vertex and its closer: H=r+2M. This is a literal
count, not a formal linear observation. Repeatedly adding ord_n(2) sinks
before closing gives infinitely many distinct counts in the same residue.

No useful bound uniform in n is hidden: the BFS bound is quadratic, not
logarithmic in a prescribed integer target. The two excluded odd primes of
this proof are structural to its specific shear extraction; they are not
claimed spectral obstructions.

## 2. The old full-state letters contain an exact affine clock

The physical F3,F5 operations append TT2 or TT3, close, and select its first
and last vertices. They add three or four vertices. For a in {3,5},

    r'=-r,
    M'=(a²+1)M+a q-E,
    q'=(2a-2)M,
    E'=2r+aM-q.                                             (3)

Put x=(M,q,E), y=r x (only for unit offsets). Because r²=1,

    y'=C_a y+c,
    C_a=-[a²+1 a -1; 2a-2 0 0; a -1 0],
    c=(0,0,-2).                                            (4)

This affine system is independent of the current offset sign; each letter
still flips that sign physically.

### Theorem 3 (exact integer clock)

For every literal unit source, define z=r(M-4q+E). Then

    z(F5(s))=z(s)-2,
    z(F5²(s))=z(s)-4.                                      (5)

Consequently a repeated fixed-offset physical operation moves this integer
observable by an exact constant, although the positive counts M,h,E grow.

*Proof.* With w=(1,-4,1), w C5=w and w c=-2. Equivalently,
M'-4q'+E'=-(M-4q+E)+2r; multiply by r'=-r. ∎

The observable is a signed combination, not itself a path count or a source
slope. Equation (5) alone cannot hit a prescribed H(T). Its constructive
value is the exact resonance used next. Starting from the marked TT2,
z=-4, so the F5^k clock is exactly -4-2k at every depth.

## 3. Full three-coordinate residue control over almost every prime

Let

    v=(1,-8,-13),             w=(1,-4,1).

Then C5 v=v, w C5=w, and wv=20. The characteristic polynomial is

    det(tI-C5)=(t-1)(t²+27t-8),                             (6)

and the quadratic discriminant is 761.

### Theorem 4 (full modular target controller)

Let p be any prime other than 2,5,761. Given ANY actual unit source,
ANY desired final sign r*, and ANY three residues (M*,q*,E*) modulo p,
a positive F3/F5 word produces an actual unit source with exactly that
sign and those residues. A shortest residue-state BFS word has at most
2p³-1 letters, adding at most 8p³-4 vertices.

*Proof.* Put d=p²-1, which is even and nonzero modulo p. The roots of the
quadratic in (6) are distinct in F_(p²), nonzero, and different from1.
The exclusions follow respectively from p!=761, p!=2, and p!=2,5 because
the quadratic at1 is20. Thus C5^d=I, and its geometric sum is dP, where
P=v w/20 is projection onto the 1-eigenspace. The d-fold physical operation
F5^d therefore fixes the offset sign and acts on y by

    y -> y+dP c = y-(d/10)v = y+(1/10)v   (mod p).          (7)

It is a PURE TRANSLATION of the complete three-coordinate state, not just
a scalar clock. Its p powers supply every multiple of v.

Conjugating a translation by any affine map replaces its vector by the
linear part of that map. Inverses are positive words in the finite group
of pairs (exact sign, affine state), as in Theorem1. Hence the available
fixed-sign translation vectors include

    v, C3 v, C3² v, C5 C3 v.

Writing L_a=-C_a, the corresponding unsigned columns are

    v          =( 1,-8,-13),
    L3 v       =(-1, 4, 11),
    L3² v      =(-9,-4, -7),
    L5 L3 v    =(-17,-8, -9).

The determinants of the first, second, third columns and of the first,
second, fourth columns are344 and632. Since gcd(344,632)=8, for EVERY odd
prime at least one triple spans the whole state space. Arbitrary translations
of F_p³ are consequently available with exact offset unchanged. An initial
F3 if needed supplies the desired offset sign, after which these translations
set all three coordinates. Multiplication by r* converts target x to target y.

The graph of the two positive generators has 2p³ residue states; their
matrices are invertible for odd p because det(F3)=-4 and det(F5)=-8.
The established transitivity therefore gives the stated BFS bound and a
fully constructive word. Every intermediate state comes from its literal
positive word, rather than a realizability test on formal tuples. ∎

**Finite exact extension at p=5.** `controls.py` gives a rooted BFS tree
containing ALL250 states (r,M,q,E) and checks every stored transition.
Both generators permute that finite set. Therefore it is one orbit, and
Theorem4's existence and size statement holds at p=5 too. This finite case
is explicitly computation-supported; it is not obtained by dividing20
modulo5. Combining the algebraic proof and this finite certificate gives
full controllability at every odd prime except possibly761.

The prime761 is an unresolved exception of this proof, not a counterexample.
Its repeated quadratic root prevents the pure-translation identity (7) with
d=p²-1. The failure is recorded exactly in the controls. No coverage at
prime powers or products of distinct primes follows just from Theorem4.
Theorem1 separately handles ALL moduli coprime to30, including those
containing761 and its powers, but controls only the q,M quotient there.

## 4. What is physically gained and what remains missing

1. Theorem1 is an initialized constructive slope-residue cover at either
   exact unit offset over infinitely many composite moduli. It is stronger
   than proving that a free source family has many words.
2. Theorem4 sets the complete source interface modulo almost every prime;
   it can prescribe endpoint residue E independently, for subsequent
   gadget operations which need that sidecar.
3. The exact clock (5) identifies a small integer displacement hidden in
   large expanding source counts. Its affine resonance, not empirical
   random-walk mixing, proves (7).
4. Positive finite-group inverse replacement is legal tournament growth.
   It never removes vertices, returns original selected vertices, realizes
   negative intermediate counts, or gives an inverse operation over Z.
5. Modular target cover does NOT provide a smaller integer predecessor.
   The graph's residue representative can overshoot a requested integer by
   an arbitrarily large amount. Even all odd moduli would not close that
   gap. Chan–Pak is still open.

A focused next question is an integer-height-controlled lifting theorem:
can one prescribe progressively higher powers of a fixed odd prime while
bounding the least physical slope at that residue tightly enough to retain
a target interval? The present BFS bounds do not imply such a statement.
The separate finite all-state experiments for odd composite n<=31 are
motivation only, not a prime-power or CRT lifting theorem.

## 5. Reproduction

Run `python3 -B 05-knowledge/results/session264_constructive_20261111/control/controls.py`.
It writes `RESULTS.json` and `MOD5_TREE.json`. Controls concern only these new
operations, shears, clock, modular orbits and constructor. Literal graphs
are independently recounted for tiny A,D/F3/F5 words; larger residue-state
orbits are exact modular computations, not graph censuses. Source provenance
and maps are recorded in `SOURCES.md`.
