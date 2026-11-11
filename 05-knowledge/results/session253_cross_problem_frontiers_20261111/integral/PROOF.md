# The exact integral and affine lattices of the coloured-core compiler

**Session253. PROVED IN NOTE from inherited physical compiler identities.**
This is an exact statement about integer spans of actual response vectors,
not realization of arbitrary vectors or target coverage. No canon promotion
or literature-priority claim.

## 1. Inheritance and question

Use the seven nonnegative coordinates

    v=(h,d,L,i,o,s,e),
    Phi_v(A,B)=h+d(g(A)+g(B)-2)+LAB+iA+oB+s(A mod2)+e(B mod2).

Here d=H(C), h=H(C+x), L=d+K_C(1); the remaining coordinates retain
endpoint/cut data. The precise definitions and physical observer are in
[the incoming positive compiler](../../chanpak_five_core_bridges_20261110/algebra/POSITIVE_COMPILER.md).
Session252 gave a rationally minimal seven-dimensional carrier; the
integral-observer lane asks what torsion this rational statement loses.

Every nonempty coloured tournament satisfies

    h≡d≡1 (mod2), L≡0 (mod2), i+s≡o+e≡0 (mod2).        (1)

These are inherited endpoint/parity identities. No old theorem is re-audited.
The new question is whether more universal LINEAR congruences survive.

## 2. The first three cores already span the entire universal linear lattice

Let Lambda be the Z-span of all actual coloured-tournament response vectors.
Then exactly

    Lambda={v∈Z^7: h≡d, L≡0, i+s≡0, o+e≡0 (mod2)}.   (2)

Its index in Z^7 is16 and its Smith invariant factors are

    (1,1,1,2,2,2,2).

Even the three cores K1,C3,S4 already generate Lambda. Neither H9 core adds
an integral linear direction.

**Proof.** Equation(1) gives containment in the lattice on the right of(2),
which has index16. The following eight ACTUAL vectors use only K1,C3,S4;
the cut is a bit mask on the literal rows (0), (2,4,1), (2,4,9,3):

    core cut    h d  L  i o s e
    K1     0    1 1  2  1 0 1 0
    K1     1    1 1  2  0 1 0 1
    C3     0    3 3 12  9 0 3 0
    C3     1    5 3 12  6 3 2 1
    C3     3    5 3 12  3 6 1 2
    S4     0    5 5 32 27 0 5 0
    S4     1    9 5 32 18 6 4 2
    S4     2    9 5 32 20 5 4 1.

Delete each column in turn from the 7×8 matrix with these columns. Its
maximal minors are

    (-288,0,176,0,0,48,0,0),

whose gcd is16. The lattice index is the gcd of maximal minors (Smith
normal form), so containment is equality. The displayed parity conditions
immediately identify the quotient as (Z/2)^4. QED.

**Consequence.** At every odd prime the integer span of the small-parent
bank is already full. No additional universal homogeneous linear congruence
in these seven coefficients can distinguish the H9 repair language. This
says nothing about the actual positive/nonlinear image of the constructor.

## 3. Affine spans do see a factor five, but not the two H9 branches

An affine integer combination has coefficients summing to1. Actual vectors
lie in the coset defined by(1). Their difference lattice is

    Delta={w∈Z^7: h,d,L even; i+s,o+e even},           (3)

of index32, with Smith factors (1,1,2,2,2,2,2).

For the first THREE cores only, replace "L even" in(3) by "L divisible by10".
The resulting difference lattice has index160. Equivalently their exact
affine lattice adds the condition

    L≡2 (mod10).                                    (4)

Adding EITHER H9 core removes this factor5 and gives the full universal
coset(1). Adding both does not further enlarge that affine lattice.

**Proof.** The first three cores have L=2,12,32, so all differences have
L divisible by10; their other conditions follow from(1). Subtract the
first displayed vector from the other seven. The determinant of this
7×7 difference matrix is−160, equal to the index of the proposed lattice.
For either H9 core the zero cut gives the same vector

    (9,9,84,75,0,9,0).

Append its difference from the first vector. The eight maximal minors of
the resulting difference matrix, in increasing seven-subset order, are

    (-160,0,0,-672,0,0,-704,0).

Their gcd is32. Thus it spans(3), which already contains every possible
tournament difference. QED.

This is a small but real integral gain supplied by H9. Its limitation is
just as useful: both H9 cores erase exactly the same affine5-torsion, yet
the inherited exact holes21901 and86641 require different branches in the
five-core grammar. Their difference cannot be explained by this lattice
invariant. We must retain the positive feasible states and bounded target
fibres, not merely the affine hull.

## 4. Seven physical measurements have a one-extra-bit congruence

Order the inherited physical measurements as

    p=(p00,p10,p30,p20,p01,p03,p11), p_ab=Phi(a,b).

The evaluation matrix E in coordinates v is

    [1 0 0 0 0 0 0]
    [1 0 0 1 0 1 0]
    [1 0 0 3 0 1 0]
    [1 2 0 2 0 0 0]
    [1 0 0 0 1 0 1]
    [1 0 0 0 3 0 1]
    [1 0 1 1 1 1 1].

It has determinant−8. The integer image E Lambda is EXACTLY the lattice

    all seven p coordinates have the same parity;
    p20−p30+p10−3p00≡0 (mod4).                       (5)

Its index is128 and its Smith factors are

    (1,2,2,2,2,2,4).

For actual AFFINE combinations, all p coordinates are odd and the same
mod4 condition holds. The difference observer has factors (2,2,2,2,2,2,4).

**Proof.** The inverse formulas from the physical observer are

    h=p00, i=(p30−p10)/2, o=(p03−p01)/2,
    d=(p20−p00−p30+p10)/2,
    s=p10−h−i, e=p01−h−o,
    L=p11−p10−p01+p00.

Common parity makes each coordinate integral and imposes L,i+s,o+e even.
The mod4 condition is precisely d≡h mod2. Conversely these conditions
recover a member of Lambda, proving(5). To read the Smith factors, change
target coordinates unimodularly to p00 and its six coordinate differences.
Write those differences as twice integers; the mod4 condition is one
primitive parity condition among them. Five factors are2 and one is4.
QED.

## 5. What transfers, and what remains open

- **Source→target:** actual tournament core/cut → seven coefficients → seven
  small physical measurements. The map preserves every two-flank response.
- **New preserved integral structure:** exact linear and affine lattices,
  including the one-extra-bit compatibility and the small-bank5-torsion.
- **Destroyed information:** realizability, positivity cone, literal core
  type, absolute target size, and the available address fibres.
- **Necessary sidecar:** actual row masks/cut and a legal nonnegative
  address pair. Integer linear combinations are NOT tournament operations.
- **Cheapest next test:** take a proposed target-preserving chart rewrite;
  check(5) and its action on the affine lattice first, then test actual
  positive states and target cost. A lattice-preserving rewrite alone does
  not supply a physical rewrite.

The strongest practical conclusion is negative but focused: prime-local
linear algebra cannot distinguish the two H9 repairs. This tells us to
spend effort on integral-height/positive-fibre geometry rather than another
linear modular obstruction search. It complements, rather than replaces,
the Laurent/jet lane where full Smith data DO distinguish hidden directions.

## Reproduction

Run `python3 integral/probe.py` from the session directory. It recomputes
only the nine small coefficient vectors, exact minors, observer determinant,
and congruence roundtrips. This is a new lattice calculation, not a replay
of inherited coverage or an audit of the inherited compiler.
