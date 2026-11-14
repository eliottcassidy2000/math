# Every finite congruence, yet an absolute hole in the same tournament family

**Theorem.** Fix any positive modulus M and any odd target N.
There is a strong tournament consisting of a transitive spine
and one apex whose Hamiltonian path count is N modulo M.
Nevertheless the same family cannot attain 13 exactly (quite
apart from the universal exclusions 7 and 21). The local
statement holds simultaneously across every factor of M, not
just separately for prime powers. In fact the address P and
count H can be prescribed jointly modulo M, subject only to
P even and H odd wherever parity is defined.

This provides a precise model of the local-to-global failure
motivating the Jacobian analogy. It does not prove Chan–Pak.
It shows why congruence solvability of a fixed construction
cannot by itself contradict absence at an absolute target.

## Two physical appends

Paint a transitive word by writing 1 when the apex beats its
spine vertex and 0 when that vertex beats the apex. Let P be
the ordinary binary address of the word, and H the full
Hamiltonian count. For a word ending in 0 the ordinal cut
transfer of [the selective-cut session](
../chanpak_selective_target_20261114/CUT_STATE.md) gives
these two physical operations, each retaining a terminal zero:

    A: append 0,       (P,H) → (2P,H+P);
    B_j: append 1^j0, (P,H) → (2^(j+1)(P+1)−2, H+2^j(P+1)).

For the second formula the ones append neutrally in count,
then the zero reads the new address plus its final bias 1.
Both formulas also follow directly by enumerating odd cycles
through the apex. Start with word 10, hence P=2, H=3; every
continuation used below remains strong because the first
letter is 1 and the last is 0, giving an explicit spanning
cycle through the entire ordinal spine.

## First set the binary part, then freeze it

Write the modulus as 2^a m with m odd, using a≥1
(by strengthening the matching modulus if necessary).
The terminal-zero state has P even. When j≥1, B_j adds
exactly 2^j times an odd integer to H, toggling bit j
without affecting lower bits. Starting at 10, choose
whether to append B_j for j=1,...,a−1 to set the target's
binary residue. Then append B_a once. It does not change
that residue and makes P≡−2 modulo 2^a.

From now on use high B_j (j≥a). They preserve the binary
residue of H and the anchor P≡−2. It remains to choose
the odd residue, without undoing the part already set.

## The odd-modulus transformations and the exceptional prime

Over Z/m with m>1, put X=P+2 and J=H−(P+1)/2. By
selecting j≥a in an appropriate period of 2, two high
blocks realize respectively 2^j≡1/2 or 1 modulo m.
Call these U and V. In these coordinates a high block
with parameter b=2^j reads

    X → 2b(X−1),    J → J+X/2.

Thus U subtracts 1 from X and adds X/2 to J. The
composition R=V U^(-1) doubles X and leaves J unchanged.
Here words in transformations apply from right to left.
Every inverse modulo m used in this proof is a nonnegative
power, since these transformations permute a finite state
space. As actual positive append words, they preserve
the already-set binary anchor and output residue.

Let S=R U R^(-1) U^(-2). Then S fixes X and adds
−3X/4−1 to J. Comparing SU against US shows that their
difference is a pure translation of J by 3/4. Thus the
positive words realize this translation with X fixed.
If 3 does not divide m, this alone reaches every H residue
at the initial X. If 3 divides m, it reaches exactly the
residue class already present modulo 3; the next small
change of word supplies the missing coordinate.

Choose L≥a with 2^L≡2 modulo m and let D be A^L followed
by U. Starting in the binary anchor, the zeros temporarily
set P≡0 modulo 2^a and increase H by 2 modulo 2^a; U
restores P≡−2 without further change to H modulo 2^a.
Thus D^q with q=2^max(1,a−1) preserves both binary
coordinates. On the odd part D acts by

    X → 2X−3,    J → J+3X/2−2.

Modulo 3, this is (X,J)→(−X,J+1). Consequently D^q
fixes X modulo 3 and shifts J by q, a nonzero residue
modulo 3. To restore the *full* original X modulo m,
append U^k with k≡X_after−X_before modulo m. This k
is divisible by 3 and U³ is identity modulo 3, so the
restoration preserves the new J modulo 3. Zero, one or
two such corrections can set J modulo 3. The pure 3/4
translation then sets J exactly modulo m with X fixed.
Since H differs from J by a function of X, this sets the
desired H. Agreement modulo 2^a and m yields agreement
modulo the prescribed M. This entire prescription consists
of physical append words. For m=1 the binary part alone
suffices.

## Strengthening: no joint congruence obstruction in the terminal-zero state

At the binary anchor, modulo the odd part one may first
choose X at will by powers of U. The subsequent correction
preserves that chosen X. Thus P can be prescribed on the
odd part while leaving P≡−2 on the binary part. To prescribe
a different even binary P as well, append a final fixed
a-bit word carrying that P residue, ending in zero. Its
effect on the terminal-zero state is triangular:

    P_final = 2^a P_initial + address_of_word,
    H_final = H_initial + c P_initial + d

for fixed integers c,d obtained by the append recurrence.
On the odd part this is invertible; on the binary part it
sets the address as requested, and its count shift is a
known constant because P_initial≡−2 there. (That shift is
even.) Choose the initial H residue and the initial odd-part
P residue by these formulas, and then run the compiler.
This proves joint saturation without finding an arbitrary
*integral* P and H simultaneously.

## Where the witnesses go as precision increases

Strip neutral leading zeros and trailing ones from any
nontransitive painted word. The remaining word starts 1
and ends 0. At active spine length n≥3 its first-to-last
odd cycles alone contribute 2^(n−2) to H: after the two
forced endpoint vertices, choose an even subset of the
n−2 internal vertices, with weight two per cycle. Thus a
witness for H=13 has active n≤5. The direct recurrence
at lengths 2 through 5 lists no 13. The replay checks this
finite exhaustion; it also checks 7 and 21 with their
own cutoffs.

For a target 13 modulo M with M>13 the witnesses must
instead have actual H≥M+13. The sum of all potential
oriented pair-cut weights on an active spine is at most
2^n (for i<j the contribution to H−1 is
2^(j−i−1)+[j=i+1]). Thus active n must tend to infinity
at least logarithmically as the modulus grows. Every
finite observation is compatible, but the lifts escape
in actual height. This is not the absence of a real
tournament witness at fixed bounded count: it reflects
unbounded witnesses to the *relaxed congruence problem*.

The missing value 13 is an obstruction only in this
grammar, not a counterexample to Chan–Pak: the five-vertex
tournament with rows (18,8,3,5,14) has H=13. In the full
problem, a scalar exact target must force a change of
source, or allow a descent that bounds the height of
the lifts. Arithmetic symmetry alone can be complete
on all finite quotient rings and still leave the
positive-integer image incomplete.
