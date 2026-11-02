# Strong one-apex local universality and the profinite blindness boundary

**PROVED, subject to independent audit.** The full integer spectrum remains
OPEN. No literature-priority claim is made.

## Scope and inherited input

For a finite tournament T, H(T) counts directed Hamiltonian vertex words.
Let S1str be the H-values of nontrivial strongly connected tournaments that
become transitive after deletion of one distinguished vertex. Let Sstr and S
be the strong and unrestricted spectra. We prove

    closure(S1str) = closure(Sstr) = closure(S) = 1 + 2 Zhat.

Equivalently, for every m>=1, the residues attained by each spectrum are all
residues when m is odd, and precisely the odd residues when m is even.
Every compatible residue contains arbitrarily large values of S1str. All
CRT requirements are simultaneous, witnessed by one actual tournament.

Inherited inputs, checked at origin/main d11ba3658601:

* `05-knowledge/results/session133_fractran_spectrum_refinement_20261101/path_counter/THEORY.md`,
  Sections 1 and 8: actual binary one-apex encoding, transfers, strongness,
  and the alternating-word formula. These are not new results here.
* `05-knowledge/results/session124_digital_constructor_inverse_20261031/PROOF.md`,
  Section 3: exact one-apex inverse and its all-length H=13 exclusion.
* `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`:
  the unrestricted omissions 7 and 21, not reproved or contradicted here.

A word b=b1...bs encodes the transitive chain v1->...->vs plus apex x,
with x->vi iff bi=1. Starting in 1 and ending in 0 supplies the spanning
cycle x,v1,...,vs,x, hence strongness. Set c=(H-1)/2. On appending a bit,
the inherited exact state is

    0: (c,t,l) -> (c+t+l, 2t+l, 0),
    1: (c,t,l) -> (c,     2t+l, 1),

starting at (0,0,0). At a block boundary with l=0 put d=c-t. Two actual
append operations are

    A=0:  (d,t) -> (d,   2t),
    B=10: (d,t) -> (d-t, 4t+1).                         (1)

Chronological words always consist of nonnegative repetitions of these
blocks; inverse symbols below occur only in a finite modular group and are
explicitly replaced by positive powers. No fabricated state is used.

## 1. Odd-modulus controllability

Fix odd n>1. Both operations (1) are permutations of (Z/nZ)^2.
Use ordinary function-composition order, rightmost first, and set

    C=A^(-2) B:          (d,t) -> (d-t,   t+1/4),
    D=A C A^(-1):        (d,t) -> (d-t/2, t+1/2),
    J=D C D^(-1) C^(-1): (d,t) -> (d+3/8,t),
    E=D C:              (d,t) -> (d-3t/2-1/8,t+3/4).    (2)

The fractions mean inverses of powers of two modulo n. Direct substitution
proves every formula. Thus J translates c=d+t by 3/8. E sends

    c -> c - 3t/2 + 5/8,

which increments c by 1 modulo 3 whenever 3 divides n.

Given any starting state and any target z modulo n, put g=gcd(3,n).
Choose j=0 if g=1; if g=3, choose j in {0,1,2} so that after E^j the
c-coordinate agrees with z modulo 3. The increment 3/8 generates exactly
the multiples of g modulo n. Consequently some 0<=k<n/g makes J^k E^j
land at c=z modulo n. This works from EVERY starting state; its use below
is not restricted to the zero state.

Positive inverse replacement is quantitative. Let e=ord_n(2). Then
A^e=I and B^(ne)=I modulo n. To prove the latter, the homogeneous matrix
of B is upper triangular with diagonal (1,4,1). Hence B^e=I+N, where N is
strictly upper triangular and N^3=0. Since n is odd,

    (I+N)^n=I+nN+binom(n,2)N^2=I modulo n.

Thus A^-1 can be replaced by A^(e-1), and B^-1 by B^(ne-1).
For a sharper word bound, C and D themselves have homogeneous matrices
I+N with N^3=0; the same odd-binomial argument gives C^n=D^n=I. Thus use
C^-1=C^(n-1) and D^-1=D^(n-1) directly. No factorial bound on a
permutation group is needed.

## 2. Simultaneous dyadic/odd control

It suffices to attain every c-residue modulo every M. Write
M=2^a n with n odd. First choose

    r >= max(1,ceil(a/2)),
    3r+1 = 9z (mod 2^a),                               (3)

where z is the desired c-residue modulo M. Such r can be chosen in any
interval of 2^a consecutive integers, since 3 is a dyadic unit.
Start with the actual word B^r=(10)^r. Its inherited alternating formula
and the recurrence (1) give

    t_r=(4^r-1)/3,
    c_r=(2*4^r+3r-2)/9,
    d_r=(3r+1-4^r)/9.                                 (4)

All three fractions are integers. By (3) and 2r>=a, d_r=z modulo 2^a.
This starts in 1 and ends in 0, so its tournament is nontrivial and strong.

If n>1 let e=ord_n(2); if n=1 set e=1. Choose L to be a positive multiple
of e with L>=a, for example the least such multiple. Define

    P=A^L,          X=P B P.                            (5)

Then P is the identity modulo n, and resets t to zero modulo 2^a while
leaving d unchanged. After the initial B^r, apply P. On the invariant
subset t=0 modulo 2^a, both A and X fix the whole (d,t) state modulo 2^a:
B sends (d,0) to (d,1), and the final P resets t. Modulo n, A and X act
exactly as A and B. Therefore perform the odd-control construction (2)
with X in place of B; every positive A/X word preserves the already
correct dyadic c=d, while attaining the desired odd c. Inverses are
replaced by A^(e-1), and by C^(n-1),D^(n-1) after constructing their
positive A/X versions. If n=1 no odd step is
needed. This establishes c=z modulo M, simultaneously, by a genuine word.

Every appended block ends in 0; no operation deletes the initial 1.
Thus the witness remains strong. This is a direct CRT synchronizer, not an
inference from separate prime-power witnesses.

For an original H-residue h modulo m, solve 2z=h-1 modulo m. This is
possible exactly when gcd(2,m) divides h-1, and the c theorem produces it.
Rédei parity excludes the other residues. The claimed closure follows.

## 3. Unbounded fibers and explicit size bound

After the construction, t=0 modulo 2^a. Appending any number of P blocks
preserves the entire state modulo M. Over the integers t>0 because the
word began with B. Hence appending P^q changes c by (2^(Lq)-1)t and gives
arbitrarily large, distinct counts in the same residue. Strongness remains.

For a numeric bound let e=ord_n(2)<=n for n>1, e=1 for n=1,
and w=2L+2, the bit length of X. Choose
r<max(1,ceil(a/2))+2^a and L<=a+e. Positive versions of C,D,E,J above
can be selected with lengths

    |C| <= w+e,
    |D| <= w+2e,
    |E| <= 2w+3e,
    |J| <= n(2w+3e).

Since j<=2 and k<n, the constructed bit word has length at most

    2r+L+(2+n(n-1))(2w+3e).                            (6)

The graph has one additional apex vertex. This is polynomial in the
NUMERIC modulus (e.g. O(2^a+n^2(a+n))), not in its bit length. It is not
a useful minimal-order bound. The compact constructor stores repetitions;
its full expansion is an explicit adjacency recipe of the stated size.

## 4. What local observers provably cannot decide

A finite collection of congruence predicates combines into a single modulus
by lcm. Every parity-compatible joint tuple actually occurs in a strong
one-apex tournament. Thus no such observer can certify nonattainment of an
odd target; in particular 7,21 cannot be excluded by any finite congruence
filter on H alone. This does NOT rule out congruence arguments using
additional graph structure or noncongruential magnitude information.

The hostile is stronger than merely the known unrestricted omissions:
S1str is already profinitely dense among odd integers, yet 13 is absent
from the entire one-apex family and is attained by an unrestricted order-five
tournament. The exact inherited inverse gives, for reduced address Q=2u,

    H=f(u)+u+(u mod2),   (H-1)/2 <= u <= H-1,
    f(0)=f(1)=1,
    f(2v+1)=f(v),
    f(2v)=f(v)+v+(v mod2).

At H=13 only u=6,...,12 need checking; none works. Thus even full local
universality of an actual strongly connected constructor does NOT imply its
integer completeness. Likewise the theorem does not settle S=odds\{7,21},
prime-lane realization, natural density, short witnesses, or uniformity in
any fixed graph order. Different moduli generally use different graphs.

## 5. Checks and provenance

`check.py` is standalone standard-library code, with no producer imports.
It builds compressed positive block words, evaluates their affine actions,
checks the complete residue universe for 1<=M<=160, and tests selected larger
mixed moduli. Independent direct Hamiltonian subset-DP checks all binary
words through length 9, including the strongness criterion and H formula.
Additional checks verify the modular commutator and odd-order bounds, compare
literal expansion with compressed execution for short constructor words,
prove the bounded H13 hostile through the exact inverse, and reject deliberate
commutator-sign and CRT-reset mutations. Normal and -O runs must agree.

These finite checks validate implementation and hostile boundaries; the
universal theorem is the argument above, not extrapolation from them.
