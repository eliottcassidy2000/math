# Arbitrary expanding banks followed by arbitrary reset banks

**Status, 2026-10-19: PROVED, conditional cylinder theorem.** Arbitrary finite
words from a bounded-center expanding bank, followed by arbitrary words from
a bounded-center contraction bank, have a uniform excluded-source-lift head
once the total multiplier contracts. Actual oriented adjacent returns admit
the same bound simultaneously at both endpoints. A concrete two-center bank
has first-descent-free lower expanding segments. No source-entry or global
Collatz conclusion is proved.

This extends the one-expander/one-reset bounded-head theorem in
[the odd-adjacent packet](collatz_odd_adjacent_return_20261019.md), section 3.
The new quantifier is arbitrary word length and order in **both** banks,
not merely repetitions of one chosen word. The elementary affine envelope
argument is not claimed as a novel general affine-systems result. The older
[finite translation-bank certificate](collatz_finite_bank_certificate_20261007.md)
concerns different state translations and does not supply this source-lift
head theorem. Searches of the Collatz result notes found the single-word
bounded-head version, not this bank extension.

## 1. The two composition envelopes

Let the expanding bank consist of lawful finite Terras parity-word branches

    U_i(x)=a_i x+b_i,  a_i>1, b_i>0,
    c_i=b_i/(a_i-1) <= C.

Let the reset bank consist of lawful finite branches

    W_j(x)=r_j x+d_j,  0<r_j<1, d_j>=0,
    p_j=d_j/(1-r_j) <= P.

Thus C>0 and P>=0. Finiteness ensures that such constants exist, though the
proof needs only the displayed uniform bounds. Here and below a word is
applied chronologically from left to right, so U_1 U_2 denotes U_2 composed
with U_1. Each block is a guarded branch, not an unrestricted transition.

For any nonempty expanding word E(x)=Ax+B,

    A>1, 0<B<=C(A-1).                            (1)

For any nonempty reset word R(x)=rx+d,

    0<r<1, 0<=d<=P(1-r).                         (2)

Indeed appending U_i changes B to a_i B+b_i, at most
C[a_i(A-1)+(a_i-1)]=C(a_i A-1). Appending W_j changes d to
r_j d+d_j, at most P[r_j(1-r)+(1-r_j)]=P(1-r_j r).
These inductive identities are also telescoping convex-center formulas:
B/(A-1) is a weighted average of the c_i, and d/(1-r) is a weighted
average of the p_j. Thus no ordering or commutation assumption is needed.

If rA<1, the final map is R(E(x))=rAx+e, where

    0<e=rB+d
       <= C r(A-1)+P(1-r)
       = C rA+P-(C+P)r
       < C+P.                                   (3)

This holds at the first reset boundary whose multiplier contracts and at
every subsequent reset boundary. If an infinite reset schedule is selected
from a finite bank, such a boundary exists since max_j r_j<1. Before it,
the full-boundary multiplier is at least one and the intercept is positive,
so no positive source descends at a full reset boundary. This says nothing
about intermediate steps inside arbitrary reset blocks.

## 2. Exact integer cylinders and the uniform excluded head

Concatenate the actual parity words of E and R. Let its exact affine form be

    T^N(x)=(Sx+H)/M,  S=3^h, M=2^N,
    Delta=M-S >=1,
    rho=-H S^(-1) modulo M.                      (4)

Here H is the integer carry, and H/M=e. The complete word contains an odd
step because E expands, so 0<rho<M. Its lawful positive source cylinder is
exactly x=rho+tM, t>=0. Full-word integrality gives every prefix guard:
reducing its carry modulo 2^k leaves the prefix carry multiplied by the
invertible remaining power of three. Induction then gives the required
parity at each step. This establishes actual trajectory admissibility;
a coarse block graph without these congruences would not suffice.

The exact descent criterion and cutoff are

    Delta x>H,
    t_min=max(0, floor((H-Delta rho)/(Delta M))+1). (5)

Equations (3) and Delta>=1 give H/Delta<(C+P)M. Hence

    t_min <= ceil(C+P).                          (6)

At most ceil(C+P) initial lifts are excluded, independently of expanding
word length, its order, and the reset choices. This bound is not a claim
that these initial lifts all fail: (5) gives the exact answer. The ceiling
wording is sharp for the general theorem: U=(3x+1)/2, W=x/2 yield the word
10 on x=1 modulo 4, for which x=1 is equality and every later lift descends.
Here C+P=1 and t_min=1.

## 3. Paired adjacent-return version, including residue indexing

Suppose every lower block has an independently verified upper parity word
of the same length, with the exact adjacent relation

    U_i^+(x+1)=U_i(x)+1,
    W_j^+(x+1)=W_j(x)+1

on the lower block's lawful cylinder. The upper branch maps are

    U_i^+(z)=a_i z+b_i+1-a_i,
    W_j^+(z)=r_j z+d_j+1-r_j.

Their expanding center sizes and contraction fixed points are respectively
c_i-1 and p_j+1. Since these are genuine Terras branches, upper intercepts
are nonnegative; in fact expanding upper branches have positive carry,
so c_i>1. Consequently one may use C^+=C-1 and P^+=P+1, with

    C^++P^+=C+P.                                (7)

The paired guard propagates through every concatenated block. The upper
full-word carry is H+M-S=H+Delta, and its residue is (rho+1) modulo M.
There is **no lift-index wrap** in this situation: if rho=M-1, the upper
word would have residue 0, which forces every parity bit to be even and
its numerator to be 1. It instead has S>1, the same odd count as the
expanding lower word followed by resets. Therefore 0<rho<M-1, and the
upper positive lifts are exactly rho+1+tM with the same t>=0.

Both endpoints have the identical descent test:

    Delta(x+1)>H+Delta  iff  Delta x>H.

Thus (5) is the same exact cutoff for the pair, and (6) is a common bound.
The preceding reasoning also handles the general least-positive-residue
convention: zero residues would be represented by M, but they do not occur
in this paired expanding-then-contracting setting.

## 4. A concrete nontrivial two-center expanding bank

Use the odd-source adjacent return K from the preceding packet:

    K_lower = 110110110110110101011011101010,
    K_upper = 001111111111111101010110000100,
    K(x)=(1162261467x+4912087879)/1073741824,
    x=815464443 modulo 1073741824,
    c_K=4912087879/88519643.

Use the adjacent contraction F, with words 10 and 01 and map (3x+1)/4.
Define V to be the genuine composite adjacent return K^4 F, with the two
words K_lower^4 10 and K_upper^4 01, each of length 122. Independent
parity-word coefficients give

    V(x)=(5474401089420219382077155933569751763 x
           +83829391903560839836419813590128166431)
          /5316911983139663491615228241121378304,

    x=1918171157796698145963164157509042171
          modulo 5316911983139663491615228241121378304.

The numerator exceeds the denominator by
157489106280555890461927692448373459. Hence V expands, and

    c_V=83829391903560839836419813590128166431
             /157489106280555890461927692448373459,
    c_V>c_K.

V is not a power of K: its center differs. The block deliberately includes
an internal contraction, but it is a macro-expander. This is fully allowed
by the theorem; it does not assume each individual Terras step expands.
Every nonempty lower prefix of V has multiplier greater than one: the K
prefixes already do, the first step of its final F multiplies by 3/2,
and its final multiplier 3^77/2^122 exceeds one. Thus any lawful expanding
word from {K,V} has no descending lower prefix. Its joint adjacent cylinder
is checked, rather than assuming arbitrary numerical sources can follow it.

Take the reset bank {F,H}, where

    F_lower=10, F_upper=01, p_F=1,
    H_lower=110101000, H_upper=001101001,
    H(x)=(81x+101)/512, x=43 modulo 512, p_H=101/431.

For every arbitrary nonempty word over {K,V} and every reset word over
{F,H} with contracting total multiplier, the paired cutoff satisfies

    t_min <= ceil(c_V+1)=534.                    (8)

This is a uniform worst-case certificate, not an optimized constant. On
all 396 sampled schedules the actual cutoff is zero; this finite fact is
not substituted for (8) or extrapolated to untested schedules.

## 5. One extra arbitrary reset removes the head in this example

There is a useful full-cylinder strengthening for this concrete bank.
Every expansion word begins with K, so every lawful positive source is
at least r_0=815464443. This satisfies

    r_0>3c_V+4.

At the first reset boundary with rA<1, append either F or H. Both have
slope at most 3/4, so the new total multiplier s is strictly below 3/4.
Apply (3) to the enlarged reset word: its intercept is at most C s+P,
therefore its final endpoint minus x is strictly less than

    -x/4+3C/4+P <0

for every x>=r_0. Hence every positive source of the **enlarged complete
cylinder** descends at both endpoints, with no excluded lifts. This does
not say that the extra reset is lawful at every source of the shorter
cylinder. It refines that cylinder by appending the actual reset guard.

## 6. Why the factorization and coverage limitations matter

The proof permits any order within the expanding bank and within the reset
bank, but all expanding macroblocks come before the final reset bank word.
There is no analogous normalized-intercept bound under unrestricted
alternation. For example U(x)=(3x+1)/2 has c=1 and W(x)=x/2 has p=0.
The lawful parity word 0^h 1^h, in the wrong order W^h U^h, has

    slope=(3/4)^h<1,
    normalized intercept=(3/2)^h-1,

which is unbounded and already exceeds C+P=1 at h=4. This explicitly
refutes extending inequality (3) without the factorization hypothesis.
It is not asserted to refute every conceivable excluded-lift bound.

Even for the correct factorization, none of these results forces an
arbitrary source into an expanding cylinder, forces an arbitrary expanding
exit into any prescribed reset, or supplies all required reset capacity.
The positive source density of a residue cylinder supplies no deterministic
frequency along one integer trajectory. First descent during arbitrary
reset blocks can occur before their final boundary. The new bank theorem
supplies explicit sufficient source classes and actual paired descent
receipts, not a proof of the Collatz conjecture.

## 7. Exact replay controls

The [self-contained checker](../experiments/collatz_bank_reset_20261019.py) uses integer coefficients, exact rational
arithmetic, and exception-raising checks that remain active under Python
optimization. It verifies 396 scheduled paired cylinders, arbitrary bank
orders of lengths 1 through 5, additional words up to 32 macroblocks,
six different reset schedules, both extra-reset choices, all lower
expanding prefixes, the exact lift cutoff, independent upper guard and
carry identities, and wrong-order controls. There are **1,234,513 explicit
checks**. Normal and optimized outputs agree byte for byte.

Run from the repository root (the checker itself only prints to stdout):

    python3 05-knowledge/experiments/collatz_bank_reset_20261019.py > /tmp/collatz-bank-normal.json
    python3 -O 05-knowledge/experiments/collatz_bank_reset_20261019.py > /tmp/collatz-bank-optimized.json
    cmp /tmp/collatz-bank-normal.json /tmp/collatz-bank-optimized.json
    cmp /tmp/collatz-bank-normal.json 05-knowledge/results/collatz_bank_reset_20261019.json
    cmp /tmp/collatz-bank-optimized.json 05-knowledge/results/collatz_bank_reset_20261019.optimized.json

Finite replay audits the formulas. The arbitrary-word theorems are proved
by (1)-(7), not by these finite controls.
