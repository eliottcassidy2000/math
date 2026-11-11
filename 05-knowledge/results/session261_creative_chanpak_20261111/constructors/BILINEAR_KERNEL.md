# A physical two-port exchange kernel and a norm of discriminant −1316

**PROVED**, with explicit finite skeleton and literal graph certificates.
This fourth construction wave enlarges the preceding dyadic family to
arbitrary tournament modules. Chan–Pak remains OPEN. No external novelty or
new numerical-coverage claim is made.

## 1. Literal quotient and universal substitution law

Let Q have outgoing-bit rows

    Q=(198,196,136,211,167,207,20,64).

Its two ports are vertices5 and7. Replace them by arbitrary nonempty
tournaments A and B, keeping all quotient orientations uniform toward each
module. Denote the resulting tournament Q[A,B]. It has |A|+|B|+6 vertices.
Let h=H(A), g=H(B), c=C_2(A), d=C_2(B), where C_2 counts ordered nonempty
spanning two-path covers.

**Theorem.** For every such A,B,

    H(Q[A,B])=149hg+72(cg+hd)+37cd.                  (1)

All coefficients are positive integers and this is an actual tournament
substitution. In particular it supplies a count-preserving module exchange

    Q[A,B]  -->  Q[B,A].

No claim that the marked response or the distinguished-port continuation
is preserved is included. The equality is a direct Hamilton-count return
across changing representative graphs.

### Derivation from the earlier all-length family

In the seven-vertex base S=(6,4,8,19,103,79,15), vertices5→6 are consecutive
unselected transitive twins. Cloning5 to length m and adjoining6 makes a
single module of length a=m+1. The t-vertex tail is the second module, of
length b=t. Taking m=0,t=1 leaves the displayed Q, with ports5 and7.
The proved tail formula from CONTINUATION.md is therefore exactly

    H(Q[TT_a,TT_b])
      =37*2^(a+b)−2^(a+1)−2^(b+1)+9,  a,b>=1.      (2)

A joint compressed-run expansion gives

    H(Q[A,B])=sum_{r,s>=1} c_rs C_r(A) C_s(B),

with finitely many nonnegative integral skeleton coefficients. Suppose a
positive c_rs had r>=3. Fix b>=s and let a grow. Its contribution on TT_a,TT_b
would grow as r^a, contradicting the O(2^a) bound in (2). Thus all r>=3
coefficients vanish, and symmetrically all s>=3 coefficients vanish.
For a nonempty transitive a-chain, C_1=1 and C_2=2^a−2. Substituting these in
(2) yields

    149+72 C_2(TT_a)+72 C_2(TT_b)+37 C_2(TT_a) C_2(TT_b),

which identifies the four remaining coefficients and proves (1).

This is an all-size argument. The finite table through r,s=3 is a certificate
of its low coefficients, not a claim that testing three run numbers alone
would exclude larger ones. The exclusion uses positivity and the proved
all-chain formula.

## 2. A quadratic norm and its actual domain

Use Q_2(A)=C_2(A)+2H(A), which allows empty path slots, and define

    x(A)=37Q_2(A)−2H(A).

Then (1) is equivalently

    37H(Q[A,B])=x(A)x(B)+329H(A)H(B).                (3)

For A=B this becomes the integral norm identity

    37H(Q[A,A])=x(A)^2+329H(A)^2.                   (4)

The discriminant of X^2+329Y^2 is−1316;329=7*47. The arithmetic identity is
physical because the two coordinates come from the same actual tournament.
An arbitrary representation by that form does not supply the needed A.
For example a singleton gives H(Q[1,1])=149; a directed triangle in both
ports gives5265. Both are counted literally in the receipt.

The inherited path-cut injection gives Q_2(A)>=(|A|+1)H(A), and consequently

    x(A)>=(37|A|+35)H(A).                           (5)

Thus the diagonal grammar has a complete finite decision procedure for a
fixed target N: enumerate positive integers x,h satisfying37N=x^2+329h^2,
then test tournaments of order at most floor((x/h−35)/37) for H=h and
Q_2=(x+2h)/37. Reject nonintegral or negative bounds. This is a decision
procedure for this grammar alone, not for arbitrary tournament counts.
Its construction is parallel to the inherited norm-interface decision
principle in Session260, with a different actual quotient and norm.

The diagonal also has the immediate residue restriction N mod7 in{0,1,2,4},
from37N=x^2 mod7. This restriction does not apply to the unrestricted
bilinear kernel A!=B, where x(A)x(B) need not be a square.

## 3. Complete classification of the nonnegative-offset rows

Take B=TT_b, b>=1. Equation (1) reduces to

    H(Q[A,TT_b])
      =(37Q_2(A)−2H(A))*2^b +9H(A)−2Q_2(A).       (6)

**Theorem.** A row of (6) has nonnegative constant offset iff A is a
singleton, TT_2, or C_3. Their offsets are5,1,3 respectively. In particular
**TT_2 is the unique unit-offset row** among all nonempty tournaments A.

Proof: a nonnegative offset gives Q_2(A)<=9H(A)/2. Combine this with the
cut bound to obtain |A|+1<=9/2, hence |A|<=3. The complete finite base is:

    singleton: H=1,Q_2=2, offset5;
    TT_2:      H=1,Q_2=4, offset1;
    TT_3:      H=1,Q_2=8, offset−7;
    C_3:       H=3,Q_2=12,offset3.

These are all tournaments on at most three vertices, up to isomorphism.
The size reduction, not a finite search cutoff, proves the classification.

This explains the earlier slices exactly. A singleton yields9*2^k+5 after
an index shift, TT_2 gives the inherited slope73 unit ray, and TT_3 gives
the direct family7(21*2^k−1). The latter is the first negative-offset row,
not an accidental factorization at one integer. Using an arbitrarily more
complicated A cannot produce a new nonnegative-offset row in this kernel.
It can still produce new direct-count rows with retained negative offsets.

## 4. What the exchange does and does not supply

For chain modules, writing the state as (a,b) gives two actual clone moves

    increment a: N -> 2N +2^(b+1)−9,
    increment b: N -> 2N +2^(a+1)−9.

Exchanging the modules preserves N and swaps which offset belongs to the
specified next port. It is therefore a count-level state return with an
explicit sidecar: the other module length or, for general modules, its
(H,C_2) profile. The offsets include−5,−1,7,23,... and cannot be treated as
an arbitrary digit alphabet. Both directions already exist on the carrier;
the exchange does not secretly create a missing +1 branch.

Transfer ledger: source = the failed unit-seed binary family; target = a
universal two-module count kernel; map = coarsen actual twins, then replace
both chain modules by arbitrary tournaments; preserved = the exact positive
run expansion and Hamilton count under exchange; lost = selected seed state
and a distinguished port; sidecar = both two-coordinate module profiles;
cheapest test = the four joint run coefficients plus arbitrary nontransitive
module controls. All of these are supplied.

The next genuinely different operation would alter the constant9 or the
coefficient37 while preserving a count across this exchange, so the two
port offsets could cover missing digits. Equations (5) and (6) explain why
simply inserting larger modules into this same kernel does not do so.

## Reproduction

    python3 bilinear_kernel_certificate.py

The program checks the explicit joint skeleton coefficients and16 literal
module pairs drawn from the singleton, TT_2, TT_3, and C_3. The receipt also
contains the complete small base justified by the all-order cut bound.
