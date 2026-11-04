# Modulo-eight ambient transfer for even-order modules

**PROVED:** Every same-size even-order module pair has its Hamilton-count difference multiplied by an odd unit modulo eight under arbitrary ambient substitution. Consequently every nonzero modulo-eight physical two-port separator of an even-order port core survives every ambient quotient. The odd-order correction is explicit and involves Hamilton-cycle parity. No full-response reconstruction theorem is claimed.

## 1. Classical input and conventions

Write H(A) for the number of directed Hamilton paths and q_k(A) for ordered spanning covers by k nonempty directed paths. Let o(A) count all directed cycles of odd length at least three, up to cyclic rotation. For odd order N, let z(A) count directed Hamilton cycles, again up to cyclic rotation. For order one set z(A)=0. For even order set z(A)=0 in the formulas below: even Hamilton cycles are not among o(A).

The sole external counting input is the classical modulo-four refinement in Grinberg and Stanley, [The Rédei–Berge symmetric function of a directed graph, detailed version](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf), Theorem 7.1, printed page 119 (June 13, 2026 version):

    H(A) = 1+2o(A) (mod 4).                              (1)

This also implies that every tournament H(A) is odd. The more general odd-cycle expansion is their Theorem 1.39; it is not being claimed as new here. All profile congruences and ambient consequences below are derived explicitly.

## 2. Exact profile information modulo eight

Let S(n,k) be the Stirling number of the second kind, zero outside its usual support. For every nonempty N-vertex tournament A,

    q_2(A) = 2H(A)+2^N-4-4z(A) (mod 8).                 (2)

For even N,

    q_3(A) = 6S(N,3) (mod 8),                            (3e)

while for odd N,

    q_3(A) = 6S(N,3)+2H(A)-2-4z(A) (mod 8).             (3o)

Finally q_k(A)=0 (mod 8) for every k>=4.

To prove these statements, first forget the ordering of the path slots. The exact partition identity gives

    q_k(A)/k! = sum_(unordered k-partitions P of V)
                                  product_(B in P) H(A[B]).

Using (1) on each nonempty part, modulo four the product is one plus twice the total number of odd cycles contained in individual parts. A fixed odd cycle of length ell is contained in one part in exactly S(N-ell+1,k) partitions: contract its entire vertex set to one atom. Therefore

    q_k(A)/k! = S(N,k)
       +2 sum_(odd directed cycles gamma) S(N-|gamma|+1,k) (mod 4).  (4)

For k=2, S(t,2) is odd for every t>=2, and S(1,2)=0. Thus exactly the proper odd cycles contribute, giving q_2=2S(N,2)+4(o-z) modulo eight. Since 2S(N,2)=2^N-2 and 4o=2H-2 modulo eight by (1), this is (2).

For k=3, S(t,3) is even at even t and odd at odd t>=3, with S(1,3)=0. This follows immediately from the Stirling recurrence modulo two, starting at S(2,3)=0 and using S(t,2) odd for t>=2. Since ell is odd, t=N-ell+1 has the parity of N. Multiplying (4) by six therefore gives (3e) or (3o), respectively. The terms ignored modulo four become multiples of 24, so the modulo-eight conclusion is justified. Finally q_k is divisible by k! because ordered cover slots can be permuted freely, and 8 divides k! for k>=4.

All boundary cases N=1,2,3 are included with the stated zero conventions. In particular the even-order formulas do not silently require N>=4.

## 3. Arbitrary ambient substitution

Fix a tournament Q with a distinguished root r. For any nonempty module A, the maximal-A-run decomposition gives exact nonnegative integers c_k depending only on (Q,r), with

    H(Q[A])=sum_k c_k q_k(A),       c_1=H(Q).

Let A,B have the same EVEN order. Write Delta H=H(A)-H(B). Equations (2)–(3e) give

    Delta q_1=Delta H,
    Delta q_2=2 Delta H,
    Delta q_k=0 for k>=3                         (all modulo 8).

Consequently

    H(Q[A])-H(Q[B]) = (c_1+2c_2) Delta H (mod 8).          (5)

The multiplier is odd. More concretely, because q(TT2)=(1,2),

    c_1+2c_2=H(Q[TT2]).

Thus the ambient difference is zero modulo eight if and only if the core difference is zero modulo eight. A core gap congruent to four remains congruent to four in every ambient quotient. This requires no sign or coordinatewise domination of the integer profile differences.

### Physical two-port lift

Suppose T[M] and U[M] differ only inside the same named module M, with the same exterior incidences. Take any actual separator on M plus x->y. Give both ports M's original incidences to every old exterior vertex. Then M together with x,y is a common module in both extensions, and the exterior rooted quotient is identical. If |M| is even, the port core has even order |M|+2, so (5) lifts every nonzero modulo-eight local query difference to the same physical full named response in the ambient tournaments.

No cubic or quartic term has been discarded: this is a statement about literal Hamilton counts of the entire two-port extension. In particular it upgrades any proved nonzero modulo-eight separator for a uniform TT2 inflation Q[TT2,...,TT2] to arbitrary ambient contexts, since the port-core order is 2|Q|+2. The scope of the local separator theorem, such as a restriction on the quotient order, is unchanged; this result supplies only the formerly missing ambient lift.

## 4. Odd order: exact correction and a necessary boundary

For same-size ODD-order modules A,B, write Delta z=z(A)-z(B). Equations (2)–(3o) instead give

    Delta q_2=Delta q_3=2 Delta H-4 Delta z (mod 8).

Hence

    H(Q[A])-H(Q[B])
      = (c_1+2c_2+2c_3) Delta H-4(c_2+c_3)Delta z (mod 8). (6)

The first multiplier is odd. If Delta z is even, the same nonzero-modulo-eight transfer conclusion follows. Without that condition, a residue-four core gap can disappear modulo eight in an ambient quotient. This is a modular obstruction, not an exact full-response collision.

For example the actual port cores from the physical query (3,12) on C3[TT2,TT2,TT3] and its quotient converse have

    (Delta q_1,Delta q_2,Delta q_3)=(-68,-396,-1116),

all congruent to four modulo eight. Substitution into the directed triangle context, whose coefficients are c_1=3,c_2=1, gives exact ambient difference -600, zero modulo eight. The analogous C3 third-child example has differences (-84,-468,-1188) and ambient gap -720. Both remain nonzero as integers; neither is an exact collision.

## 5. Controls

The verifier checks all labelled tournaments through order five, additional fixed-seed orders six through nine, and both the even and odd formulas using independently counted Hamilton cycles. It checks literal even-core substitutions in every rooted quotient through order four, including residue-two and residue-four core pairs, and retains the full profiles and actual strong-context counts for four cyclic port-core controls.

Run `python verify.py > receipt.json` and `python -O verify.py > receipt.optimized.json`. The sibling `covered_response` and `response` code must be present. All arithmetic is integral; explicit exceptions remain enabled with -O. Finite checks support the conventions; the all-order proofs are (1)–(6).
