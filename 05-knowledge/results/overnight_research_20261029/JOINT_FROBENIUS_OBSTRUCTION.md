# A native norm polynomial excludes every fixed-name linear Frobenius lift

2026-10-29. **PROVED from FINITE-EXACT matrix/polynomial identities.**
Same-agent alternative-path verification. This is a universal linear
obstruction for the specified eighteen-dimensional ternary observer,
not a statement about every nonlinear map or every native operation.

Incoming9c876584a1 supplies
[JOINT_PHASE](../overnight_rank_return_20261029/JOINT_PHASE.md): its full
six-selector module has an explicit chart K, bulk clock T, native name
matrices M_pi, and an order-five map S with ST=T^487 S. That particular
S moves the name action. Its question whether a different LINEAR lift
could keep the native names fixed is now closed negatively.

Let R=M_120 be the supplied native three-cycle and let N=I+R+R^2.
The exact source matrices satisfy

    N=f(T),
    f(X)=2+2X^3+2X^4+X^5+X^6+X^7+X^8+2X^10+X^11.     (J1)

Although R commutes with T, the norm polynomial is not invariant under
the proposed clock exponent:

    rank_F3(f(T^487)-f(T))=5.                          (J2)

**Theorem.** No invertible F3-linear map A satisfies BOTH

    A T=T^487 A,       A R=R A.                        (J3)

**Proof.** The second equality and linearity imply A N=N A. The first
implies A f(T)=f(T^487) A. By J1 these give
(f(T^487)-f(T)) A=0. Invertibility contradicts J2. QED.

Thus even the single native three-cycle cannot remain fixed; requiring
the full six-element name action cannot repair it. The proof does not
enumerate candidate lifts or infer impossibility from failure of one
chosen splitting. It preserves the native relation N=f(T).

The mechanism is also visible in the incoming decomposition

    F3[X]/h^2 + F3[X]/h,    h=(X+1)g,
    g=X^5+X^4-X^3-X^2-1.

The polynomial f is divisible by h. In the h^2 quotient,

    gcd(f(X^487)-f(X), h^2)=(X+1)^2 g.                 (J4)

Consequently the discrepancy occupies precisely the five-dimensional
nilpotent F243 layer. The one-dimensional nilpotent factor and the
extra semisimple amplitudes do not remove it. The norm N itself has
rank6 and satisfies N^2=0, exactly as the characteristic-three group
norm predicts. Matrix and polynomial calculations verify J1--J4 by
different paths.

The successful inherited lift remains valid with transported names
A M_pi A^-1 and the appropriate source functional. It is the demand
to KEEP the original native action that fails. The source chart,
original selectors, and matrix direction are retained in the receipts.

There is a second, operation-level boundary. Incoming
[TERMINAL_SELECTOR S1](../overnight_rank_return_20261029/TERMINAL_SELECTOR.md)
proves that every native path-name permutation commutes with the bulk
clock. Every word in these generators therefore centralizes the clock.
It cannot realize a map that sends T to T^487, since T^487!=T. If such
a word descends to the old scalar observer, the same commutation
constraint descends with it. This excludes the declared bulk/name
language, not new complements or arbitrary vertex operations outside it.

This distinction connects the earlier centralizer results to a precise
new limitation: the Frobenius map is in the clock normalizer, and a
language contained in the centralizer cannot generate its nontrivial
clock exponent. It is not an identification with the separate Paley7
126-element action or the older LRC Frobenius no-go.

At higher integer precision, the independently proved
[image-cardinality obstruction](PERIODIC_PRECISION_OBSTRUCTION.md)
excludes even nonlinear conjugacies on the old full weighted carrier
at modulus243. J3 is instead a residue-field LINEAR obstruction with
the native names fixed. Neither statement proves an arithmetic Selmer
compatibility theorem or Sha finiteness.

Reproduce with `python -B joint_frobenius_obstruction.py` and -O.
The receipt stores the incoming canonical hash, norm polynomial, literal
matrices, rank-five witness and quotient-polynomial check. The incoming
positive lift is rechecked, so the obstruction does not accidentally
forbid transported-name lifts.
