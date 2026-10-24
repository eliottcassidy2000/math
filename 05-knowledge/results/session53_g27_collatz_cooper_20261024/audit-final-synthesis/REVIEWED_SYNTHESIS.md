# Session 53: exact returns without forgetting the state

24 October 2026. This synthesis distinguishes proved deductions, inherited
mechanisms, exact finite experiments, and remaining questions. Detailed
proofs, checks and provenance are linked from the [session report](README.md).
No external novelty or peer-review claim is made.

## What 27 does in each problem

**Square-sum graph.** Adding vertex 27 permits 3 extensions of Hamiltonian
paths and 32 joins of two-path covers. Those covers expose 9 and 22 on distinct
components. G27 therefore has 35 paths up to reversal, but it is neither the
first traceable graph nor a new degree threshold. The meaningful state is
the exposed endpoint pairing, not the vertex count alone. The general
linear-forest polynomial transfer retains this state.

**Collatz.** The first 59 shortcut steps of 27 form an exact 37-odd-step
descent word. Its entire dyadic cylinder descends at shortcut time 59,
ordinary time 96. Repeating the word gives nested cylinders with an exact
valuation budget, not a collection covering every integer. Modulo 7 the
word induces 1-x, so two repeats restore the residue; the actual integer map
still contracts and is not an involution.

**Cooper coefficients.** The branch polynomial is
s^2=(1+t)(1-27t). Modulo 3 its degeneration gives a two-state Cartier swap;
modulo 9 it gives a full base 9 digit law. An exact prime-3 modular polynomial
also proves the sharp excellent-Frobenius improvement modulo 27. These are
literal operator and polynomial identities, not consequences of Collatz's
long orbit or G27's path count.

## The strongest new deduction in this session

For the normalized root G(t)=sqrt(sum s7(n)t^n)=sum c(n)t^n,

    c(p^2 n) = c(n)  (mod p^2)

for every prime p and n>=0. The inherited prime iff gives one-step return
at 7 and square residues modulo 7. The new inert calculation uses a
weight-three companion and an explicitly killed coordinate-lift error.
Independent audits cover both the new reduction and the remaining modular
inputs. The earlier conditional packet is retained, with its dependencies
discharged by a later audit rather than silently rewritten.

Put K=C_pG over Z/p^2. At inert primes G,K form a free rank-two module on
which Cartier swaps them; at split primes and 7 the module has rank one.
One unit coefficient of G-K proves freeness at inert primes. Thus the prime
dichotomy becomes an exact rank and minimal-period classification. G itself
is not generally the negative eigenvector; the eigenchannels are G+K and
G-K. The prime 2 belongs to rank one and is treated without dividing by 2.

For p different from 7, the determinant of this Cartier action is exactly
the sign (p/7): +1 on the rank-one branch and -1 on the rank-two branch.
This is also the orientation sign of multiplication by p on Paley 7.
The ramified prime 7 is excluded from that comparison because multiplication
by 7 is not a permutation of its residue field.

This is a genuine shared invariant, not an identification of operators.
The Paley 7 orientation character in a fixed Collatz word also takes values
in this sign group, but no map from these coefficient modules to native
tournament wiring or unrestricted integer trajectories has been proved.
CRT shows that every finite parity prefix occurs in every mod 7 residue.

## Exact digit budgets and generating functions

Modulo 9,

    G(t)=(1+2t+4t^2+3t^3+3t^4+6t^5)G(t^9).

Each base 9 digit contributes its own weight. Digits 6,7,8 force zero, and two
digits from 3,4,5 supply two factors of 3 and also force zero. Below 9^r there
are exactly 3^r units and r*3^r nonzero multiples of 3. The residue-zero count
has rational generating function

    1/(1-9z) - 1/(1-3z)^2.

All nine residue counts have explicit rational generating functions. The
six-phase unit oscillation comes from the order of 2 modulo 9, with exact
eigenvalues 3,1,0,0 and 2exp(plus or minus i*pi/3). The initial place count
r=0 is treated separately where necessary.

Modulo 4, the only positive indices with nonzero coefficients are powers
of 2, with residue 2. Hence 4 and 9, and therefore 36, divide c(n) on
natural-density-one sets. This conclusion needs no independence of the
binary and base 9 exceptional sets. Their exact intersection remains a
separate question, not a consequence of density.

## A second sharp theorem, with the correct coordinate

The exact modular polynomial

    27-18x^2-8s x^3+(24t-1)x^4

vanishes at x=G/G^sigma. Formal unit-root algebra gives

    F(t)=(1+4t+48t^2)F(t^sigma)  (mod 27),
    t^sigma=t(q(t)^3).

This closes the named local modulo 27 question in the inherited note,
under its audited modular framework. It cannot be raised uniformly to 81:
the t^3 discrepancy is756. Nor may t^sigma be replaced by t^3: that false
replacement has residue 18 already at t^5. Exact coordinates matter as much
as the modulus.

## What the repeated Collatz certificate does and does not buy

For q=2^59, a=3^37, b=23q-27a and d=q-a, the maximum number of complete
copies of the specified word is floor(v_2(d*n-b)/59). The nonintegral fixed
point b/d is between 8 and 9. The depth-k cylinder has density q^(-k), while
all positive-depth cylinders together still cover only the first cylinder.
Arbitrarily deep certified repetition is therefore compatible with very
small coverage. It supplies neither universal descent nor uniform 153N
coverage, and it does not establish any unrestricted stopping-time law.

This is the useful shared lesson: an endpoint pairing, a digit address,
or a valuation clock can make a local transition exact. Forgetting that
state can create a spurious global claim. The repeated number 27 is a prompt
to inspect these operations, not a proof that the objects are equivalent.

## Precise next obligations

1. Prove or refute the higher-depth two-step bound
   c(p^r n)=c(p^(r-2)n) modulo p^(2r-2). The r=2 base is now proved;
   15,407 checked comparisons do not close r>2.
2. Explain the calibrated quadratic next-digit operator
   (C_p^2G-G)/p^2 = a_p*(t*d/dt)^2G modulo p at odd primes. Its 7,995
   comparisons outside n=1 coefficient calibration pass, but a proof must
   control second-order lift errors and
   Eisenstein corrections, without confusing Euler differentiation with
   the formal derivative in a different Gessel-Lucas statement.
3. Keep graph growth endpoint-sensitive. The exact polynomial response
   offers a reusable local counting object; it is not yet an efficient
   global square-sum Hamiltonian algorithm or tournament-class prediction.
4. Keep Collatz coverage separate from repeatability. Any proposed uniform
   argument must account for the shrinking cylinders and actual clocks,
   rather than recycling the same local certificate at multiple scales.
