# Level eleven: exact bridges between an eta product, finite geometry, and guarded memory

**PROVED elementary representation results; CITED modular-form and elliptic-curve
inputs; FINITE-EXACT independent controls. Universal Collatz coverage OPEN.**

This session extends the inherited [level-eleven recurrence](level11_short_20260922.md),
[golden/tournament transfer audit](five_eight_nine_transfer_20261005.md),
[paid controller](paid_guard_budget_20261005.md), and
[recursive trace construction](recursive_small_structures_20261004.md).
The closest proved mechanisms are Hecke recurrence and exact guarded affine
composition. The hostile examples are equal target representatives with opposite
payment and rational subdivision points that cycle. The corrected near miss is
confusing a Frobenius power trace with a prime-power Fourier coefficient. The
least-used sidecar is the integral tangent lattice, which matters even when two
rational matrices have the same characteristic polynomial.

Anchor: lossless Collatz certificate storage. Niche: the level-11 eta curve.
Wildcard: rational subdivision membership and marked modular dessins. The live
board is **history, terminal marker, source guard, lattice, local operator, graded
count**. The following connections explicitly retain different subsets of it.

## 1. Identify the product, without inferring a theorem from matching digits

For q=exp(2 pi i tau), the intended product and its corrected initial series are

```text
f(q)=q product_(n>=1)(1-q^n)^2(1-q^(11n))^2
    =eta(tau)^2 eta(11tau)^2
    =q-2q^2-q^3+2q^4+q^5+2q^6-2q^7+0q^8-2q^9-2q^10+q^11+... .
```

**CITED:** this is the normalized weight-two level-11 eigenform associated with
`E: y^2+y=x^3-x^2-10x-20`, a model of X0(11). See
[Elkies, Elliptic Curves in Nature, conductor11](https://people.math.harvard.edu/~elkies/nature.html)
and [LMFDB11.a2 = Cremona11a1](https://www.lmfdb.org/EllipticCurve/Q/11/a/2).
The same sources identify rational5-torsion; our script separately checks that
(5,5) has successive multiples (16,-61),(16,60),(5,-6),O. The curve has no CM.
The quadratic fields below belong to particular local operators, not to a
claimed global complex-multiplication structure.

For good primes p!=11, the cited relation is
`a_p=p+1-#E(F_p)`, with local factor `1-a_p T+p T^2`.
At11 split multiplicative reduction gives `1-T`. Our finite checks give

| place | local factor | exact consequence |
|---|---|---|
| 2 | 1+2T+2T^2 | recurrence roots -1+i,-1-i |
| 3 | 1+T+3T^2 | roots (-1+sqrt(-11))/2,(-1-sqrt(-11))/2 |
| 11 | 1-T | bad prime; not a good-prime quadratic factor |

Direct extension-field enumeration gives
`#E(F2)=#E(F3)=#E(F4)=#E(F8)=5`, `#E(F9)=15`, `#E(F16)=25`.
The field F9 has8 nonzero scalars, whereas this curve over F9 has15 points.
These are different sets; an eight-state field clock does not turn its curve
into an eight-point orbit.

## 2. A lossless formal-series coordinate, and what it forgets (PROVED)

Every F(q) in 1+q Z[[q]] has a unique expansion

```text
F(q)=product_(n>=1)(1-q^n)^(b_n),                 b_n in Z.
```

Proof: remove factors of degrees below n. The remaining coefficient of q^n is
`-b_n`; multiplying by `(1-q^n)^(-b_n)` kills it without changing any earlier
coefficient. This inductively gives existence and uniqueness. Negative exponents
use the ordinary integral binomial series. The proof is formal and needs no
analytic convergence.

An explicit inverse coordinate is

```text
-q F'(q)/F(q)=sum_(m>=1)c_m q^m,
c_m=sum_(d|m)d b_d,
b_n=(1/n) sum_(d|n)mu(n/d)c_d.
```

For this eta product divided by q,
`b_n=2+2*[11 divides n]` and
`c_m=2 sigma_1(m)+22 sigma_1(m/11)` with the second summand zero unless11|m.
Thus the complete coefficient sequence retains every indexed factor
multiplicity, including the extra multiplicity at multiples of11.

**Loss boundary:** commutative multiplication discards the order of the factors.
A coefficient zero such as a8=0 is cancellation, not a missing factor. A finite
coefficient prefix through degree N determines b1,...,bN only: multiplying by
`1-q^(N+1)` gives a distinct series with the same prefix. Additional global
structure must be supplied before a finite modular-form comparison proves an
infinite identity. This is a factor-multiplicity codec, not a controller-word
codec and not a proof of termination.

## 3. Two independently checkable bridges to the same local arithmetic

The [dessin and torsion package](level11_dessin_golden_torsion_20261005.md)
recovers two actual maps, rather than matching the set {2,3,11} by appearance.

* Modular generators on P1(F11) have cycle types `2^6`, `3^4`, `11.1`.
  Their degree12 dessin has genus1 and is the cited X0(11) construction.
  Marking infinity leaves its Borel stabilizer, of order55, acting on F11 by
  `x -> a^2 x+b`. It preserves the Paley tournament `x->y iff y-x is a nonzero
  square`. Forgetting the marked point loses this tournament: the full action
  is2-transitive. A graph relation needs this sidecar.
* The dyadic Hecke companion `A=[-2,-2;1,0]` reduces modulo3 to
  `[1,1;1,0]`, swap-conjugate to multiplication by phi in
  `F3[phi]/(phi^2-phi-1)=F9`. The irreducible characteristic polynomial
  determines Frobenius-at2 on E[3] up to conjugacy. Its order is8; projectivizing
  loses sign and gives order4. This F9 is an endomorphism coefficient field:
  full3-torsion point coordinates require F256, as independently enumerated.
  Naively identifying the same two integer matrices modulo9 fails their
  trace/determinant test.

The [centroid membership package](centroid_membership_20261005.md) discovers a
second local bridge. Its denominator9 decoding cycle has transverse return

```text
L=[-7,4;-5,-1],             charpoly(L)=t^2+8t+27.
```

On the index2 tangent sublattice with basis `(2(e0-e2), e1-e2)`, its matrix is
`L'=[-7,2;-10,-1]`. Put

```text
D=[-2,1;-5,1],    S=[0,1;1,2],    B=[-1,-3;1,0].
D S=S B,         det(S)=-1,       D^3=-L'.
```

Hence `det(I-TD)=1+T+3T^2`, exactly the eta curve's local factor at3.
This is an integral conjugacy **after** retaining the index2 lattice. It is
stronger than a shared discriminant, but still local: it does not establish
the global conductor from the geometric cycle. Also det(D)=3, so D mod3 is
singular; it is not the invertible golden clock from Frobenius-at2.

The earlier odd-square theorem still says that the primes p sharing their
odd-square bracket with2p are exactly {2,3,11}; see
[the original proof](procgen_numerology_20260923_odd_square_brackets.md).
The separate prime step `2^3+3=11` is recovered in
[the Mills seam](seam_mills_20260925.md). Neither construction yet maps its
defining predicate to the eta curve's conductor. The bridges above do not
retroactively establish that missing implication.

## 4. Trace doubling retains its observable (PROVED from the cited local factor)

Let alpha,beta be the local roots, so alpha+beta=a_p and alpha beta=p.
There are two useful sequences with the same recurrence but different seeds:

```text
power traces: s0=2, s1=a_p, s_(r+2)=a_p s_(r+1)-p s_r;
Hecke coefficients: h0=1, h1=a_p, h_(r+2)=a_p h_(r+1)-p h_r.
```

Here `s_r=alpha^r+beta^r` and `h_r=a_(p^r)`. Consequently
`s_(2r)=s_r^2-2p^r`; normalized by p^(r/2), this is `J -> J^2-2`.
It is the same reciprocal-eigenvalue identity inherited in the Collatz
word-doubling note. The matrix identity transfers; an integer trajectory,
source guard, or payment inequality does not follow from the trace.

The cheapest hostile is p3: `a9=-2`, whereas `s2=-5` and `#E(F9)=15`.
Likewise `a27=5`, `s3=8`, and `tr(L)=-8`. Comparing the period3 geometric
return directly with a27 would choose the wrong sequence. At2 the inherited
law `a_(2^(4j+r))=(-4)^j*(1,-2,2,0)_r` remains unchanged.

## 5. The controller supplies another exact series (PROVED)

The [translation decoder](translation_phase_decoder_20261005.md) proves that
the five letters H,G,A,B,L have disjoint last-translation residues modulo9.
Exact translation alone recovers the formal word. Native realizability is
equivalent to avoiding the adjacent pair LB. Thus the word and its complete
native guard can be reconstructed, including L's three extra source bits.
Adding K as a sixth independent letter breaks injectivity, since K=LG as a
guarded affine map. This is an alphabet-dependent theorem.

Weight H,G,A,B,L by their reduced dyadic denominator costs10,3,7,7,4.
Let W(z) count native words by total cost, including the empty word. The two
states are whether the previous letter is L. With

```text
X=z^3+2z^7+z^10, Y=z^3+z^7+z^10,
T(z)=[X,z^4;Y,z^4],
W(z)=[1,0](I-T(z))^(-1)[1;1]
    =1/(1-z^3-z^4-2z^7-z^10+z^11).
```

This derives the positive z11 term from the forbidden pair's cost4+7.
It supplies a rigorous history-counting counterpart to the eta product;
it does not equate their series (already their degree1 coefficients differ).
Euler inversion applies to W too, with initial nonzero exponents
`b3=-1,b4=-1,b7=-3,b10=-4,b11=-2`. The exponents encode the aggregate series.
They do not specify an individual controller history: A and B alone already
have equal cost7. Unweighted counts satisfy `N0=1,N1=5,Nr=5N_(r-1)-N_(r-2)`.

The incoming phase-tower work provides a further representation. For actual
valuation words with a **known number m of odd steps**, the real fractional part
of the exact translation recovers the whole word. If that fraction is N/2^A,
peel the first valuation from `v2((N-3^(m-1)) mod2^A)` and recurse on the tail.
The [decoder package](translation_phase_decoder_20261005.md) proves the general
positive odd-multiplier version, including finite membership and the absence of
antipodal pairs at fixed m. These facts explain the inherited exact Haar-start
second-moment diagonal. Here the compact group is R/Z and the denominator is
retained by the exact fraction. It is not the target's residue modulo3^m.
Dropping m breaks injectivity: words(2) and(1,1) have translations1/4 and5/4.

## 6. Transfer audit and next useful target

| source -> target | preserved predicate | loss / necessary sidecar | decisive control |
|---|---|---|---|
| Euler factors -> full coefficients | indexed multiplicities | order; finite-tail assumptions | signed exponents and invisible degree65 factor |
| finite controller word -> exact translation | ordered history | alphabet; native grammar; actual source for payment | K=LG; forbidden LB |
| finite fan word -> image of centroid | word with canonical ports | terminal port permutation | exact denominator3^(depth+1) |
| arbitrary rational point -> geometric decoder | rational coordinate | terminal membership may fail | denominator9 period3 |
| local companion -> golden F9 operator | mod3 linear action | prime/operator label; higher precision | mod9 nonconjugacy |
| denominator9 return -> local3 factor | integral linear recurrence on a sublattice | lattice index; local versus global; observable seeds | a9 != s2 |
| Collatz Fourier target phase -> target representative | canonical d when range/parity retained | source c, clock, slope and payment | 3->5 grows;13->5 shrinks |

The incoming [Pascal Fourier phase identities](collatz_fourier_experiments_20261004.md)
also preserve exact phases while changing depth. They complement the last row;
they do not restore its discarded source or certify a pointwise paid exit.

The geometric lane also produces a sharp limit warning with a positive result.
Let W=A001 and p=(1,2,6)/9. On the tangent plane the positive quadratic form
`Q=[5,-3;-3,4]` satisfies `L^T Q L=27Q`. Thus W contracts squared error by1/27,
and `W^j(center)` converges to p with squared error `4/27^(j+1)`.
Nevertheless each term is a different finite address with exact denominator
`3^(3j+1)`, and p itself is periodic under stripping, not a finite centroid.
This gives an explicit nonclosed finite-code language. Real contraction can
coexist with unbounded exact-address depth. This geometric contraction alone
supplies no well-founded arithmetic rank; any such use needs an additional
discrete stopping argument, possibly from a different coordinate. No Collatz
payment is inferred from this geometric norm.

**Next target:** use the exact translation decoder as a verifier front end:
decode, validate noLB and the native cylinder, recover the source-dependent
credit ledger, and only then emit a common-future receipt. A proposed extension
must supply its last-letter separator or record the new relation such as K=LG.
The unresolved mathematical obligation remains production of a paid guard for
every uncovered integer. Neither finite point counts nor better history storage
settles that coverage obligation. No coverage percentage is increased here.

## Reproduction and audit scope

Run `python -X utf8 -B 04-computation/experiments/eta11_lossless_coordinates_20261005.py`.
Use `--write` to refresh the [saved output](eta11_lossless_coordinates_20261005.out).
The stdlib-only companion compares direct product expansion through q257 with
an independent divisor recurrence; reverses its Euler exponents; tests signed
exponents and a finite-prefix hostile; checks coprime multiplicativity through
256 and good-prime point counts below100 by two paths; enumerates the displayed
extension fields directly; checks rational order5; and compares two weighted
word counts through cost64. Explicit checks also run under `python -O`.

The three linked peer packages contain the dessin/torsion, geometric membership,
and guarded translation proofs and their separate exact computations. Computed
agreement is FINITE-EXACT; universal formal identities have proofs above;
modularity and general Frobenius facts remain CITED inputs. No Lean claim.
