# Independent audit of the Hamiltonian transport and carrier criteria

**Status: INDEPENDENT ANALYTIC/SOURCE AUDIT PASS + EXACT REPLAY.** This
audit accepts the finite-jet identities and the fixed/universal source-form
criteria in [the primary note](planar_jc48_sep06_hamiltonian.md), including
the stronger non-local-nilpotence statement for every nonconstant fixed
source carrier. No planar Keller, global flow-termination, or arbitrary
source-entry conclusion is inferred.

## 1. Finite jets and the terminal representative

The source uses the fixed coefficients of `A,C` in rows zero through
three. Those coefficients are not zero. The draft phrase that could
suggest otherwise was repaired before this acceptance; the current scope
and the nonzero literal arrays in the verifier agree. Later background
rows start at four and cannot affect the displayed correction through
row fifteen, since `S_t` starts at twelve and `S_x` at thirteen.

The five extra terms of `Sstar-S0` all begin in generator row sixteen.
They can therefore change coordinate row fifteen while retaining every
earlier correction row. The stated terminal difference is a common
multiple of `(A0_x,C0_x)`: its scalar is
`16x R(x)/419175`. This explains the tangent cancellation, while the
literal `P_2/P_3` witnesses independently establish depth. The source
checks each witness's depth and valuation before substituting it.

The complete linear source packet and bracket identities have the
correct cutoffs: source equation modulo `t^16`, bracket modulo `t^15`.
The quadratic correction terms have valuation at least twenty-four in
the source equation and twenty-three in the bracket. Thus the finite
parameter statement is valid for every scalar parameter. It does not
silently require that parameter to be infinitesimal.

The Hamiltonian derivation for `S0` raises `t`-adic order by at least
twelve, so its exponential exists in the stated completed source ring
and has the same finite jet. Formal integration and a polynomial
additive-group action remain separate statements.

The old representative's `P_0` obstruction is complete. The background
Jacobian has constant term one; hence equality of both coordinate
corrections modulo `t^16` forces both derivatives of the difference of
their Hamiltonians to vanish modulo `t^16`. The difference is therefore
a scalar modulo `t^17`. At each valuation, the leading monomials of
`p^c y^e` have distinct `x` exponents. Recursive elimination leaves
`-52451x^9 t^16/3105`; a new `P_0` generator monomial in that row has
`e<=8`. Earlier cancellation cannot supply the missing exponent.
The terminal adjustment, rather than an identification of unknown depth
with generator depth, is what makes `S0` possible.

## 2. Independent Poisson reconstruction and denominator separation

I differentiated the actual source substitution independently and obtained

```
{u,p}=2xp,       {u,y}=3up,       {p,y}=-tp,
D=p^3-y^2=tp^2,
t=D/p^2,        x=yp/D,          u=y^2/D.
```

Thus the primary rational formula is exactly

```
{-u/2+H,S}=-py E(S)/(2D)-(D/p) J_(p,y)(H,S),
E=2p partial_p+3y partial_y.
```

The discarded `2y` formula survives only as an explicitly rejected
hostile and a correction explanation. It is not used by the finite-jet
derivatives or either carrier proof.

The cusp residue for `S0` is nonzero. At the generic point of the prime
`D`, the second term is regular because `p` is a unit. It cannot cancel
the first term's pole. This gives the asserted failure for every `H`.
These poles concern polynomial membership in the compressed `(p,y)`
ring; they are not poles of the original polynomial `(x,t)` response.

## 3. Both exact carrier criteria

For fixed `H`, regularity at `D` forces `D|E(S)`. In the cusp quotient
`k[p,y]/D=k[z^2,z^3]`, the operator is `z partial_z`; its kernel in
characteristic zero consists exactly of constants. Hence
`S=c+D R`. The first term then becomes polynomial, and coprimality of
`p,D` makes the remaining condition precisely `p|J(H,S)`.
This independently proves both directions of the fixed-`H` criterion.

For the universal criterion it is enough to test `H=0,p,y`. Subtracting
the first response from the other two gives `p|S_y,S_p`, respectively.
Writing `S` by powers of `p` shows `S=c+p^2 R`: the degree-zero
coefficient is constant in `y`, and the degree-one coefficient vanishes.
The cusp condition gives `S=c'+D R'`; evaluating at the origin identifies
the constants. Coprimality then gives `S-c` divisible by `p^2D`.
The explicit polynomial expression in the primary note verifies the
converse for all `H`, including both derivative slots of the arbitrary
polynomial `R`.

The positive control `H=0,S=D` satisfies the fixed criterion and has
response `-3py`, while lying outside the universal carrier. Thus the
two quantifiers are genuinely different. Only the universal carrier is
claimed to preserve the whole ring `k[p,y]` under its derivation, which
justifies the further coefficientwise formal source-form statement.

## 4. Local nilpotence and the exact stopping scope

The general obstruction is valid already for every nonconstant
`S=c+D R`, without the additional fixed-`H` divisibility. Its source
factorization is

```
S-c=t^3(1+x^2t)^2 R(p,y).
```

If its Hamiltonian derivation were locally nilpotent, the usual derivative
degree is additive on products by the highest surviving Leibniz term.
Its kernel is therefore factorially closed. Since `S-c` is invariant,
both `t` and `1+x^2t` would be invariant. Subtracting one and factoring
would also make `x` invariant. The derivation would be zero, making `S`
constant; the birational inverse proves injectivity of the substitution.
This is a contradiction for every claimed nonconstant carrier.

The separate nontermination proof for `S0` is also sound. Its unique
highest total-degree term is `-2848(xt)^25/6075`. Bracketing a leading
monomial increases both exponents by twenty-four and multiplies its
coefficient by `(-2848/243)(a-b)`. Starting from `x`, the exponent
difference stays one. Lower-degree Hamiltonian terms cannot reach this
highest degree, so the displayed nonzero leading term holds for every
iterate. This rules out a polynomial additive-group action generated by
that derivation. It does not rule out isolated specializations, other
generators with the same finite jet, or maps from a different mechanism.

## 5. Source review and separate exact calculation

The complete primary source was read. Its only inherited mathematical
input is the two literal old depth tables, extracted by `ast.literal_eval`
from the hash-pinned weight14 producer. That producer is neither imported
nor executed. Ordinary differentiation and full coefficient comparisons
check every declared finite jet. The symbolic formal derivative slots
in the repaired carrier identity retain both denominator contributions.

In addition to the full source replay, I reconstructed the substituted
polynomials using a separate standard-library `Fraction` dictionary
implementation of multiplication and differentiation, with no SymPy or
producer import. Thirty exact checks recover the source packet modulo
`t^16`, bracket response modulo `t^15`, background bracket `1+O(t^3)`,
both leading corrections and all eight degree caps, the complete row-sixteen
triangular obstruction, all three source Poisson brackets, both carrier
factorizations, and five leading-iterate controls. Its row-fifteen
corrections are independently

```
dA0: -6400x^4/729-934504x^6/139725-18928x^8/1215
      +289198x^10/3105+56x^12/3,
dC0: 31232x^5/1215+247744x^7/5175-6706382x^9/46575
      -133882x^11/1035-14x^13.
```

This separate calculation is an audit of the displayed identities;
the all-order conclusions rest on Sections 3–4's analytic arguments.

The normal, optimized and frozen primary outputs are byte-identical,
with **241 always-active exact gates**. Reproduce from the repository root:

```bash
python3 -B 04-computation/planar_jc48_sep06_hamiltonian.py > /tmp/planar_jc48_hamiltonian_audit_normal.out
python3 -B -O 04-computation/planar_jc48_sep06_hamiltonian.py > /tmp/planar_jc48_hamiltonian_audit_optimized.out
cmp /tmp/planar_jc48_hamiltonian_audit_normal.out /tmp/planar_jc48_hamiltonian_audit_optimized.out
cmp /tmp/planar_jc48_hamiltonian_audit_normal.out 05-knowledge/results/planar_jc48_sep06_hamiltonian.out
```

```
source SHA256
ba65f911e64907cd41d7a950263e281a8da370045c9da40668796f1245a8df6b
frozen and both replay outputs SHA256
842d3d0b6c8300900c44cb34317dde4f7f64adac72ac18067cb8eafc736952a3
semantic SHA256
270f7657c66819230a5990639875d1faedb15e453b8ac6b22e22d6220120ea57
inherited depth source SHA256
a4a140ab538620e7885d8b77758c02e16a5c3a8de7e53e3daac51f01bda321f7
```

The primary source and output were not edited by this audit. Its accepted
scope is finite transport, exact infinitesimal source-form carriers and
the stated formal/additive-action boundaries.
