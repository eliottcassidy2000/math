# Independent audit: all-exponent binomial algebraic exactness

**Status: PROVED / INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.**
Root independently read the full [binomial proof](planar_jc48_sep08_binomial_exactness.md), its standalone source and the cited classical primary. Both source modes reproduce the frozen 1,654-gate output byte for byte. A separate calculation in the original u-coordinate passes 1,170 full Laurent coefficient systems, including odd total degrees. The all-exponent result is analytic; these finite systems independently check its mechanism.

## Exact statement and field

For A>=0, r>=1 and nonzero complex a,b, let N=u^A(a u^r+b) and K=C(u)(y), y²=N. The accepted iff is

    du/y exact in K
       iff (A,r)=(0,1), or A=r+2+2r ell for an integer ell>=0.

The residual binomial is squarefree and has r nonzero roots. Thus N is not square, even when its u-power contains a large square factor, and K is the connected radical field. The root translation is an ordinary curve coordinate; it is not being extended to an automorphism of the original graph-complement surface.

Allowing a primitive in a further finite algebraic extension makes no difference. In characteristic zero the derivation extends uniquely and commutes with field trace. Dividing the trace of a primitive by the extension degree returns a primitive in K. This descent does not assume the larger extension is Galois. Conversely a primitive in K is already algebraic.

The low-degree exceptions are complete. Degree one has the primitive 2y/a. Degree two has two simple roots and nonzero infinity residues on its two branches; all other possible degree-two binomial exponent pairs are explicitly included. For odd degree at least three, the inherited pole-degree argument forces a pure root, which this binomial cannot be. No numerical range justifies that unbounded statement.

## Complete primitive space and the exponent obstruction

For even n=A+r=2k>=4, inversion gives X=1/u, Y=y/u^k, Y²=a+bX^r and the actual differential -X^(k-2)dX/Y. The numerator exponent k-2 is nonnegative. This differential is regular on the entire smooth affine model: at a branch point Y is a local parameter, and X=0 is not a branch point because a is nonzero.

A primitive has no poles at any affine point, since differentiating a pole increases its order. Normality therefore places it in the full affine coordinate ring C[X,Y]/(Y²-a-bX^r), not just a guessed finite subspace. Taking its odd part leaves YB(X), with B polynomial, and preserves the desired derivative. The complete equation is

    (a+bX^r)B' + (rb/2)X^(r-1)B = -X^j,   j=k-2.

For an input monomial X^q the two output degrees are q-1 and q+r-1, with coefficients aq and b(q+r/2). The latter is nonzero for every q>=0. It proves injectivity and forces deg B=j-r+1. Different input residue classes modulo r have disjoint output residue classes, so injectivity eliminates every class except j+1 modulo r. The least nonzero exponent of B must be zero: otherwise its lower derivative term is nonzero, cannot be cancelled by any smaller input, and lies strictly below the target exponent. Thus the only remaining class is zero and j+1=r(ell+1), with ell>=0. This is precisely the accepted staircase. The r=1 case obeys the same argument without a fictitious decomposition into several classes.

The recurrence in the primary cancels all lower coefficients and sets the top one to -1. Every denominator is a nonzero integer times a power of b. It consequently supplies the primitive over the original coefficient field, even for parameter-dependent coefficients, without adjoining square roots of parameters. I independently checked its sign after returning to the u-coordinate: for a primitive yR, the equation is NR'+N'R/2=1. The pure quadratic mate -R/(2t) has the same positive Jacobian sign.

## Genus, poles, degenerations and claims not implied

For general inputs the branch count is r+(A mod2)+((A+r) mod2). On the accepted staircase A and r have the same parity and n is even, giving genus floor((r-1)/2). This genus is not assigned to rejected exponent pairs. In particular A=1,r=2 has genus one and detects the omitted infinity branch that a careless formula would lose.

The primitive has poles only above u=0. If A is odd the normalized pole order is A-2 at one point; if A is even there are two poles of order A/2-1. In either case the total is A-2=r(2ell+1). The statement is confined to the staircase; the separate degree-one primitive has its pole at infinity.

At fixed r,a,b, changing ell retains the same radical curve while changing the differential. The primary correctly records the same-field hostile at A=5,7,11 and r=3: exact, nonexact, exact. The missing datum in a field-only argument is the exponent of the differential after inversion.

Every nonzero member of the stated sparse pencils is exact. When a coefficient vanishes, the primary switches to a pure-power primitive and does not evaluate the generic formula across its 1/b denominator. Neither exponent endpoint is two. The linear pencil's constant endpoint is included. Through degree eight the six spaces agree with the already proved exact-pencil classification. No exhaustion of all higher-degree pencils is inferred.

For F=N(u)t², the leading-field necessity and the explicit primitive give a rational-mate iff. Lower t-coefficients would change the original fibre and are not covered by that sufficiency. Polynomiality or global regularity of the displayed rational mate is not asserted. None of these conclusions closes the planar Jacobian conjecture.

## Classical antecedent checked directly

I opened the original [Tchebichef 1853 paper](https://www.numdam.org/item/JMPA_1853_1_18__87_0.pdf). Section I, printed pp87–88, separates algebraic from logarithmic contributions. Section VIII, printed pp106–108, treats rational-exponent binomial integration. The primary correctly gives this classical context while proving its stronger algebraic-exactness predicate directly. Elementary integrability, which permits logarithms, is not substituted for exactness in a finite algebraic function field. No external priority claim is accepted or made here.

## Independent original-coordinate calculation

The producer's finite matrix is formed after inversion. I instead retained the original polynomial N and looked for the odd primitive yR(u). For n>=3 a primitive has no infinity pole and no finite pole except possibly above zero. Consequently R is a Laurent polynomial with all exponents in

    min(0,1-A) <= e <= -ceil(n/2).

This full interval follows from the actual local orders. At a simple residual root R cannot have a pole, because multiplying it by y would still leave a pole in the normalized coordinate. At zero the primitive pole bound yields e>=1-A; the A=0 case is ordinary regularity. At infinity regularity of yR yields e<=-ceil(n/2). A negative-length interval means there is no nonzero candidate, not an omitted case. No predicted congruence is imposed on the interval.

For each e, the original operator gives

    N(u)(u^e)' + N'(u)u^e/2
       = b(e+A/2)u^(A+e-1)+a(e+n/2)u^(n+e-1).

I solved these complete rational coefficient systems by a separate full-row Gaussian elimination. Every A=0..48, r=1..12 with n>=3 was checked at (a,b)=(3,-2) and (7,11). All 1,170 systems agree with the analytic iff, including the odd-degree failures. This uses neither the producer's matrix code nor its recurrence. The accepted values remain a finite control, not an extension of the theorem by sampling.

The exact independent calculation can be reproduced with:

```python
from fractions import Fraction as F
from hashlib import sha256
import json

def solvable(rows,cols):
    pivot=0
    for col in range(cols):
        i=next((i for i in range(pivot,len(rows)) if rows[i][col]),None)
        if i is None:continue
        rows[pivot],rows[i]=rows[i],rows[pivot]
        z=rows[pivot][col];rows[pivot]=[v/z for v in rows[pivot]]
        for k in range(len(rows)):
            if k!=pivot and rows[k][col]:
                z=rows[k][col];rows[k]=[a-z*b for a,b in zip(rows[k],rows[pivot])]
        pivot+=1
    return all(any(row[:cols]) or not row[-1] for row in rows)
records=[]
for A in range(49):
 for r in range(1,13):
  n=A+r
  if n<3:continue
  powers=list(range(min(0,1-A),-((n+1)//2)+1))
  exponents=sorted({0}|{A+e-1 for e in powers}|{n+e-1 for e in powers})
  for a,b in [(F(3),F(-2)),(F(7),F(11))]:
   rows=[]
   for j in exponents:
    row=[]
    for e in powers:
     row.append((b*(e+F(A,2)) if j==A+e-1 else F(0))+(a*(e+F(n,2)) if j==n+e-1 else F(0)))
    row.append(F(j==0));rows.append(row)
   actual=solvable(rows,len(powers))
   expected=A>=r+2 and (A-r-2)%(2*r)==0
   if actual!=expected:raise RuntimeError((A,r,a,b,actual,expected,powers))
   records.append([A,r,str(a),str(b),actual])
print('Independent original-coordinate Laurent systems:',len(records),'PASS')
print('record sha256:',sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
```

The independent record SHA256 is `edfe15c63b1971db1b30a2f31e56d664dc34255e3be796814aed384df9fa61e4`.

## Frozen source and replay audit

I read every source gate. The producer tests all 588 declared exponent pairs, 584 unrestricted exact systems at two coefficient choices, all 48 symbolic primitive cases in both coordinates, the complete displayed monomial operator range, parameter rationality, pure-power endpoints and the named hostile examples. Always-active checks use explicit exceptions and remain enabled under optimization. Some bookkeeping gates only confirm the declared scope; they are not counted as separate analytic proofs.

Reproduce both modes from the worktree root:

    python3 04-computation/planar_jc48_sep08_binomial_exactness.py
    python3 -O 04-computation/planar_jc48_sep08_binomial_exactness.py

Both produce 628 identical bytes and pass 1,654 gates. Pins independently checked:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source | 6562 | `ebc2df4e64b7abf29f9a1caf939cbe29329c94dfc7361590d0e61c50873d15bb` |
| Frozen output and both replays | 628 | `e70e6daab742bfa00fc7cc26daaacc9fb62f20df90395c01a99fc30550618683` |
| Primary before promotion | 15819 | `f1c881eb13258493b9137c618f055c4fd8a4f642c7a7fb950291f231d2d232c1` |

No mathematical or producer-source repair was needed. The source/output stay frozen, including the source's historical pre-audit wording. This audit authorizes promotion of the stated theorem and its limited consumers, with the precise boundaries above.
