# Independent audit of the shared-boundary-root jet theorem

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE + REPLAY AUDIT PASS.**
I read the complete
[primary proof](planar_jc48_sep08_shared_roots.md) and its standalone
source. The three conclusions are accepted with their declared scope:
local regularity at simple restriction roots or a nonzero normal
N-derivative; a finite-point logarithmic obstruction from any of the
three indicated first jets; and exclusion of rational mates when every
complete octic multiplicity is at most two, for arbitrary M. No blanket
shared-root exclusion, infinity version of the finite logarithmic gate,
or global-regularity conclusion for a rational mate is inferred.

## 1. The actual surface and relative differential

The input is the fixed surface W=(P1×P1) minus S, not an assumed Keller
finite envelope. I recovered the full chart and divisor calculation in
[the DG filtration](planar_jc48_sep08_dg_quadratic.md) and
[the quartic differential theorem](planar_jc48_sep08_quartic_differential.md).
The numerator sections are in O(4,2) and O(2,1); their restrictions have
boundary degrees eight and four. The complete generic equation is

    E=N²+s³M−cs⁴=0.

The relative form is a unit times s²dw/E_s or s²ds/E_w at a finite
boundary point. At the unique intersection of S with the closure of D,
its canonical numerator has an extra w² factor. In particular this
factor can remove poles, and cannot be replaced by a unit in the
infinity calculation. All local coordinate changes in the proof have
invertible Jacobian; their extra factors do not introduce poles.

## 2. Independent local regularity check

For m=1, replacing the boundary coordinate by v=N(s,w) is invertible.
Since the point is shared, M has no constant term. With weights
wt(s)=1 and wt(v)=2, the full lowest equation is
v²+(alpha−c)s⁴. Terms s³v, s⁵, and higher have larger weight.
Outside one exceptional c-value its two leading roots are nonzero and
simple. The two actual branches have v of order two in s and E_v of
order two; s²ds/E_v is regular. These are the two Weierstrass branches,
not roots of an unrelated truncated equation.

For N_s a unit, the implicit centre N=0 is s=psi(w), of order m.
Writing ell for the order of M(psi(w),w)−c psi(w) gives
1<=ell<=m generically: terms below m are independent of c, and at
order m the coefficient has the nonzero slope supplied by psi.
The actual displacement from that centre has order (3m+ell)/2>m.
It therefore does not change the leading value of M−cs, whose order
is ell. N_s remains a unit. In E_s, the term 2N N_s has order
(3m+ell)/2, strictly below both 2m+ell and 3m, the possible orders
of the remaining terms. Consequently the relative form has order
(m−ell)/2 before the possible quadratic normalization. This is
nonnegative; normalization adds the nonnegative order of dw. The
Weierstrass degree in s is two, so no further boundary branches are
omitted.

This verifies the unbounded local assertion, including arbitrarily
high tangential M order or identically zero M|S. It also verifies why
(m,n) alone is insufficient: the two actual global m=2,n=1 controls
have distinct normal derivatives and fall into different cases.

## 3. The finite logarithmic gate and the full small-root consumer

At dN=0 and M=0, I independently substituted s=wZ in the full
quadratic and linear jets. The degree-four tangent polynomial is

    P_c(Z)=(a+bZ+eZ²)²+dZ³+(f−c)Z⁴.

For a nonzero, its constant term is a², and the equation for a repeated
nonzero root is ZP0'−4P0=0. That polynomial has constant term −4a²,
so only finitely many critical roots and c-values occur; the leading
quartic coefficient is nonzero generically. There are four simple
nonzero roots. For a=0,b nonzero, factoring Z² leaves a quadratic
with nonzero constant b² and the same generic-simplicity argument.
For a=b=0,d nonzero, factoring Z³ leaves a linear polynomial with a
nonzero generic root. No branch at Z=0 needs to be classified for the
necessary obstruction: one simple nonzero root already suffices.

The implicit function theorem gives an actual smooth branch
s=wZ(w). Along it, E_s has order three with nonzero leading
coefficient P_c'(Z0), while s²dw has order two. Its residue is
Z0²/P_c'(Z0), times the nonzero canonical finite-chart unit. This
is nonzero. A rational primitive cannot produce that simple pole.
The three coefficients in the primary are exactly a,b,d; the stated
necessary disjunction is therefore correct.

At infinity the extra w² factor removes these logarithmic poles.
Thus the finite-point obstruction cannot be transported there merely
from a shared root or an identical tangent polynomial.

For the global multiplicity-at-most-two theorem, the branch exhaustion
is paid: m=1 is regular; m=2 with N_s nonzero is regular; m=2 with
N_s=0 has a nonzero, hence precisely the four generic simple nonzero
tangent branches. The latter exhaust Weierstrass degree four. Nonshared
roots of these multiplicities are regular by the prior local theorem.
All remaining possible finite poles are simple, and every infinity
branch is regular.

Choose c outside the finite critical-value set and, if F|D is constant,
its one value. Every compact normalized generic component meets W away
from D; neither S nor the closure of D is a generic component. The
relative form is nonzero there and regular in W. A proposed rational
mate restricts meromorphically to each component for generic c. Its
derivative cannot have a simple pole, nor can a pole of that function
have a derivative of pole order at most one. Therefore its restriction
would have no poles anywhere and would be constant on the compact
component, contradicting the nonzero form. This proof neither assumes
irreducible generic fibres nor needs a transverse generic point on D.

## 4. Independent hostile and stopping-boundary checks

I independently derived and checked the nonzero-M rational-mate control
before its inclusion in the producer:

    h=x²+x⁴t=−b, H=h², L=h, F=h⁴+h,
    G=1/[3x³(4h³+1)].

Its Jacobian is exactly one. Its actual global numerator sections are
N=x⁴z² and M=x²z; boundary restrictions are x⁸ and x⁴, so both
sections are nonzero and share only the finite zero. Their first-normal
and tangential orders are (m,j,n)=(8,6,4), inside the surviving class.
The complete other chart gives F=b⁴−b and
G=r³/[3(1−4b³)]. Thus F|D is nonconstant, while G has real rational
poles and is not a regular mate.

For generic c, each root alpha of alpha⁴+alpha=c gives one compact
component h=alpha with parameter x and

    t=(alpha−x²)/x⁴,
    G=x^(−3)/[3(4alpha³+1)].

Its pole degree is exactly three, entirely at x=0. At D, the parameter
r=1/x gives local degree three. This is a same-component equality
control, not a sum over four components. The common-root hypothesis
cannot be removed even with nonzero M and nonconstant F|D.

I also checked the root-supplied stopping object: with beta nonzero,
alpha=−4beta³/27 and
N=beta x²z−beta x⁴+alpha x²z², M=−1, the finite root has order
six and is M-unit, while infinity is a shared order-two root. The full
infinity numerators are alpha r²+beta r²q−beta q² and r²q.
The actual restrictions H|D=−beta,L|D=0 are constant. This correctly
falls outside the all-small-multiplicity theorem, and the primary makes
no rational-mate existence or nonexistence claim for it.

## 5. Bounded antecedent and correction recovery

A targeted current-canon and MISTAKES search, supplemented by a tracked
HEAD search to cover sparse omissions, excluded this wave's own notes.
It found no earlier exact statement about a common zero of these two DG
boundary sections. This is a bounded recovery result, not a literature
priority claim.

The primitive-degree versus local-multiplicity mechanism is explicitly
inherited: [THM-2071](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
uses map degree d−h versus local multiplicity d−1 for a rational
primitive of 1/U. The same argument is recovered in
[THM-2214 §7.3](../../01-canon/theorems/THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten.md)
and [THM-2723 §4](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md).
The current same-component degree-three hostile pays that classical
inequality exactly; the new information is the actual DG boundary
valuation and jet calculation.

[THM-2129](../../01-canon/theorems/THM-2129-quartic-faber-three-coefficient-boundary-classification.md)
uses the phrase common boundary zero for a triple of Faber observables,
not for N|S and M|S. It is not the same predicate.
[THM-2778](../../01-canon/theorems/THM-2778-all-degree-complete-chosen-sheet-split-exact-prefix-closure.md)
retains its complete chosen-sheet split polynomial/Faber entry, as
active guardrail 61 requires. Neither removes the hypotheses of this
fixed-surface theorem. MISTAKE-248 repairs the false independence of two
quartic pole congruences; no independent pair of those divisibility
conditions is used here. The close polar-family precedents
[THM-3598](../../01-canon/theorems/THM-3598-danielewski-rational-exact-polar-graph-family-and-classification.md)
and [THM-3755](../../01-canon/theorems/THM-3755-composite-monomial-generic-fibre-residue-obstruction.md)
already distinguish rational exactness from regular mates, on their own
Danielewski and monomial objects. No unproved map to those objects is
assumed.

## 6. Source and frozen replay

The source has an explicit bounded control universe: normal-unit
m=1..8 with tangential orders n=1..9 or infinity; the three nonzero
tangent regimes and the first unsupported regime; the actual low-order
global controls; and the two displayed rational/constant-D hostiles.
The proof of arbitrary higher jets comes from the valuation arguments,
not from extrapolating this finite bank. Every gate is always active;
there are no assertion-only proof checks.

Independent normal and optimized runs reproduce every frozen output
byte: **454 gates**, **418 bytes**.

```bash
python3 04-computation/planar_jc48_sep08_shared_roots.py
python3 -O 04-computation/planar_jc48_sep08_shared_roots.py
```

* Source SHA256, 6,255 bytes:
  `695b0da5d8fb5dec95b85d0d32086432ac6aaa8cb1fd8bc00fa2eef27f4da2ad`.
* Frozen output and both independent replays:
  `d156179f74f97998556c7c60f6c2f27166543abe36321cf229a5a38c6a24d085`.

No mathematical or source correction is requested. The primary may be
promoted with exactly its three statements and the surviving finite
shared-root jet class m>=3,j>=2,n>=2. That class contains the actual
rational-mate hostile and is not excluded by this bundle.
