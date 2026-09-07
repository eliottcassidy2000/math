# Independent audit of the completed-carrier finite response

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS.**
The complete theorem in [the primary note](planar_jc48_sep07_completed_response.md)
is accepted, including necessity, sufficiency, the earlier-prefix obstruction,
and the actual completed-flow interpretation. No mathematical correction is
requested. This audit does not assert later compatibility or termination.

## 1. The all-polynomial cutoff and actual source are complete

I read the full producer and the incoming valuation-eleven response theorem.
Its literal A5,C5 agree with the incoming arbitrary-Xi background. Deriving
P(A,C) through row five removes Xi and gives exactly the displayed H5;
Xi has not been specialized away in the coordinate corrections.

The actual valuation of R is n=min(a+2b), because the distinct monomials
at a fixed valuation have distinct leading powers x^b. The Hamiltonian
p²Delta R starts at t^(n+5)r_n(x). Since A0_x=x/2, its action on A starts
at (n+5)x r_n(x)t^(n+4)/2, which is nonzero in characteristic zero.
Thus preserving rows through nine is equivalent to n>=6. Every carrier
of valuation at least twelve is invisible through row fifteen; precisely
the thirty columns of valuations6..11 suffice, without a hidden degree
cutoff or generic-coefficient assumption.

The source starts one row later because G0=0. Its leading operator is
((x²/2)+3)r'_n-(n+5)x r_n. Its highest coefficient is
(b/2-n-5)[x^b]r_n, which is nonzero for b<=floor(n/2).
This proves injectivity at each successive source valuation. The six
valuation-eleven columns are exactly the source-invisible kernel, and
their two displayed terminal coordinate coefficients follow directly
from A0_x and C0_x. Kernel dimensions are therefore tied to actual
coordinates, not inferred from a source projection alone.

Only rows zero through five enter at the horizon, since a row of valuation
six bracketed with a Hamiltonian of valuation eleven begins at sixteen.
I independently rederived the cusp-coordinate operator from
{p,y}=-Delta/p and delta u=p³y(10R+2pR_p+3yR_y). Its signs agree with
the literal source action. The source checks those two routes on every
column, including the full chain derivative of the original polynomial P.

## 2. Exact complete image, supplier, and obstruction

The primary source retains every coefficient of the original response,
the actual projected-depth annihilator, and both caps. Its bracket loop
indexed by n checks power n-1: n=10..15 is precisely bracket powers9..14.
The complete91 depth rows are inherited from the audited actual P2/P3
annihilator, rather than replaced by a degree bound. All checks are
polynomial identities with Xi retained.

I also independently reconstructed the source-image matrix using
[a standalone Fraction dictionary engine](../../04-computation/planar_jc48_sep07_completed_response_audit.py).
It imports no mathematical producer and uses no SymPy. It freshly
generates the carrier and target monomials, constructs the cusp operator
from H5, and substitutes p=t(1+x²t),y=xtp by binomial convolution.
Its50 nonzero coefficient rows correspond to the primary51-row matrix
after removal of a zero row. Independent rational Gaussian elimination
gives rank36 and exactly the same eleven free columns. The carrier
source map has rank24, so its kernel is exactly the six terminal columns.

The primary output is used only as proposed certificate data: all five
literal source columns and their thirty-entry carrier lifts are
substituted into the independently assembled matrix. Each vanishes,
and the five declared source coordinates give the identity matrix.
Together with rank36 this proves completeness of the source dimension5
and coordinate fibre dimension6. The primary verifier independently
checks that the incoming seven equations and the five displayed new
equations have joint rank12 and annihilate this entire basis. Hence the
written conditions are an iff, not merely necessary equations.

The Fraction engine separately verifies the six-term Rstar and all eight
terms of Tstar, the failure of the earlier compensated packet, and the
exact nonzero leading minor1751787/413887398592. Vanishing of coordinate
row ten forces every valuation-six coefficient to vanish; in particular
the two displayed coefficients then force rho=k=0. This correctly proves
that every nonzero high payer in the declared image must change row ten.
It does not deny the older finite Hamiltonian realization, whose generator
is outside the universal carrier.

The constant rational matrix and its explicit inverse lifts contain no
denominator in Xi. The source/image statements therefore hold for every
Xi in every characteristic-zero field, not only generically. Multiplying
the supplier by the actual scalar j(s) introduces no division by that
function, so the high coefficient is cancelled also at its zeros, where
the selected scalar-time map is the identity.

## 3. Completed flow and preserved scope

The literal derivation of a Hamiltonian of valuation at least eleven
raises t-valuation by at least ten. Its second iterates on the background
coordinates begin at twenty. The nonlinear finite-coordinate source
terms begin at twenty and the quadratic bracket terms at nineteen.
Thus arbitrary scalar time has the displayed linear response through the
required horizon, with no infinitesimal-only inference.

The completed source automorphism is supplied by the existing all-order
universal-carrier theorem. Finite depth representatives give the required
rows; no assertion identifies arbitrary Delta-adic elements with rational
Laurent expressions or claims that the displayed truncated background
already solves the later Keller equations. The result preserves complete
finite compatibility and the completed source structure, while allowing
the prefix motion from row ten. Neither rationality of a nonzero scalar
time nor polynomial termination follows. In particular the new genus
obstruction for y-linear invariants and this positive finite response are
compatible statements on different parts of the proof graph.

## 4. Exact replays and pins

I replayed the primary source normally and under python -O. Both outputs
equal the frozen3152-byte report, passing4025 always-active gates. The
independent Fraction audit passes20 gates and its normal/optimized outputs
are byte-identical240-byte reports.

```
python3 -B 04-computation/planar_jc48_sep07_completed_response.py
python3 -B -O 04-computation/planar_jc48_sep07_completed_response.py
python3 -B 04-computation/planar_jc48_sep07_completed_response_audit.py
python3 -B -O 04-computation/planar_jc48_sep07_completed_response_audit.py
```

```
primary source 669567f67d3febc2c05089b4215a656d7bbfb1315527df319019b3b1c71d0b47
primary output 2427ecb11e54c8a9e386e13c061f7914a0616af5867f8907f5049751885d73da
audit source 5a0a3495f8a14aeecaa46b97362b67bdedb32e1ea20753cc20859d4dbfd7b682
audit output 378fd8139f6603cef904bc1e93006caf8af307caaba51408776b44de59ca7355
audit RREF 6641791bc96ae9a2dce07bb47e98f400060d66861bc1078ba41514152dfe4d50
```

All primary source/output pins match those announced at freeze. The
three independent audit files are frozen for parent integration.
