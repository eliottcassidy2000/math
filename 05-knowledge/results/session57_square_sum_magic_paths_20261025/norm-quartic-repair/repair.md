# Norm matching versus fourth-harmonic matching

Independent audit of Ajai Choudhry, arXiv:2508.07806v1, section 5.1, printed page 12.
Sources: https://arxiv.org/pdf/2508.07806v1 and https://arxiv.org/html/2508.07806v1#S5.SS1.
PDF SHA256: 964fdfc1046e0f34e652f97ce507e803bf50581e7a6042bb0c461d7a5bd3ad66.
Both HTML and rendered PDF were checked; the discrepancy is present in the PDF, not merely a conversion artifact. No full source paper is included in the public deliverable.

## Exact correction

Let A_i,B_j be the expressions in (5.4), over Q with m nonzero. Write P=p²+q², R=r²+s², and H(a,b)=a⁴−6a²b²+b⁴=Re((a+ib)⁴). Define D_cross=A3+B4−A4−B3 and D_diag=A3+B3−A4−B4.

Direct expansion gives

D_cross=(m²−1)(H(r,s)−m²H(p,q))/m²,
D_diag=(m²−1)(R²−m²P²)/m².

In fact A3+B3=(P−R)² and A4+B4=(m²P−R)²/m²; their subtraction immediately proves the second formula. Expanding the first formula verifies the coefficients independently. The printed (5.5) is −m²D_diag=0, not D_cross=0. Thus the premise preceding (5.5) should be diagonal matching, not cross matching. The auxiliary diagonal condition is a permissible ansatz, not a sufficient condition for all sixteen sums to be squares.

For (p,q,r,s,m)=(1,2,3,4,5), the four values (25A3,25A4,25B3,25B4)=(2500,6724,7500,3276). Both diagonal sums scaled by 25 equal 10000; the cross sums are 5776 and 14224. This is a counterexample to the displayed reduction, not to Theorem 3.

If both diagonal and cross equalities hold in characteristic other than 2, adding and subtracting them gives A3=A4 and B3=B4. Therefore imposing both would force repeated entries in each corresponding quadruple, rather than establish a nondegenerate four-by-four construction.

## Repair without refuting the construction

Use u=u1,v=u2,x=v1,y=v2, and set
D=x²−2xy−y², N=x²+2xy−y², t=N/D, m=t²,
p=ux+vy, q=uy−vx, r=t(ux−vy), s=t(uy+vx).
Assume ND≠0, so the original rational construction is defined.
Set K=8xy(x²−y²)(x²+y²)/D².
Then direct polynomial calculation gives

A3+B4=[K(u²−2uv−v²)]²,
A4+B3=[K(u²+2uv−v²)]².

The other fourteen square roots, in row order, are

2(pr−qs), 2(ps−qr), 2(pq−rs), 2(m²pq−rs)/m;
2(pr+qs), 2(ps+qr), 2(pq+rs), 2(m²pq+rs)/m;
p²−q²+r²−s², p²−q²−r²+s², P−R, [the first exceptional root];
(m²p²−m²q²+r²−s²)/m, (m²p²−m²q²−r²+s²)/m, [the second exceptional root], (m²P−R)/m.

Squaring each root verifies its corresponding A_i+B_j. These fourteen are generic identities before substitution; the two exceptional roots are identities after substitution. This is a direct repair of the rational construction underlying Theorem 3 on its nonsingular domain. It does not by itself claim verification of every cleared-denominator boundary specialization of the separately printed polynomial theorem.

The published numerical example is also independently checked: the root rows are (6692,8400,35592,42000), (18445,19131,39525,45381), (41692,42000,54408,58800), (52700,52944,63240,67056).

## Structural lesson

For complex z of modulus rho and argument theta, Re(z⁴)=rho⁴ cos(4theta), whereas |z|⁴=rho⁴. The diagonal condition records only squared norm; the actual cross condition retains a fourth angular harmonic. Thus a norm-balance argument cannot replace the fourth-harmonic condition. This is a mathematical connection, not a priority claim or a claim about the author's intent.

The certificate checks both exact identities, the control, all fourteen generic square roots, both repaired exceptional square roots, a bounded rational sample grid, and the published numerical example. Bounded checks are supplementary; the symbolic identities provide the proof on ND≠0.
