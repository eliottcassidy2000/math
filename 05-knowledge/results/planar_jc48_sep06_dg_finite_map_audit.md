# Independent audit: finite quartic and the prescribed relative primitive

**Status: PASS / PROVED in the stated scope.** Root independently read
the complete [proof](planar_jc48_sep06_dg_finite_map.md) and
[exact source](../../04-computation/planar_jc48_sep06_dg_finite_map.py).
The supplier surface was already independently audited. No classification
of actual Keller envelopes enters this argument.

The two affine charts of the relative projective quartic recover W
scheme-theoretically by solving for B and T respectively. No fibre
polynomial vanishes identically, because its middle coefficient is one.
Thus proper plus quasi-finite really proves finite; the affine degree
drop at T=0 loses no projective point. The quartic has length four on
every projective fibre. The three connecting Laurent classes are the
whole H1(O(-4)) basis, so the four displayed functions form a saturated
free module, not merely a generic vector-space basis. The exact sequence
also supplies a direct check of finite flatness.

I checked the multiplication table and its trace pairing, which has
determinant 16TB(1-4TB)². The two chart Jacobians give every ramification
component: D, Z and C, each with coefficient one. The asserted residue
degrees follow directly from the restrictions -x², -r² and -x²/2.
Thus both axes have generic length ledger 2·1+1·2=4, while the hyperbola
has 2·2=4. Its double pullback and the two double origin points agree
with the projective quartic and discriminant. The generic source is
ramified; this is correctly retained as a positive finite-map control
with a failed Keller condition.

The divisor calculation t=-r²(1+r²b) has exact order two on D and
the other component is precisely the disjoint affine line L. For the
unbounded theorem, restriction of any global g to the source is a
polynomial. The bracket equation forces g=-cx+h(t) in characteristic
zero. The pole -c/r cannot cancel against the regular polynomial in t.
The argument proves nonexistence for every global g, not just the
quadratic h sampled by the verifier. For (F(t),g), the product
-F'(t)g_x being a nonzero constant forces both factors to be units;
this reduces to the same contradiction. A finite such map is dominant
and generically separable, so its nonconstant nonzero Jacobian must
vanish somewhere in the affine-plane chart.

Finally a finite map cannot contract D; its image is the whole line with
first coordinate F(0). Pulling back that line's local uniformizer gives
e_D=2 ord_0(F-F(0)). This is the valuation of a target prime, not a
confusion between the ramification index and total degree. Its parity
really excludes the required index three within this coordinate class.

Normal and optimized runs both match the 41-gate frozen output exactly.

* Source SHA256: `e1f4474f580b9cd076f1c00ebdb34aa9451226ed8716a74e044bf5775e31315f`.
* Output SHA256: `f12e3a3f7d0b558fef894aab4d3addf3284b25ac2b0559e0014e08c4409c6b62`.

The result retains both distinctions needed by the consumer: a boundary
separator can supply a finite map with interior ramification, and a
globally exact two-form need not have a regular prescribed relative
primitive. No arbitrary pair of global coordinates is excluded.
