# Retain the old Phi6 identity, not its retracted extremality claim

The June30 reflection `07-reflections/the-covering-min-as-a-function-of-n-one-transition-at-7-razor-margin-1-over-n-squared-opus-20260630.md` already records Phi6(19)=7³. Session108 credits that identity. The reflection still displayed a covering-minimum claim later retracted by MISTAKE-087 and HYP-3730's status warning. A historical notice was inserted below its title without rewriting the old prose.

For a speed set S define M(S)=max_(0<=t<=1) min_(v in S)||vt||. Exact critical-time enumeration gives

* M({1,2,5,6,7,8})=2/13, first attained at3/13.
* M({1,2,3,4,5,42})=7/43, first attained at7/43.

Since2/13<7/43, the latter construction is not the minimum at seven runners. Both values remain above1/7; this is not a counterexample to the lonely runner conjecture. It is an inherited correction, independently replayed here, not a new bound.

Completeness of the finite check: every distance ||vt|| is piecewise affine, with its cusps at times k/(2v). An extremum of their lower envelope occurs at a cusp, an endpoint, or an intersection of two affine branches. Such intersections have t=k/(u+v) or k/|u-v|. The checker enumerates all these times in[0,1] using rational arithmetic and evaluates the original minimum directly. No floating-point grid or extrapolation is used.

`python3 check_covering_boundary.py` and `python3 -O check_covering_boundary.py` produce the same receipt. Only the arithmetic Phi6 identity is imported into the divisor-sum theorem; no old geometric optimality inference is revived.
