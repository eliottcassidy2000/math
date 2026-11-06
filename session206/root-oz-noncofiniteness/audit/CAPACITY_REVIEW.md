# Independent parameter-capacity review

Verdict: PASS for `session206-arithmetic-bank/PARAMETER_CAPACITY.md`.

The digit identity, uniform approximation on the inverse-bounded parameter range, residue-class counting, and Riemann-sum passage establish

A_s(X) = X log(R/(R-d))/d + O_s(sqrt(X))

for each fixed seed with R>d>0. The count is of inputs with multiplicity, not distinct outputs. Summation over a fixed finite bank is legitimate, but no uniform assertion for growing banks follows. The hypothetical d=0 case is separately A_s(X)=X/R+O_s(1), provided R>0. All physical seeds under discussion have d>0.

Independent derivation: for M=2^k and 0<=r<M,

Q(Mv+r)=Q(v)+Q(r)+v rev_k(r)+topbit_k(r)(v mod2).

Here rev_k is reversal of exactly k binary digits. Splitting the defining pair sum into low-low, high-high, and cross pairs proves this identity. Since 0<=Q(v)<=v, each residue-class count differs from X/(MR-d rev_k(r)) by O_s(1+X/M^2). Summing gives O_s(M+X/M); M comparable with sqrt(X) proves the same asymptotic without any equidistribution assumption.

The strict capacity inequality is a sufficient obstruction only. Equality is not excluded. Physical attainability is required when calling an arithmetic source a tournament seed.
