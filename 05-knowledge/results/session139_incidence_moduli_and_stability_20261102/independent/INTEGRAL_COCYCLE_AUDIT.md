# Integral four-face cocycle: independent audit

Status: VERIFIED / PROVED by the following algebra, with separate exhaustive
checks. This strengthens the recovered F2 boundary equation; it is not a
claim of literature priority or a Riemann-hypothesis consequence.

Let S be a binary tournament skew matrix with entries s_ij in {+1,-1}, and
interpret every expression as an alternating cochain on ordered distinct
vertices. Set

    t_ijk=s_ij s_jk s_ki,
    c_ijk=(s_ij+s_jk+s_ki+t_ijk)/4.

The latter is exactly +1 for the cyclic orientation i->j->k->i, -1 for its
reverse, and zero for a transitive triangle. Thus c is integral. Using the
usual alternating simplicial coboundary, put eta=delta c. The linear term
s_ij+s_jk+s_ki is delta s, so delta^2s=0 gives

    eta=delta t/4.

For ordered vertices 0,1,2,3 write a=s01,b=s02,c=s03,d=s12,e=s13,f=s23.
Then

    4 eta=abd-ace+bcf-def.

These four signed monomials are ±1 and their product is +1; the sum is
therefore -4,0,or4. In particular eta belongs to {-1,0,1}. Reduction of
signed c modulo2 is the indicator of a cyclic triangle. Hence eta modulo2
is the old odd-fourset indicator epsilon, and because eta has those three
integer values, eta^2=epsilon. The elementary four-vertex Pfaffian identity
therefore reads

    eta^2=(Pf(S[0,1,2,3])^2-1)/8.

All formulas were checked independently by direct oriented triangle tests,
not by assuming this Pfaffian identity.

Each t_ijk is invariant under switching S->DSD, since each diagonal sign
occurs twice. Consequently eta is switching-invariant and alternating under
relabeling. Global reversal S->-S sends eta->-eta. Most importantly,

    delta eta=delta^2 c=0 over the integers.

This is an oriented five-face conservation law modulo every M after reduction
of the integer lift. At odd M one may equivalently divide delta t by4 using
its inverse. At even M one must form the integral lift first: S mod2 is the
same all-one off-diagonal matrix for every tournament, but eta is not constant.
The theorem does not say that the unsigned support size is divisible by every
M. Dropping the alternating boundary signs gives a false identity and is an
explicit rejected mutation in the checker.

The support is exactly the old bad-fourset support. A nonzero witness face X
forces another support face in X union {v} for every outside v, giving support
at least n-3. An ordinal cyclic-triangle block with singleton blocks attains
this floor. No classification of all equality cocycles or all tournament
realizations follows from this counting argument.

## Weighted firewall

The binary square identity must not be applied to a weighted skew matrix K.
At K=0, eta=0 and Pf(K)=0, so the claimed weighted equality would be
0=-1/8. The polynomial eta(K) is the mean of the binary signed lift under
independent edge rounding, but squaring it is not the mean of eta^2. In
particular a variance/second-moment sidecar is needed for weighted support
probabilities; the integral conservation law alone does not prove weighted
five-vertex stability.

## Exact controls

`audit_integral_cocycle.py` exhausts all64 labelled four-vertex tournaments
and all1024 labelled five-vertex tournaments. It checks the signed triangle
formula, the cubic eta formula, the Pfaffian square identity, all1024 binary
switchings at order4, all1536 alternating relabelings, and all five-boundary
integer conservation laws. The signed four-face distribution is
`eta=-1:8, eta=0:48, eta=1:8`.

Checks modulo2 through33 are redundant controls; the proof over Z supplies
the all-modulus quantifier. Normal and -O receipts agree byte-for-byte. The
unsigned-five-boundary mutation is rejected by explicit RuntimeError in both
modes. The K=0 weighted hostile and mod2 information-loss hostile are retained.
