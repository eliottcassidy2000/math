#!/usr/bin/env python3
"""Primary exact certificate for THM-4397.

Compare Christopher D. Long's arXiv:2608.23777v1 rank-two Poisson map with
the earlier THM-2044 suspension.  The certificate constructs the source
Hamiltonian translation and the target symplectic quarter-turn explicitly.
"""

import sympy as sp


x, q, p, z, e, w = sp.symbols("x q p z e w")
VARS = (x, q, p, z)
s = x * q


def check(condition, label):
    if not condition:
        raise AssertionError(label)


def check_zero(value, label):
    check(sp.expand(value) == 0, label)


def bracket(f, h):
    """Canonical bracket with {p,x}={z,q}=1."""
    return sp.expand(
        sp.diff(f, p) * sp.diff(h, x)
        - sp.diff(f, x) * sp.diff(h, p)
        + sp.diff(f, z) * sp.diff(h, q)
        - sp.diff(f, q) * sp.diff(h, z)
    )


# Shared linear source coordinates.
R = x * (2 - 3 * s)
D0 = (1 + 3 * s) * p / 2 - 3 * q**2 * z
L = 3 * x**2 * p + (2 - 6 * s) * z


# THM-2044 gauge.
G = 252 * s**3 + 1008 * s**2 + 1379 * s + 659
g = -q**2 * G / 140
ell = sp.expand(L + g)
y_repo_e = q - x * e / 3
v_repo_e = x * y_repo_e
T_repo_e = sp.expand(
    y_repo_e
    + 3 * x * (1 + v_repo_e) ** 2 * e
    + 3 * x * y_repo_e**2 * (4 + 3 * v_repo_e)
)
S_repo_e = sp.expand(
    ((1 + v_repo_e) ** 3 * e
     + y_repo_e**2 * (1 + v_repo_e) * (4 + 3 * v_repo_e)) / 2
)
B_repo_e = (
    2 * e**4 * x**6 * (3 * s - 2)
    + e**3 * x**4 * (-90 * s**2 - 30 * s + 55)
    + e**2 * x**2 * (540 * s**3 + 720 * s**2 - 120 * s - 270)
    + e * (-1620 * s**4 - 3780 * s**3 - 1890 * s**2 + 810 * s + 540)
    + q**2 * (2430 * s**3 + 8100 * s**2 + 8640 * s + 2430)
)
H_repo_e = sp.expand(-e * B_repo_e / 1620)
T_repo = sp.expand(T_repo_e.subs(e, ell))
S_repo = sp.expand(S_repo_e.subs(e, ell))
D_repo = sp.expand(D0 + H_repo_e.subs(e, ell))
PHI_REPO = (R, T_repo, D_repo, S_repo)


# Long's arXiv:2608.23777v1 gauge.
a = 1 - 3 * s
beta = sp.expand(L - 9 * q**2)
y_paper = sp.expand(q - x * beta / 3)
v_paper = x * y_paper
R_paper = sp.expand(2 * x - 3 * x**2 * y_paper - x**3 * beta)
S_paper = sp.expand(
    y_paper
    + 3 * x * (1 + v_paper) ** 2 * beta
    + 3 * x * y_paper**2 * (4 + 3 * v_paper)
)
T_paper = sp.expand(
    -((1 + v_paper) ** 3 * beta
      + y_paper**2 * (1 + v_paper) * (4 + 3 * v_paper)) / 2
)
H_paper = sp.expand(
    y_paper**4 * (18 * v_paper**2 + 78 * v_paper + 125) / 20
    + 3 * beta * y_paper**2
      * (v_paper**3 + 5 * v_paper**2 + 10 * v_paper - 5) / 10
    - beta**2 * (9 * v_paper + 2) / 6
    - x**2 * beta**3 / 6
)
D_paper = sp.expand(D0 + H_paper)
PHI_PAPER = (R_paper, T_paper, D_paper, S_paper)


# The source Hamiltonian translation U_K and target quarter-turn tau.
k = (63 * w**2 + 318 * w + 601) / 840
K = sp.expand(q**3 * k.subs(w, s))
Kx = sp.diff(K, x)
Kq = sp.diff(K, q)
U_SUBS = {x: x, q: q, p: p + Kx, z: z + Kq}


def through_U(f):
    return sp.expand(f.subs(U_SUBS, simultaneous=True))


TAU_PHI_PAPER_AFTER_U = tuple(
    through_U(f) for f in (R_paper, S_paper, D_paper, -T_paper)
)


# 1. The one-variable equation identifies the generator without guessing.
h = (601 - 1379 * w - 1008 * w**2 - 252 * w**3) / 140
check_zero(
    w * (2 - 3 * w) * sp.diff(k, w) + 6 * (1 - 3 * w) * k - h,
    "Hamiltonian coefficient ODE",
)
check_zero(
    through_U(L) - L - (g + 9 * q**2),
    "linear coordinate shift",
)


# 2. U_K is a polynomial symplectomorphism, with inverse U_{-K}.
pU = p + Kx
zU = z + Kq
check_zero(bracket(pU, x) - 1, "{pU,x}=1")
check_zero(bracket(zU, q) - 1, "{zU,q}=1")
for label, f, h0 in (
    ("{x,q}", x, q),
    ("{x,zU}", x, zU),
    ("{q,pU}", q, pU),
    ("{pU,zU}", pU, zU),
):
    check_zero(bracket(f, h0), label)
check_zero(through_U(p - Kx) - p, "U inverse p")
check_zero(through_U(z - Kq) - z, "U inverse z")
# The identical mixed-partial identity is the only nontrivial generator
# relation for the exact Weyl lift p |-> p+K_x, z |-> z+K_q.
check_zero(sp.diff(Kq, x) - sp.diff(Kx, q), "exact Weyl lift mixed commutator")


# 3. The adapted coordinates agree after U_K.
check_zero(through_U(beta) - ell, "beta after U equals ell")
check_zero(through_U(y_paper) - y_repo_e.subs(e, ell), "adapted y equality")
check_zero(through_U(D_paper) - D_repo, "adapted D equality")


# 4. Full point-map identity: Phi_repo = tau o Phi_paper o U_K, where
#    tau(r,t,d,s)=(r,s,d,-t).
for index, (left, right) in enumerate(zip(PHI_REPO, TAU_PHI_PAPER_AFTER_U)):
    check_zero(left - right, f"full map coordinate {index}")


# The target quarter-turn preserves dr^dD+dt^dS.  Its matrix is symplectic.
Omega = sp.Matrix(
    [[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]
)
Jtau = sp.Matrix(
    [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0]]
)
check(Jtau.T * Omega * Jtau == Omega, "target quarter-turn symplectic")


# 5. The published and repository exact fibres transport point by point.
paper_fibre = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(1, 24), sp.Rational(-1, 8)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
    (sp.Rational(-1), sp.Rational(-2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
)
repo_fibre = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(1, 24), sp.Rational(-1, 8)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(224839, 90720), sp.Rational(-173417, 60480)),
    (sp.Rational(-1), sp.Rational(-2, 3), sp.Rational(224839, 90720), sp.Rational(-173417, 60480)),
)
transported = []
for xx, qq, pp, zz in paper_fibre:
    base = {x: xx, q: qq}
    transported.append(
        (xx, qq, sp.factor(pp - Kx.subs(base)), sp.factor(zz - Kq.subs(base)))
    )
check(tuple(transported) == repo_fibre, "inverse-U fibre transport")
for point in repo_fibre:
    image = tuple(sp.factor(f.subs(dict(zip(VARS, point)))) for f in PHI_REPO)
    check(image == (0, 0, 0, sp.Rational(-1, 8)), "repo fibre image")


# 6. Complexity is gauge-dependent, not a new noninvertibility mechanism.
def terms_and_degree(poly):
    P = sp.Poly(poly, *VARS)
    return len(P.terms()), P.total_degree()


paper_profile = tuple(terms_and_degree(f) for f in PHI_PAPER)
repo_profile = tuple(terms_and_degree(f) for f in PHI_REPO)
check(paper_profile == ((2, 3), (47, 15), (139, 23), (22, 11)), "paper profile")
check(repo_profile == ((2, 3), (35, 21), (246, 48), (78, 30)), "repo profile")


print("THM-4397 RANK-TWO POISSON GAUGE EQUIVALENCE")
print("source_generator=K=q^3*(63*(xq)^2+318*xq+601)/840")
print("source_map=U_K:(x,q,p,z)->(x,q,p+K_x,z+K_q)")
print("target_map=tau:(r,t,d,s)->(r,s,d,-t)")
print("identity=Phi_repo=tau o Phi_paper o U_K")
print("K_x=", sp.factor(Kx))
print("K_q=", sp.factor(Kq))
print("paper_profile_(terms,degree)=", paper_profile)
print("repo_profile_(terms,degree)=", repo_profile)
print("fibre_transport=paper fibre under U_K^{-1} equals THM-2044 fibre")
print("scope=polynomial symplectic right-equivalence over Q")
print("weyl_gauge=U_K and tau lift to exact A2 Weyl automorphisms")
print("consequence=finite A2 quantizability and termination are gauge-invariant")
print("NO_CLAIM=DC2_or_JC2_or_existence_of_an_A2_quantization")
print("RESULT=PASS")
