#!/usr/bin/env python3
"""Exact coefficient-height anatomy and JC(2) rational-search controls.

This companion has four deliberately separate scopes.

1. VERIFIED-EXACT: expand THM-1300, put it in tangent-to-identity form, and
   minimize its rational coefficient height over the residual diagonal gauge.
2. VERIFIED-EXACT: in a seven-parameter torus-equivariant factored support,
   eliminate the Keller equations and recover only the two-parameter diagonal
   orbit of THM-1300 on the nondegenerate branch.
3. VERIFIED-EXACT / FINITE-EXACT: run small planar Gauss-chart controls and
   print exact rational-box universe sizes.
4. VERIFIED-EXACT POINT INCIDENCE ONLY: check the 30 rational points currently
   stored for ICARM curve #273.  Independence and the database's rank lower
   bound are NOT re-proved here.

Reproduce with both

    python3 04-computation/jc2_small_rational_height_and_elliptic_slice_prior.py
    python3 -O 04-computation/jc2_small_rational_height_and_elliptic_slice_prior.py
"""

from fractions import Fraction
from itertools import product
from math import gcd, isqrt, log10

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def qheight(value):
    """Multiplicative height of one rational number, with H(0)=1."""
    value = Fraction(value)
    return max(abs(value.numerator), value.denominator)


def as_fraction(value):
    value = sp.Rational(value)
    return Fraction(int(value.p), int(value.q))


def coefficient_profile(expressions, variables):
    coeffs = []
    term_counts = []
    degrees = []
    for expression in expressions:
        poly = sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
        cs = [as_fraction(c) for c in poly.coeffs()]
        coeffs.extend(cs)
        term_counts.append(len(cs))
        degrees.append(poly.total_degree())
    return {
        "terms": term_counts,
        "degrees": degrees,
        "height": max(qheight(c) for c in coeffs),
        "max_denominator": max(c.denominator for c in coeffs),
        "coefficients": coeffs,
    }


def rational_box_size(bound):
    # #{a/b in Q: gcd(a,b)=1, b>0, max(|a|,b)<=B} = 4 sum phi - 1.
    phi_sum = 0
    for n in range(1, bound + 1):
        phi_sum += sum(gcd(a, n) == 1 for a in range(1, n + 1))
    return 4 * phi_sum - 1


def exact_decimal_digits_of_power(base, exponent):
    estimate = int(exponent * log10(base)) + 1
    value = base**exponent
    while value >= 10**estimate:
        estimate += 1
    while value < 10 ** (estimate - 1):
        estimate -= 1
    return estimate


print("JC2 SMALL-RATIONAL HEIGHT AND ELLIPTIC-SLICE PRIOR")
print("scope: exact algebra/finite controls; no JC(2) proof or counterexample")


# ---------------------------------------------------------------------------
# I. THM-1300, natural normalization, and the exact diagonal height minimum.
# ---------------------------------------------------------------------------

x, y, z = sp.symbols("x y z")
u = 1 + x * y
F1 = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)
F = sp.Matrix([F1, F2, F3])
det_F = sp.expand(F.jacobian([x, y, z]).det())
require(det_F == -2, "THM-1300 determinant mismatch")

original_profile = coefficient_profile(F, (x, y, z))
G = sp.Matrix([F3 / 2, F2, F1])
require(sp.expand(G.jacobian([x, y, z]).det()) == 1, "normalized determinant")
require(G.jacobian([x, y, z]).subs({x: 0, y: 0, z: 0}) == sp.eye(3), "normalized jet")
natural_profile = coefficient_profile(G, (x, y, z))

ss, tt = sp.symbols("s t", nonzero=True)
sub_diag = {x: x, y: ss * y, z: tt * z}
G_diag = sp.Matrix(
    [
        sp.expand(G[0].subs(sub_diag)),
        sp.expand(G[1].subs(sub_diag) / ss),
        sp.expand(G[2].subs(sub_diag) / tt),
    ]
)
G_diag_expected = sp.Matrix(
    [
        x - sp.Rational(3, 2) * ss * x**2 * y - sp.Rational(1, 2) * tt * x**3 * z,
        y
        + 3 * tt / ss * x * z
        + 6 * tt * x**2 * y * z
        + 3 * ss * tt * x**3 * y**2 * z
        + 12 * ss * x * y**2
        + 9 * ss**2 * x**2 * y**3,
        z
        + 3 * ss * x * y * z
        + 3 * ss**2 * x**2 * y**2 * z
        + ss**3 * x**3 * y**3 * z
        + 4 * ss**2 / tt * y**2
        + 7 * ss**3 / tt * x * y**3
        + 3 * ss**4 / tt * x**2 * y**4,
    ]
)
require(all(sp.cancel(a - b) == 0 for a, b in zip(G_diag, G_diag_expected)), "diagonal formula")
require(sp.cancel(G_diag.jacobian([x, y, z]).det()) == 1, "diagonal determinant")


def diagonal_coefficients(s_value, t_value):
    s_value, t_value = Fraction(s_value), Fraction(t_value)
    return [
        Fraction(1),
        -Fraction(3, 2) * s_value,
        -Fraction(1, 2) * t_value,
        Fraction(1),
        3 * t_value / s_value,
        6 * t_value,
        3 * s_value * t_value,
        12 * s_value,
        9 * s_value**2,
        Fraction(1),
        3 * s_value,
        3 * s_value**2,
        s_value**3,
        4 * s_value**2 / t_value,
        7 * s_value**3 / t_value,
        3 * s_value**4 / t_value,
    ]


def complete_height_gauges(bound):
    # H(s^3)<=B bounds the reduced numerator/denominator of s by cubert(B).
    s_candidates = set()
    cube_cap = 0
    while (cube_cap + 1) ** 3 <= bound:
        cube_cap += 1
    for den in range(1, cube_cap + 1):
        for num in range(-cube_cap, cube_cap + 1):
            if num and gcd(abs(num), den) == 1:
                s_candidates.add(Fraction(num, den))

    # If t=m/n is reduced and H(6t)<=B, then |m|<=B and n<=6B.
    t_candidates = set()
    for den in range(1, 6 * bound + 1):
        for num in range(-bound, bound + 1):
            if num and gcd(abs(num), den) == 1:
                t_candidates.add(Fraction(num, den))

    answers = []
    for s_value in sorted(s_candidates):
        for t_value in sorted(t_candidates):
            cs = diagonal_coefficients(s_value, t_value)
            if max(qheight(c) for c in cs) <= bound:
                answers.append((s_value, t_value))
    return answers


height8 = complete_height_gauges(8)
height9 = complete_height_gauges(9)
expected_height9 = [
    (s_value, t_value)
    for s_value in (Fraction(-1, 2), Fraction(1, 2))
    for t_value in (
        Fraction(-3, 4),
        Fraction(-1, 2),
        Fraction(-1, 4),
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(3, 4),
    )
]
require(height8 == [], "height-eight diagonal hostile")
require(sorted(height9) == sorted(expected_height9), "height-nine gauge classification")

optimal_s = Fraction(1, 2)
optimal_t = Fraction(1, 2)
optimal_map = G_diag.subs({ss: sp.Rational(1, 2), tt: sp.Rational(1, 2)})
optimal_profile = coefficient_profile(optimal_map, (x, y, z))
require(optimal_profile["height"] == 9, "optimal coefficient height")

collision_sources = [
    (Fraction(0), Fraction(0), Fraction(-1, 2)),
    (Fraction(1), Fraction(-3), Fraction(13)),
    (Fraction(-1), Fraction(3), Fraction(13)),
]
collision_target = (Fraction(0), Fraction(0), Fraction(-1, 2))
for point in collision_sources:
    image = tuple(
        as_fraction(expression.subs({x: point[0], y: point[1], z: point[2]}))
        for expression in optimal_map
    )
    require(image == collision_target, "optimal-gauge collision")

print("\n[I] THM-1300 coefficient profiles")
print(
    "original expanded: terms=%s degrees=%s H=%d maxden=%d"
    % (
        original_profile["terms"],
        original_profile["degrees"],
        original_profile["height"],
        original_profile["max_denominator"],
    )
)
print(
    "tangent identity G=(F3/2,F2,F1): terms=%s degrees=%s H=%d maxden=%d"
    % (
        natural_profile["terms"],
        natural_profile["degrees"],
        natural_profile["height"],
        natural_profile["max_denominator"],
    )
)
print("diagonal characters collapse D=diag(a,b,c) to s=ab, t=a^2 c")
print("H<=8 gauges:", height8)
print("H<=9 gauges:", [(str(a), str(b)) for a, b in height9])
print("proved diagonal-orbit minimum: H=9")
print(
    "secondary representative s=t=1/2: terms=%s H=%d maxden=%d"
    % (optimal_profile["terms"], optimal_profile["height"], optimal_profile["max_denominator"])
)
print("secondary collision: (0,0,-1/2),(+/-1,-/+3,13) -> (0,0,-1/2)")


# ---------------------------------------------------------------------------
# II. The seven-parameter equivariant support is only the diagonal orbit.
# ---------------------------------------------------------------------------

p, q, r, aa, bb, cc, dd = sp.symbols("p q r A B C D")
u_general = 1 + p * x * y
w_general = q + r * x * y
G_general = sp.Matrix(
    [
        x + cc * x**2 * y + dd * x**3 * z,
        y + aa * x * u_general**2 * z + bb * x * y**2 * w_general,
        u_general**3 * z + y**2 * u_general * w_general,
    ]
)
det_error = sp.Poly(sp.expand(G_general.jacobian([x, y, z]).det() - 1), x, y, z)
det_equations = [coefficient for _, coefficient in det_error.terms()]
require(len(det_equations) == 12, "unexpected determinant coefficient count")

groebner_basis = sp.groebner(det_equations, dd, cc, bb, aa, r, q, p, order="lex")
certificates = [
    p**4 * (aa - bb),
    p**4 * (aa * q - 12 * p),
    p**3 * (aa * r - 9 * p**2),
    -2 * aa * q + 2 * bb * q + 2 * cc + 3 * p,
    -2 * aa**2 * q + 2 * aa * bb * q + aa * p + 6 * dd,
]
for certificate in certificates:
    require(groebner_basis.reduce(certificate)[1] == 0, "Groebner certificate")

solution = {
    bb: aa,
    q: 12 * p / aa,
    r: 9 * p**2 / aa,
    cc: -sp.Rational(3, 2) * p,
    dd: -aa * p / 6,
}
G_solution = G_general.subs(solution).applyfunc(sp.cancel)
require(sp.cancel(G_solution.jacobian([x, y, z]).det()) == 1, "family sufficiency")
orbit_substitution = {ss: p, tt: aa * p / 3}
G_orbit = G_diag_expected.subs(orbit_substitution).applyfunc(sp.cancel)
require(all(sp.cancel(a - b) == 0 for a, b in zip(G_solution, G_orbit)), "orbit identity")

general_collision_sources = [
    (0, 0, -3 / (4 * aa * p)),
    (1, -3 / (2 * p), 39 / (2 * aa * p)),
    (-1, 3 / (2 * p), 39 / (2 * aa * p)),
]
general_collision_target = (0, 0, -3 / (4 * aa * p))
for point in general_collision_sources:
    image = tuple(
        sp.cancel(expression.subs({x: point[0], y: point[1], z: point[2]}))
        for expression in G_solution
    )
    require(all(sp.cancel(a - b) == 0 for a, b in zip(image, general_collision_target)), "family collision")

print("\n[II] seven-parameter factored-support Keller elimination")
print("universe: parameters (p,q,r,A,B,C,D), branch p!=0")
print("det(DG)-1 has 12 coefficient equations:")
for index, ((monomial, _), equation) in enumerate(zip(det_error.terms(), det_equations), start=1):
    print("  e%02d monomial=%s: %s" % (index, monomial, sp.factor(equation)))
print("lex Groebner basis size:", len(groebner_basis.polys))
print("selected ideal certificates:")
for certificate in certificates:
    print(" ", sp.factor(certificate))
print("on p!=0 these force A!=0 and")
print("  B=A, q=12p/A, r=9p^2/A, C=-3p/2, D=-Ap/6")
print("sufficiency: det(DG)=1 identically")
print("orbit identity: s=p, t=Ap/3; no nongauge modulus remains in this support")


# ---------------------------------------------------------------------------
# III. Height-box sizes and planar Gauss-chart controls.
# ---------------------------------------------------------------------------

print("\n[III] exact rational universes")
for bound in (1, 2, 4, 9, 12):
    size = rational_box_size(bound)
    print("B=%d: |R_B|=%d" % (bound, size))
size9 = rational_box_size(9)
for parameter_count in (2, 3, 4, 5, 7, 13):
    print("B=9 k=%d: %d" % (parameter_count, size9**parameter_count))

monomial_count_72 = (72 + 1) * (72 + 2) // 2
monomial_count_108 = (108 + 1) * (108 + 2) // 2
dense_slots_72_108 = (monomial_count_72 - 3) + (monomial_count_108 - 3)
require(dense_slots_72_108 == 8690, "dense slot arithmetic")
print("typed dense (72,108) tangent-identity slots:", dense_slots_72_108)
for bound in (1, 2, 9):
    size = rational_box_size(bound)
    print(
        "dense box B=%d has %d decimal digits"
        % (bound, exact_decimal_digits_of_power(size, dense_slots_72_108))
    )

# A finite hostile: every affine U,V,W coefficient lies in {-1,0,1}.
u0, ux, uy, v0, vx, vy, w0, wx, wy = sp.symbols("u0 ux uy v0 vx vy w0 wx wy")
gauss_parameters = (u0, ux, uy, v0, vx, vy, w0, wx, wy)
U_aff = u0 + ux * x + uy * y
V_aff = v0 + vx * x + vy * y
W_aff = w0 + wx * x + wy * y
a_aff = sp.expand(1 + U_aff * V_aff)
b_aff = sp.expand(V_aff + W_aff * (1 + U_aff * V_aff))
c_aff = U_aff
d_aff = sp.expand(1 + U_aff * W_aff)
curl1_aff = sp.Poly(sp.diff(a_aff, y) - sp.diff(b_aff, x), x, y)
curl2_aff = sp.Poly(sp.diff(c_aff, y) - sp.diff(d_aff, x), x, y)
affine_curl_equations = list(dict.fromkeys(curl1_aff.coeffs() + curl2_aff.coeffs()))
require(len(affine_curl_equations) == 9, "affine curl equation count")
curl_function = sp.lambdify(gauss_parameters, affine_curl_equations, "math")

affine_solutions = []
for values in product((-1, 0, 1), repeat=9):
    if all(value == 0 for value in curl_function(*values)):
        affine_solutions.append(values)
require(len(affine_solutions) == 375, "affine Gauss census")

# Independent identity test: a bidegree-at-most-two polynomial vanishing on
# the 3x3 integer grid vanishes identically.
grid_count = 0
for values in product((-1, 0, 1), repeat=9):
    U0, Ux, Uy, V0, Vx, Vy, W0, Wx, Wy = values
    good = True
    for x0 in (0, 1, 2):
        for y0 in (0, 1, 2):
            Uv = U0 + Ux * x0 + Uy * y0
            Vv = V0 + Vx * x0 + Vy * y0
            Wv = W0 + Wx * x0 + Wy * y0
            av = 1 + Uv * Vv
            curl2_value = Uy - (Ux * Wv + Uv * Wx)
            curl1_value = (Uy * Vv + Uv * Vy) - (
                Vx + Wx * av + Wv * (Ux * Vv + Uv * Vx)
            )
            if curl1_value or curl2_value:
                good = False
                break
        if not good:
            break
    if good:
        grid_count += 1
require(grid_count == 375, "independent affine Gauss census")

affine_histogram = {}
for values in affine_solutions:
    substitution = dict(zip(gauss_parameters, values))
    entries = [sp.expand(entry.subs(substitution)) for entry in (a_aff, b_aff, c_aff, d_aff)]
    monomials = set()
    for entry in entries:
        monomials.update(sp.Poly(entry, x, y).monoms())
    coefficient_vectors = []
    for monomial in sorted(monomials):
        if monomial != (0, 0):
            coefficient_vectors.append(
                [sp.Poly(entry, x, y).coeff_monomial(monomial) for entry in entries]
            )
    span = sp.Matrix(coefficient_vectors).rank() if coefficient_vectors else 0
    degree_p = max(sp.Poly(entries[0], x, y).total_degree(), sp.Poly(entries[1], x, y).total_degree()) + 1
    degree_q = max(sp.Poly(entries[2], x, y).total_degree(), sp.Poly(entries[3], x, y).total_degree()) + 1
    key = (degree_p, degree_q, span)
    affine_histogram[key] = affine_histogram.get(key, 0) + 1
require(max(key[2] for key in affine_histogram) == 1, "affine B1 span hostile")

print("affine Gauss chart B=1 universe:", 3**9)
print("curl-integrable (coefficient and grid implementations):", len(affine_solutions), grid_count)
print("(deg P,deg Q,coefficient-span) histogram:")
for key in sorted(affine_histogram):
    print(" ", key, affine_histogram[key])

# Positive tame control in the first balanced/full-span Gauss cell.
U_control = x + y**2 / 2
V_control = y + U_control**2 / 2
W_control = y
a_control = sp.expand(1 + U_control * V_control)
b_control = sp.expand(V_control + W_control * (1 + U_control * V_control))
c_control = sp.expand(U_control)
d_control = sp.expand(1 + U_control * W_control)
M_control = sp.Matrix([[a_control, b_control], [c_control, d_control]])
require(sp.expand(M_control.det()) == 1, "Gauss control determinant")
require(sp.expand(sp.diff(a_control, y) - sp.diff(b_control, x)) == 0, "Gauss control curl 1")
require(sp.expand(sp.diff(c_control, y) - sp.diff(d_control, x)) == 0, "Gauss control curl 2")
P_control = sp.expand(U_control + V_control**2 / 2)
Q_control = sp.expand(V_control)
require(sp.Matrix([P_control, Q_control]).jacobian([x, y]) == M_control, "Gauss control integration")

control_monomials = set()
for entry in (a_control, b_control, c_control, d_control):
    control_monomials.update(sp.Poly(entry, x, y).monoms())
control_vectors = []
for monomial in sorted(control_monomials):
    if monomial != (0, 0):
        control_vectors.append(
            [sp.Poly(entry, x, y).coeff_monomial(monomial) for entry in (a_control, b_control, c_control, d_control)]
        )
control_span = sp.Matrix(control_vectors).rank()
require(control_span == 4, "Gauss control full span")
T_control = sp.expand(c_control * b_control)
require(sp.expand(a_control * d_control - T_control) == 1, "Gauss consecutive fibres")

fake = sp.Matrix([[x * y + 1, x * y + 2], [x * y, x * y + 1]])
require(sp.expand(fake.det()) == 1, "determinant-only hostile")
fake_curls = (
    sp.expand(sp.diff(fake[0, 0], y) - sp.diff(fake[0, 1], x)),
    sp.expand(sp.diff(fake[1, 0], y) - sp.diff(fake[1, 1], x)),
)
require(fake_curls == (x - y, x - y), "determinant-only curl defects")
print("positive Gauss control: determinant=1, both curls=0, coefficient span=4, tame")
print("determinant-only hostile: both curl defects = x-y")


# ---------------------------------------------------------------------------
# IV. ICARM curve #273: exact point incidence, not an independent rank proof.
# ---------------------------------------------------------------------------

curve_a4 = int("-201769035260418549083594900060734240952308696994802735114305555")
curve_a6 = int(
    "1151107939141058565733479426024323225135665982951300586808823640527729578307228357301072889377"
)
curve_discriminant = int(
    "-46714661255308767314567688733841531918983356002159772613256840842851650254036518701100578342601553513579222272710220496887616034526983492843954090554197033638137245037791044053017600000000"
)
curve_points_raw = [
    ("-4761204159891138283979053265906", "44764265461782973805868732003346421827415264953"),
    ("-14158422539541566469588779426546", "34199834254251713784176619895082644508395077433"),
    ("-11522667358396562420423130332066", "-44115070023357103726405378140637465204943359607"),
    ("-204839531927226269712122049566", "-34531574232452693997231136031282772551453427107"),
    ("3899324051227528532535432912094", "20582352852872417675268569815574934013539218953"),
    ("149851368287976334870008075289384", "-1826442728148288630645637436047625928557963231657"),
    ("240440240734591134232325971191694", "3721941824016160691689265341458606456425791434553"),
    ("58446054919170749975942104376446/9", "-289145377197241504032247540119122580900747897469/27"),
    ("25642661602146479458845459929344", "-113306861325798987289137854129016658652160209297"),
    ("25720885078613923889202869994094", "-113918565504468051036791617945007239588074855047"),
    ("4956414590296956229584100339596814", "348939117745197814060339374812186839746231405619513"),
    ("725964821994104294477684670330094", "-19556488133953913131900670560396205869420775943047"),
    ("20802191136944676997135829374", "33866070189878993817062821320678356972094522793"),
    ("79052318332408565020526148386446/9", "-202982221452031541387733280916787231176177841869/27"),
    ("-11232245340662775388535509780886", "-44725045659073489550941507272219743508825024527"),
    ("8362456338772815315335239525614", "-6972475614865802969141741730862843401376795527"),
    ("2011658715643038193607509024534", "27447371869432010931671648582500375378543228993"),
    ("5027695440284894460797358334207726/529", "-116678851641395817353208818767411586490893208148849/12167"),
    ("24649144267565165528439068441554", "105612540318783792474731275264940719325335867213"),
    ("-12211389420609043025008816968566", "42356225616159991618318584560811156010370207173"),
    ("7798692390172953821075781106768126/1369", "-691870045568822811690292896396241871567834072004011/50653"),
    ("87157992815740534253438806045216/9", "277139378791840529410298740802253693668472375061/27"),
    ("30786757706172245427369935940751/4", "58841476683002984849182029306774218124047405249/8"),
    ("245309280348041323814668746104926/25", "1346501028820415725958868015485008289981037919061/125"),
    ("-343878076324392159036619356326", "-34934957027779219869199839566344035316624147307"),
    ("544211807917340289404451270094", "-32271721754226832038590040491036826507284103047"),
    ("20286216384652039303944170492166526/9409", "-24593234902246769413006506020777691495865223432164871/912673"),
    ("-25558163204018019740775243468600589874/1760929", "74707049582033426178338768659390679551818201954095350999/2336752783"),
    ("4546264873863829383537112534021848799606/398521369", "-145386763829319577901520209264368135012688669041886277075149/7955682089347"),
    ("1709164065046406773620054102684586/169", "26450264171408287955631955124255640794301463854841/2197"),
]
curve_points = [(Fraction(x_value), Fraction(y_value)) for x_value, y_value in curve_points_raw]
require(len(curve_points) == 30, "ICARM point count")
for x_value, y_value in curve_points:
    lhs = y_value**2 + x_value * y_value
    rhs = x_value**3 + curve_a4 * x_value + curve_a6
    require(lhs == rhs, "ICARM point off curve")

# Integral generalized Weierstrass coordinates have denominators d^2,d^3.
denominator_bases = []
for x_value, y_value in curve_points:
    d = isqrt(x_value.denominator)
    require(d * d == x_value.denominator, "x denominator not square")
    require(d**3 == y_value.denominator, "y denominator not cube")
    denominator_bases.append(d)

b2 = 1
b4 = 2 * curve_a4
b6 = 4 * curve_a6
b8 = curve_a6 - curve_a4**2
computed_discriminant = -(b2**2) * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
require(computed_discriminant == curve_discriminant, "ICARM discriminant")

point_heights = [max(qheight(x_value), qheight(y_value)) for x_value, y_value in curve_points]
integral_points = sum(x_value.denominator == y_value.denominator == 1 for x_value, y_value in curve_points)
print("\n[IV] ICARM elliptic-rank leaderboard curve #273")
print("model: y^2+xy=x^3+a4*x+a6")
print("a-invariant digit lengths: a4=%d a6=%d" % (len(str(abs(curve_a4))), len(str(abs(curve_a6)))))
print("displayed coefficient height digits:", len(str(max(abs(curve_a4), abs(curve_a6)))))
print("stored points checked exactly on curve: %d/%d" % (len(curve_points), len(curve_points)))
print("integral/rational-nonintegral points: %d/%d" % (integral_points, len(curve_points) - integral_points))
print("denominator bases d (den x=d^2, den y=d^3):", sorted(set(denominator_bases)))
print("maximum denominator base:", max(denominator_bases))
print("affine point-height digit range: %d..%d" % (len(str(min(point_heights))), len(str(max(point_heights)))))
print("minimal discriminant matches the frozen primary-data value")
print("CITED/CURRENT DATABASE FIELD: rank_lower_bound=30")
print("NOT REPROVED HERE: independence, rank>=30, analytic rank, or exact rank")


# A collision-normalization hostile for height: normalizing one rational
# THM-1300 collision to 0 and the other to (1,1,1) expands both support/height.
p0 = sp.Matrix([0, 0, -sp.Rational(1, 4)])
p1 = sp.Matrix([1, -sp.Rational(3, 2), sp.Rational(13, 2)])
Jp0 = F.jacobian([x, y, z]).subs({x: p0[0], y: p0[1], z: p0[2]})
collision_vector = Jp0 * (p1 - p0)
collision_scale = sp.diag(*list(collision_vector))
source_linear = Jp0.inv() * collision_scale
new_variables = sp.Matrix(sp.symbols("X Y Z"))
X, Y, Z = list(new_variables)
source_point = p0 + source_linear * new_variables
translated = F.subs({x: source_point[0], y: source_point[1], z: source_point[2]}) - F.subs(
    {x: p0[0], y: p0[1], z: p0[2]}
)
collision_normalized = (collision_scale.inv() * translated).applyfunc(sp.expand)
require(collision_normalized.subs({X: 0, Y: 0, Z: 0}) == sp.zeros(3, 1), "collision gauge origin")
require(
    collision_normalized.jacobian([X, Y, Z]).subs({X: 0, Y: 0, Z: 0}) == sp.eye(3),
    "collision gauge jet",
)
require(collision_normalized.subs({X: 1, Y: 1, Z: 1}) == sp.zeros(3, 1), "collision gauge second point")
collision_profile = coefficient_profile(collision_normalized, (X, Y, Z))
print("\n[V] height hostile")
print(
    "forcing collision 0,(1,1,1) and identity jet: terms=%s H=%d maxden=%d"
    % (collision_profile["terms"], collision_profile["height"], collision_profile["max_denominator"])
)
print("lesson: collision normalization can destroy small displayed height; search orbit/generator height first")

print("\nVERDICT: ALL LOCAL EXACT AND FINITE CONTROLS PASS")
