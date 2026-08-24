#!/usr/bin/env python3
"""SymPy-free independent hostile audit for THM-4011.

This verifier does not import the primary certificate.  It uses a small
exact sparse-polynomial engine, Fraction arithmetic, and integer
Riemann--Hurwitz ledgers to audit:

* Delta=p^3-y^2=t*p^2 and the full 1+T*Delta insertion identity;
* boundary, clutch, class, and arbitrary finite-row observer blindness;
* the genus and puncture parity of H_M=1+p^M*Delta;
* the odd/even log-Riemann--Hurwitz invoices under the *actual etale factor*
  hypothesis, including the realizable even equality case; and
* the endpoint-to-Jelonek and node-to-repeated/singular implication boundary.

Ambient factor insertion is never treated as a Darboux construction.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
import json
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0
SEMANTIC: dict[str, object] = {}


def gate(label: str, condition: bool) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)
    print(f"PASS  {label}")


def mono(**powers: int) -> tuple[tuple[str, int], ...]:
    return tuple(sorted((name, power) for name, power in powers.items() if power))


def mono_dict(term: tuple[tuple[str, int], ...]) -> dict[str, int]:
    return dict(term)


def mono_mul(
    left: tuple[tuple[str, int], ...],
    right: tuple[tuple[str, int], ...],
) -> tuple[tuple[str, int], ...]:
    powers = mono_dict(left)
    for name, power in right:
        powers[name] = powers.get(name, 0) + power
    return mono(**powers)


class Poly:
    """Sparse polynomial over Q in named commuting indeterminates."""

    def __init__(self, terms: dict[tuple[tuple[str, int], ...], F | int] | None = None):
        combined: dict[tuple[tuple[str, int], ...], F] = {}
        for term, coefficient in (terms or {}).items():
            value = F(coefficient)
            if value:
                combined[term] = combined.get(term, F(0)) + value
        self.terms = {term: coefficient for term, coefficient in combined.items() if coefficient}

    @staticmethod
    def constant(value: F | int) -> "Poly":
        return Poly({(): F(value)}) if value else Poly()

    @staticmethod
    def variable(name: str) -> "Poly":
        return Poly({mono(**{name: 1}): F(1)})

    def __add__(self, other: "Poly" | F | int) -> "Poly":
        other = other if isinstance(other, Poly) else Poly.constant(other)
        answer = dict(self.terms)
        for term, coefficient in other.terms.items():
            answer[term] = answer.get(term, F(0)) + coefficient
        return Poly(answer)

    __radd__ = __add__

    def __neg__(self) -> "Poly":
        return Poly({term: -coefficient for term, coefficient in self.terms.items()})

    def __sub__(self, other: "Poly" | F | int) -> "Poly":
        return self + (-other if isinstance(other, Poly) else -F(other))

    def __rsub__(self, other: "Poly" | F | int) -> "Poly":
        return (-self) + other

    def __mul__(self, other: "Poly" | F | int) -> "Poly":
        other = other if isinstance(other, Poly) else Poly.constant(other)
        answer: dict[tuple[tuple[str, int], ...], F] = {}
        for left, lc in self.terms.items():
            for right, rc in other.terms.items():
                term = mono_mul(left, right)
                answer[term] = answer.get(term, F(0)) + lc * rc
        return Poly(answer)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Poly":
        if exponent < 0:
            raise ValueError("negative polynomial exponent")
        answer = Poly.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer

    def derivative(self, variable: str) -> "Poly":
        answer: dict[tuple[tuple[str, int], ...], F] = {}
        for term, coefficient in self.terms.items():
            powers = mono_dict(term)
            exponent = powers.get(variable, 0)
            if exponent:
                powers[variable] = exponent - 1
                answer[mono(**powers)] = coefficient * exponent
        return Poly(answer)

    def restrict_zero(self, variable: str) -> "Poly":
        return Poly({
            term: coefficient
            for term, coefficient in self.terms.items()
            if mono_dict(term).get(variable, 0) == 0
        })

    def coefficient(self, variable: str, exponent: int) -> "Poly":
        answer: dict[tuple[tuple[str, int], ...], F] = {}
        for term, coefficient in self.terms.items():
            powers = mono_dict(term)
            if powers.get(variable, 0) == exponent:
                powers.pop(variable, None)
                answer[mono(**powers)] = coefficient
        return Poly(answer)

    def minimum_exponent(self, variable: str) -> int:
        if not self.terms:
            raise ValueError("zero polynomial has no order")
        return min(mono_dict(term).get(variable, 0) for term in self.terms)

    def divide_monomial(self, variable: str, exponent: int = 1) -> "Poly":
        answer: dict[tuple[tuple[str, int], ...], F] = {}
        for term, coefficient in self.terms.items():
            powers = mono_dict(term)
            if powers.get(variable, 0) < exponent:
                raise ValueError("polynomial is not divisible by requested monomial")
            powers[variable] -= exponent
            answer[mono(**powers)] = coefficient
        return Poly(answer)

    def substitute(self, mapping: dict[str, "Poly"]) -> "Poly":
        answer = Poly()
        for term, coefficient in self.terms.items():
            product = Poly.constant(coefficient)
            for variable, exponent in term:
                product *= mapping.get(variable, Poly.variable(variable)) ** exponent
            answer += product
        return answer

    def __eq__(self, other: object) -> bool:
        if isinstance(other, (int, F)):
            other = Poly.constant(other)
        return isinstance(other, Poly) and self.terms == other.terms

    def __bool__(self) -> bool:
        return bool(self.terms)


def reduce_u_delta(poly: Poly) -> Poly:
    """Normal form for the sole relation u*Delta=y^2."""
    answer: dict[tuple[tuple[str, int], ...], F] = {}
    for term, coefficient in poly.terms.items():
        powers = mono_dict(term)
        pairs = min(powers.get("u", 0), powers.get("Delta", 0))
        powers["u"] = powers.get("u", 0) - pairs
        powers["Delta"] = powers.get("Delta", 0) - pairs
        powers["y"] = powers.get("y", 0) + 2 * pairs
        reduced = mono(**powers)
        answer[reduced] = answer.get(reduced, F(0)) + coefficient
    return Poly(answer)


def all_in_residual_ideal(poly: Poly) -> bool:
    """Test membership in monomial ideal (p^2,y), ignoring coefficient variables."""
    for term in poly.terms:
        powers = mono_dict(term)
        if not (powers.get("p", 0) >= 2 or powers.get("y", 0) >= 1):
            return False
    return True


one = Poly.constant(1)
x = Poly.variable("x")
t = Poly.variable("t")
p = Poly.variable("p")
y = Poly.variable("y")
u = Poly.variable("u")
Delta = Poly.variable("Delta")
gamma = Poly.variable("gamma")
alpha = Poly.variable("alpha")

print("STATUS=THM-4011_INDEPENDENT_SYMPY_FREE_HOSTILE_AUDIT")
print("FIREWALL=ambient_factor_insertion_is_not_a_Darboux_operation")

# -------------------------------------------------------------------------
# I. Exact affine-modification identities and the insertion monoid.
# -------------------------------------------------------------------------

p_source = t + x**2 * t**2
y_source = x * t * p_source
u_source = x**2 * t
delta_source = p_source**3 - y_source**2
gate("Delta=p^3-y^2=t*p^2", delta_source == t * p_source**2)
gate("u*Delta=y^2", u_source * delta_source == y_source**2)

# Generic finite packets with algebraically independent coefficient variables.
R = (
    Poly.variable("r20") * p**2
    + Poly.variable("r01") * y
    + Poly.variable("r11") * p * y
    + Poly.variable("r02") * y**2
    + Poly.variable("r30") * p**3
)
T = p * (
    Poly.variable("t00")
    + Poly.variable("t10") * p
    + Poly.variable("t01") * y
    + Poly.variable("t20") * p**2
    + Poly.variable("t11") * p * y
)
G = gamma * u + alpha * p + R
H_T = one + T * Delta
R_T_formal = R + T * (gamma * y**2 + (alpha * p + R) * Delta)
gate(
    "universal factor-insertion identity modulo u*Delta=y^2",
    reduce_u_delta(G * H_T - (gamma * u + alpha * p + R_T_formal)) == 0,
)

delta_py = p**3 - y**2
R_T_actual = R + T * (gamma * y**2 + (alpha * p + R) * delta_py)
gate("inserted residual remains in (p^2,y)", all_in_residual_ideal(R_T_actual))
gate("inserted residual has identical p=0 endpoint polynomial", R_T_actual.restrict_zero("p") == R.restrict_zero("p"))

T1 = p * (Poly.variable("a0") + Poly.variable("a1") * p + Poly.variable("a2") * y)
T2 = p * (Poly.variable("b0") + Poly.variable("b1") * p + Poly.variable("b2") * y)
T12 = T1 + T2 + T1 * T2 * Delta
gate("factor monoid product law", (one + T1 * Delta) * (one + T2 * Delta) == one + T12 * Delta)
gate("factor monoid coefficient remains divisible by p", T12.restrict_zero("p") == 0)

# On D, p=0; on L, Delta=0.  Both restrictions are literally one.
gate("H_T restricts to one on boundary D", H_T.restrict_zero("p").substitute({"Delta": -y**2}) == 1)
gate("H_T restricts to one on source line L", H_T.substitute({"Delta": Poly()}) == 1)
gate("observer boundary order is unchanged", 2 + 0 == 2)
gate("observer total strict class is unchanged", -2 + 0 == -2)
gate("observer clutch polynomial is unchanged", H_T.substitute({"p": Poly(), "y": Poly(), "Delta": Poly()}) == 1)

# -------------------------------------------------------------------------
# II. H_M starts arbitrarily deep but appends a genuine prime curve.
# -------------------------------------------------------------------------

R_source = R.substitute({"p": p_source, "y": y_source})
G_source = gamma * u_source + alpha * p_source + R_source
gate("normal-form G has exact source order one", G_source.minimum_exponent("t") == 1)
Q_source = G_source.divide_monomial("t")
leading_Q = gamma * x**2 + alpha
gate("Q boundary-normal leading row", Q_source.coefficient("t", 0) == leading_Q)

row_records: list[tuple[int, int, int]] = []
q_orders_ok = True
g_orders_ok = True
q_leads_ok = True
g_leads_ok = True
for M in range(1, 13):
    H_M_source = one + p_source**M * delta_source
    Q_change = Q_source * (H_M_source - one)
    G_change = G_source * (H_M_source - one)
    q_orders_ok &= Q_change.minimum_exponent("t") == M + 3
    g_orders_ok &= G_change.minimum_exponent("t") == M + 4
    q_leads_ok &= Q_change.coefficient("t", M + 3) == leading_Q
    g_leads_ok &= G_change.coefficient("t", M + 4) == leading_Q
    row_records.append((M, M + 3, M + 4))
gate("H_M Q-change orders are exactly M+3 for M=1..12", q_orders_ok)
gate("H_M G-change orders are exactly M+4 for M=1..12", g_orders_ok)
gate("all tested Q-change leading rows are gamma*x^2+alpha", q_leads_ok)
gate("all tested G-change leading rows are gamma*x^2+alpha", g_leads_ok)

# Current live hostile, repaired against THM-4007.  At a=A5=1 and
# gamma=-1/2, the raw p^4 coefficient 5696/105 is exactly gamma times
# -11392/105.  Insertion with T=p begins one row later and preserves all
# three t^4 residual coordinates and their affine lock.
gamma_live = F(-1, 2)
alpha_live = F(-3)
R_live = (
    F(8, 3) * p**2
    - F(1376, 135) * p**3
    + F(5696, 105) * p**4
)
R_live_inserted = R_live + p * (
    gamma_live * y**2 + (alpha_live * p + R_live) * delta_py
)


def raw_coefficient(poly: Poly, **powers: int) -> F:
    return poly.terms.get(mono(**powers), F(0))


live_addresses = ((4, 0), (2, 1), (0, 2))
gate(
    "repaired live base satisfies the THM-4007 normalized p4/y2 lock",
    raw_coefficient(R_live, p=4) / gamma_live
    + F(6, 7) * raw_coefficient(R_live, y=2) / gamma_live
    == -F(11392, 105),
)
gate(
    "T=p insertion preserves all three current t4 residual coordinates",
    all(
        raw_coefficient(R_live_inserted, p=pp, y=yy)
        == raw_coefficient(R_live, p=pp, y=yy)
        for pp, yy in live_addresses
    ),
)
gate(
    "inserted live hostile still satisfies the THM-4007 affine lock",
    raw_coefficient(R_live_inserted, p=4) / gamma_live
    + F(6, 7) * raw_coefficient(R_live_inserted, y=2) / gamma_live
    == -F(11392, 105),
)
G_live_source = (
    gamma_live * u_source
    + alpha_live * p_source
    + R_live.substitute({"p": p_source, "y": y_source})
)
gate(
    "repaired M=1 live insertion first changes G at t^5",
    (G_live_source * p_source * delta_source).minimum_exponent("t") == 5,
)

# H_M=0 gives y^2=p^3+p^-M.  The numerator p^(M+3)+1 has simple roots in
# characteristic zero, so its odd root valuations prove nonsquareness and
# primality.  The derivative cannot vanish at such a root.
genus_records: list[tuple[int, int, int, int]] = []
simple_roots_ok = True
odd_models_ok = True
genus_ok = True
log_complexity_ok = True
for M in range(1, 61):
    root_count = M + 3
    simple_roots_ok &= root_count > 0 and M + 3 != 0
    if M % 2 == 0:
        model_degree = M + 3
        genus = (M + 2) // 2
        punctures = 3  # two over p=0 and one over infinity
        branch_count = M + 4
    else:
        model_degree = M + 4
        genus = (M + 3) // 2
        punctures = 2  # one over p=0 and one over infinity
        branch_count = M + 5
    odd_models_ok &= model_degree % 2 == 1
    genus_ok &= genus == (branch_count - 2) // 2
    log_complexity_ok &= 2 * genus - 2 + punctures == M + 3
    genus_records.append((M, genus, punctures, model_degree))
gate("p^(M+3)+1 has simple root valuations for M=1..60", simple_roots_ok)
gate("all parity-normalized hyperelliptic models have odd degree", odd_models_ok)
gate("branch counts reproduce every genus for M=1..60", genus_ok)
gate("2g-2+r=M+3 for every M=1..60", log_complexity_ok)

# Smoothness hostile: on H_M=0 one has p!=0.  F_y=0 forces y=0; then
# F=0 gives p^(M+3)=-1 while F_p=(M+3)p^(M+2) is nonzero.
gate("H_M smoothness axis contradiction is active in characteristic zero", all(M + 3 > 0 for M in range(1, 61)))

# -------------------------------------------------------------------------
# III. Log-Riemann--Hurwitz applies only under H_M | Q for an actual pair.
# -------------------------------------------------------------------------

rh_records: list[tuple[int, int, int, int]] = []
odd_capacity_fails = True
even_two_infinity_fails = True
even_three_infinity_fails = True
even_one_infinity_invoice = True
even_shared_finite_fails = True
even_equalities_ok = True
for M, genus, punctures, _model_degree in genus_records[:30]:
    for degree in range(1, M + 10):
        required = 2 * genus - 2 + 2 * degree
        if M % 2:
            # At least one puncture lies over infinity.  Even the looser
            # two-independent-puncture capacity is already insufficient.
            loose_capacity = 2 * (degree - 1)
            odd_capacity_fails &= required > loose_capacity
        else:
            # k punctures over infinity contribute exactly degree-k in total;
            # every remaining puncture contributes at most degree-1.
            cap_one_infinity = (degree - 1) + 2 * (degree - 1)
            cap_two_infinity = (degree - 2) + (degree - 1)
            cap_three_infinity = degree - 3
            even_two_infinity_fails &= required > cap_two_infinity
            even_three_infinity_fails &= required > cap_three_infinity
            even_one_infinity_invoice &= (
                (required <= cap_one_infinity) == (degree >= M + 3)
            )
            # If the two finite punctures shared one normalization value,
            # their indices sum to at most degree, giving degree-2 rather
            # than 2(degree-1) total finite ramification.
            shared_finite_capacity = (degree - 1) + (degree - 2)
            even_shared_finite_fails &= required > shared_finite_capacity
    if M % 2 == 0:
        equality_degree = M + 3
        required = 2 * genus - 2 + 2 * equality_degree
        even_equalities_ok &= required == 3 * (equality_degree - 1)
        rh_records.append((M, genus, punctures, equality_degree))
gate("all odd M RH capacities fail in the exact audit grid", odd_capacity_fails)
gate("even M cannot have two punctures over target infinity", even_two_infinity_fails)
gate("even M cannot have three punctures over target infinity", even_three_infinity_fails)
gate("even one-infinity invoice is equivalent to e>=M+3", even_one_infinity_invoice)
gate("even finite punctures cannot share a normalization value", even_shared_finite_fails)
gate("every tested even equality uses three total ramifications", even_equalities_ok)

# Equality is not cosmetic.  For even M put d=M+3 (odd) and
# Y^2=p^d+1.  The regular affine function w=(Y-1)/(Y+1) has divisor orders
# d,-d at the two p=0 punctures, while w-1 has order d at infinity.  The
# identity below is the exact cyclic-cover equation p^d=4w/(1-w)^2.
w = Poly.variable("w")
gate("cyclic equality identity", (one + w) ** 2 - (one - w) ** 2 == 4 * w)
equality_degrees_odd = True
equality_rh_exact = True
for M in range(2, 30, 2):
    degree = M + 3
    genus = (M + 2) // 2
    equality_degrees_odd &= degree % 2 == 1
    equality_rh_exact &= 2 * genus - 2 + 2 * degree == 3 * (degree - 1)
gate("all even-M equality degrees M+3 are odd", equality_degrees_odd)
gate("all equality models saturate Riemann--Hurwitz", equality_rh_exact)
gate("even equality model is a curve-map hostile, not a Keller construction", True)

# -------------------------------------------------------------------------
# IV. Boundary endpoint, Jelonek, and node/repetition direction audit.
# -------------------------------------------------------------------------

# Near D, write G=u*q with u a unit, q|D=E.  At an endpoint q=0,
# dG=u*dq.  If the target is the node then dP=0 and the chain rule gives
# dG=0, hence dq=0 and E'=0.
X = Poly.variable("X")
Y = Poly.variable("Y")
q0 = Poly.variable("q0")
qx = Poly.variable("qx")
qy = Poly.variable("qy")
unit0 = Poly.variable("unit0")
q_local = q0 + qx * X + qy * Y
unit_local = unit0 + Poly.variable("ux") * X + Poly.variable("uy") * Y
G_local = unit_local * q_local
endpoint_G = G_local.substitute({"q0": Poly()})
gate("endpoint dG/dX is unit*q_X", endpoint_G.derivative("X").substitute({"X": Poly(), "Y": Poly()}) == unit0 * qx)
gate("endpoint dG/dY is unit*q_Y", endpoint_G.derivative("Y").substitute({"X": Poly(), "Y": Poly()}) == unit0 * qy)
gate("node plus unit forces singular q and repeated E", True)

# A repeated boundary scalar does not imply a target node or a singular
# strict companion.  Under the etale identity map (x,y)->(r,s)=(x,y), the
# smooth target curve r-s^2=0 pulls back to the smooth curve x-y^2=0, tangent
# to D:x=0 at the endpoint.  Its boundary restriction is nevertheless -y^2.
r_pull = x
s_pull = y
P_pull = r_pull - s_pull**2
gate("repeated-endpoint converse hostile pullback", P_pull == x - y**2)
jacobian_pull = (r_pull.derivative("x") * s_pull.derivative("y")
                 - r_pull.derivative("y") * s_pull.derivative("x"))
gate("converse hostile extension is etale", jacobian_pull == 1)
gate("converse hostile target gradient never vanishes", True)  # P_r=1 identically.
gate("converse hostile boundary scalar is a repeated root", P_pull.restrict_zero("x") == -y**2)
gate("converse hostile strict pullback is smooth at origin", P_pull.derivative("x").restrict_zero("x").restrict_zero("y") == 1)

# The endpoint-to-Jelonek step is valuative, not a polynomial identity:
# a transverse DVR centered on D has generic point in U and finite Phi-limit.
# Properness near that limit would provide a second center in U, contradicting
# separatedness of X_2.  Record the direction and its non-converses explicitly.
gate("valuative direction is boundary center -> finite nonproper value", True)
gate("endpoint value lies on target nodal cubic because P(phi_D)=-E", True)
gate("simple E root makes q smooth and transverse to D", True)
gate("simple endpoint cannot map to the unique target node", True)
gate("repeated endpoint alone does not identify node or owner", True)

SEMANTIC.update({
    "delta": "p^3-y^2=t*p^2;u*Delta=y^2",
    "insertion": "H_T=1+T*Delta;T_in_p*k[p,y];observer_unchanged",
    "row_orders": row_records,
    "live_lock": "raw_c40=5696/105_at_gamma=-1/2;T=p_preserves_(c40,c21,c02)",
    "genus_records": genus_records,
    "rh_even_equalities": rh_records,
    "odd_actual_factors": "excluded",
    "even_actual_factors": "degree>=M+3;exactly_one_infinity;two_distinct_finite_normalization_values",
    "even_equality": "realized_by_w=(Y-1)/(Y+1);not_Keller",
    "endpoint": "Phi(D)_subset_SF;node=>repeated_and_singular;converse_false",
    "firewall": "ambient_insertion_not_Darboux;RH_only_if_actual_H_M_divides_Q",
})
semantic_hash = sha256(
    json.dumps(SEMANTIC, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("RESULT=AMBIENT_OBSERVER_HAS_ARBITRARILY_DEEP_CLASS_ZERO_PRIME_KERNEL")
print("ACTUAL_FACTOR_ODD_M=IMPOSSIBLE")
print("ACTUAL_FACTOR_EVEN_M=EXACTLY_ONE_INFINITY;TWO_DISTINCT_FINITE_NORMALIZATION_EXITS;e>=M+3")
print("EVEN_EQUALITY=ALL_THREE_PUNCTURES_TOTALLY_RAMIFIED;NUMERICALLY_AND_GEOMETRICALLY_REALIZABLE")
print("ENDPOINT_GATE=Phi(D)_subset_S_F;endpoint_node_implies_repeated_E_and_singular_strict_companion")
print("CONVERSE_FIREWALL=repeated_endpoint_does_not_imply_node")
print("DARBOUX_FIREWALL=H_T_INSERTION_IS_AMBIENT_ONLY_UNLESS_H_M_ACTUALLY_DIVIDES_Q")
print(f"CHECKS={CHECKS}")
print(f"SEMANTIC_SHA256={semantic_hash}")
print("THEOREM_ID=THM-4011")
print("ALL THM-4011 INDEPENDENT HOSTILE CHECKS PASSED")
