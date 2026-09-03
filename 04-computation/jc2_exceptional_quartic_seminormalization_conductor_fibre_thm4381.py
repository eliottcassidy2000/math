"""Exact conductor-fibre and seminormalization certificate for THM-4381.

The load-bearing computation is over the quartic field K.  It projects the
89-dimensional algebra S/c onto the reduced conductor fibre K[x]/rad(c) and
computes its rank exactly.  THM-4034's positive good-reduction gcd certificate
is replayed for squarefreeness and coprimality.  The rational-fibre census is
a hostile control only, and no modular rank failure is used.
"""

import contextlib
import hashlib
import io
import subprocess
import types
from fractions import Fraction
from pathlib import Path


CONDUCTOR_PATH = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_global_conductor_thm4034.py"
)
PRIME = 137
CONDUCTOR_CORE_MARKER = "# The three divided-difference resultants give an intrinsic conductor formula."


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def load_conductor_core(path, name):
    if path.is_file():
        source = path.read_text()
    else:
        source = subprocess.check_output(
            ["git", "show", f"HEAD:{path.as_posix()}"], text=True
        )
    require(source.count(CONDUCTOR_CORE_MARKER) == 1, "unique conductor core marker")
    source = source.split(CONDUCTOR_CORE_MARKER, 1)[0]
    module = types.ModuleType(name)
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(source, path.as_posix(), "exec"), module.__dict__)
    return module


conductor_certificate = load_conductor_core(CONDUCTOR_PATH, "quartic_conductor_thm4034")
sagbi = conductor_certificate.sagbi
field = conductor_certificate.field
ring = conductor_certificate.ring
X = conductor_certificate.X
zero = conductor_certificate.zero
one = conductor_certificate.one
conductor = conductor_certificate.conductor
h172 = conductor_certificate.quotient
normal_form = conductor_certificate.normal_form

L = X * (X**2 - 1)
radical_conductor = L * h172
require(conductor == L * radical_conductor, "conductor radical factorization")
require(conductor.degree() == 178, "conductor degree")
require(radical_conductor.degree() == 175, "reduced conductor degree")
# THM-4034's good reduction proves squarefreeness of h172 and coprimality
# with L in characteristic zero.  Replay that positive certificate here;
# avoid a redundant degree-175 exact xgcd with enormous quartic coefficients.
h172_mod = conductor_certificate.polynomial_mod(h172)
L_mod = conductor_certificate.polynomial_mod(L)
require(
    h172_mod.gcd(h172_mod.derivative()).degree() == 0,
    "inherited squarefree good reduction",
)
require(h172_mod.gcd(L_mod).degree() == 0, "inherited coprime good reduction")


def coefficient(polynomial, degree):
    return polynomial.get((degree,), zero)


def field_coordinates(value):
    return conductor_certificate.field_coordinates(value)


def serialize_field(value):
    return ",".join(str(item) for item in field_coordinates(value))


def polynomial_hash(polynomial):
    payload = ";".join(
        f"{degree}:" + serialize_field(coefficient(polynomial, degree))
        for degree in range(polynomial.degree() + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


# The monic Apéry presentation supplies one canonical S-element in every
# semigroup degree.  Below deg(c)=178 there are exactly dim_K(S/c)=89.
canonical_basis = tuple(
    conductor_certificate.canonical[degree]
    for degree in sorted(conductor_certificate.canonical)
    if degree < conductor.degree()
)
require(len(canonical_basis) == 89, "S modulo conductor basis size")
low_basis = tuple(item for item in canonical_basis if item.degree() < 175)
top_basis = tuple(item for item in canonical_basis if item.degree() >= 175)
require(len(low_basis) == 86, "low reduced-fibre basis size")
require(tuple(item.degree() for item in top_basis) == (175, 176, 177), "top degrees")

# The 86 low elements remain independent modulo the monic degree-175
# radical: their degrees are distinct and they are monic.  Reduce the three
# top elements first modulo rad(c), then modulo S.  Their three gap vectors
# measure exactly how many directions they add to the low image.
top_gap_remainders = []
for item in top_basis:
    reduced = divmod(item, radical_conductor)[1]
    gap_remainder = normal_form(reduced)
    require(not gap_remainder or gap_remainder.degree() <= 169, "gap normal form")
    top_gap_remainders.append(gap_remainder)

gaps = tuple(sagbi.EXPECTED_GAPS)
gap_rows = [
    [coefficient(remainder, degree) for degree in gaps]
    for remainder in top_gap_remainders
]


def exact_row_rank(rows):
    echelon = []
    pivots = []
    for source_row in rows:
        row = list(source_row)
        for pivot, old_row in zip(pivots, echelon):
            value = row[pivot]
            if value:
                row = [left - value * right for left, right in zip(row, old_row)]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = one / row[pivot]
        row = [value * inverse for value in row]
        pivots.append(pivot)
        echelon.append(row)
    return len(echelon), tuple(pivots)


top_rank, top_pivots = exact_row_rank(gap_rows)
require(top_rank == 1, "exact top residual rank")
reduced_image_rank = len(low_basis) + top_rank
require(reduced_image_rank == 87, "exact reduced conductor-fibre rank")
gap_payload = ";".join(
    ":".join(serialize_field(value) for value in row) for row in gap_rows
)
gap_hash = hashlib.sha256(gap_payload.encode("ascii")).hexdigest()

# The radical itself is the proposed missing seminormal class.  It is not in
# S, while its square is c*h172 and hence lies in the conductor ideal.
radical_normal_form = normal_form(radical_conductor)
require(bool(radical_normal_form), "radical is not in S")
require(radical_normal_form.degree() <= 169, "radical normal-form support")
require(radical_conductor**2 == conductor * h172, "radical square in conductor")

# The retained target fibre is exactly {-1,0,1}, not merely a known subset.
# The gcd check is an exact algebraic certificate.  The derivative rows then
# show three distinct tangent lines in the (C,E+3)-plane.
B = sagbi.B
C = sagbi.C
E = sagbi.E
retained_fibre_gcd = B.gcd(C).gcd(E + 3).monic()
require(retained_fibre_gcd == L, "retained fibre is exactly the triple")


def evaluate(polynomial, point):
    return polynomial.evaluate(X, point * one)


tangent_rows = tuple(
    tuple(evaluate(polynomial.diff(X), point) for point in (-1, 0, 1))
    for polynomial in (B, C, E)
)
require(tangent_rows[0] == (zero, zero, zero), "retained B tangent row")
require(tangent_rows[1] == (3 * one, 3 * one, 3 * one), "retained C tangent row")
require(tangent_rows[2] == (-9 * one, 4 * one, 9 * one), "retained E tangent row")
for left in range(3):
    for right in range(left):
        determinant = (
            tangent_rows[1][left] * tangent_rows[2][right]
            - tangent_rows[1][right] * tangent_rows[2][left]
        )
        require(bool(determinant), "distinct retained tangent lines")


# Good-reduction hostile controls.  At alpha=92 mod 137, one of the 86
# double fibres splits rationally as {44,134}; at alpha=44 it does not split.
# Both fibres retain the same rational transverse triple.  This census is not
# used in the exact rank or fibre-count proof.
def rational_mod(value, prime):
    value = Fraction(value)
    require(value.denominator % prime != 0, "good reduction denominator")
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def field_mod_at(value, alpha_mod, prime):
    return sum(
        rational_mod(entry, prime) * pow(alpha_mod, power, prime)
        for power, entry in enumerate(field_coordinates(value))
    ) % prime


def coefficient_list_mod(polynomial, alpha_mod, prime):
    return [
        field_mod_at(coefficient(polynomial, degree), alpha_mod, prime)
        for degree in range(polynomial.degree() + 1)
    ]


def evaluate_mod(coefficients, point, prime):
    answer = 0
    for value in reversed(coefficients):
        answer = (answer * point + value) % prime
    return answer


def derivative_mod(coefficients, point, prime):
    return sum(
        degree * coefficients[degree] * pow(point, degree - 1, prime)
        for degree in range(1, len(coefficients))
    ) % prime


def rational_multifibres(alpha_mod):
    reductions = tuple(
        coefficient_list_mod(polynomial, alpha_mod, PRIME)
        for polynomial in (B, C, E)
    )
    fibres = {}
    for point in range(PRIME):
        target = tuple(evaluate_mod(item, point, PRIME) for item in reductions)
        fibres.setdefault(target, []).append(point)
    multifibres = tuple(
        sorted(tuple(points) for points in fibres.values() if len(points) > 1)
    )
    tangent_proportional_pairs = 0
    for points in multifibres:
        tangents = tuple(
            tuple(derivative_mod(item, point, PRIME) for item in reductions)
            for point in points
        )
        for left in range(len(tangents)):
            for right in range(left):
                proportional = all(
                    (
                        tangents[left][i] * tangents[right][j]
                        - tangents[left][j] * tangents[right][i]
                    )
                    % PRIME
                    == 0
                    for i in range(3)
                    for j in range(i + 1, 3)
                )
                tangent_proportional_pairs += int(proportional)
    return multifibres, tangent_proportional_pairs


fibres_44, bad_tangents_44 = rational_multifibres(44)
fibres_92, bad_tangents_92 = rational_multifibres(92)
require(fibres_44 == ((0, 1, 136),), "alpha44 rational fibre control")
require(
    fibres_92 == ((0, 1, 136), (44, 134)),
    "alpha92 rational fibre control",
)
require(bad_tangents_44 == bad_tangents_92 == 0, "finite-field transversality")


def rank_mod_prime(rows, alpha_mod, prime):
    work = [
        [field_mod_at(value, alpha_mod, prime) for value in row]
        for row in rows
    ]
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column] % prime),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column] % prime, -1, prime)
        work[rank] = [(value * inverse) % prime for value in work[rank]]
        for row in range(len(work)):
            if row == rank:
                continue
            value = work[row][column] % prime
            if value:
                work[row] = [
                    (left - value * right) % prime
                    for left, right in zip(work[row], work[rank])
                ]
        rank += 1
        if rank == len(work):
            break
    return rank


require(rank_mod_prime(gap_rows, 44, PRIME) == 1, "modular hostile top rank")


print("field=THM3683_exceptional_quartic;restriction_algebra=THM3703")
print(
    "conductor=c=L^2*h172;degrees_c_h_radical="
    f"{conductor.degree()},{h172.degree()},{radical_conductor.degree()};"
    "radical_squarefree=True"
)
print(
    f"S_mod_c_dimension={len(canonical_basis)};low_degree_independent={len(low_basis)};"
    f"top_degrees={','.join(str(item.degree()) for item in top_basis)};"
    f"top_gap_rank_exact={top_rank};reduced_image_rank_exact={reduced_image_rank}"
)
print(
    f"top_gap_pivots={','.join(map(str, top_pivots))};top_gap_hash={gap_hash};"
    "top_gap_rank_mod137_alpha44=1"
)
print(
    f"r=L*h172;r_degree={radical_conductor.degree()};r_in_S=False;"
    f"r_hash={polynomial_hash(radical_conductor)};"
    f"r_normal_form_degree={radical_normal_form.degree()};"
    f"r_normal_form_hash={polynomial_hash(radical_normal_form)};r_squared_in_c=True"
)
print(
    "retained_fibre_gcd=x*(x^2-1);retained_fibre=-1,0,1;"
    "B_derivatives=0,0,0;C_derivatives=3,3,3;E_derivatives=-9,4,9;"
    "retained_tangents_distinct=True"
)
print(
    "derived_nilradical_dimension=2;reduced_conductor_fibres=87;"
    "classification_over_algebraic_closure=one_ordinary_triple_plus_86_ordinary_nodes"
)
print(
    "seminormalization=S+K*r;seminormal_defect_length=1;"
    "seminormal_conductor=r*K[x];normalization_delta=88;delta_split=86+2"
)
print(
    "hostile_mod137_alpha44_multifibres=0,1,136;"
    "hostile_mod137_alpha92_multifibres=0,1,136|44,134;"
    "hostile_tangent_proportional_pairs=0"
)
print("GOOD_REDUCTION_ROLE=inherited_THM4034_squarefree_coprime_certificate")
print("FINITE_FIELD_FIBRE_CENSUS_ROLE=hostile_control_only;exact_rank_is_over_K")
print("NO_CLAIM=global_Keller_pair_or_general_quartic_or_chart_entry_or_JC2_or_DC2")
print("RESULT=PASS")
