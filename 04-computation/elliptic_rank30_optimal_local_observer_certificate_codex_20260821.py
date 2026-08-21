#!/usr/bin/env python3
"""Independent finite-reduction certificate for ICARM curve 273.

The public curve and its thirty public rational points are frozen below.  A
bank of good-prime maps to E(F_p)/2E(F_p) has rank thirty.  Two explicit
finite-group realizations (the displayed long model and an integral
odd-prime change of variables) are checked.  The result is an unconditional
rank lower bound, never an exact-rank or construction claim.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import itertools
import json
from math import gcd, isqrt


A4 = -201769035260418549083594900060734240952308696994802735114305555
A6 = int(
    "1151107939141058565733479426024323225135665982951300586808823640527729578307228357301072889377"
)
DISCRIMINANT = int(
    "-46714661255308767314567688733841531918983356002159772613256840842851650254036518701100578342601553513579222272710220496887616034526983492843954090554197033638137245037791044053017600000000"
)

POINTS_RAW = (
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
)
POINTS = tuple((Fraction(x), Fraction(y)) for x, y in POINTS_RAW)

EXPECTED_PRIMES = (43, 61, 101, 211, 223, 241, 263, 271, 283, 311,
                   313, 457, 521, 569, 571)
EXPECTED_ORDERS = (56, 72, 120, 212, 248, 260, 292, 296, 304, 344,
                   336, 464, 516, 600, 604)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    return all(value % divisor for divisor in range(3, isqrt(value) + 1, 2))


def on_curve(point, modulus: int, coefficients: tuple[int, ...]) -> bool:
    if point is None:
        return True
    a1, a2, a3, a4, a6 = coefficients
    x_value, y_value = point
    return (
        y_value * y_value + a1 * x_value * y_value + a3 * y_value
        - x_value**3 - a2 * x_value**2 - a4 * x_value - a6
    ) % modulus == 0


def add(point, other, modulus: int, coefficients: tuple[int, ...]):
    if point is None:
        return other
    if other is None:
        return point
    a1, a2, a3, a4, a6 = coefficients
    x1, y1 = point
    x2, y2 = other
    if x1 == x2 and (y1 + y2 + a1 * x1 + a3) % modulus == 0:
        return None
    if point != other:
        denominator = (x2 - x1) % modulus
        inverse = pow(denominator, -1, modulus)
        slope = (y2 - y1) * inverse % modulus
        intercept = (y1 * x2 - y2 * x1) * inverse % modulus
    else:
        denominator = (2 * y1 + a1 * x1 + a3) % modulus
        if denominator == 0:
            return None
        inverse = pow(denominator, -1, modulus)
        slope = (3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1) * inverse % modulus
        intercept = (-x1**3 + a4 * x1 + 2 * a6 - a3 * y1) * inverse % modulus
    x3 = (slope * slope + a1 * slope - a2 - x1 - x2) % modulus
    y3 = (-(slope + a1) * x3 - intercept - a3) % modulus
    require(on_curve((x3, y3), modulus, coefficients), (modulus, "group law"))
    return x3, y3


def finite_points(modulus: int, coefficients: tuple[int, ...]) -> tuple[object, ...]:
    values = [None]
    for x_value in range(modulus):
        for y_value in range(modulus):
            if on_curve((x_value, y_value), modulus, coefficients):
                values.append((x_value, y_value))
    return tuple(values)


def reduce_rational(point: tuple[Fraction, Fraction], modulus: int):
    x_value, y_value = point
    # At good reduction, a point with negative x-valuation reduces to O.
    if x_value.denominator % modulus == 0:
        require(y_value.denominator % modulus == 0,
                (modulus, "compatible nonintegral reduction"))
        return None
    require(y_value.denominator % modulus != 0,
            (modulus, "affine denominator"))
    return (
        x_value.numerator * pow(x_value.denominator, -1, modulus) % modulus,
        y_value.numerator * pow(y_value.denominator, -1, modulus) % modulus,
    )


def quotient_data(modulus: int, coefficients: tuple[int, ...],
                  rational_points=POINTS):
    group = finite_points(modulus, coefficients)
    doubles = {add(item, item, modulus, coefficients) for item in group}
    require(len(group) % len(doubles) == 0, (modulus, "double subgroup"))
    quotient_size = len(group) // len(doubles)
    require(quotient_size in (1, 2, 4), (modulus, "elliptic 2-quotient size"))
    if quotient_size != 4:
        return group, doubles, quotient_size, (), (), ()

    first = next(item for item in group if item not in doubles)
    first_coset = {add(first, item, modulus, coefficients) for item in doubles}
    second = next(item for item in group
                  if item not in doubles and item not in first_coset)
    second_coset = {add(second, item, modulus, coefficients) for item in doubles}
    rows = [0, 0]
    point_bits = []
    for index, rational_point in enumerate(rational_points):
        reduced = reduce_rational(rational_point, modulus)
        require(reduced in group, (modulus, index, "point reduction"))
        if reduced in doubles:
            bits = (0, 0)
        elif reduced in first_coset:
            bits = (1, 0)
        elif reduced in second_coset:
            bits = (0, 1)
        else:
            bits = (1, 1)
        point_bits.append(bits)
        rows[0] |= bits[0] << index
        rows[1] |= bits[1] << index
    return group, doubles, quotient_size, (first, second), tuple(rows), tuple(point_bits)


def rank_two(rows: tuple[int, ...] | list[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = row
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def transformed(point, modulus: int):
    if point is None:
        return None
    x_value, y_value = point
    return 4 * x_value % modulus, (8 * y_value + 4 * x_value) % modulus


def transformed_crosscheck(modulus: int, original) -> None:
    coefficients_b = (0, 1, 0, 16 * A4, 64 * A6)
    group_a, doubles_a, _, basis_a, rows_a, point_bits_a = original
    group_b = finite_points(modulus, coefficients_b)
    doubles_b = {add(item, item, modulus, coefficients_b) for item in group_b}
    require({transformed(item, modulus) for item in group_a} == set(group_b),
            (modulus, "integral model isomorphism"))
    require({transformed(item, modulus) for item in doubles_a} == doubles_b,
            (modulus, "double image isomorphism"))
    first_b, second_b = (transformed(item, modulus) for item in basis_a)
    first_coset_b = {add(first_b, item, modulus, coefficients_b) for item in doubles_b}
    second_coset_b = {add(second_b, item, modulus, coefficients_b) for item in doubles_b}
    rows_b = [0, 0]
    bits_b = []
    for index, rational_point in enumerate(POINTS):
        reduced_a = reduce_rational(rational_point, modulus)
        reduced_b = transformed(reduced_a, modulus)
        if reduced_b in doubles_b:
            bits = (0, 0)
        elif reduced_b in first_coset_b:
            bits = (1, 0)
        elif reduced_b in second_coset_b:
            bits = (0, 1)
        else:
            bits = (1, 1)
        bits_b.append(bits)
        rows_b[0] |= bits[0] << index
        rows_b[1] |= bits[1] << index
    require(tuple(rows_b) == rows_a and tuple(bits_b) == point_bits_a,
            (modulus, "two finite-group paths"))


def main() -> None:
    require(len(POINTS) == 30, "point count")
    for index, (x_value, y_value) in enumerate(POINTS):
        require(y_value**2 + x_value * y_value == x_value**3 + A4 * x_value + A6,
                (index, "point incidence"))

    b2 = 1
    b4 = 2 * A4
    b6 = 4 * A6
    b8 = A6 - A4**2
    delta = -(b2**2) * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    # The public ICARM record labels this the minimal discriminant.  Here we
    # independently check the invariant formula, not global minimality.
    require(delta == DISCRIMINANT, "stored discriminant invariant")

    coefficients_a = (1, 0, 0, A4, A6)
    selected = []
    selected_rows: list[int] = []
    selected_records = []
    all_prime_records = []
    for modulus in range(3, 572, 2):
        if not prime(modulus) or DISCRIMINANT % modulus == 0:
            continue
        data = quotient_data(modulus, coefficients_a)
        group, _, quotient_size, _, rows, _ = data
        before = rank_two(selected_rows)
        after = rank_two(selected_rows + list(rows)) if rows else before
        gain = after - before
        all_prime_records.append((modulus, len(group), quotient_size, gain))
        if quotient_size == 4 and gain == 2:
            selected.append(modulus)
            selected_rows.extend(rows)
            selected_records.append((modulus, len(group), before, after,
                                     tuple(f"{row:08x}" for row in rows)))
            transformed_crosscheck(modulus, data)
            if len(selected_rows) == 30:
                break

    require(tuple(selected) == EXPECTED_PRIMES, "optimal prime bank")
    require(tuple(item[1] for item in selected_records) == EXPECTED_ORDERS,
            "finite group orders")
    require(rank_two(selected_rows) == 30, "stacked quotient rank")
    require(len(selected) == 15, "certificate cardinality")

    # Every E(F_p)/2E(F_p) has F_2-dimension at most two.  Thus fifteen local
    # blocks are an information-theoretic lower bound for thirty columns.
    # Since the increasing scan accepts a rank-two block whenever possible
    # and reaches the lower bound, this is the lexicographically least
    # cardinality-minimal certificate in the scanned odd-good-prime universe.
    require(all(item[3] - item[2] == 2 for item in selected_records),
            "two dimensions per selected prime")

    torsion_group_23 = finite_points(23, coefficients_a)
    require(DISCRIMINANT % 23 != 0 and len(torsion_group_23) == 33,
            "odd good-reduction torsion control")

    column_masks = tuple(
        sum(((row >> column) & 1) << row_index
            for row_index, row in enumerate(selected_rows))
        for column in range(30)
    )
    require(rank_two(list(column_masks)) == 30, "column rank transpose")

    semantic = (
        A4,
        A6,
        DISCRIMINANT,
        POINTS_RAW,
        tuple(all_prime_records),
        tuple(selected_records),
        tuple(f"{mask:08x}" for mask in column_masks),
        len(torsion_group_23),
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("ELLIPTIC_RANK30_OPTIMAL_LOCAL_OBSERVER_CERTIFICATE_20260821")
    print("status=PROVED_RANK_LOWER_BOUND_30+FINITE_EXACT_CERTIFICATE;exact_rank_open_unconditionally")
    print("model=y^2+x*y=x^3+A4*x+A6;public_points=30;incidence=PASS")
    print(f"stored_minimal_discriminant={DISCRIMINANT};"
          "invariant_formula=PASS;minimality_not_checked")
    print("prime_bank=" + repr(EXPECTED_PRIMES).replace(" ", ""))
    print("group_orders=" + repr(EXPECTED_ORDERS).replace(" ", ""))
    print("rank_trace=" + repr(tuple(selected_records)).replace(" ", ""))
    print("local_dimension_cap=2;cardinality_lower_bound=15;achieved=15")
    print("lexicographically_least_minimal_bank_in_odd_good_primes_through_571=PASS")
    print("integral_model_crosscheck=W^2=X^3+X^2+16*A4*X+64*A6;X=4x;W=8y+4x;PASS")
    print("torsion_sidecar=good_p23_group_order_33;E(Q)[2]=0")
    print("column_masks=" + repr(tuple(f"{mask:08x}" for mask in column_masks)).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print("scope=rank_at_least_30_only;no_rank_upper_bound;no_curve_generation_mechanism")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
