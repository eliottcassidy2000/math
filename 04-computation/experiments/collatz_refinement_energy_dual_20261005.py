"""Arithmetic-shell moment duals; no orbit, ROOT receipt, or weight bank input.

Moment intervals are independent conditional premises, not authenticated by
this finite checker.  The law must live on the stated integer-distance shell
set.  Normal and optimized execution run the same exact checks.
"""
from fractions import Fraction as F
from math import prod
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


def natural(value):
    need(type(value) is int and value >= 0, "exact natural required")


def rational(value):
    need(type(value) in (int, F), "exact rational required")
    return F(value)


def shell(distance):
    natural(distance)
    t = 1 << distance
    return F(4*t, (t+1)**2)


def kernel(source_index, atom_index):
    natural(source_index); natural(atom_index)
    return shell(abs(source_index-atom_index))


def odd_degree(degree):
    need(type(degree) is int and degree > 0 and degree % 2 == 1,
         "positive odd exact degree required")


def multiply(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)


def minorant(degree):
    odd_degree(degree)
    q = (F(1),)
    for d in range(1, degree+1):
        h = shell(d)
        q = multiply(q, (-h/(1-h), 1/(1-h)))
    return q


def evaluate(poly, x):
    x = rational(x)
    out = F(0)
    for c in reversed(poly):
        out = out*x+c
    return out


def error_bound(degree):
    odd_degree(degree)
    return prod((shell(d)/(1-shell(d)) for d in range(1, degree+1)), start=F(1))


def coefficient_norm(degree):
    return sum(map(abs, minorant(degree)), F(0))


def moment_packet(intervals, degree):
    need(type(intervals) is tuple and len(intervals) == degree+1,
         "one exact interval per moment required")
    out = []
    for pair in intervals:
        need(type(pair) is tuple and len(pair) == 2, "exact interval pair required")
        lo, hi = map(rational, pair)
        need(0 <= lo <= hi <= 1, "bounded ordered moment interval required")
        out.append((lo, hi))
    need(out[0] == (F(1), F(1)), "moment zero must equal one")
    return tuple(out)


def certify(source_index, degree, intervals, tail_cap=F(1)):
    """Conditional atom interval; passing does not validate moment provenance.

    tail_cap bounds the mass at shell distances strictly larger than degree.
    This is an additional independently justified premise if less than one.
    """
    natural(source_index); odd_degree(degree)
    bounds = moment_packet(intervals, degree)
    tail_cap = rational(tail_cap)
    need(0 <= tail_cap <= 1, "tail cap must lie in [0,1]")
    q = minorant(degree)
    lo = sum((c*(a if c >= 0 else b) for c, (a, b) in zip(q, bounds)), F(0))
    hi = sum((c*(b if c >= 0 else a) for c, (a, b) in zip(q, bounds)), F(0))
    lower, upper = max(F(0), lo), min(F(1), hi+error_bound(degree)*tail_cap)
    need(lower <= upper, "packet fails necessary arithmetic-shell condition")
    return {"source_index": source_index, "degree": degree,
            "raw_lower": lo, "raw_upper": hi, "lower": lower, "upper": upper,
            "coefficient_norm": sum(map(abs, q), F(0)), "tail_cap": tail_cap}


def refinement(source_index, packets):
    """Intersect repeated moment intervals and retain all previous atom bounds.

    Input entries are (odd_degree, intervals, independently justified tail_cap).
    Degrees increase from at least three.  This checks necessary consistency,
    not the existence of one actual law or the truth of any measurement.
    """
    natural(source_index)
    need(type(packets) is tuple and len(packets) > 0, "nonempty tuple required")
    retained = ()
    previous_degree = 1
    lower, upper = F(0), F(1)
    history = []
    for entry in packets:
        need(type(entry) is tuple and len(entry) == 3, "three-field packet required")
        degree, intervals, cap = entry
        odd_degree(degree)
        need(degree > previous_degree, "strictly increasing degrees from three required")
        bounds = list(moment_packet(intervals, degree))
        for k, (old_lo, old_hi) in enumerate(retained):
            lo, hi = bounds[k]
            bounds[k] = (max(lo, old_lo), min(hi, old_hi))
            need(bounds[k][0] <= bounds[k][1], "inconsistent repeated moment interval")
        retained = tuple(bounds)
        result = certify(source_index, degree, retained, cap)
        lower, upper = max(lower, result["lower"]), min(upper, result["upper"])
        need(lower <= upper, "inconsistent retained atom bounds")
        history.append((degree, lower, upper))
        previous_degree = degree
    return tuple(history)


def exact_moments(law, degree):
    return tuple(sum((mass*x**k for x, mass in law), F(0)) for k in range(degree+1))


def exact_intervals(moments):
    return tuple((x, x) for x in moments)


def conditional_degree(floor):
    """Least odd degree >=3 with the uniform tail error <= floor/4.

    This does not assert that the selected source has the supplied floor.
    """
    floor = rational(floor)
    need(0 < floor <= 1, "positive conditional floor at most one required")
    degree = 3
    while error_bound(degree) > floor/4:
        degree += 2
    return degree


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for x in range(total+1):
            for rest in compositions(total-x, parts-1):
                yield (x,)+rest


def main():
    checks = 0

    def check(condition, label):
        nonlocal checks
        need(condition, label)
        checks += 1

    # Independent direct products versus coefficient expansion, with all
    # source indices 0..6 and atom indices 0..40, odd degrees 1..17.
    for degree in range(1, 18, 2):
        q = minorant(degree)
        eps = error_bound(degree)
        check(evaluate(q, F(1)) == 1, "target normalization")
        check(evaluate(q, F(0)) == -eps, "sharp support-closure tail")
        check(coefficient_norm(degree) == prod(
            ((1+shell(d))/(1-shell(d)) for d in range(1, degree+1)), start=F(1)),
            "alternating coefficient norm")
        for m in range(7):
            for j in range(41):
                x = kernel(m, j)
                direct = prod(((x-shell(d))/(1-shell(d))
                               for d in range(1, degree+1)), start=F(1))
                value = evaluate(q, x)
                check(value == direct, "independent product evaluation")
                check(value <= int(m == j) <= value+eps, "pointwise atom bracket")
                if degree >= 3:
                    check(value <= evaluate(minorant(degree+2), x), "pointwise refinement")

    c5 = coefficient_norm(5)
    check(c5 == F(33835804361, 95355225), "first-five norm")
    check(c5*F(15, 11) == F(33835804361, 69927165) < 512, "uniform norm ceiling")
    for degree in range(1, 42, 2):
        exponent = -degree*(degree-3)//2
        power = F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))
        check(error_bound(degree) <= F(1024, 81)*power, "quadratic exponent error")
        check(coefficient_norm(degree) < 512, "finite norm controls")

    # Every denominator-four law on target and first four shells: 70 laws.
    nodes = (F(1),)+tuple(shell(d) for d in range(1, 5))
    laws = 0
    for numerators in compositions(4, 5):
        laws += 1
        law = tuple((x, F(n, 4)) for x, n in zip(nodes, numerators))
        moments = exact_moments(law, 9)
        previous = None
        for degree in (3, 5, 7, 9):
            q = minorant(degree)
            value = sum((c*h for c, h in zip(q, moments)), F(0))
            check(value <= law[0][1] <= value+error_bound(degree), "law bracket")
            if previous is not None:
                check(previous <= value, "law monotone refinement")
            previous = value
        packets = tuple((d, exact_intervals(moments[:d+1]), F(1)) for d in (3, 5, 7, 9))
        history = refinement(0, packets)
        check(all(lo <= law[0][1] <= hi for _, lo, hi in history), "retained law bounds")

    # The arithmetic support, not an extra measured moment, supplies the gain.
    law = ((F(1), F(1, 1000)),)+tuple((shell(d), F(333, 1000)) for d in range(1, 4))
    moments = exact_moments(law, 3)
    result = certify(0, 3, exact_intervals(moments))
    check(result["raw_lower"] == F(1, 1000), "strict arithmetic-support gain")
    c = shell(1)
    matrix = tuple(tuple(c*moments[i+j]-moments[i+j+1] for j in range(2)) for i in range(2))
    det = matrix[0][0]*matrix[1][1]-matrix[0][1]**2
    check(matrix[0][0] > 0 and det > 0, "old interval localizer positive definite")
    tangent = matrix[1][1]-2*matrix[0][1]+matrix[0][0]
    linear = matrix[0][1]-matrix[0][0]
    old_floor = -(matrix[0][0]-linear**2/tangent)/(1-c)
    check(old_floor == -F(811338319, 7894961000), "old normalized optimum")
    ghost = tuple(zip((shell(3), F(1, 2), shell(2), shell(1)),
                      (F(111782799, 337280000), F(9, 2125),
                       F(1305719, 3968000), F(107289, 320000))))
    check(all(p > 0 for _, p in ghost) and sum(p for _, p in ghost) == 1,
          "zero-target continuous completion")
    check(exact_moments(ghost, 3) == moments, "same four moments")
    check(shell(3) < F(1, 2) < shell(2), "impossible integer shell")
    check(evaluate(minorant(3), F(1, 2)) > 0, "support hypothesis is essential")

    # Exact interval errors, and monotone retained bounds under independent
    # nested interval refinements, for a law with a genuine off-head tail.
    tail_law = ((F(1), F(1, 10)), (shell(1), F(2, 5)), (shell(4), F(1, 2)))
    tail_moments = exact_moments(tail_law, 13)
    history_packets = []
    for degree in (3, 5, 7, 9, 11, 13):
        width = F(1, 10**degree)
        intervals = ((F(1), F(1)),)+tuple((max(F(0), h-width/2), min(F(1), h+width/2))
                                                    for h in tail_moments[1:degree+1])
        result_i = certify(0, degree, intervals)
        true_readout = sum((a*b for a, b in zip(minorant(degree), tail_moments)), F(0))
        check(true_readout-result_i["raw_lower"] <= coefficient_norm(degree)*width,
              "joint interval cost")
        check(result_i["lower"] <= F(1, 10) <= result_i["upper"], "interval soundness")
        history_packets.append((degree, intervals, F(1)))
    history = refinement(0, tuple(history_packets))
    check(all(history[i][1] <= history[i+1][1] and history[i][2] >= history[i+1][2]
              for i in range(len(history)-1)), "monotone retained intervals")

    # Energy of minorants cannot certify an atom: a law at shell four has
    # target mass zero, finite total squared increments, and zero thereafter.
    missing = shell(4)
    q3_missing = evaluate(minorant(3), missing)
    check(q3_missing < 0 and evaluate(minorant(5), missing) == 0, "missing atom refinement")
    energy = sum((evaluate(minorant(d+2), missing)-evaluate(minorant(d), missing))**2
                 for d in range(3, 16, 2))
    check(energy == q3_missing**2 <= error_bound(3)**2, "bounded energy missing-atom hostile")
    check(evaluate(minorant(3), shell(10)) < evaluate(minorant(1), shell(10)),
          "degree-one monotonicity hostile")

    conditional_degrees = []
    for k in (1, 4, 10, 20, 40, 80, 160):
        eta = F(1, 2**k)
        d = conditional_degree(eta)
        check(error_bound(d) <= eta/4, "conditional degree guarantee")
        if d > 3:
            check(error_bound(d-2) > eta/4, "least degree convention")
        conditional_degrees.append((k, d))

    # Malformed inputs are rejected even under -O; no empty-word alias.
    hostiles = (
        lambda: shell(True), lambda: shell(-1), lambda: kernel(0, 1.0),
        lambda: minorant(0), lambda: minorant(2), lambda: minorant(True),
        lambda: certify(True, 3, exact_intervals(moments)),
        lambda: certify(0, 3, exact_intervals(moments), 1.0),
        lambda: certify(0, 3, exact_intervals(moments), -1),
        lambda: certify(0, 3, ((F(0), F(1)),)+exact_intervals(moments)[1:]),
        lambda: conditional_degree(0), lambda: conditional_degree(F(2)),
        lambda: refinement(0, ((3, exact_intervals(moments), 1),
                              (3, exact_intervals(moments), 1))),
        lambda: refinement(0, ((3, exact_intervals(moments), 1),
                              (5, ((F(1), F(1)), (F(0), F(0)))+
                               ((F(0), F(1)),)*4, 1))),
    )
    for hostile in hostiles:
        try:
            hostile()
        except ValueError:
            check(True, "invalid input rejected")
        else:
            check(False, "invalid input accepted")

    print(json.dumps({
        "status": "PASS; conditional arithmetic-shell certificates, no ROOT inputs",
        "checks": checks, "finite_laws": laws,
        "source_indices": "0..6", "atom_indices": "0..40", "degree_controls": "odd1..17",
        "uniform_direct_moment_norm_bound": "<512",
        "sharper_norm_ceiling": str(c5*F(15, 11)),
        "same_moment_arithmetic_floor": str(result["raw_lower"]),
        "same_moment_interval_optimizer_floor": str(old_floor),
        "continuous_zero_target_completion": [[str(x), str(p)] for x, p in ghost],
        "conditional_degrees_for_eta_2_minus_k": conditional_degrees,
        "refinement_history": [[d, str(lo), str(hi)] for d, lo, hi in history],
        "missing_atom_energy": str(energy), "type_hostiles": len(hostiles),
        "root_orbit_or_certificate_inputs": 0,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
