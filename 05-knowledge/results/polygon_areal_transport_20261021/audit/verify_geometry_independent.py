#!/usr/bin/env python3
"""Independent exact, assertion-free audit controls; standard library only.

Does not import the producer. The loneliness certificates use exact interval
intersection/cover, not the producer's lower-envelope candidate enumeration.
Hexagon checks use coordinates (a, b*sqrt(3)) over rational coefficients.
General transport theorems are proved in the accompanying proof audit.
"""
from fractions import Fraction as Q
from math import gcd
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def norm(t):
    t %= 1
    return min(t, 1-t)


def normalize(v, starts):
    require(bool(v) and all(type(x) is int and x > 0 for x in v), "positive integer speeds required")
    require(len(v) == len(starts), "speed/start lengths differ")
    return [Q(s) % 1 for s in starts]


def intersect_unions(a, b):
    intersections = sorted((max(x,u), min(y,w))
        for x,y in a for u,w in b if max(x,u) <= min(y,w))
    result = []
    for left, right in intersections:
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1],right))
        else:
            result.append((left,right))
    return result


def safe_times(v, starts, radius):
    """All t in [0,1] where each circular distance >= radius, endpoints kept."""
    starts = normalize(v, starts)
    require(0 <= radius <= Q(1,2), "radius out of range")
    result = [(Q(0),Q(1))]
    for speed, s in zip(v, starts):
        # s+speed*t lies in [s,s+speed] subset [0,speed+1).
        individual = []
        for integer in range(speed+1):
            left = max(Q(0), (integer+radius-s)/speed)
            right = min(Q(1), (integer+1-radius-s)/speed)
            if left <= right:
                individual.append((left,right))
        result = intersect_unions(result, individual)
    return result


def upper_cover(v, starts, radius):
    """A certificate that some distance <= radius at every t in [0,1]."""
    starts = normalize(v, starts)
    intervals = []
    for runner, (speed,s) in enumerate(zip(v, starts)):
        for integer in range(speed+2):
            left = max(Q(0), (integer-radius-s)/speed)
            right = min(Q(1), (integer+radius-s)/speed)
            if left <= right:
                intervals.append((left,right,runner,integer))
    reach = Q(0)
    certificate = []
    for left,right,runner,integer in sorted(intervals):
        if left > reach:
            return None
        if right > reach:
            certificate.append({"left":str(left),"right":str(right),
                                "runner_index":runner,"integer":integer})
            reach = right
        if reach == 1:
            return certificate
    return None


def certify_maximum(v, starts, radius, expected_times):
    cover = upper_cover(v, starts, radius)
    require(cover is not None, "closed bad intervals do not cover [0,1]")
    allowed = safe_times(v, starts, radius)
    target = [(t,t) for t in expected_times]
    require(allowed == target, f"maximizer set mismatch: {allowed} != {target}")
    require(all(min(norm(speed*t+s) for speed,s in zip(v,starts)) == radius
                for t in expected_times), "claimed maximizer has wrong value")
    translated_starts = [s+(i+1)*137 for i,s in enumerate(starts)]
    require(safe_times(v, translated_starts, radius) == target,
            "integer shifts of starts changed safe set")
    return {"speeds":v,"starts":list(map(str,starts)),"maximum":str(radius),
            "all_maximizers":list(map(str,expected_times)),
            "closed_bad_cover":cover,"endpoint_convention":"closed"}


def hexagon(x):
    """Coordinates are (x, y/sqrt(3)); vertices lie on the unit circle."""
    vertices = [(Q(1),Q(0)),(Q(1,2),Q(1,2)),(Q(-1,2),Q(1,2)),
                (Q(-1),Q(0)),(Q(-1,2),Q(-1,2)),(Q(1,2),Q(-1,2))]
    position = 6*(x % 1)
    edge = int(position)
    amount = position-edge
    p, q = vertices[edge], vertices[(edge+1)%6]
    return tuple((1-amount)*a+amount*b for a,b in zip(p,q))


def det_coeff(p,q):
    return p[0]*q[1]-p[1]*q[0]


def dot(p,q):
    return p[0]*q[0]+3*p[1]*q[1]


def sign(x):
    return (x > 0)-(x < 0)


def phase_orientation(x,y):
    gap = (y-x) % 1
    return 0 if gap in (0,Q(1,2)) else (1 if gap < Q(1,2) else -1)


def main():
    results = {}
    cases = [
        ([1,2,3,4,5],[Q(a,94) for a in (0,46,38,47,72)],Q(15,94),
         [Q(63,188),Q(39,94),Q(45,94),Q(55,94)]),
        ([1,2,3,4,5],[Q(0)]*5,Q(1,6),[Q(1,6),Q(5,6)]),
        ([1,3],[Q(0)]*2,Q(1,2),[Q(1,2)]),
        ([1,5],[Q(0)]*2,Q(1,2),[Q(1,2)]),
        ([2,4],[Q(0)]*2,Q(1,3),[Q(1,6),Q(1,3),Q(2,3),Q(5,6)]),
        ([1,2,6],[Q(0)]*3,Q(2,7),[Q(2,7),Q(5,7)]),
        ([1,3,5],[Q(0)]*3,Q(1,2),[Q(1,2)]),
    ]
    results["exact_loneliness_certificates"] = [certify_maximum(*case) for case in cases]
    require(Q(15,94) < Q(1,6), "shifted hostile is not strict")
    require(safe_times(*cases[0][:2],Q(1,6)) == [], "shifted hostile has an LRC witness")
    # Arithmetic-progression boundary controls, including N=2 and composite N.
    ap_counts = {}
    for n in range(2,31):
        wanted = [(Q(a,n),Q(a,n)) for a in range(1,n) if gcd(a,n) == 1]
        require(safe_times(list(range(1,n)),[Q(0)]*(n-1),Q(1,n)) == wanted,
                f"primitive AP endpoint control failed for N={n}")
        ap_counts[str(n)] = len(wanted)
    results["AP_closed_lonely_set_counts_N_2_to_30"] = ap_counts
    # Exact rational/radical-coordinate hexagon orientation and radial unit edges.
    phases = [Q(a,48) for a in range(48)]
    checks = 0
    for x in phases:
        for y in phases:
            p,q = hexagon(x),hexagon(y)
            require(sign(det_coeff(p,q)) == phase_orientation(x,y), "hexagon orientation failure")
            radial_unit_edge = dot(p,q) > 0 and 4*dot(p,q)**2 == dot(p,p)*dot(q,q)
            require(radial_unit_edge == (norm(y-x) == Q(1,6)), "radial unit-distance failure")
            # N=2 angular safe iff opposite; N=3 angular safe iff cosine <= -1/2.
            polar_antipodes = det_coeff(p,q) == 0 and dot(p,q) < 0
            require(polar_antipodes == (norm(y-x) >= Q(1,2)), "delta=1/2 endpoint failure")
            polar_third_safe = dot(p,q) < 0 and 4*dot(p,q)**2 >= dot(p,p)*dot(q,q)
            require(polar_third_safe == (norm(y-x) >= Q(1,3)), "hexagon third-threshold failure")
            polar_sixth_safe = dot(p,q) <= 0 or 4*dot(p,q)**2 <= dot(p,p)*dot(q,q)
            require(polar_sixth_safe == (norm(y-x) >= Q(1,6)), "hexagon sixth-threshold failure")
            checks += 1
    results["hexagon_exact_rational_pair_checks"] = checks
    # Triangle hostile uses independent perimeter interpolation in Q(sqrt(3)).
    triangle_vertices = [(Q(1),Q(0)),(Q(-1,2),Q(1,2)),(Q(-1,2),Q(-1,2))]
    def triangle(t):
        pos = 3*(t % 1); edge = int(pos); amount = pos-edge
        return tuple((1-amount)*a+amount*b
            for a,b in zip(triangle_vertices[edge],triangle_vertices[(edge+1)%3]))
    x,y = Q(1,12),Q(9,16)
    p,q = triangle(x),triangle(y)
    require(p == (Q(5,8),Q(1,8)) and q == (Q(-1,2),Q(-3,16)), "triangle coordinates")
    require(phase_orientation(x,y) == 1 and det_coeff(p,q) == Q(-7,128), "triangle hostile")
    results["triangle_hostile"] = {"phase_gap":"23/48","det_sqrt3_coefficient":"-7/128"}
    # Same labelled tournament with different observer loneliness.
    first,second = [Q(0),Q(2,5),Q(3,5)],[Q(0),Q(1,5),Q(3,5)]
    word = lambda p: [phase_orientation(p[i],p[j]) for i in range(3) for j in range(i+1,3)]
    require(word(first) == word(second) == [1,-1,1], "labelled tournament mismatch")
    require(min(map(norm,first[1:])) == Q(2,5) >= Q(1,3) > min(map(norm,second[1:])), "observer hostile")
    results["same_tournament_hostile"] = {"orientation_word":[1,-1,1],"observer_slacks":["1/15","-2/15"]}
    # Unit-side square: same 1/4 perimeter phase, different physical chord lengths.
    length2 = lambda p,q: sum((a-b)**2 for a,b in zip(p,q))
    a = length2((Q(1,2),Q(1,2)),(Q(-1,2),Q(1,2)))
    b = length2((Q(0),Q(1,2)),(Q(-1,2),Q(0)))
    require((a,b) == (Q(1),Q(1,2)), "square chord hostile")
    results["square_chord_squared"] = list(map(str,(a,b)))
    # Exact focus-ellipse witness a=5,e=3/5,b=4 has positive angular determinant.
    require((-3)*(-4)-4*(-3) == 24, "ellipse orientation determinant")
    results["focus_ellipse_literal_control"] = {"a":5,"b":4,"e":"3/5","focus_vectors":[[-3,4],[-3,-4]],"determinant":24,
        "mean_forward_gap":"pi+6/5 > pi","true_forward_gap":"2*acos(3/5) < pi"}
    # Explicitly test that validation cannot be disabled by Python optimization.
    try:
        require(False,"sentinel")
    except RuntimeError as exc:
        require(str(exc) == "sentinel", "sentinel exception changed")
    else:
        raise RuntimeError("require(False) was not enforced")
    results["validation"] = "PASS; explicit checks remain active under python -O"
    results["scope"] = "Exact literal controls only; no general LRC proof, no shifted-start optimization, no preprint computation audit."
    print(json.dumps(results,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
