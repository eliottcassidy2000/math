#!/usr/bin/env python3
"""Small hostile controls. Rational results exact; polygon diagnostics float only."""
from fractions import Fraction as Q
from math import floor, cos, sin, atan2, pi
import json


def require(condition, label):
    if not condition:
        raise RuntimeError("verification failed: " + label)


def validated(v, s=None):
    v = list(v)
    if not v or any(not isinstance(x, int) or x <= 0 for x in v):
        raise ValueError("speeds must be nonempty positive integers")
    s = [Q(0)] * len(v) if s is None else [Q(x) % 1 for x in s]
    if len(v) != len(s):
        raise ValueError("one rational start is required per speed")
    return v, s


def norm(x):
    r = x % 1
    return min(r, 1-r)


def maximum_loneliness(v, s=None):
    """Exact PL lower envelope on [0,1] for positive integer v, rational s.

    Between half-integer crossings every distance is affine. Its lower
    envelope can change slope only at an endpoint or a pair intersection.
    Thus this finite candidate list is exhaustive, not a time-grid sample.
    """
    v, s = validated(v, s)
    cuts = {Q(0), Q(1)}
    for vi, si in zip(v, s):
        for k in range(-2, 2*(vi+2)+1):
            t = (Q(k, 2)-si)/vi
            if 0 <= t <= 1:
                cuts.add(t)
    cuts = sorted(cuts)
    candidates = set(cuts)
    for left, right in zip(cuts, cuts[1:]):
        mid = (left+right)/2
        lines = []
        for vi, si in zip(v, s):
            integer = floor(si+vi*mid)
            residue = si+vi*mid-integer
            lines.append((Q(vi), si-integer) if residue < Q(1, 2)
                         else (Q(-vi), integer+1-si))
        for i, (a, b) in enumerate(lines):
            for c, d in lines[i+1:]:
                if a != c:
                    t = (d-b)/(a-c)
                    if left <= t <= right:
                        candidates.add(t)
    evaluated = [(min(norm(si+vi*t) for vi, si in zip(v,s)), t)
                 for t in candidates]
    maximum = max(value for value, _ in evaluated)
    return maximum, sorted(t for value, t in evaluated if value == maximum)


def bad_cover(v, s, radius):
    """Independent interval-union proof of a global upper bound."""
    v, s = validated(v, s)
    radius = Q(radius)
    if radius < 0:
        raise ValueError("radius must be nonnegative")
    intervals = []
    for vi, si in zip(v, s):
        for k in range(-1, vi+2):
            left, right = (k-si-radius)/vi, (k-si+radius)/vi
            if left <= 1 and right >= 0:
                intervals.append((max(Q(0), left), min(Q(1), right)))
    rightmost = Q(0)
    for left, right in sorted(intervals):
        if left > rightmost:
            return False
        rightmost = max(rightmost, right)
    return rightmost == 1


def polygon(m, x):
    x %= 1
    j = floor(m*x)
    a = m*x-j
    u, v = 2*pi*j/m, 2*pi*(j+1)/m
    return ((1-a)*cos(u)+a*cos(v), (1-a)*sin(u)+a*sin(v))


def angle(m, x):
    u, v = polygon(m, x)
    return (atan2(v, u)/(2*pi)) % 1


def main():
    result = {"exact": {}, "floating_diagnostics": {}}
    v = [1,2,3,4,5]
    s = [Q(a,94) for a in [0,46,38,47,72]]
    value, times = maximum_loneliness(v,s)
    require(value == Q(15,94) < Q(1,6), 'value == Q(15,94) < Q(1,6)')
    require(times == [Q(63,188), Q(39,94), Q(45,94), Q(55,94)], "all shifted maximizing times")
    shifted = [a + (i-2)*13 for i,a in enumerate(s)]
    require(maximum_loneliness(v, shifted) == (value, times), "integer start shifts preserve full result")
    require(bad_cover(v, shifted, value), "integer start shifts preserve cover")
    require(bad_cover(v,s,value), 'bad_cover(v,s,value)')
    result["exact"]["shifted_hostile"] = {"value":str(value),
        "maximizing_times":list(map(str,times)), "independent_bad_cover":True}
    value, times = maximum_loneliness(v)
    require(value == Q(1,6), 'value == Q(1,6)')
    result["exact"]["unshifted_AP_control"] = {"value":str(value),
        "maximizing_times":list(map(str,times))}
    for velocities in ([1,3],[1,5],[2,4],[1,2,6],[1,3,5]):
        value, times = maximum_loneliness(velocities)
        result["exact"]["speeds_"+"_".join(map(str,velocities))] = {
            "value":str(value),"maximizing_times":list(map(str,times))}
    require(maximum_loneliness([1,3])[0] == Q(1,2), 'maximum_loneliness([1,3])[0] == Q(1,2)')
    require(maximum_loneliness([1,2,6])[0] == Q(2,7), "same-area first maximum")
    require(maximum_loneliness([1,3,5])[0] == Q(1,2), "same-area second maximum")
    # Triangle P(1/12)=(5/8,sqrt(3)/8), P(9/16)=(-1/2,-3sqrt(3)/16).
    det_sqrt3_coefficient = Q(5,8)*Q(-3,16)-Q(1,8)*Q(-1,2)
    require(Q(9,16)-Q(1,12) == Q(23,48) < Q(1,2), 'Q(9,16)-Q(1,12) == Q(23,48) < Q(1,2)')
    require(det_sqrt3_coefficient == Q(-7,128), 'det_sqrt3_coefficient == Q(-7,128)')
    result["exact"]["odd_triangle_orientation_hostile"] = {
        "phase_difference":"23/48", "determinant":"-7*sqrt(3)/128"}
    # A centered unit square: vertex-to-vertex and adjacent midpoint-to-midpoint
    # both have normalized perimeter distance 1/4 but distinct Euclidean chords.
    def squared_distance(p, q):
        return sum((a-b)**2 for a,b in zip(p,q))
    side_chord = squared_distance((Q(1,2),Q(1,2)),(Q(-1,2),Q(1,2)))
    midpoint_chord = squared_distance((Q(0),Q(1,2)),(Q(-1,2),Q(0)))
    require(side_chord == 1 and midpoint_chord == Q(1,2), 'side_chord == 1 and midpoint_chord == Q(1,2)')
    result["exact"]["square_chord_squared_at_same_phase_distance"] = [
        str(side_chord),str(midpoint_chord)]
    def edge_word(phases):
        return [(phases[j]-phases[i]) % 1 < Q(1,2)
                for i in range(len(phases)) for j in range(i+1,len(phases))]
    lonely = [Q(0),Q(2,5),Q(3,5)]
    not_lonely = [Q(0),Q(1,5),Q(3,5)]
    require(edge_word(lonely) == edge_word(not_lonely), 'edge_word(lonely) == edge_word(not_lonely)')
    require(min(norm(x) for x in lonely[1:]) == Q(2,5) >= Q(1,3), 'min(norm(x) for x in lonely[1:]) == Q(2,5) >= Q(1,3)')
    require(min(norm(x) for x in not_lonely[1:]) == Q(1,5) < Q(1,3), 'min(norm(x) for x in not_lonely[1:]) == Q(1,5) < Q(1,3)')
    result["exact"]["same_labelled_tournament_different_loneliness"] = {
        "lonely_phases":list(map(str,lonely)),
        "not_lonely_phases":list(map(str,not_lonely)),
        "edge_word":edge_word(lonely)}
    count = boundary_count = orientation_count = 0
    for N in range(3,21):
        for m in [N,2*N]:
            for a in range(101):
                x = (a+0.173)/101
                d = (angle(m,x+1/N)-angle(m,x)) % 1
                require(abs(d-1/N) < 1e-12, 'abs(d-1/N) < 1e-12')
                boundary_count += 1
                for b in [1,17,39,72,99]:
                    y = (b+0.419)/101
                    dx, dy = float(norm(x-y)), float(norm(angle(m,x)-angle(m,y)))
                    if min(abs(dx-1/N),abs(dy-1/N)) > 1e-10:
                        require((dx >= 1/N) == (dy >= 1/N), '(dx >= 1/N) == (dy >= 1/N)')
                        count += 1
                    if m % 2 == 0:
                        p, q = polygon(m,x), polygon(m,y)
                        det = p[0]*q[1]-p[1]*q[0]
                        phase_sine = sin(2*pi*(y-x))
                        if abs(det) > 1e-10 and abs(phase_sine) > 1e-10:
                            require((det > 0) == (phase_sine > 0), '(det > 0) == (phase_sine > 0)')
                            orientation_count += 1
    result["floating_diagnostics"] = {"N_range":[3,20],"m_choices":["N","2N"],
        "boundary_checks":boundary_count,"threshold_checks":count,
        "even_polygon_orientation_checks":orientation_count,
        "status":"sanity checks only; general statements proved in report"}
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
