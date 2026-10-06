"""Exact controls for Collatz Poisson duals with controlled unbounded growth.

The all-source operator estimates are proved in the companion note. Finite
rows below independently check arithmetic, rational radical enclosures, and
the declared growth boundaries; they are not a coverage census.
"""
from fractions import Fraction as F
import json

CHECKS = 0
KP = F(3280, 3367)
KL = F(641, 686)


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)


def v2(n):
    return (n & -n).bit_length() - 1


def is_base(n):
    return type(n) is int and n > 0 and n % 2 == 1 and v2(3*n+1) in (1, 2)


def decompose(n):
    depth = 0
    while n % 8 == 5:
        n = (n-1)//4
        depth += 1
    return n, depth


def parent(b):
    if not is_base(b) or b == 1:
        raise ValueError("nonroot base required")
    u = (3*b+1) >> v2(3*b+1)
    return decompose(u)


def children(c, count):
    if not is_base(c):
        raise ValueError("base required")
    y = c
    out = []
    for k in range(count):
        if y % 3:
            a = 2 if y % 3 == 1 else 1
            b = ((1 << a)*y-1)//3
            if b != 1:
                out.append((b, k, a, F(1, 1 << (k+1))))
        y = 4*y+1
    return out


def radical_interval(x, degree, bits=36):
    """Dyadic lower/upper enclosure, certified by integer powers only."""
    x = F(x)
    if x < 0 or type(degree) is not int or degree < 1:
        raise ValueError("nonnegative rational and positive integral degree")
    scale = 1 << bits
    target = x.numerator*scale**degree
    lo, hi = 0, scale
    while hi**degree*x.denominator <= target:
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**degree*x.denominator <= target:
            lo = mid
        else:
            hi = mid
    lower = F(lo, scale)
    upper = lower if lo**degree*x.denominator == target else F(hi, scale)
    require(lower**degree <= x <= upper**degree, "radical interval checked by powers")
    return lower, upper


def log_weight(b):
    return F(25+(3*b+1).bit_length(), 28)


def polynomial_row_interval(c, count, degree=12):
    lo, hi = F(0), F(0)
    for b, _, _, mass in children(c, count):
        lower, upper = radical_interval(F(3*b+1, 3*c+1), degree)
        lo += mass*lower
        hi += mass*upper
    if degree == 12:
        # All depths, even the forbidden class, dominate the omitted tail.
        hi += F(41, 80)*F(9, 16)**count/(1-F(9, 16))
    return lo, hi


def allowed_cost_sum(residue):
    """Exact sum of 2^(-k-1)(2k+1), before root-loop deletion."""
    t = F(1, 8)
    total = F(0)
    for s in range(3):
        if (residue+s) % 3:
            total += F(1, 1 << (s+1))*((2*s+1)/(1-t)+6*t/(1-t)**2)
    return total


def finite_path_packet(target):
    chain, factors = [target], []
    while chain[-1] != 1:
        p, k = parent(chain[-1])
        factors.append(F(1, 1 << (k+1)))
        chain.append(p)
        if len(chain) > 128:
            raise ValueError("finite audit path cap")
    packet = {target: F(1)}
    for i, factor in enumerate(factors):
        packet[chain[i+1]] = factor*packet[chain[i]]
    return packet


def defect(packet):
    out = dict(packet)
    for b, value in packet.items():
        if b != 1:
            p, k = parent(b)
            out[p] = out.get(p, F(0))-F(1, 1 << (k+1))*value
    return {b: value for b, value in out.items() if value}


def main():
    report = {}
    require(3*41**12 >= 4*40**12, "(4/3)^(1/12) <= 41/40")
    require(9**6 >= 1 << 19, "2^(-5/6) <= 9/16")
    require(F(41,80)*(1+F(9,16))/(1-F(9,16)**3) == KP < 1,
            "exact polynomial contraction")
    costs = [allowed_cost_sum(c) for c in range(3)]
    require(costs == [F(95,49), F(106,49), F(93,49)], "three phase costs")
    require(F(6,7)+F(106,49*28) == KL < 1, "exact logarithmic contraction")
    report["constants"] = {
        "polynomial_exponent":"1/12", "polynomial_contraction":str(KP),
        "polynomial_normalized_mass_bound":str(1/(1-KP)),
        "log_weight":"(25+bitlength(3b+1))/28", "log_contraction":str(KL),
        "log_normalized_mass_bound":str(1/(1-KL)),
        "phase_additive_costs":list(map(str,costs))}

    bases = [b for b in range(1, 512, 2) if is_base(b)]
    max_poly_upper, max_log_upper = F(0), F(0)
    for c in bases:
        mass = F(0)
        log_head = F(0)
        for b,k,a,edge in children(c, 24):
            require(is_base(b) and parent(b) == (c,k), "inverse child and forward parent agree")
            require(3*b+1 == (1 << a)*((4**k*(3*c+1)-1)//3),
                    "shifted-height identity")
            require((3*b+1).bit_length() <= (3*c+1).bit_length()+2*k+1,
                    "bitlength growth bound")
            require(F(3*b+1, 3*c+1) < F(4,3)*4**k,
                    "polynomial height ratio")
            mass += edge
            log_head += edge*log_weight(b)/log_weight(c)
        require(mass <= F(6,7), "unweighted column bound")
        lower, upper = polynomial_row_interval(c,24)
        require(lower <= upper <= KP, "finite polynomial row with certified infinite tail")
        max_poly_upper = max(max_poly_upper,upper)
        # Exact arithmetic all-depth tail: sum_{k>=K} 2^-k-1(V(c)+2k+1).
        k = 24
        raw = 25+(3*c+1).bit_length()
        log_tail = F(raw+2*k+3, raw*(1 << k))
        require(log_head+log_tail <= KL, "finite logarithmic row with infinite tail")
        max_log_upper = max(max_log_upper,log_head+log_tail)
    report["finite_rows"] = {"bases_below_512":len(bases),"sibling_cut":24,
        "largest_certified_poly_upper":str(max_poly_upper),
        "largest_certified_log_upper":str(max_log_upper)}

    # A genuinely unbounded dual: a checked path packet minus eps*V_log.
    packet = finite_path_packet(3)
    require(packet == {3:F(1),1:F(1,4)}, "source3 path packet")
    require(defect(packet) == {3:F(1)}, "exact global packet defect including boundary")
    eps = F(1,8)
    require(log_weight(1) == 1 and packet[1]-eps == F(1,8), "positive unbounded-dual root value")
    require(1/(1-KL) == F(686,45), "weighted residual bill")
    # For every base, (I-P)V>= (1-KL)V; therefore phi=packet-eps V
    # has defect <=delta3. The finite rows are controls, not that global proof.
    report["unbounded_positive_control"] = {
        "target":3,"function":"path_packet_3 - V_log/8",
        "root_value":"1/8","implied_g3_lower":"1/8",
        "actual_g3":"1/4", "scope":"A strictly larger witness class, not a new convergence basin."}

    low_quarter, _ = polynomial_row_interval(25,18,degree=4)
    require(low_quarter > 1, "the exponent1/4 does not have one-step contraction in this norm")
    # At exponent1/2, every admissible large-k term has nonvanishing size:
    # 2^(-k-1)*sqrt((3b+1)/(3c+1)) tends to sqrt(2^a/3)/2.
    for b,k,a,edge in children(25,24):
        if k >= 6:
            require(edge**2*F(3*b+1,76) > F(1,9),
                    "half-power row terms stay bounded below")
    report["growth_hostiles"] = {
        "quarter_power_parent":25,"quarter_power_partial_row_lower":str(low_quarter),
        "half_power":"infinite row sum diverges; terms do not tend to zero"}

    # Abstract killed-root inverse ray: g_j=2^-j, P phi(j)=phi(j+1)/2.
    # Critical phi_j=2^j has zero defect and positive root but boundary never dies.
    for n in range(1,25):
        require(F(1,1 << n)*(1 << n) == 1, "critical harmonic boundary survives")
        require(F(1,2)*F(3,2) == F(3,4) < 1, "subcritical normalized contraction")
        require(F(1,1 << n)*F(3,2)**n == F(3,4)**n, "subcritical boundary decays")
    report["boundary_hostile"] = {
        "critical_weight":"2^j", "weighted_contraction":"1", "boundary":"1",
        "subcritical_weight":"(3/2)^j", "weighted_contraction":"3/4",
        "scope":"Dropping the strict growth/tail hypothesis gives the false readout1<=0."}
    report["checks"] = CHECKS
    report["status"] = "PROVED operator/dual domain; FINITE-EXACT controls; all-source positivity OPEN"
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
