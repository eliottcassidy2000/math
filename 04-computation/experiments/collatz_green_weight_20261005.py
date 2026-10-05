"""Rooted Green weights for the Collatz sibling-base operator.

All arithmetic is exact. Infinite statements have proofs in the companion note.
No arbitrary source is declared rooted by a small approximation error.
"""
from fractions import Fraction as F
import json

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def natural(x, zero=True):
    if type(x) is not int or x < (0 if zero else 1):
        raise ValueError("exact nonnegative/positive integer required")


def odd(x):
    natural(x, False)
    if x % 2 != 1:
        raise ValueError("odd integer required")


def v2(x):
    natural(x, False)
    return (x & -x).bit_length()-1


def U(n):
    odd(n)
    return (3*n+1) >> v2(3*n+1)


def S(n, k):
    odd(n)
    natural(k)
    return 4**k*n+(4**k-1)//3


def base(n):
    odd(n)
    k = 0
    while n % 8 == 5:
        n = (n-1)//4
        k += 1
    return n, k


def primitive(b):
    odd(b)
    if v2(3*b+1) not in (1, 2):
        raise ValueError("primitive sibling base required")


def parameters(r, rho):
    if type(r) not in (int, F) or type(rho) not in (int, F):
        raise ValueError("exact rational parameters required")
    r, rho = F(r), F(rho)
    if not (0 < r < 1 and 0 < rho <= 1):
        raise ValueError("require 0<r<1 and 0<rho<=1")
    return r, rho


def edge(b):
    primitive(b)
    if b == 1:
        raise ValueError("ROOT row is killed")
    return base(U(b))


def inverse_child(c, k):
    """The unique possible base at target S^k(c); ROOT source is deleted."""
    primitive(c)
    natural(k)
    y = S(c, k)
    if y % 3 == 0:
        return None
    b = ((4 if y % 3 == 1 else 2)*y-1)//3
    if b == 1:
        return None
    return b


def d(k, r, rho):
    natural(k)
    r, rho = parameters(r, rho)
    return rho*(1-r)*r**k


def constants(r, rho):
    r, rho = parameters(r, rho)
    den = 1+r+r*r
    kappa = rho*(1+r)/den
    root_column = rho*r*(1+r*r)/den
    return rho*(1-r), kappa, root_column


def column_mass(c, r, rho):
    primitive(c)
    r, rho = parameters(r, rho)
    j = (-c) % 3
    result = rho*(1-r)*(1/(1-r)-r**j/(1-r**3))
    if c == 1:
        result -= rho*(1-r)
    return result


def column_tail(c, cutoff, r, rho):
    """Exact column mass at depths k>cutoff (cutoff>=0)."""
    primitive(c)
    natural(cutoff)
    r, rho = parameters(r, rho)
    start = cutoff+1
    excluded = start+((-c-start) % 3)
    return rho*(1-r)*(r**start/(1-r)-r**excluded/(1-r**3))


def green_value(b, depth, r=F(1, 16), rho=F(1, 2)):
    """Return a certified interval [value,value+error] for the global weight.

    A positive value is returned only with an actual finite G-path to ROOT.
    An unresolved zero approximation is not an assertion that the limit is zero.
    """
    primitive(b)
    natural(depth)
    r, rho = parameters(r, rho)
    c0, _, _ = constants(r, rho)
    current, weight = b, F(1)
    for length in range(depth+1):
        if current == 1:
            return {'value': weight, 'error': F(0), 'root_depth': length}
        if length == depth:
            break
        current, k = edge(current)
        weight *= d(k, r, rho)
    return {'value': F(0), 'error': c0**(depth+1), 'root_depth': None}


def green_to_error(b, epsilon, r=F(1, 16), rho=F(1, 2)):
    primitive(b)
    if type(epsilon) not in (int, F) or epsilon <= 0:
        raise ValueError("positive exact error tolerance required")
    r, rho = parameters(r, rho)
    c0, _, _ = constants(r, rho)
    depth = 0
    while c0**(depth+1) > epsilon:
        depth += 1
    return green_value(b, depth, r, rho), depth


def finite_kernel(depth, cutoff, r=F(1, 16), rho=F(1, 2)):
    """Finite support for sum_{j=0}^depth A_cutoff^j delta_ROOT."""
    natural(depth)
    natural(cutoff)
    r, rho = parameters(r, rho)
    layer = {1: F(1)}
    weights = dict(layer)
    rows = {1: None}
    for _ in range(depth):
        nxt = {}
        for parent, weight in layer.items():
            for k in range(cutoff+1):
                b = inverse_child(parent, k)
                if b is None:
                    continue
                need(b not in weights and b not in nxt, "unique first-hit inverse ancestry")
                nxt[b] = d(k, r, rho)*weight
                rows[b] = (parent, k)
        weights.update(nxt)
        layer = nxt
    for b, weight in weights.items():
        primitive(b)
        if b == 1:
            need(weight == 1 and rows[b] is None, "root normalization and killed row")
        else:
            parent, k = rows[b]
            need(edge(b) == (parent, k), "independent actual base edge")
            need(U(b) == S(parent, k), "literal common-future equation")
            need(weight == d(k, r, rho)*weights[parent], "exact supported equality")
    return weights, rows


def approximation_bound(depth, cutoff, r, rho):
    """An L1 bound for global Green weight minus the finite kernel."""
    natural(depth)
    natural(cutoff)
    _, kappa, a0 = constants(r, rho)
    r, rho = parameters(r, rho)
    depth_error = a0*kappa**depth/(1-kappa)
    operator_error = rho*r**(cutoff+1)
    branch_error = operator_error/(1-kappa)**2
    return depth_error+branch_error


def finite_to_error(epsilon, r=F(1, 16), rho=F(1, 2)):
    if type(epsilon) not in (int, F) or epsilon <= 0:
        raise ValueError("positive exact error tolerance required")
    epsilon = F(epsilon)
    r, rho = parameters(r, rho)
    _, kappa, a0 = constants(r, rho)
    depth = 0
    while a0*kappa**depth/(1-kappa) > epsilon/2:
        depth += 1
    cutoff = 0
    while rho*r**(cutoff+1)/(1-kappa)**2 > epsilon/2:
        cutoff += 1
    weights, rows = finite_kernel(depth, cutoff, r, rho)
    return weights, rows, depth, cutoff, approximation_bound(depth, cutoff, r, rho)


def nu(n):
    odd(n)
    return F(8, 3*4**n.bit_length())


def main():
    report = {}
    bases = [b for b in range(1, 512, 2) if v2(3*b+1) in (1, 2)]
    param_grid = [(r, rho) for r in (F(1,16), F(1,4), F(1,2))
                  for rho in (F(1,2), F(2,3), F(1))]
    column_controls = 0
    for r, rho in param_grid:
        c0, kappa, a0 = constants(r, rho)
        need(column_mass(7, r, rho) == kappa, "operator norm attained")
        need(column_mass(1, r, rho) == a0, "killed root column")
        need(0 < a0 < kappa < rho <= 1, "strict contraction")
        for c in bases:
            seen = set()
            head = F(0)
            for k in range(9):
                b = inverse_child(c, k)
                if b is not None:
                    need(b not in seen and b != 1, "distinct child and no root self-loop")
                    seen.add(b)
                    need(edge(b) == (c, k), "inverse child has exact decoded parent")
                    head += d(k, r, rho)
                need(head+column_tail(c, k, r, rho) == column_mass(c, r, rho),
                     "exact column tail completion")
                need(column_mass(c, r, rho) <= kappa, "all three phases bounded")
                column_controls += 1
    report['column_controls_bases_lt512_cutoffs0_to8_nine_parameter_pairs'] = column_controls

    critical_bounds = {}
    for radius in (F(1,16), F(1,4), F(1,2)):
        _, critical_kappa, critical_a0 = constants(radius, F(1))
        critical_bound = 1+critical_a0/(1-critical_kappa)
        need(critical_bound == 1+radius+1/radius,
             "critical root-column improvement")
        critical_weights, _ = finite_kernel(3,3,radius,F(1))
        need(sum(critical_weights.values()) <= critical_bound,
             "critical finite inverse kernel under sharp root-start bound")
        critical_bounds[str(radius)] = {
            'root_start_bound': str(critical_bound),
            'inherited_general_column_bound': str(1/(1-critical_kappa)),
            'finite_depth3_cutoff3_mass': str(sum(critical_weights.values()))}
    report['critical_rho1_bounds'] = critical_bounds

    r, rho = F(1,16), F(1,2)
    c0, kappa, a0 = constants(r, rho)
    need(kappa == F(136,273) and a0 == F(257,8736), "default exact constants")
    mass_bound = 1+a0/(1-kappa)
    need(mass_bound == F(4641,4384), "root-specific global mass bound")
    report['default_constants'] = {
        'r': str(r), 'rho': str(rho), 'row_bound': str(c0),
        'operator_L1_norm': str(kappa), 'root_column': str(a0),
        'global_total_mass_upper_bound': str(mass_bound)}

    precision_controls = 0
    for b in bases:
        truth = green_value(b, 256, r, rho)
        need(truth['root_depth'] is not None, "bounded control sources have finite checked routes")
        for depth in (0,1,2,3,4,8,16,32,64):
            x = green_value(b, depth, r, rho)
            need(x['value'] <= truth['value'] <= x['value']+x['error'],
                 "global precision interval against independent longer control")
            if x['value'] > 0:
                need(x['value'] == truth['value'], "grounded value is exact")
            precision_controls += 1
        # The depth variable is a finite base path, not an assumed root rank.
        target, product_weight = b, F(1)
        for _ in range(truth['root_depth']):
            target, k = edge(target)
            product_weight *= d(k, r, rho)
        need(target == 1 and product_weight == truth['value'], "literal product control")
    report['point_precision_controls'] = precision_controls
    report['selected_exact_weights'] = {
        str(b): {'weight': str(green_value(b,256)['value']),
                 'base_root_depth': green_value(b,256)['root_depth']}
        for b in (1,3,7,27,127)}

    need(green_value(27, 0)['value'] == 0
         and green_value(27, 0)['error'] > 0
         and green_value(27,256)['value'] > 0,
         "zero finite approximation is not zero limiting weight")
    for b in bases[:48]:
        current, k = edge(b) if b != 1 else (1,0)
        if b != 1:
            multiplier = d(k,r,rho)*nu(current)/nu(b)
            need(multiplier == rho*(1-r)*nu(U(b))/nu(b), "sibling depth cancels at r1/16")
            need(multiplier in (F(15,128),F(15,32),F(15,8)), "three factored multipliers")
    report['factored_multipliers'] = ['15/128','15/32','15/8']

    kernel_rows = []
    for depth, cutoff in ((0,0),(1,1),(2,2),(4,3),(5,5)):
        weights, rows = finite_kernel(depth,cutoff,r,rho)
        larger, _ = finite_kernel(depth+1,cutoff+1,r,rho)
        need(all(larger.get(b,0) == value for b,value in weights.items()),
             "nested kernel values are retained exactly")
        difference = sum(larger.values())-sum(weights.values())
        bound = approximation_bound(depth,cutoff,r,rho)
        need(0 <= difference <= bound, "finite enlargement obeys global bound")
        need(sum(weights.values()) <= mass_bound, "kernel mass under global bound")
        for b,value in weights.items():
            need(green_value(b,depth,r,rho)['value'] == value, "independent forward decode of every kernel base")
        kernel_rows.append({'depth': depth, 'branch_cutoff': cutoff, 'bases': len(weights),
                            'weight_mass': str(sum(weights.values())),
                            'global_L1_error_bound': str(bound),
                            'certified_nu_sibling_mass': str(F(16,15)*sum(nu(b) for b in weights))})
    report['finite_kernels'] = kernel_rows

    finite_precision = []
    for epsilon in (F(1,100),F(1,1000)):
        weights, _, depth, cutoff, error = finite_to_error(epsilon,r,rho)
        need(error <= epsilon, "terminating finite global precision compiler")
        finite_precision.append({'epsilon': str(epsilon), 'depth': depth,
                                 'branch_cutoff': cutoff, 'bases': len(weights),
                                 'proven_L1_error': str(error),
                                 'contains_base27': 27 in weights,
                                 'certified_nu_sibling_mass': str(F(16,15)*sum(nu(b) for b in weights))})
    report['global_precision_compiler'] = finite_precision
    huge_weights, _, _, _, huge_error = finite_to_error(10**400)
    need(huge_weights[1] == 1 and type(huge_error) is F,
         "large exact integer tolerance never coerces to floating point")

    hostile = [
        lambda: green_value(True, 1), lambda: green_value(2,1),
        lambda: green_value(5,1), lambda: green_value(3,False),
        lambda: green_value(3,1,0.0625,F(1,2)), lambda: green_to_error(3,0),
        lambda: finite_kernel(1,False), lambda: edge(1),
        lambda: approximation_bound(True,1,r,rho),
        lambda: approximation_bound(1.0,1,r,rho),
        lambda: approximation_bound(-1,1,r,rho),
        lambda: approximation_bound(1,False,r,rho),
        lambda: approximation_bound(1,-1,r,rho),
        lambda: d(-1,r,rho), lambda: d(1,0.0625,rho),
    ]
    for job in hostile:
        try:
            job()
        except ValueError:
            need(True, "malformed or ungrounded type rejected")
        else:
            raise ValueError("accepted hostile")
    report['invalid_input_controls'] = len(hostile)

    # Abstract hostile: G(n)=n+1 for n>=2, ROOT isolated, d=1/2.
    # A has L1 norm1/2 but A delta_ROOT=0, so its Green vector is delta_ROOT.
    need(all((n+1) != 1 for n in range(2,101)), "toy root has no incoming nonroot edge")
    report['abstract_contraction_does_not_imply_positive_support'] = {
        'map': 'G(n)=n+1 for n>=2; ROOT1 isolated',
        'operator_norm': '1/2', 'Green_at_1': 1, 'Green_at_every_n_ge2': 0}

    report['exact_checks'] = CHECKS
    print(json.dumps(report,indent=2,sort_keys=True))
    print('PASS: unconditional computable summable Green weight; positivity everywhere remains OPEN.')


if __name__ == '__main__':
    main()
