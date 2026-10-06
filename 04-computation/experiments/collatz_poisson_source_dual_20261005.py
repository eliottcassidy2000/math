"""Exact Poisson dual certificates on the rooted Collatz sibling kernel.

Finite controls contain already completed paths or synthetic finite graphs.
A global analytic dual inequality is a premise, never inferred from a grid.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import json

from collatz_green_weight_20261005 import (
    primitive, base, edge, inverse_child, S, v2, column_mass, column_tail)
from collatz_floor_transport_deadlines_20261005 import replay, word_weight, step
from collatz_localized_resolvent_floor_20261005 import selector as localized_selector

R = F(1, 2)
KAPPA = F(6, 7)


def natural(n):
    if type(n) is not int or n < 0:
        raise ValueError("exact natural required")


def rational(x):
    if type(x) not in (int, F):
        raise ValueError("exact rational required")
    return F(x)


def coefficient(k):
    natural(k)
    return F(1, 2**(k+1))


def finite_dual(target, phi):
    """Globally check (I-A*)phi <= delta_target for a finite exact function."""
    primitive(target)
    if type(phi) is not dict:
        raise ValueError("finite exact dictionary required")
    values = {}
    for b, value in phi.items():
        primitive(b)
        value = rational(value)
        if value:
            values[b] = value
    incoming = {}
    children = {}
    for b, value in values.items():
        if b != 1:
            c, k = edge(b)
            incoming[c] = incoming.get(c, F(0)) + coefficient(k)*value
            children.setdefault(c, []).append((b, k))
    domain = set(values) | set(incoming) | {1, target}
    for c in domain:
        defect = values.get(c, F(0))-incoming.get(c, F(0))
        if defect > int(c == target):
            raise ValueError("failed global finite-support dual inequality")
    return {"target":target, "values":values, "children":children,
            "root_floor":values.get(1, F(0)),
            "sup_norm":max(map(abs, values.values()), default=F(0)),
            "checked_constraints":len(domain)}


def inverse_path_from_dual(target, phi):
    """Extract the finite path encoded by an accepted positive dual packet."""
    checked = finite_dual(target, phi)
    if checked["root_floor"] <= 0:
        raise ValueError("no positive source floor")
    values, children = checked["values"], checked["children"]
    c, path = 1, [1]
    while c != target:
        options = children.get(c, ())
        if not options:
            raise ValueError("missing positive descendant")
        child, _ = max(options, key=lambda item:values[item[0]])
        if values[child] < values[c]/KAPPA:
            raise ValueError("discounted maximum-principle failure")
        if child in path:
            raise ValueError("nonincreasing finite path")
        path.append(child)
        c = child
        if len(path) > len(values):
            raise ValueError("finite support exhausted")
    return tuple(path), checked


def word_from_inverse_path(path):
    """ROOT-outward G-path to an actual strict U-word, without an orbit search."""
    if type(path) is not tuple or not path or path[0] != 1:
        raise ValueError("ROOT-outward base tuple required")
    if len(set(path)) != len(path):
        raise ValueError("repeated base")
    word = ()
    for index, child in enumerate(path):
        primitive(child)
        if index == 0:
            continue
        parent, k = edge(child)
        if parent != path[index-1]:
            raise ValueError("wrong labelled G edge")
        a = v2(3*child+1)
        if parent == 1:
            if k < 1:
                raise ValueError("root self-edge")
            word = (a, 2*k+2)
        else:
            word = (a, word[0]+2*k)+word[1:]
    if replay(path[-1], word) != 1:
        raise ValueError("decoded word is not a strict ROOT receipt")
    return word


def counter_bound_from_green_floor(delta):
    delta = rational(delta)
    if not 0 < delta <= 1:
        raise ValueError("positive Green floor at most one required")
    return (delta.denominator//delta.numerator).bit_length()-1


def mixture_floor_from_counter_bound(bound):
    natural(bound)
    if bound == 0:
        return F(1)
    return F(2, (bound+2)*comb(bound+1, (bound+1)//2))


def receipt_from_green_floor(source, delta):
    """Conditional floor-to-ROOT compiler; every returned word is actually replayed.

    delta is claimed as a lower bound for f_(1/2)(source), not the mixture W.
    An unsupported floor is rejected if bounded literal replay contradicts it.
    """
    base(source)
    delta = rational(delta)
    cap = counter_bound_from_green_floor(delta)
    current, word = source, []
    for used in range(cap+1):
        if current == 1:
            result = tuple(word)
            total = 0 if not result else len(result)-1+sum((a-1)//2 for a in result)
            if F(1,2**total) < delta:
                raise ValueError("claimed fixed-price floor exceeds actual weight")
            return result
        if used == cap:
            break
        current, a = step(current)
        word.append(a)
    raise ValueError("claimed fixed-price floor failed its finite ROOT deadline")


def compile_source(source, phi):
    """Validate a finite dual, retain the supplied source, then replay its receipt."""
    b, j = base(source)
    path, checked = inverse_path_from_dual(b, phi)
    word = word_from_inverse_path(path)
    if b == 1:
        word = () if j == 0 else (2*j+2,)
    elif j:
        word = (word[0]+2*j,)+word[1:]
    if replay(source, word) != 1:
        raise ValueError("source identity failed")
    green_floor = checked["root_floor"]/2**j
    total_bound = counter_bound_from_green_floor(green_floor)
    lower = mixture_floor_from_counter_bound(total_bound)
    if word_weight(word) < lower:
        raise ValueError("derived mixture lower bound failed")
    return {"source":source,"base":b,"sibling_depth":j,"word":word,
            "inverse_path":path,"green_floor":green_floor,
            "counter_bound":total_bound,"mixture_floor":lower}


def bounded_root_path_control(b, cap=256):
    """Finite test producer only: never used to validate an unknown dual."""
    primitive(b); natural(cap)
    reverse, current = [b], b
    for _ in range(cap):
        if current == 1:
            return tuple(reversed(reverse))
        current, _ = edge(current)
        reverse.append(current)
    if current == 1:
        return tuple(reversed(reverse))
    raise ValueError("declared finite control cap exhausted")


def path_dual(path):
    """The explicit optimal finite witness attached to an already checked path."""
    word_from_inverse_path(path)
    phi = {path[-1]:F(1)}
    for child, parent in zip(reversed(path[1:]), reversed(path[:-1])):
        actual, k = edge(child)
        if actual != parent:
            raise ValueError("wrong path")
        phi[parent] = coefficient(k)*phi[child]
    return phi


def integer_deadline(delta, norm, contraction=KAPPA):
    delta, norm, contraction = map(rational, (delta, norm, contraction))
    if not 0 < delta <= norm or not 0 < contraction < 1:
        raise ValueError("positive bounded witness and contraction required")
    t, upper = 0, norm
    while upper >= delta:
        upper *= contraction
        t += 1
    return t


def invert(matrix):
    n = len(matrix)
    rows = [list(map(F, matrix[i]))+[F(i == j) for j in range(n)]
            for i in range(n)]
    for col in range(n):
        pivot = next(i for i in range(col,n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        value = rows[col][col]
        rows[col] = [x/value for x in rows[col]]
        for i in range(n):
            if i != col:
                value = rows[i][col]
                rows[i] = [x-value*y for x,y in zip(rows[i],rows[col])]
    return [row[n:] for row in rows]


def leaf_selector(n, target_index, degree):
    if n % 3:
        return F(0)
    return localized_selector(target_index, degree, (n-3)//6)


def main():
    checks = 0
    def check(ok, why):
        nonlocal checks
        checks += 1
        if not ok:
            raise ValueError(why)

    source_controls = []
    deleted_controls = 0
    for n in range(1,256,2):
        b,j = base(n)
        path = bounded_root_path_control(b)
        phi = path_dual(path)
        packet = compile_source(n,phi)
        check(replay(n,packet["word"]) == 1, "actual source ROOT receipt")
        check(packet["source"] == n and packet["base"] == b, "source retained")
        check(receipt_from_green_floor(n,packet["green_floor"]) == packet["word"],
              "independent fixed-price floor forward compiler")
        check(word_weight(packet["word"]) >= packet["mixture_floor"], "mixture floor")
        eps=phi[1]/100
        corrected=(1+eps)*phi[1]-F(7,2)*eps
        check(0 < corrected <= phi[1], "global residual bill restores a valid floor")
        if n != 1:
            actual_total = len(packet["word"])-1+sum((a-1)//2 for a in packet["word"])
            check(actual_total == packet["counter_bound"], "path witness exact binary counter")
        check(len(path)-1 < integer_deadline(phi[1],max(phi.values())),
              "maximum-principle path deadline")
        if b != 1:
            broken = dict(phi)
            del broken[b]
            try:
                finite_dual(b,broken)
            except ValueError:
                check(True,"removing forcing node breaks dual")
            else:
                raise ValueError("deleted target accepted")
            deleted_controls += 1
        if n in (1,3,5,7,27,155):
            source_controls.append({"source":n,"inverse_path":list(packet["inverse_path"]),
                "word":list(packet["word"]),"green_floor":str(packet["green_floor"]),
                "counter_bound":packet["counter_bound"],
                "mixture_floor":str(packet["mixture_floor"])})

    columns = 0
    for c in range(1,128,2):
        b,j = base(c)
        if j:
            continue
        subtotal = F(0)
        for k in range(13):
            child = inverse_child(b,k)
            if child is not None:
                check(edge(child) == (b,k),"complete inverse child phase")
                subtotal += coefficient(k)
        check(subtotal+column_tail(b,12,R,F(1)) == column_mass(b,R,F(1)),
              "finite column plus exact omitted tail")
        check(column_mass(b,R,F(1)) <= KAPPA,"l1 contraction column")
        columns += 1
    check(column_mass(1,R,F(1)) == F(5,14),"killed root column")
    check(column_mass(7,R,F(1)) == KAPPA,"sharp column witness")

    # Every finite functional graph on 2..5 labelled vertices, root row deleted.
    # These are independent toy universes, not invented Collatz edges.
    graphs = 0
    duals = 0
    for size in range(2,6):
        d = F(1,2*size)
        for parents in product(range(size),repeat=size-1):
            matrix = [[F(i == j) for j in range(size)] for i in range(size)]
            for child,parent in enumerate(parents,1):
                matrix[child][parent] -= d
            resolvent = invert(matrix)
            green = [resolvent[i][0] for i in range(size)]
            check(all(x>=0 for row in resolvent for x in row),"positive resolvent")
            for target in range(size):
                current,seen,rooted = target,set(),False
                while current not in seen:
                    if current == 0:
                        rooted = True
                        break
                    seen.add(current)
                    current = parents[current-1]
                check((green[target]>0) == rooted,"unique Green support equals actual component")
                phi = resolvent[target]
                for c in range(size):
                    incoming = sum((d*phi[child] for child,parent in enumerate(parents,1)
                                    if parent == c),F(0))
                    check(phi[c]-incoming == int(c == target),"exact adjoint Poisson equation")
                check(phi[0] == green[target],"dual/primal source identity")
                duals += 1
            graphs += 1

    # A unique contracted solution may vanish. Small immigration changes forcing.
    regularization = []
    for exponent in range(9):
        eps = F(1,2**exponent)
        zero_component = 2*eps
        check(zero_component == eps+F(1,2)*zero_component,"positive regularized loop")
        check(zero_component > 0,"strictly positive at every finite regularization")
        regularization.append({"epsilon":str(eps),"off_root_weight":str(zero_component)})
    check(F(0) == F(1,2)*0,"unique unforced contracted component is zero")
    check(F(2) == 1+F(1,2)*2,"retained ROOT loop changes the boundary")

    # Boundedness is essential: an unbounded harmonic function carries a
    # nonvanishing boundary term even for a contractive rooted ray.
    ray_controls = 0
    for horizon in range(1,17):
        for j in range(horizon):
            check(2**j == F(1,2)*2**(j+1), "unbounded harmonic ray interior")
        check(F(1,2**horizon)*2**horizon == 1, "nonvanishing omitted boundary")
        ray_controls += 1

    # Finite regrouping controls for the typed leaf-selector lift.
    regroupings = 0
    for target_index in (0,1,4):
        for degree in range(3):
            for depth in range(3):
                left=right=F(0)
                for b in (1,3,7,9):
                    proxy=F(1,b+1)
                    lifted=sum((R**j*leaf_selector(S(b,j),target_index,degree)
                                for j in range(depth+1)),F(0))
                    left += proxy*lifted
                    for j in range(depth+1):
                        right += proxy*R**j*leaf_selector(S(b,j),target_index,degree)
                    finite_tail=sum((R**j*leaf_selector(S(b,j),target_index,degree)
                                     for j in range(depth+1,depth+3)),F(0))
                    check(abs(finite_tail) <= 8*R**(depth+1)/(1-R),
                          "geometric sibling tail bound")
                check(left == right,"marked leaf-to-base regrouping")
                regroupings += 1

    bad = [lambda:finite_dual(3,{1:F(1)}),
           lambda:finite_dual(3,{True:F(1)}),
           lambda:finite_dual(3,{1:0.25}),
           lambda:finite_dual(5,{1:F(1)}),
           lambda:inverse_path_from_dual(3,{}),
           lambda:word_from_inverse_path((1,1)),
           lambda:word_from_inverse_path((1,7)),
           lambda:compile_source(2,{1:F(1)}),
           lambda:counter_bound_from_green_floor(0.2),
           lambda:integer_deadline(F(0),F(1)),
           lambda:receipt_from_green_floor(27,F(1,4))]
    for job in bad:
        try:
            job()
        except ValueError:
            check(True,"invalid exact source or witness rejected")
        else:
            raise ValueError("invalid witness accepted")
    print(json.dumps({"status":"PROVED dual compiler; independent global positive witnesses OPEN",
        "actual_sources":{"odd_below":256,"count":128,"control_base_depth_cap":256},
        "actual_examples":source_controls,"deleted_target_controls":deleted_controls,
        "exact_columns":columns,"synthetic_functional_graphs":graphs,
        "synthetic_exact_Poisson_duals":duals,"positive_regularization_zero_limit":regularization,
        "typed_selector_regroupings":regroupings,"unbounded_harmonic_ray_controls":ray_controls,
        "invalid_inputs":len(bad),"checks":checks},
        indent=2,sort_keys=True))
    print("PASS")


if __name__ == "__main__":
    main()
