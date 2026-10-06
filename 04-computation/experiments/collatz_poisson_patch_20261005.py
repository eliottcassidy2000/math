"""Exact Schur patches for the bounded Collatz Poisson dual.

Finite matrices are test kernels, not inferred quotients of the infinite graph.
Actual Collatz controls use explicitly replayed, parent-closed finite sets.
"""
from fractions import Fraction as F
from itertools import product
import json

from collatz_green_weight_20261005 import primitive, edge

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def rational(value):
    if type(value) not in (int, F):
        raise ValueError("exact rational required")
    return F(value)


def matrix(rows, square=False):
    if type(rows) not in (list, tuple):
        raise ValueError("matrix sequence required")
    result = tuple(tuple(rational(v) for v in row) for row in rows)
    if result and any(len(row) != len(result[0]) for row in result):
        raise ValueError("rectangular matrix required")
    if square and any(len(row) != len(result) for row in result):
        raise ValueError("square matrix required")
    return result


def vector(values, size):
    if type(values) not in (list, tuple) or len(values) != size:
        raise ValueError("vector dimension")
    return tuple(rational(v) for v in values)


def ident(n):
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def transpose(a):
    return tuple(zip(*a))


def add(a, b):
    return tuple(tuple(x+y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def mul(a, b):
    if not a:
        return ()
    if not b:
        return tuple(() for _ in a)
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


def mv(a, x):
    return tuple(sum((v*w for v, w in zip(row, x)), F(0)) for row in a)


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), F(0))


def inverse(a):
    a = matrix(a, True)
    n = len(a)
    rows = [list(a[i])+list(ident(n)[i]) for i in range(n)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [x/scale for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[j])]
    return tuple(tuple(row[n:]) for row in rows)


def defect_matrix(p):
    return tuple(tuple(F(i == j)-v for j, v in enumerate(row))
                 for i, row in enumerate(p))


def kernel(rows):
    p = matrix(rows, True)
    if not p or any(v < 0 for row in p for v in row):
        raise ValueError("nonempty nonnegative kernel required")
    if max(map(sum, p)) >= 1:
        raise ValueError("strict row contraction required")
    return p


def submatrix(a, rows, cols):
    return tuple(tuple(a[i][j] for j in cols) for i in rows)


def schur_patch(rows, forcing, ports):
    """Exact finite patch: eliminate the complement of retained labelled ports."""
    p = kernel(rows)
    n = len(p)
    q = vector(forcing, n)
    if (type(ports) not in (tuple, list) or not ports or
            any(type(i) is not int or not 0 <= i < n for i in ports) or
            len(set(ports)) != len(ports)):
        raise ValueError("distinct nonempty exact port indices required")
    j = tuple(ports)
    interior = tuple(i for i in range(n) if i not in j)
    pjj = submatrix(p, j, j)
    if not interior:
        return {"P":p, "q":q, "ports":j, "interior":(), "R":(),
                "K":pjj, "b":tuple(q[i] for i in j)}
    pii = submatrix(p, interior, interior)
    pji = submatrix(p, j, interior)
    pij = submatrix(p, interior, j)
    resolvent = inverse(defect_matrix(pii))
    k = add(pjj, mul(mul(pji, resolvent), pij))
    charge = mv(mul(pji, resolvent), tuple(q[i] for i in interior))
    b = tuple(q[v]+charge[t] for t, v in enumerate(j))
    return {"P":p, "q":q, "ports":j, "interior":interior,
            "R":resolvent, "K":k, "b":b}


def harmonic_extension(patch, port_values):
    p, q = patch["P"], patch["q"]
    j, interior = patch["ports"], patch["interior"]
    h = vector(port_values, len(j))
    phi = [F(0)]*len(p)
    for v, value in zip(j, h):
        phi[v] = value
    if interior:
        coupling = mv(submatrix(p, interior, j), h)
        rhs = tuple(q[v]+coupling[t] for t, v in enumerate(interior))
        for v, value in zip(interior, mv(patch["R"], rhs)):
            phi[v] = value
    return tuple(phi)


def reduced_readout(patch, root_port=0):
    """Return exact maximal root value and its nonnegative dual occupation row."""
    j = patch["ports"]
    if type(root_port) is not int or root_port not in j:
        raise ValueError("root must be a retained port")
    r = inverse(defect_matrix(patch["K"]))
    z = r[j.index(root_port)]
    return dot(z, patch["b"]), z, mv(r, patch["b"])


def lumping(rows, blocks):
    """Check exact row lumpability; averaging unequal rows is rejected."""
    p = kernel(rows)
    if type(blocks) not in (tuple, list) or not blocks:
        raise ValueError("partition required")
    if any(type(block) not in (tuple, list) or not block for block in blocks):
        raise ValueError("nonempty blocks required")
    flat = [i for block in blocks for i in block]
    if any(type(i) is not int for i in flat) or sorted(flat) != list(range(len(p))):
        raise ValueError("partition indices")
    result = []
    for block in blocks:
        profiles = [tuple(sum((p[i][j] for j in target), F(0))
                          for target in blocks) for i in block]
        if any(profile != profiles[0] for profile in profiles):
            raise ValueError("not an exact lumping")
        result.append(profiles[0])
    return tuple(result)


def actual_parent_closed(sources):
    """Finite controls only: explicitly resolve these selected G-paths."""
    vertices = {1}
    for source in sources:
        primitive(source)
        current = source
        seen = set()
        while current != 1:
            if current in seen or len(seen) >= 1000:
                raise ValueError("control path did not close within declared bound")
            seen.add(current)
            vertices.add(current)
            current, _ = edge(current)
    vertices = tuple(sorted(vertices))
    index = {v:i for i, v in enumerate(vertices)}
    p = [[F(0)]*len(vertices) for _ in vertices]
    for source in vertices:
        if source != 1:
            parent, k = edge(source)
            p[index[parent]][index[source]] = F(1, 2**(k+1))
    return vertices, kernel(p)


def partitions(n):
    """All set partitions of range(n), in restricted-growth order."""
    def visit(word):
        if len(word) == n:
            yield tuple(tuple(i for i,a in enumerate(word) if a==j)
                        for j in range(max(word)+1))
            return
        for value in range(max(word)+2):
            yield from visit(word+(value,))
    yield from visit((0,))


def run():
    global CHECKS
    CHECKS = 0
    universes = {"all_binary_quarter_kernels_3x3":0, "patches":0,
                 "nested_patches":0, "actual_sources":0,
                 "functional_tree_partitions":0, "exact_lumpings":0}
    for bits in product((0, 1), repeat=9):
        p = tuple(tuple(F(bits[3*i+j],4) for j in range(3)) for i in range(3))
        q = (F(1), F(-1,3), F(2,5))
        full_r = inverse(defect_matrix(p))
        exact = mv(full_r, q)
        universes["all_binary_quarter_kernels_3x3"] += 1
        for mask in range(4):
            ports = (0,)+tuple(i for i in (1,2) if mask & (1 << (i-1)))
            patch = schur_patch(p, q, ports)
            value, z, h = reduced_readout(patch)
            phi = harmonic_extension(patch, h)
            need(phi == exact, "Schur global solution")
            need(value == exact[0] == dot(z, patch["b"]), "dual sign/readout")
            need(all(v >= 0 for v in z), "occupation nonnegative")
            need(max(map(sum, patch["K"])) <= max(map(sum,p)), "return contraction")
            need(mv(transpose(defect_matrix(patch["K"])), z) ==
                 tuple(F(i==0) for i in range(len(ports))), "dual feasibility")
            # A global subsolution with an explicit nonnegative residual slack.
            phi0 = tuple(a-b for a,b in zip(exact, mv(full_r,(F(0),F(1),F(0)))))
            h0 = tuple(phi0[i] for i in ports)
            lifted = harmonic_extension(patch, h0)
            need(all(a >= b for a,b in zip(lifted,phi0)), "harmonic domination")
            need(all(a <= b for a,b in zip(mv(defect_matrix(p),lifted),q)),
                 "subsolution gluing")
            need(dot(z, mv(defect_matrix(patch["K"]),h0)) == h0[0],
                 "localized residual bill")
            universes["patches"] += 1
    for seed in range(128):
        p = tuple(tuple(F(((seed+1)*(i+3)*(j+5)+i+j)%3,16)
                        for j in range(4)) for i in range(4))
        q = (F(1),F(-2),F(3),F(-4))
        direct = schur_patch(p,q,(0,3))
        first = schur_patch(p,q,(0,2,3))
        second = schur_patch(first["K"],first["b"],(0,2))
        need(direct["K"] == second["K"] and direct["b"] == second["b"],
             "associative labelled elimination")
        universes["nested_patches"] += 1
    actual = []
    for source in (3,7,11,17,27,31,47,63,95,127):
        vertices,p = actual_parent_closed((source,))
        root, target = vertices.index(1),vertices.index(source)
        q = tuple(F(i==target) for i in range(len(vertices)))
        patch = schur_patch(p,q,(root,target))
        value,z,h = reduced_readout(patch)
        need(patch["K"][0][1] == value > 0, "actual source segment")
        need(all(sum(v != 0 for v in col) <= 1 for col in transpose(patch["K"])),
             "one forward parent after compression")
        need(mv(defect_matrix(p),harmonic_extension(patch,h)) == q,
             "actual extension defect")
        actual.append({"source":source,"vertices":len(vertices),"root_floor":str(value)})
        universes["actual_sources"] += 1
    vertices,p = actual_parent_closed((7,))
    q = tuple(F(i==vertices.index(7))-F(1,32)*F(i==vertices.index(3))
              for i in range(len(vertices)))
    patch = schur_patch(p,q,tuple(vertices.index(v) for v in (1,11,7)))
    signed,z,_ = reduced_readout(patch)
    need(signed == F(1,128), "signed source margin")
    need(z == (F(1),F(1,32),F(1,64)), "localized error prices")
    # A non-Collatz countable renewal patch: root->v0, vi->v(i+1), vi->target.
    a,b,c = F(1,2),F(1,4),F(1,4)
    infinite_readout = a*c/(1-b)
    for depth in range(20):
        partial = a*c*sum((b**j for j in range(depth)),F(0))
        need(infinite_readout-partial == a*c*b**depth/(1-b), "infinite renewal tail")
        need(infinite_readout-partial <= F(1,2)**(depth+2), "first-return tail bound")
    # Missing port charge: local interior equation does not validate ROOT.
    ptoy = ((F(0),F(1,2),F(0)),(F(0),F(0),F(1,2)),(F(0),F(0),F(0)))
    false_phi = (F(1),F(0),F(0))
    false_defect = mv(defect_matrix(ptoy),false_phi)
    need(false_defect[1] <= 0 and false_defect[0] == 1, "omitted root boundary hostile")
    try:
        lumping(ptoy,((0,1),(2,)))
    except ValueError:
        need(True,"averaging is not lumping")
    else:
        raise ValueError("unfaithful quotient accepted")
    coarse = lumping(((F(0),F(1,4)),(F(1,4),F(0))),((0,1),))
    need(coarse == ((F(1,4),),), "non-singleton target loses singleton forcing")
    # A disconnected contracted cycle has a unique solution, but no ROOT readout.
    cycle = ((F(0),F(0),F(0)),(F(0),F(0),F(1,2)),(F(0),F(1,2),F(0)))
    cp = schur_patch(cycle,(0,1,0),(0,1))
    value,_,h = reduced_readout(cp)
    need(cp["K"][1][1] == F(1,4) and h[1] == F(4,3) and value == 0,
         "contracted cycle is not ROOT forcing")
    for n in range(2,6):
        for parents in product(*(range(i) for i in range(1,n))):
            p = [[F(0)]*n for _ in range(n)]
            for child,parent in enumerate(parents,1):
                p[parent][child] = F(1,2*n)
            for blocks in partitions(n):
                universes["functional_tree_partitions"] += 1
                try:
                    lumping(p,blocks)
                except ValueError:
                    continue
                universes["exact_lumpings"] += 1
                singletons = {block[0] for block in blocks if len(block)==1}
                for target in singletons:
                    current = target
                    while current != 0:
                        current = parents[current-1]
                        need(current in singletons, "singleton ancestry is forced")
    bad = [lambda:schur_patch(((1,),),(1,),(0,)),
           lambda:schur_patch(((0,),),(1,),(True,)),
           lambda:schur_patch(((0.0,),),(1,),(0,)),
           lambda:schur_patch(((0,),),(1,),(0,0)),
           lambda:reduced_readout(schur_patch(((0,),),(1,),(0,)),True)]
    for call in bad:
        try:
            call()
        except (ValueError,TypeError):
            need(True,"exact API guard")
        else:
            raise ValueError("hostile accepted")
    result = {"status":"FINITE-EXACT controls; global claims proved in note",
              "universe":universes,"actual_parent_closed_controls":actual,
              "signed_margin_7_minus_3_over_32":str(signed),
              "localized_error_prices_at_1_11_7":[str(v) for v in z],
              "synthetic_infinite_patch_readout":str(infinite_readout),
              "checks":CHECKS,
              "coverage":"No new Collatz ROOT coverage or global positivity claim."}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == "__main__":
    run()
