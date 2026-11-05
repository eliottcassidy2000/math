"""Exact controls for primitive-source branching obstructions and rebasing."""
from fractions import Fraction
from itertools import product
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "dependencies"))
import positive_source as s
from core import inventories, cycles, literal_h


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def family(i, j, k):
    require(i >= 0 and j >= 1 and k >= 1, "primitive family domain")
    z, y, q = 2 ** (i+1) + 1, 2 ** j, 2 ** k
    D = (3*y+1)*z-2
    response = (D, 0, q*(D+1), (y+1)*z)
    word = ("O",)*i + ("J",) + ("O",)*j + ("Z", "J") + ("O",)*k
    return response, word, (z,y,q)


def bilinear_value(P, x, y):
    return s.dot(x, s.mv(P, y))


def complete_binary_menu():
    words = [w for length in range(3) for w in product("OZJ", repeat=length)]
    output = []
    for left in words:
        for right in words:
            if len(left) + len(right) <= 2:
                mat = s.mm(s.mm(s.tr(s.word_matrix(left)), s.P["N"]), s.word_matrix(right))
                output.append(("N", left, right, mat))
    output.append(("M", (), (), s.P["M"]))
    for kind in ("J0", "J1", "J2", "J3"):
        for a,b in product(range(3), repeat=2):
            left, right = ("J",)*a, ("J",)*b
            mat = s.mm(s.mm(s.tr(s.word_matrix(left)), s.P[kind]), s.word_matrix(right))
            output.append((kind,left,right,mat))
    return output


def primitive_grid_witness(menu):
    """Exact bounded search promised by the finite-menu theorem.

    This routine is only exercised on small menus. The theorem proves
    termination of the full bound; no runtime claim is made for large m.
    """
    m = len(menu)
    indices = list(product(range(m+1), range(1,m+2), range(1,m+2)))
    family_bank = [(ijk, s.new(family(*ijk)[0])) for ijk in indices]
    for i,x in family_bank:
        for j,y in family_bank:
            n = bilinear_value(s.P["N"],x,y)
            if all(bilinear_value(P,x,y) != b*n+c for P,b,c in menu):
                return {"left_indices": i, "right_indices": j,
                        "left_response": s.old(x), "right_response": s.old(y),
                        "parent": n, "candidate_values": [bilinear_value(P,x,y) for P,_,_ in menu]}
    raise RuntimeError("finite interpolation-grid witness missing")


def main():
    result = {"status": "PASS"}
    basis = [tuple(int(i==j) for i in range(4)) for j in range(4)]
    for name, P in s.P.items():
        fn = s.bilinear.value if name in ("N","M") else s.c3.value
        recovered = tuple(tuple(fn(name,s.old(a),s.old(b)) for b in basis) for a in basis)
        require(recovered == P, "independent positive-matrix reconstruction")

    checks = []
    for i,j,k in product(range(5), range(1,6), range(1,6)):
        response, word, parameters = family(i,j,k)
        actual = s.old(s.mv(s.word_matrix(word),s.BASE))
        require(actual == response, "primitive family exact word")
        D,E,U,d = actual
        require(gcd(D,U) == 1 and gcd(*actual) == 1, "primitive physical family")
        z = Fraction(3*d-D-2,2)
        y = Fraction(d,z)-1
        q = Fraction(U,D+1)
        require((z,y,q) == parameters, "exact rational parameter inverse")
        require(len(word)+2 == i+j+k+5, "physical source order bound")
        checks.append((i,j,k))
    result["primitive_parameter_checks"] = len(checks)
    literal_family = []
    for indices in ((0,1,1),(1,1,1),(0,2,1)):
        response, word, _ = family(*indices)
        V,c = s.build(word)
        inv = inventories(V)
        require(s.profile(V,c) == response and inv["s"][c] == 1, "literal primitive source")
        require(literal_h(V) == response[0], "independent primitive Hamilton words")
        literal_family.append({"indices":indices,"response":response,"root":c,"rows":V})
    result["literal_primitive_family"] = literal_family

    # All matrix bounds required to make the 71-case search all-order.
    certificate = {}
    for kind, P in s.P.items():
        certificate[kind] = {"relative_N_minimum": str(s.minimum_ratio(P,s.P["N"]))}
        for letter,G in s.GEN.items():
            left = s.mm(s.tr(G),P)
            right = s.mm(P,G)
            q = min(s.minimum_ratio(left,P),s.minimum_ratio(right,P))
            certificate[kind][letter] = str(q)
            require(q >= 1, "every source letter nondecreasing")
            if kind in ("N","M"):
                require(q >= Fraction(4,3), "C2 one-letter growth")
            if kind.startswith("J") and letter in ("O","Z"):
                require(q >= Fraction(3,2), "Type-II O/Z growth")
        if kind in ("A","B","C","T"):
            delta = [[P[i][j]-2*s.P["N"][i][j] for j in range(4)] for i in range(4)]
            require(min(sum(delta,[])) >= 0, "Type-I excess matrix")
            ep_min = min(delta[i][j] for i,j in product((2,3),repeat=2))
            require(ep_min >= 3, "Type-I physical strict excess")
            certificate[kind]["ep_excess_minimum"] = ep_min
        if kind.startswith("J"):
            require(s.minimum_ratio(P,s.P["N"]) >= Fraction(3,2), "Type-II relative floor")
            J3 = s.word_matrix(("J",)*3)
            q = min(s.minimum_ratio(s.mm(s.tr(J3),P),P),s.minimum_ratio(s.mm(P,J3),P))
            require(q >= 2, "Type-II JJJ growth")
            certificate[kind]["JJJ"] = str(q)
    require(s.minimum_ratio(s.P["M"],s.P["N"]) == Fraction(11,6), "M relative floor")
    result["all_word_reduction_certificate"] = certificate

    menu = complete_binary_menu()
    require(len(menu) == 71, "complete binary candidate count")
    require(len(set(P for _,_,_,P in menu)) == 65, "distinct candidate matrix count")
    v = s.new((5,1,5,1))
    t = s.BASE
    n = bilinear_value(s.P["N"],v,t)
    require(n == 41 and gcd(*s.old(v)) == gcd(*s.old(t)) == 1, "primitive bad pair")
    values = sorted(set(bilinear_value(P,v,t) for _,_,_,P in menu))
    expected = [41,69,73,75,77,79,85,93,99,115,123,125,129,131,133,135,137,
                141,143,147,149,153,155,167,225,249,257,263,267,273,283,287,
                289,303,435,441,467,469,479,499,501,819,855,885,897]
    require(values == expected and 81 not in values and 83 not in values, "all-length primitive bad pair")
    result["binary_menu_size"] = len(menu)
    result["binary_menu_distinct_matrices"] = len(set(P for _,_,_,P in menu))
    result["bad_pair"] = {"left":s.old(v),"right":s.old(t),"parent":n,
                          "excluded_children":[81,83],"complete_candidate_values":values}
    result["binary_menu"] = [{"family":k,"left_word":a,"right_word":b,"matrix":P}
                             for k,a,b,P in menu]

    # Small constructive hostile grids are positive controls for the
    # interpolation theorem, not the premise for its all-menu quantifier.
    menus = [
        [(s.P["N"],2,1)],
        [(s.P["M"],2,1),(s.P["J1"],2,-1)],
        [(menu[i][3],2,1 if i % 2 else -1) for i in (1,4,15,40)],
    ]
    result["finite_grid_controls"] = [primitive_grid_witness(M) for M in menus]

    # Directly verify the scalar-preserving rebase and the two children.
    edge, eroot = s.build(())
    triangle, hroot = s.build(("Z",))
    V, vr = s.append(triangle,hroot,"Z")
    W, wr = s.append(triangle,hroot,"J")
    OW, owr = s.append(W,wr,"O")
    graph_cases = [
        ("parent_Z","N",V,vr,edge,eroot,41,2),
        ("parent_J","N",W,wr,edge,eroot,41,2),
        ("parent_repartition","N",triangle,hroot,triangle,hroot,41,2),
        ("child_81","N",OW,owr,edge,eroot,81,2),
        ("child_83","J1",triangle,hroot,triangle,hroot,83,3),
    ]
    literal_graphs = []
    for label,kind,A,a,B,b,H,C in graph_cases:
        make = s.bilinear.construct if kind in ("N","M") else s.c3.construct
        G,root = make(kind,A,a,B,b)
        require(inventories(G)["H"] == H and cycles(G) == C, "physical rebase graph count")
        require(literal_h(G) == H, "independent rebase Hamilton words")
        literal_graphs.append({"label":label,"family":kind,"H":H,"C":C,"rows":G,"root":root})
    result["literal_rebase_graphs"] = literal_graphs

    bank = s.source_bank(2000)
    zero = []
    primitive_zero = []
    for x,word in bank.items():
        a,f,e,p = x
        difference = bilinear_value(s.P["N"],s.mv(s.GEN["Z"],x),t)-bilinear_value(s.P["N"],s.mv(s.GEN["J"],x),t)
        require(difference == 4*a+6*f-4*e, "exact twin defect")
        if difference == 0:
            require(x == t or (a == e and f == p == 0 and a >= 1), "physical twin zero classification")
            zero.append(x)
            if gcd(*s.old(x)) == 1:
                primitive_zero.append(x)
    require(set(primitive_zero) == {t,(1,0,1,0)}, "primitive twin zero set")
    result["twin_bank_cap"] = 2000
    result["twin_bank_size"] = len(bank)
    result["twin_zero_profiles"] = len(zero)
    result["primitive_twin_zero_profiles"] = primitive_zero

    result["hashes"] = {str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted(HERE.rglob("*")) if p.is_file() and p.suffix in (".py",".md")}
    output = json.dumps(result,indent=2,sort_keys=True)+"\n"
    (HERE/"VALIDATION.json").write_text(output)
    print(output)


if __name__ == "__main__":
    main()
