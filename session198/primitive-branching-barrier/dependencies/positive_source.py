"""Exact, standalone checks for the positive physical source grammar."""
from functools import lru_cache
from fractions import Fraction
from itertools import product, permutations
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "dependencies"))
import bilinear
import c3
from core import inventories, cycles


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


GEN = {
    "O": ((1, 0, 1, 0), (1, 2, 1, 2), (0, 0, 0, 0), (0, 0, 1, 1)),
    "Z": ((2, 1, 1, 1), (0, 1, 0, 0), (0, 0, 1, 1), (0, 0, 0, 0)),
    "J": ((0, 1, 0, 1), (0, 0, 0, 0), (1, 0, 1, 1), (0, 0, 1, 0)),
}
P = {
    "N": ((6,4,10,9),(4,3,7,6),(10,7,15,14),(9,6,14,13)),
    "M": ((11,8,20,18),(8,6,14,13),(20,14,34,31),(18,13,31,29)),
    "A": ((14,10,21,20),(10,7,16,15),(21,16,33,31),(20,15,31,29)),
    "B": ((14,10,26,23),(11,8,20,18),(20,15,35,32),(19,14,34,31)),
    "C": ((23,18,44,40),(18,14,34,31),(44,34,80,73),(40,31,73,67)),
    "T": ((14,11,20,19),(10,8,15,14),(26,20,35,34),(23,18,32,31)),
    "J0": ((12,12,24,24),(12,12,21,21),(15,12,30,27),(15,12,30,27)),
    "J1": ((12,12,21,21),(12,12,18,18),(18,15,32,30),(18,15,31,29)),
    "J2": ((12,12,18,18),(12,12,15,15),(21,18,32,31),(21,18,30,29)),
    "J3": ((12,12,15,15),(12,12,12,12),(24,21,30,30),(24,21,27,27)),
}
UROW = (1, 1, 2, 2)
BASE = (0, 0, 0, 1)
ONE = (1, 1, 1, 1)
I = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))


def mv(A, v):
    return tuple(sum(a * b for a, b in zip(row, v)) for row in A)


def tr(A):
    return tuple(zip(*A))


def mm(A, B):
    return tuple(tuple(sum(a * b for a, b in zip(row, col)) for col in tr(B)) for row in A)


def rm(row, A):
    return tuple(sum(a * b for a, b in zip(row, col)) for col in tr(A))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def old(y):
    a, f, e, p = y
    return a + 2 * e + p, e, a + f + 2 * e + 2 * p, e + p


def new(s):
    D, E, U, d = s
    return D - E - d, U - D + E - d, E, d - E


def source_bank(cap):
    if cap < 2:
        return {}
    bank = {BASE: ()}
    queue = [BASE]
    for y in queue:
        for label, A in GEN.items():
            z = mv(A, y)
            if dot(UROW, z) <= cap and z not in bank:
                bank[z] = bank[y] + (label,)
                queue.append(z)
    return bank


def predecessors(y):
    a, f, e, p = y
    if e == 0 and (f - a) % 2 == 0:
        for q in range(min(a, p) + 1):
            yield "O", (a-q, (f-a-2*p)//2+q, q, p-q)
    if p == 0 and a-f-e >= 0 and (a-f-e) % 2 == 0:
        for q in range(e + 1):
            yield "Z", ((a-f-e)//2, f, q, e-q)
    if f == 0:
        for q in range(min(a, e-p) + 1):
            yield "J", (e-p-q, a-q, p, q)


@lru_cache(None)
def inverse(y):
    if min(y) < 0:
        return None
    D, E, U, d = old(y)
    if d < 1 or D % 2 == 0 or U < 2:
        return None
    if y == BASE:
        return ()
    for label, z in predecessors(y):
        if min(z) < 0 or old(z)[3] < 1 or dot(UROW, z) < 2:
            continue
        require(mv(GEN[label], z) == y, "inverse formula")
        require(dot(UROW, z) < U, "strict inverse decrease")
        path = inverse(z)
        if path is not None:
            return path + (label,)
    return None


def append(T, c, label):
    n = len(T)
    out = list(T) + [0]
    for v in range(n):
        if label in ("Z", "J") and v == c:
            out[n] |= 1 << v
        else:
            out[v] |= 1 << n
    return out, n if label == "J" else c


def profile(T, c):
    q = inventories(T)
    return q["H"], q["e"][c], q["R"][c], q["allH"][((1 << len(T)) - 1) ^ (1 << c)]


def build(word):
    T, c = [2, 0], 0
    for label in word:
        T, c = append(T, c, label)
    return T, c


def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(edges)):
        rows = [0] * n
        for k, (a, b) in enumerate(edges):
            if mask >> k & 1:
                rows[a] |= 1 << b
            else:
                rows[b] |= 1 << a
        yield rows


def direct_counts(T, root):
    """An independent literal-word and literal-cover control."""
    n = len(T)
    paths = {}
    for mask in range(1 << n):
        V = [v for v in range(n) if mask >> v & 1]
        paths[mask] = [w for w in permutations(V)
                       if all(T[a] >> b & 1 for a, b in zip(w, w[1:]))]
    full = (1 << n) - 1
    H = len(paths[full])
    E = sum(w[-1] == root for w in paths[full])
    starts = sum(w[0] == root for w in paths[full])
    U = sum(len(paths[mask]) * sum(bool(w) and w[0] == root for w in paths[full ^ mask])
            for mask in range(1 << n))
    return (H, E, U, len(paths[full ^ (1 << root)])), starts


def word_matrix(word):
    A = I
    for letter in word:
        A = mm(GEN[letter], A)
    return A


def minimum_ratio(A, B):
    require(all(x > 0 for row in B for x in row), "positive denominator matrix")
    return min(Fraction(A[i][j], B[i][j]) for i in range(4) for j in range(4))


def main():
    receipt = {"status": "PASS"}
    full = source_bank(2000)
    inherited = bilinear.source_bank(2000)
    inherited.pop(ONE)
    require({old(y) for y in full} == set(inherited), "complete capped source banks")
    require(len(full) == 32548, "frozen bank count")
    for y, word in full.items():
        require(new(old(y)) == y, "integral coordinate inverse")
        require(dot(UROW, y) >= 2 ** (1 + len(word)//2), "two-letter cover growth")
        require(mv(word_matrix(word), BASE) == y, "source recipe replay")
    receipt["complete_source_bank_U_cap"] = 2000
    receipt["complete_nontrivial_source_profiles"] = len(full)
    receipt["longest_retained_word"] = max(map(len, full.values()))

    # Finite exact membership tests include the whole 0..8 coordinate cube,
    # not merely positive source controls.
    tiny = source_bank(100)
    cases = set(product(range(9), repeat=4)) | set(tiny)
    accepted = 0
    for y in cases:
        found = inverse(y)
        require((found is not None) == (y in tiny), "exact inverse physical membership")
        if found is not None:
            accepted += 1
            require(mv(word_matrix(found), BASE) == y, "inverse witness replay")
    for U in range(2, 101):
        y = new((1, 0, U, 1))
        require((inverse(y) is not None) == (U & (U-1) == 0), "D=1 power-two hostile")
    require(inverse(new((1,0,3,1))) is None, "formal 19 source rejected")
    require(bilinear.value("N", (1,0,3,1), (1,0,2,1)) == 19, "formal scalar hostile")
    require(bilinear.inverse(19) is None, "physical 19 C2 omission")
    receipt["inverse_candidates"] = len(cases)
    receipt["inverse_accepted"] = accepted

    graph_controls = literal_controls = 0
    for n in range(2, 6):
        for T in all_tournaments(n):
            inv = inventories(T)
            for c in range(n):
                if inv["s"][c] != 1:
                    continue
                before = new(profile(T, c))
                require(min(before) >= 0, "source cone")
                for label in GEN:
                    A, r = append(T, c, label)
                    after = profile(A, r)
                    require(old(mv(GEN[label], before)) == after, "physical birth response")
                    require(inventories(A)["s"][r] == 1, "physical birth start one")
                    graph_controls += 1
                    if n <= 4:
                        direct, starts = direct_counts(A, r)
                        require(direct == after and starts == 1, "literal birth response")
                        literal_controls += 1
    receipt["physical_birth_controls"] = graph_controls
    receipt["independent_literal_birth_controls"] = literal_controls

    basis = [tuple(int(i == j) for i in range(4)) for j in range(4)]
    for name, mat in P.items():
        fn = bilinear.value if name in ("N", "M") else c3.value
        computed = tuple(tuple(fn(name, old(x), old(y)) for y in basis) for x in basis)
        require(computed == mat, "bilinear matrix reconstruction " + name)
        require(all(x > 0 for row in mat for x in row), "positive bilinear matrix")
    receipt["positive_bilinear_matrices"] = P

    u_differences = {}
    for letters in product(GEN, repeat=2):
        W = word_matrix(letters)
        row = rm(UROW, W)
        delta = tuple(row[i]-2*UROW[i] for i in range(4))
        require(min(delta) >= 0, "two-letter cover certificate")
        u_differences["".join(letters)] = delta
    receipt["two_letter_cover_differences"] = u_differences

    one_step = {}
    for name in ("N", "M"):
        for label, G in GEN.items():
            A = mm(tr(G), P[name])
            delta = tuple(tuple(3*A[i][j]-4*P[name][i][j] for j in range(4)) for i in range(4))
            require(min(x for row in delta for x in row) >= 0, "sharp C2 one-step certificate")
            one_step[name+label] = delta
    receipt["C2_one_step_3WTP_minus_4P"] = one_step

    two_min = {}
    for name in ("N", "M"):
        ratios = []
        for letters in product(GEN, repeat=2):
            W = word_matrix(letters)
            ratios.append((minimum_ratio(mm(tr(W), P[name]), P[name]), "".join(letters)))
        q, word = min(ratios)
        require(q > 2, "C2 two-letter halving")
        two_min[name] = {"ratio": str(q), "word": word}
    require(two_min["N"]["ratio"] == "15/7" and two_min["M"]["ratio"] == "7/3", "C2 frozen minima")
    receipt["C2_two_letter_minimum_ratios"] = two_min

    three_min = {}
    for name, mat in P.items():
        ratios = []
        for letters in product(GEN, repeat=3):
            W = word_matrix(letters)
            for side, A in (("left", mm(tr(W), mat)), ("right", mm(mat, W))):
                q = minimum_ratio(A, mat)
                require(q >= 2, "three-letter halving certificate")
                ratios.append((q, "".join(letters), side))
        q, word, side = min(ratios)
        three_min[name] = {"ratio": str(q), "word": word, "side": side}
    receipt["three_letter_minimum_ratios"] = three_min

    boundary_rows = {"C_singleton": (21,16,36,33), "T_singleton": (12,9,15,15), "J0_singleton": (12,9,15,15)}
    for name, row in boundary_rows.items():
        if name == "C_singleton":
            computed = tuple(c3.value("C", ONE, old(b)) for b in basis)
        elif name == "T_singleton":
            computed = tuple(c3.value("T", ONE, old(b)) for b in basis)
        else:
            computed = tuple(c3.value("J0", old(b), ONE) for b in basis)
        require(computed == row, "singleton boundary row")
        ratios = []
        for letters in product(GEN, repeat=3):
            A = rm(row, word_matrix(letters))
            q = min(Fraction(A[i], row[i]) for i in range(4))
            require(q >= 2, "singleton boundary halving")
            ratios.append(q)
        boundary_rows[name] = {"row": row, "minimum_ratio": str(min(ratios))}
    receipt["singleton_boundaries"] = boundary_rows

    sharpness = []
    for u, v in ((2,2), (4,8), (8,8)):
        left = (1,0,u,1)
        right = (1,0,v,1)
        Jleft = old(mv(GEN["J"], new(left)))
        V, c = build(("O",) * (u.bit_length()-2))
        W, r = build(("O",) * (v.bit_length()-2))
        JV, jc = append(V, c, "J")
        expected = {
            "N": (3*u*v+1, 4*u*v+u+3*v-1),
            "M": (6*u*v+u+v+1, 8*u*v+2*u+6*v+1),
        }
        for kind in ("N", "M"):
            before, after = expected[kind]
            require(bilinear.value(kind, left, right) == before, "sharpness base identity")
            require(bilinear.value(kind, Jleft, right) == after, "sharpness birth identity")
            for A, a, target in ((V,c,before),(JV,jc,after)):
                G, _ = bilinear.construct(kind,A,a,W,r)
                require(inventories(G)["H"] == target and cycles(G) == 2, "literal C2 sharpness graph")
            sharpness.append({"u":u,"v":v,"family":kind,"before":before,"after":after})
    receipt["physical_sharpness_checks"] = sharpness
    receipt["hashes"] = {
        str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted([HERE / "PROOF.md", HERE / "verify.py"] + list((HERE / "dependencies").glob("*.py")) + list((HERE / "dependencies").glob("*.md")))
    }
    output = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    (HERE / "VALIDATION.json").write_text(output)
    print(output)


if __name__ == "__main__":
    main()
