"""Actual R9 block with a kernel nonnormal in G but normal in its projection."""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    return tuple(p.index(i) for i in range(len(p)))


def conjugate(p, q):
    return compose(compose(p, q), inverse(p))


def preserves(A, p):
    return all(A[i][j] == A[p[i]][p[j]] for i in range(9) for j in range(9))


def main():
    A = tuple(tuple(int((i//3 != j//3 and (j//3-i//3) % 3 == 1)
                        or (i//3 == j//3 and (j-i) % 3 == 1))
                    for j in range(9)) for i in range(9))
    G = tuple(p for p in permutations(range(9)) if preserves(A, p))
    need(len(G) == 81, "full S9 automorphism enumeration")
    identity = tuple(range(9))
    B = tuple(tuple(int((i < 3 and j < 3)
                        or (i >= 3 and j >= 3 and (j-i) % 6 in (0, 1, 3)))
                    for j in range(9)) for i in range(9))
    need({sum(row) for row in B} == {3}, "actual block row degree")
    need({sum(B[i][j] for i in range(9)) for j in range(9)} == {3},
         "actual block column degree")
    need(len(set(B[:3])) == 1 and len(set(B[3:])) == 6,
         "one repeated row triple and six distinct outside rows")
    H = {(g, h) for g in G for h in G
         if all(B[i][j] == B[g[i]][h[j]] for i in range(9) for j in range(9))}
    P = {g for g, h in H}
    Q = {h for g, h in H}
    K = {g for g, h in H if h == identity}
    L = {h for g, h in H if g == identity}
    need((len(H), len(P), len(Q), len(K), len(L)) == (9, 3, 3, 3, 3),
         "actual projection and kernel sizes")
    need(H == set(product(K, L)), "quotient has one point with nine lifts")
    need(all(conjugate(p, k) in K for p in P for k in K), "normal inside P")
    outer = tuple((i+3) % 9 for i in range(9))
    need(outer in G and outer not in P, "outer rotation is real but inadmissible")
    nontrivial = next(k for k in sorted(K) if k != identity)
    conjugated = conjugate(outer, nontrivial)
    need(conjugated not in K, "kernel not normal in full G")
    need({i for i in range(9) if nontrivial[i] != i} == {0, 1, 2},
         "source kernel support")
    need({i for i in range(9) if conjugated[i] != i} == {3, 4, 5},
         "conjugated kernel support")

    # Literal target system, not a formal re-labeling of a relation table.
    row_map = outer
    col_map = compose(outer, outer)
    ri, ci = inverse(row_map), inverse(col_map)
    target = tuple(tuple(B[ri[i]][ci[j]] for j in range(9)) for i in range(9))
    R = {(g, h) for g in G for h in G
         if all(B[i][j] == target[g[i]][h[j]] for i in range(9) for j in range(9))}
    expected = {(compose(row_map, g), compose(col_map, h)) for g, h in H}
    need(R == expected and len(R) == 9, "source-target relation is r times source H")
    need({(compose(ri, g), compose(ci, h)) for g, h in R} == H,
         "actual projected representatives translate back to P domains")
    need(all((compose(g, k), compose(h, ell)) in R
             for g, h in R for k in K for ell in L), "source kernels act on right")
    wrong = (compose(nontrivial, row_map), col_map)
    need(wrong not in R, "source kernel on left is wrong after transport")
    target_K = {conjugate(row_map, k) for k in K}
    target_L = {conjugate(col_map, k) for k in L}
    need(all((compose(k, g), compose(ell, h)) in R
             for g, h in R for k in target_K for ell in target_L),
         "conjugated target kernels act on left")
    observations = Counter(g[0] for g, h in R)
    need(dict(observations) == {3: 3, 4: 3, 5: 3},
         "nine unweighted lifts do not determine actual owner image")
    out = dict(status="FINITE-EXACT; actual regular nine-cell nonnormal-kernel hostile",
               internal_adjacency=A, cross_block=B,
               full_group=len(G), stabilizer=len(H),
               projections=[len(P), len(Q)], kernels=[len(K), len(L)],
               quotient_frames=1, lifts_per_quotient=9,
               outer_rotation=outer, source_kernel_generator=nontrivial,
               conjugated_generator=conjugated,
               target_cross_block=target, row_reference=row_map, column_reference=col_map,
               source_target_frames=len(R), wrong_left_action=wrong,
               owner_image_counts=dict(sorted(observations.items())), checks=CHECKS)
    text = json.dumps(out, indent=2, sort_keys=True)+"\n"
    frozen = HERE / "non_normal_kernel_receipt.json"
    if "--write" in sys.argv:
        frozen.write_text(text, encoding="utf-8")
    else:
        need(json.loads(frozen.read_text(encoding="utf-8")) == json.loads(text),
             "frozen exact receipt")
    print(json.dumps(dict(status="PASS", checks=CHECKS, full_group=len(G),
                          projections=[len(P), len(Q)], kernels=[len(K), len(L)],
                          lifts=len(R), owner_image_counts=dict(sorted(observations.items())))))


if __name__ == "__main__":
    main()
