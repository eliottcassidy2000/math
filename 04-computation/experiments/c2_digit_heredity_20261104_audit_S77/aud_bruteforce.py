"""aud_bruteforce.py -- permutation-level checks for the S77 audit (no dynamic programming at all).

Definitions are applied literally:
  h(T)   = number of vertex permutations that are directed paths;
  B_a(T) = number of pairs (P1, P2) of vertex-disjoint directed paths covering V, P1 possibly empty, P2 starting at a,
           enumerated as (permutation, split point);
  digit sigma_a = B_a - h;
  C2(B, a): x -> a only (every other old vertex beats x); y beaten by every old vertex; y -> x.

Checks:
  (1) the C2 law H(C2(B, a)) = h(B) + 2 B_a(B) for every labelled tournament B of order <= 4 and every a;
  (2) the native digits sigma_x(C2) = 0 and sigma_y(C2) = h(B) - e_a(B) on the same set;
  (3) the 'broken-cover' reading sigma_a = #{(P1, P2): P1 nonempty ending at an out-neighbour of a, P2 starting at a};
  (4) the three digit-2 grandchildren: B transitive of order 1-3, a = its sink, u in the 3-cycle {a, x, y} of the child:
      which marks of the grandchild carry digit 2 (the note says only the native y; the audit finds u as well).
"""
from itertools import permutations


def tour_from_code(n, code):
    arcs = set()
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> k) & 1:
                arcs.add((i, j))
            else:
                arcs.add((j, i))
            k += 1
    return arcs


def is_path(seq, arcs):
    return all((seq[i], seq[i + 1]) in arcs for i in range(len(seq) - 1))


def h_count(V, arcs):
    return sum(1 for p in permutations(V) if is_path(p, arcs))


def ends(V, arcs):
    e = {v: 0 for v in V}
    for p in permutations(V):
        if is_path(p, arcs):
            e[p[-1]] += 1
    return e


def B_and_broken(V, arcs, a):
    """B_a via (permutation, split) pairs; also the broken-cover count."""
    V = list(V)
    Bc = 0
    broken = 0
    for p in permutations(V):
        k = p.index(a)          # P1 = p[:k], P2 = p[k:]
        P1, P2 = p[:k], p[k:]
        if is_path(P1, arcs) and is_path(P2, arcs):
            Bc += 1
            if k > 0 and (a, P1[-1]) in arcs:
                broken += 1
    return Bc, broken


def c2(V, arcs, a):
    n = max(V) + 1 if V else 0
    x, y = n, n + 1
    A = set(arcs)
    A.add((x, a))
    for v in V:
        if v != a:
            A.add((v, x))
        A.add((v, y))
    A.add((y, x))
    return list(V) + [x, y], A, x, y


def main():
    cases = law_bad = nat_bad = broken_bad = 0
    for n in range(1, 5):
        for code in range(2 ** (n * (n - 1) // 2)):
            arcs = tour_from_code(n, code)
            V = list(range(n))
            h = h_count(V, arcs)
            e = ends(V, arcs)
            for a in V:
                Ba, br = B_and_broken(V, arcs, a)
                if br != Ba - h:
                    broken_bad += 1
                V1, A1, x, y = c2(V, arcs, a)
                H1 = h_count(V1, A1)
                cases += 1
                if H1 != h + 2 * Ba:
                    law_bad += 1
                Bx, _ = B_and_broken(V1, A1, x)
                By, _ = B_and_broken(V1, A1, y)
                if Bx - H1 != 0 or By - H1 != h - e[a]:
                    nat_bad += 1
    print(f"labelled (B, a), |B| <= 4: {cases} cases; C2 law failures {law_bad}; native digit failures {nat_bad}; "
          f"broken-cover reading failures {broken_bad}")

    # (4) the C3-bottom grandchildren
    for n in range(1, 4):
        # transitive B on 0..n-1 with i -> j for i < j; sink = n-1
        arcs = set((i, j) for i in range(n) for j in range(i + 1, n))
        V = list(range(n))
        a = n - 1
        V1, A1, x1, y1 = c2(V, arcs, a)
        H1 = h_count(V1, A1)
        for u, name in ((a, "a"), (x1, "x1"), (y1, "y1")):
            V2, A2, x2, y2 = c2(V1, A1, u)
            H2 = h_count(V2, A2)
            small = []
            for w in V2:
                Bw, _ = B_and_broken(V2, A2, w)
                d = Bw - H2
                if d in (1, 2):
                    role = "x2" if w == x2 else "y2" if w == y2 else ("u" if w == u else "old")
                    small.append(f"{role}:{d}")
            print(f"|B|={n} transitive, a = sink, u = {name}: H(child) = {H1}, H(grandchild) = {H2}, marks with digit 1/2: {small}")


if __name__ == "__main__":
    main()
