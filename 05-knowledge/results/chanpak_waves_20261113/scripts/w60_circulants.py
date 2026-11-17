"""Wave 60: odd-cycle data of the two Z7-circulants: interval I7 = {1,2,3} and Fano/Paley P7 = {1,2,4} (quadratic residues mod 7)."""
import itertools
def circ(S): return [sum(1 << ((i + d) % 7) for d in S) for i in range(7)]
def cycles(out, k):
    cnt = 0; sets = []
    for comb in itertools.combinations(range(7), k):
        for perm in itertools.permutations(comb[1:]):
            seq = (comb[0],) + perm
            if all(out[seq[i]] >> seq[(i + 1) % k] & 1 for i in range(k)): cnt += 1; sets.append(frozenset(comb))
    return cnt, sets
for name, S in (("interval I7 {1,2,3}", (1, 2, 3)), ("Fano/Paley P7 {1,2,4}", (1, 2, 4))):
    out = circ(S); c3, s3 = cycles(out, 3); c5, _ = cycles(out, 5); c7, _ = cycles(out, 7)
    D2 = sum(1 for a, b in itertools.combinations(s3, 2) if not (a & b))
    G = lambda x, y: y ** 7 + c3 * x * y ** 4 + c5 * x * y ** 2 + D2 * x * x * y + c7 * x
    print(f"{name}: triangles {c3}, 5-cycles {c5}, Hamiltonian 7-cycles {c7}, disjoint triangle pairs {D2}; "
          f"readings P_q = G(2q, q) for q = 1, 2, 3, 4: {[G(2 * q, q) for q in (1, 2, 3, 4)]}")
