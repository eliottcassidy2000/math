"""b_states.py - auditor's own enumeration of physical source responses (THM-4888 (B) classification, seed lists).

Moves (Session 198 bilinear-descent (4)): O(D,E,U,d)=(D,0,2U,d), Z=(D-E+U,d,2U-d,d), J=(D-E+U,D-E,D+U,D), from
t=(1,0,2,1).  Completeness of a U-bounded search: every move raises U (increments U, U-d, D, all >= 1 because
U >= 2d >= 2 by Session 198 inequalities (1)), so every ancestor of a state with U <= UMAX has smaller U.  We
check U >= 2d and inequalities (1) on every enumerated state, and that the strict increase holds on every edge.

Also: Session 198's exact inverse (section 3) is re-derived here by brute-force solving of the linear moves (5) in
positive coordinates y=(a,f,e,p)=(D-E-d, U-D+E-d, E, d-E) and compared with the printed predecessor formulas on all
enumerated states; the BFS parents must be among them.
Outputs JSON to stdout.
"""
import json, sys
from fractions import Fraction

def O(s): D, E, U, d = s; return (D, 0, 2 * U, d)
def Z(s): D, E, U, d = s; return (D - E + U, d, 2 * U - d, d)
def J(s): D, E, U, d = s; return (D - E + U, D - E, D + U, D)
MOVES = (('O', O), ('Z', Z), ('J', J))
T0 = (1, 0, 2, 1)

def bfs(umax):
    parent = {T0: None}; order = [T0]; jborn = set(); i = 0
    edges_ok = True
    while i < len(order):
        s = order[i]; i += 1
        for ch, f in MOVES:
            c = f(s)
            if c[2] <= s[2]:
                edges_ok = False
            if c[2] > umax:
                continue
            if ch == 'J':
                jborn.add(c)
            if c not in parent:
                parent[c] = (s, ch); order.append(c)
    return parent, jborn, edges_ok

def y_of(s):
    D, E, U, d = s; return (D - E - d, U - D + E - d, E, d - E)

def s_of(y):
    a, f, e, p = y; return (a + 2 * e + p, e, a + f + 2 * e + 2 * p, e + p)

def printed_predecessors(y):
    """Session 198 section 3, verbatim formulas."""
    a, f, e, p = y; out = []
    if e == 0 and (f - a) % 2 == 0:
        for q in range(0, min(a, p) + 1):
            out.append(('O', (a - q, (f - a - 2 * p) // 2 + q, q, p - q)))
    if p == 0 and a - f - e >= 0 and (a - f - e) % 2 == 0:
        for q in range(0, e + 1):
            out.append(('Z', ((a - f - e) // 2, f, q, e - q)))
    if f == 0:
        for q in range(0, min(a, e - p) + 1):
            out.append(('J', (e - p - q, a - q, p, q)))
    return out

def solved_predecessors(s, bound):
    """brute force: all integer x=(D,E,U,d) with nonnegative y(x), d>=1, U>=2, U < U(s), and move(x)=s."""
    out = []
    D1, E1, U1, d1 = s
    # O: D=D1, U=U1/2, d=d1, E free with y(x)>=0
    if E1 == 0 and U1 % 2 == 0:
        for E in range(0, d1 + 1):
            x = (D1, E, U1 // 2, d1)
            if min(y_of(x)) >= 0 and x[3] >= 1 and x[2] >= 2: out.append(('O', y_of(x)))
    # Z: d=E1=d1, U=(U1+d1)/2, D-E = D1-U
    if E1 == d1 and (U1 + d1) % 2 == 0:
        U = (U1 + d1) // 2
        for E in range(0, d1 + 1):
            x = (D1 - U + E, E, U, d1)
            if min(y_of(x)) >= 0 and x[3] >= 1 and x[2] >= 2: out.append(('Z', y_of(x)))
    # J: d1 = D, E1 = D - E, U1 = D + U, D1 = D - E + U
    D = d1; E = D - E1; U = U1 - D
    if D1 == D - E + U:
        for dd in range(1, bound + 1):
            x = (D, E, U, dd)
            if min(y_of(x)) >= 0 and x[2] >= 2: out.append(('J', y_of(x)))
    return out

def main():
    UMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4096
    parent, jborn, edges_ok = bfs(UMAX)
    st = list(parent)
    res = {'UMAX': UMAX, 'states': len(st)}
    # inequalities (1) and U >= 2d, oddness of D, every move strictly raises U
    bad = [s for s in st if not (s[1] <= s[3] and s[0] - s[1] >= s[3] and s[2] - s[0] >= s[3] - s[1] and s[2] >= 2 * s[3] and s[0] % 2 == 1)]
    res['violations_of_ineq1_or_U_ge_2d_or_D_odd'] = len(bad)
    inc_bad = 0
    for s in st:
        for ch, f in MOVES:
            if f(s)[2] <= s[2]: inc_bad += 1
    res['moves_not_raising_U'] = inc_bad
    res['states_U_eq_2d_E_eq_0'] = sorted(s for s in st if s[2] == 2 * s[3] and s[1] == 0)
    res['states_R_eq_2d_any_E'] = len([s for s in st if s[2] == 2 * s[3]])
    res['states_U_le_64'] = len([s for s in st if s[2] <= 64])
    # inverse formulas: printed vs brute-force solved, on states with U <= 512 (all states, both ways)
    mism = 0; parent_missing = 0; checked = 0
    for s in st:
        if s[2] > 512: continue
        checked += 1
        y = y_of(s)
        pr = sorted(printed_predecessors(y))
        sv = sorted(solved_predecessors(s, s[2]))
        pr_ok = sorted((ch, x) for ch, x in pr if min(x) >= 0 and (x[2] + x[3]) >= 1 and s_of(x)[2] >= 2)
        if pr_ok != sv: mism += 1
        if parent[s] is not None:
            ps, ch = parent[s]
            if (ch, y_of(ps)) not in pr_ok: parent_missing += 1
    res['inverse_check_states'] = checked
    res['printed_vs_solved_predecessor_mismatches'] = mism
    res['bfs_parent_not_among_printed_predecessors'] = parent_missing
    # the classification step of THM-4888 (B): y=(0,0,0,d) has no nonnegative predecessor for every d>=1
    res['y_000d_predecessors_d_1_to_50'] = {d: [(ch, x) for ch, x in printed_predecessors((0, 0, 0, d)) if min(x) >= 0] for d in range(1, 51) if [(ch, x) for ch, x in printed_predecessors((0, 0, 0, d)) if min(x) >= 0]}
    # seed lists (P1, bank A) and comparison with the author's files
    p1 = sorted(s for s in jborn if 64 < s[2] <= 160 and s[3] > 1 and s[2] != 2 * s[3] and s[2] != s[0])
    bankA = [T0] + sorted(s for s in jborn if s[2] <= 64 and s[3] > 1)
    res['P1_count'] = len(p1); res['bankA_count'] = len(bankA)
    import pathlib
    here = pathlib.Path(__file__).resolve().parent.parent
    def load(name):
        return [tuple(map(int, l.split())) for l in (here / name).read_text().split('\n') if l.strip()]
    p1_file = load('p1_seeds.txt'); bank_file = load('bankA.txt'); p6_file = load('p6_seeds.txt')
    as_hRed = lambda s: (s[0], s[2], s[1], s[3])
    res['P1_matches_file'] = [as_hRed(s) for s in p1] == p1_file
    res['bankA_matches_file'] = [as_hRed(s) for s in bankA] == bank_file
    def gcd(a, b):
        while b: a, b = b, a % b
        return a
    def has_prime_ge5(n):
        for q in range(5, n + 1, 2):
            if n % q == 0 and all(q % r for r in range(3, int(q ** 0.5) + 1, 2)):
                return True
        return False
    p6 = [as_hRed(s) for s in p1 if not has_prime_ge5(gcd(s[3], s[2]))]
    res['P6_count'] = len(p6); res['P6_matches_file'] = p6 == p6_file
    res['P1_excluded_from_P6'] = [as_hRed(s) for s in p1 if has_prime_ge5(gcd(s[3], s[2]))]
    # bank A anatomy
    bal = [as_hRed(s) for s in bankA if s[0] == s[2] and s[1] == s[3]]
    r2d = [as_hRed(s) for s in bankA if s[2] == 2 * s[3]]
    hR_not_bal = [as_hRed(s) for s in bankA if s[0] == s[2] and s[1] != s[3]]
    generic = [as_hRed(s) for s in bankA[1:] if not (s[0] == s[2] and s[1] == s[3]) and s[2] != 2 * s[3]]
    res['bankA_balanced'] = len(bal); res['bankA_R_eq_2d_incl_root'] = len(r2d)
    res['bankA_h_eq_R_not_balanced'] = hR_not_bal; res['bankA_generic'] = len(generic)
    # each balanced bank-A seed = J(O(y)) for a physical y, and J(y) is a non-balanced bank-A seed or chain to one
    stset = set(st); nest = {}
    bank_set = set(bankA)
    for s in bankA:
        if not (s[0] == s[2] and s[1] == s[3]): continue
        found = None
        for x in st:
            if J(O(x)) == s: found = x; break
        if found is None: nest[str(as_hRed(s))] = None; continue
        tgt = J(found); chain = [as_hRed(tgt)]
        while tgt[0] == tgt[2] and tgt[1] == tgt[3]:
            nxt = None
            for x in st:
                if J(O(x)) == tgt: nxt = x; break
            tgt = J(nxt); chain.append(as_hRed(tgt))
        nest[str(as_hRed(s))] = {'y': found, 'J(y)_chain': chain, 'end_in_bankA_nonbalanced': tgt in bank_set}
    res['balanced_nesting'] = nest
    res['balanced_all_nested_into_nonbalanced_bankA'] = all(v and v['end_in_bankA_nonbalanced'] for v in nest.values())
    res['bankA_lists'] = {'balanced': bal, 'R_eq_2d': r2d, 'generic': generic}
    print(json.dumps(res, indent=1, default=list))

if __name__ == '__main__':
    main()
