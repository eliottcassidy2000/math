"""Lift Frobenius on an inherited native response, retaining its nilpotent clock."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) % 3 for col in zip(*b)] for row in a]


def ident(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def mpow(a, e):
    v = ident(len(a))
    while e:
        if e & 1:
            v = mm(v, a)
        a = mm(a, a)
        e //= 2
    return v


def inverse(a):
    n = len(a)
    rows = [list(row)+unit for row, unit in zip(a, ident(n))]
    for j in range(n):
        r = next((i for i in range(j, n) if rows[i][j] % 3), None)
        need(r is not None, 'invertible observation chart')
        rows[j], rows[r] = rows[r], rows[j]
        inv = pow(rows[j][j] % 3, -1, 3)
        rows[j] = [x*inv % 3 for x in rows[j]]
        for i in range(n):
            if i != j:
                c = rows[i][j]
                rows[i] = [(x-c*y) % 3 for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def run():
    source_path = HERE.parent/'overnight_rank_return_20261029'/'periodic_defect.json'
    src = json.loads(source_path.read_text())
    q = src['ternary_recurrence_low_first']
    n = len(q)-1
    need(n == 12 and q[-1] == 1, 'inherited monic degree-twelve response polynomial')
    def times(a, b):
        v = [0]*(2*n-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                v[i+j] = (v[i+j]+x*y) % 3
        for j in range(2*n-2, n-1, -1):
            c = v[j]
            for i in range(n+1):
                v[j-n+i] = (v[j-n+i]-c*q[i]) % 3
        return v[:n]
    one = [1]+[0]*(n-1)
    x = [0, 1]+[0]*(n-2)
    def power(a, e):
        v = one
        while e:
            if e & 1:
                v = times(v, a)
            a = times(a, a)
            e //= 2
        return v
    exponent = next(a for a in range(726) if a % 242 == 3 and a % 3 == 1)
    need(exponent == 487 and pow(exponent, 5, 726) == 1, 'compatible CRT exponent and order')
    T = list(map(list, zip(*[times(x, [int(i == j) for i in range(n)]) for j in range(n)])))
    image = power(x, exponent)
    phi = list(map(list, zip(*[power(image, j) for j in range(n)])))
    cube = list(map(list, zip(*[power(x, 3*j) for j in range(n)])))
    I = ident(n)
    need(mpow(T, 726) == I and all(mpow(T, 726//p) != I for p in (2, 3, 11)), 'actual clock exact order726')
    need(mpow(phi, 5) == I and phi != I, 'lifted Frobenius exact order five')
    need(mm(phi, T) == mm(mpow(T, exponent), phi), 'full clock conjugacy with nilpotent directions retained')
    need(mm(inverse(phi), phi) == I, 'invertible full twelve-coordinate lift')
    # Literal algebra homomorphism on every basis pair, not only matrix conjugacy.
    for i in range(n):
        for j in range(n):
            a, b = [int(k == i) for k in range(n)], [int(k == j) for k in range(n)]
            lhs = [row[0] for row in mm(phi, [[v] for v in times(a, b)])]
            rhs = times([row[i] for row in phi], [row[j] for row in phi])
            need(lhs == rhs, 'literal algebra product preserved')
    from ternary_obstruction_verify import kernel
    need(kernel(cube, 3)[1] == 6, 'naive cubing collapses six of twelve dimensions')
    vals = src['values']
    H = [[vals[i+j] % 3 for j in range(n)] for i in range(n)]
    Hi = inverse(H)
    obs = mm(mm(H, phi), Hi)
    initial = [[vals[i] % 3] for i in range(n)]
    need(mm(obs, initial) == initial, 'actual initial response window fixed')
    old_selector = [1]+[0]*(n-1)
    pulled_back_selector = obs[0]
    transported_selector = mpow(obs, 4)[0]
    need(pulled_back_selector != old_selector, 'terminal observation needs its own transported selector')
    need(mm([transported_selector], obs) == [old_selector], 'contragredient selector preserves the original observation')
    values = []
    state = one
    for _ in range(726):
        values.append(sum(vals[j]*state[j] for j in range(n)) % 3)
        state = times(x, state)
    mismatch = next(r for r in range(726) if values[r] != values[exponent*r % 726])
    need(all([row[0] for row in mm(phi, [[v] for v in power(x, r)])] == power(x, exponent*r)
             for r in range(726)), 'every actual orbit state transported by exponent487')
    return {'status': 'PROVED scoped lift / FINITE-EXACT native scalar response application',
            'source': str(source_path.relative_to(HERE.parent)),
            'source_sha256_canonical': hashlib.sha256(json.dumps(src,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'modulus_polynomial': q, 'clock': T, 'clock_order': 726,
            'crt_exponent': exponent, 'field_jet_only_exponent_mod363': 124,
            'frobenius': phi, 'frobenius_order': 5, 'naive_cube_rank': 6,
            'observation_chart': H, 'observation_chart_inverse': Hi,
            'frobenius_in_observation_coordinates': obs, 'initial_window': initial,
            'original_selector': old_selector, 'pulled_back_selector': pulled_back_selector,
            'transported_selector': transported_selector,
            'first_scalar_mismatch': {'r': mismatch, 'b_r_mod3': values[mismatch],
                                      'reindexed_r': exponent*mismatch % 726,
                                      'b_reindexed_mod3': values[exponent*mismatch % 726]},
            'gates': GATES}


if __name__ == '__main__':
    receipt = run()
    (HERE/'periodic_frobenius_lift.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ('gates','crt_exponent','frobenius_order','naive_cube_rank','first_scalar_mismatch')}))
