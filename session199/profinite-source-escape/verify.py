"""Independent exact controls for the profinite source escape theorem."""
from pathlib import Path
from math import gcd
import hashlib
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def g(n):
    """Original upward-index recursion, evaluated without recursive stack."""
    acc = 1
    while n > 1:
        if n % 2 == 0:
            n //= 2
        else:
            n = (n + 1) // 2
            acc += n - n % 2
    return acc


def eval_bits(word):
    n, f = 0, 1
    for bit in word:
        if bit == '0':
            f += n + n % 2
            n *= 2
        else:
            n = 2 * n + 1
    return n, f


def modular_word(M, u, z, extra=0):
    need(M >= 2 and M % 2 == 0 and z % 2 == 1, 'compiler inputs')
    q, k = M, 0
    while q % 2 == 0:
        q //= 2
        k += 1
    a = 2
    while pow(2, a, q) != 1 % q:
        a += 1
    delta = ((1 << a) - 1) * ((1 << k) - 1) + 1
    need(gcd(delta, M) == 2, 'modular loop content')
    r, half = (u - 1) % M, M // 2
    j = 0 if half == 1 else (((z - g(r + 1)) // 2) * pow(delta // 2, -1, half)) % half
    j += (1 + extra) * half
    word = ('1' * a + '0' * k) * j + format(r, f'0{(M-1).bit_length()}b')
    n, gg = eval_bits(word)
    need(gg == g(n + 1), 'independent binary evaluator')
    need((n + 1 - u) % M == 0 and (gg - z) % M == 0, 'joint residues')
    return n + 1, gg, word


def interface(M, D0, U0, extra=0):
    D0, U0 = D0 % M, U0 % M
    need(D0 % 2 == U0 % 2 == 1, 'odd interface inputs')
    u = (U0 + 1) // 2
    z = (D0 - u + u % 2) % M
    L, gg, word = modular_word(M, u, z, extra)
    s = (gg + L - L % 2, 1, 2 * L - 1, 1)
    need(L >= 2 and s[0] % M == D0 and s[2] % M == U0, 'interface residues')
    need(s[0] == g(s[2]), 'unit-deletion indexed source equality')
    return s, word


def move(s, letter):
    D, E, U, d = s
    if letter == 'O':
        return D, 0, 2 * U, d
    if letter == 'Z':
        return D - E + U, d, 2 * U - d, d
    if letter == 'J':
        return D - E + U, D - E, D + U, D
    raise RuntimeError('letter')


def bank(cap):
    seed = (1, 0, 2, 1)
    out = {seed: ''}
    queue = [seed]
    for s in queue:
        for letter in 'OZJ':
            t = move(s, letter)
            need(t[2] > s[2], 'strict cover increase')
            if t[2] <= cap and t not in out:
                out[t] = out[s] + letter
                queue.append(t)
    return out


def form(kind, v, w):
    D, E, U, d = v
    g0, s, L, f = w
    if kind == 'N':
        h, e, r, k = D-E+U, D-E+d, D+3*U-d, U+E
    else:
        h, e, r, k = D+2*U, U+d, 2*D-E+6*U-d, U+E
    return h*g0-e*s+r*L-k*f


def source_word(U):
    need(U >= 2, 'unit source cover bound')
    reverse = []
    while U > 2:
        if U % 2 == 0:
            reverse.append('O')
            U //= 2
        else:
            reverse.append('Z')
            U = (U+1)//2
    return ''.join(reversed(reverse))


def append(T, root, letter):
    n = len(T)
    out = list(T) + [0]
    for i in range(n):
        if letter != 'O' and i == root:
            out[n] |= 1 << i
        else:
            out[i] |= 1 << n
    return out, n if letter == 'J' else root


def build(word):
    T, root = [2, 0], 0
    for letter in word:
        T, root = append(T, root, letter)
    return T, root


def literal(T, root=0):
    n, full = len(T), (1 << len(T))-1
    end = [[0]*n for _ in range(full+1)]
    start = [[0]*n for _ in range(full+1)]
    H = [1] + [0]*full
    for mask in range(1, full+1):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            rest = mask ^ (1 << v)
            if rest == 0:
                end[mask][v] = start[mask][v] = 1
            else:
                end[mask][v] = sum(end[rest][u] for u in range(n) if (rest >> u) & 1 and (T[u] >> v) & 1)
                start[mask][v] = sum(start[rest][u] for u in range(n) if (rest >> u) & 1 and (T[v] >> u) & 1)
        H[mask] = sum(end[mask])
    U = sum(H[mask] * start[full ^ mask][root] for mask in range(full+1))
    # Every directed cycle is counted once by fixing vertex 0 first.
    fixed = [[0]*n for _ in range(full+1)]
    fixed[1][0] = 1
    for mask in range(1, full+1, 2):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            value = fixed[mask][v]
            for u in range(n):
                if not (mask >> u) & 1 and (T[v] >> u) & 1:
                    fixed[mask | (1 << u)][u] += value
    cycle = sum(fixed[full][v] for v in range(1, n) if T[v] & 1)
    return (H[full], end[full][root], U, H[full ^ (1 << root)]), start[full][root], cycle


def triangle_graft(W, wr):
    B, br = build('Z')
    left, new = append(B, br, 'J')
    entries = {br, new}
    n, m = len(left), len(W)
    G = [((1 << m)-1) ^ (1 << v) ^ row for v, row in enumerate(W)]
    T = list(left) + [row << n for row in G]
    for x in range(n):
        for y in range(m):
            if y == wr and x in entries:
                T[n+y] |= 1 << x
            else:
                T[x] |= 1 << (n+y)
    return T


def main():
    b = (3,1,3,1)
    interfaces = 0
    small = []
    for M in range(2, 17, 2):
        for D0 in range(1, M, 2):
            for U0 in range(1, M, 2):
                s, bits = interface(M, D0, U0)
                word = source_word(s[2])
                q = (1,0,2,1)
                for letter in word:
                    q = move(q, letter)
                need(q == s, 'physical unit-interface recipe')
                need(form('N', b, s) == 5*s[0]+11*s[2]-7, 'fixed triangle row')
                interfaces += 1
                if len(word) <= 5:
                    small.append((s, word))
    scalar_tests = 0
    for m in range(1, 33):
        for r in range(1, 2*m, 2):
            z = (r-9)//2
            y = z % 5
            x = (z-11*y)//5
            need(5*x+11*y == z, 'Bezout target choice')
            q, _ = interface(2*m, 1+2*x, 1+2*y)
            h = form('N', b, q)
            need((h-r) % (2*m) == 0 and h % 2 == 1, 'actual C2 residue')
            scalar_tests += 1
    ghost = []
    for m in [1,2,3,4,6,8,12,16,24,32,48,64,96]:
        s, bits = interface(2*m, 1, 3)
        q = move(s, 'O')
        need(all((a-c) % m == 0 for a,c in zip(q,(1,0,6,1))), 'ghost residue')
        h = form('N', b, q)
        need((h-67) % m == 0 and q[1] == 0 and q[3] == 1, 'ghost exact sidecars')
        if m % 24 == 0:
            need(q[0] >= m+1 and q[2] >= 2*m+2 and h >= 27*m+23, 'height escape')
        ghost.append(dict(modulus=m, source_bits=q[2].bit_length(), count_bits=h.bit_length(), recipe_hash=hashlib.sha256(bits.encode()).hexdigest()))
    all_sources = bank(11)
    need(len(all_sources) == 14, 'complete source bank size')
    support = {15}
    comparisons = 0
    for v in all_sources:
        for w in all_sources:
            for kind in ['N','M']:
                h = form(kind,v,w)
                if h <= 67:
                    support.add(h)
                comparisons += 1
    expected = [13,15,23,25,29,41,43,45,49,53,55]
    need(sorted(support) == expected and 67 not in support, 'complete C2 inverse negative')
    need(form('N',b,(1,0,6,1)) == 67, 'false formal target')
    literal_controls = []
    small = sorted(set(small), key=lambda pair: (len(pair[1]), pair[1]))[:4]
    for s, word in small:
        W, wr = build(word)
        measured, start, cycles = literal(W, wr)
        need(measured == s and start == cycles == 1, 'literal strong unit source')
        T = triangle_graft(W,wr)
        result, _, C = literal(T)
        need(result[0] == form('N',b,s) and C == 2, 'literal fixed-left C2 graph')
        literal_controls.append(dict(word=word, profile=s, graph_order=len(T), H=result[0], C=C, rows=T))
    need(len(literal_controls) >= 1, 'nonvacuous literal controls')
    here = Path(__file__).resolve().parent
    output = dict(status='PASS', interface_residue_tests=interfaces, C2_scalar_residue_tests=scalar_tests,
                  ghost_controls=ghost, source_cap=11, source_profiles=len(all_sources),
                  all_order_cutoff_comparisons=comparisons, C2_support_through_67=expected,
                  literal_controls=literal_controls,
                  dependency_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((here/'dependencies').glob('*.md'))},
                  scope='All-moduli physical congruence density and nonphysical profinite limit; no universal odd coverage.')
    print(json.dumps(output, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
