"""Controls for the high-bit shell theorem; not a spectrum-cover census."""
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


@lru_cache(None)
def g(n):
    if n < 2:
        return 1
    q, bit = divmod(n, 2)
    return g(q) + (q + (q & 1) if not bit else 0)


def direct_g(n):
    bits = [(n >> j) & 1 for j in range(n.bit_length())]
    return 1 + sum(
        bits[i] * (1 - bits[j]) * (2 if i-j == 1 else 2**(i-j-1))
        for i in range(len(bits)) for j in range(i)
    )


def reverse(n, width):
    return int(f'{n:0{width}b}'[::-1], 2)


def reduced(sig):
    d, kp, km, ip, im, op, om, jp, jm = sig
    return d+jp+jm, d, d+kp, ip, op, -im, -om


def phi(z, a, b):
    h, d, L, i, o, s, e = z
    return h+d*(g(a)+g(b)-2)+L*a*b+i*a+o*b+s*(a&1)+e*(b&1)


def construct_padded(core, cut, a, b, left_width):
    """Literal compiler with a prescribed left width (including zeroes)."""
    left = [int(x) for x in f'{a:0{left_width}b}']
    right = [1-int(x) for x in bin(b)[2:][::-1]] if b else []
    csize = len(core)
    old = left_width+csize+len(right)
    rows = [0]*(old+1)
    for u in range(old):
        for v in range(u+1, old):
            forward = True
            if left_width <= u < v < left_width+csize:
                forward = bool(core[u-left_width] >> (v-left_width) & 1)
            winner, loser = (u, v) if forward else (v, u)
            rows[winner] |= 1 << loser
    colours = left+[(cut >> v) & 1 for v in range(csize)]+right
    for v, colour in enumerate(colours):
        winner, loser = (old, v) if colour else (v, old)
        rows[winner] |= 1 << loser
    return rows


def literal_h(rows):
    """End-subset recurrence, separate from the inherited compiler."""
    size = len(rows)
    full = (1 << size)-1
    @lru_cache(None)
    def endings(mask, v):
        if mask == 1 << v:
            return 1
        old = mask ^ (1 << v)
        return sum(endings(old, u) for u in range(size)
                   if old >> u & 1 and rows[u] >> v & 1)
    return sum(endings(full, v) for v in range(size))


def main():
    shell_count = 0
    for n in range(1 << 15):
        require(g(n) == direct_g(n), ('pair expansion', n))
    histograms = {}
    for k in range(1, 15):
        hist = Counter()
        for n in range(1 << k):
            r = reverse(n, k)
            gap = g((1 << k)+n)-g(n)
            require(gap == (1 << k)-r-(n >> (k-1)), ('shell', k, n))
            hist[gap] += 1
            shell_count += 1
        expected = {j: (1 if j in (0, 1 << k) else 2)
                    for j in range(0, (1 << k)+1, 2)}
        require(dict(hist) == expected, ('multiplicity', k))
        if k <= 4:
            histograms[k] = dict(sorted(hist.items()))
    rectangle_count = 0
    for k in range(1, 10):
        for a in range(k+1):
            for b in range(k-a+1):
                counts = Counter((n >> (k-a), reverse(n, k) >> (k-b))
                                 for n in range(1 << k))
                require(len(counts) == 1 << (a+b), ('rectangle support', k,a,b))
                require(set(counts.values()) == {1 << (k-a-b)}, ('rectangle count',k,a,b))
                rectangle_count += len(counts)

    data = json.loads((ROOT/'05-knowledge/results/chanpak_five_core_bridges_20261110/cover/signatures.json').read_text())
    signature_checks = 0
    for row in data:
        z = reduced(row['signature'])
        _, d, L, i, _, _, _ = z
        for B in range(8):
            for k in range(1, 9):
                for n in range(1 << k):
                    expected = (L*B+i)*(1 << k)+d*((1 << k)-reverse(n,k)-(n >> (k-1)))
                    require(phi(z,(1 << k)+n,B)-phi(z,n,B) == expected,
                            ('physical signature shell', row['core_index'], B,k,n))
                    signature_checks += 1

    literal_checks = 0
    for row in data:
        z = reduced(row['signature'])
        for k in (1, 2, 3):
            for n in {0, (1 << k)-1, (1 << k)//2}:
                for B in (0, 1):
                    old = construct_padded(row['core'], row['cut'], n,B,k)
                    new = construct_padded(row['core'], row['cut'], (1 << k)+n,B,k+1)
                    require(literal_h(old) == phi(z,n,B), ('literal old',row,n,B))
                    require(literal_h(new) == phi(z,(1 << k)+n,B), ('literal new',row,n,B))
                    # Deleting new vertex0 returns old exactly, with labels shifted.
                    restored = [(r >> 1) for r in new[1:]]
                    require(restored == old, ('literal deletion',row,n,B))
                    literal_checks += 2

    hostile = {'k':1, 'n':1, 'actual':g(3)-g(1), 'without_top_bit':2-reverse(1,1)}
    require(hostile['actual'] != hostile['without_top_bit'], 'hostile must fail')
    report = {
        'status':'FINITE-EXACT controls for proofs in DYADIC_SHELLS.md; target cover OPEN',
        'g_pair_formula_universe':'0<=n<32768',
        'shell_universe':'all1<=k<=14 and all0<=n<2^k',
        'shell_identity_checks':shell_count,
        'small_jump_histograms':histograms,
        'rectangle_universe':'all1<=k<=9, a,b>=0, a+b<=k, all cells',
        'rectangle_cells_checked':rectangle_count,
        'inherited_signature_count':len(data),
        'signature_shell_universe':'all44 inherited rows; B0..7; k1..8; all low words',
        'signature_shell_checks':signature_checks,
        'literal_H_checks':literal_checks,
        'literal_universe':'44 rows; k1,2,3; distinct n in{0,2^k-1,2^(k-1)}; B0,1; old/new',
        'missing_top_bit_hostile':hostile,
    }
    Path(__file__).with_name('shell_checks.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'small_jump_histograms'},sort_keys=True))


if __name__ == '__main__':
    main()
