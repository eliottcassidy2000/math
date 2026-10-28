#!/usr/bin/env python3
"""Exact finite induction certificate, plus independent native-cutoff controls.

The universal implication uses the event-order proof in ROOT_BUDGET.md.
The finite numerical samples are controls, not its replacement.
"""
from itertools import product
from pathlib import Path
import sys
import json
import hashlib

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'root_ports_closed_carriers_20261028'))
import check_root_rays as R

W = ((1, 0), (0, 1), (-1, -1))
U = ((0, 0), (1, 0), (1, 1), (0, 0))
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def add(*zs):
    if not zs:
        return 0, 0
    return tuple(map(sum, zip(*zs)))


def neg(z):
    return -z[0], -z[1]


def mul(z, w):
    a, b = z; c, d = w
    return a*c-b*d, a*d+b*c-b*d


def Q(z):
    a, b = z
    return a*a-a*b+b*b


def coefficients(mode, r):
    ell = (1, 2)  # omega - omega^2
    if mode == 'A':
        eb = add(mul(ell, U[r]), W[(r+1) % 3])
        ef = add(eb, neg(W[(r+2) % 3]))
        ob = mul(ell, U[r]); of = add(ob, W[(r+1) % 3])
    else:
        eb = mul(ell, U[r+1]); ef = add(eb, W[(r+2) % 3])
        ob = add(mul(ell, U[r]), W[(r+1) % 3])
        of = add(ob, neg(W[(r+2) % 3]))
    return ef, eb, of, ob


def words(mode, r, alpha, m, n, eps, head):
    """head=0 absent, 1 passed its event, 2 not yet passed (m=eps=0)."""
    ef, eb, of, ob = coefficients(mode, r)
    beta = (2*alpha+1) % 3
    odd = []; even = []
    if head:
        odd.append(neg(mul(W[(beta+1) % 3], of if head == 1 else ob)))
    for k in range(m+n+eps):
        even.append(mul(W[(alpha+k) % 3], ef if k < m+eps else eb))
        odd.append(neg(mul(W[(beta-k) % 3], of if k < m else ob)))
    return odd, even


def increments(zs):
    tail = (0, 0); ds = []
    for z in reversed(zs):
        total = add(z, tail)
        ds.append(Q(total)-Q(z)-Q(tail)); tail = total
    return ds[::-1], tail


def score(odd, even, head):
    od, ot = increments(odd); ed, et = increments(even)
    cross = Q(add(ot, et))-Q(ot)-Q(et)
    D = max(0, cross)+(max(0, od[0]) if head and od else 0)
    off = int(bool(head))
    for k in range(max(len(odd)-off, len(even))):
        D += max(0, (od[k+off] if k+off < len(od) else 0)+(ed[k] if k < len(ed) else 0))
    S = sum(map(Q, odd+even))
    return 2*D-S, D, S


def certificate():
    rows = []; slopes = set(); peak = None; stronger_peak = None
    for mode, r, alpha, m, n, eps, head in product('AB', range(3), range(3), range(3), range(3), range(2), range(3)):
        if head == 2 and (m or eps):
            continue
        args = mode, r, alpha, m, n, eps, head
        base = score(*words(*args), head)
        need(base[0] <= 6, ('base', args, base))
        stronger = 3*base[1]-base[2]
        need(stronger <= 12, ('one-third base', args, base))
        ds = []
        for dm, dn in ((3, 0), (0, 3)):
            if head == 2 and dm:
                continue
            extended = score(*words(mode, r, alpha, m+dm, n+dn, eps, head), head)
            delta = extended[0]-base[0]
            need(delta <= 0, ('induction slope', args, dm, delta))
            stronger_delta = 3*(extended[1]-base[1])-(extended[2]-base[2])
            need(stronger_delta <= 0, ('one-third induction slope', args, dm, stronger_delta))
            # A new complete phase triple has zero sum in EACH ray. Earlier
            # suffixes, initial cross term, and orphan-head term are unchanged.
            old = words(*args)
            new = words(mode, r, alpha, m+dm, n+dn, eps, head)
            need(all(add(*a) == add(*b) for a, b in zip(old, new)), 'separate ray charge preserved')
            twice = score(*words(mode, r, alpha, m+2*dm, n+2*dn, eps, head), head)
            need(twice[0] == base[0]+2*delta, 'two-step recurrence control')
            ds.append((dm, dn, delta)); slopes.add((mode, dm, delta))
        rows.append((args, base, ds))
        if peak is None or base[0] > peak[0]:
            peak = base[0], args, base
        if stronger_peak is None or stronger > stronger_peak[0]:
            stronger_peak = stronger, args, base
    raw = json.dumps(rows, separators=(',', ':')).encode()
    return dict(base_rows=len(rows), maximum_twice_excess=peak, slopes=sorted(slopes),
                maximum_thrice_excess=stronger_peak,
                certificate_sha256=hashlib.sha256(raw).hexdigest())


def trim(zs):
    zs = list(zs)
    while zs and zs[-1] == (0, 0):
        zs.pop()
    return zs


def native_words(u, N):
    data = []
    for a in (1, 2):
        C, D, tau = R.ray(u, a)
        zs = []; k = 0
        while R.source(C, D, k) < 4*N:
            z = R.literal_charge(R.source(C, D, k), N)
            zs.append(z if tau % 2 == 0 else neg(z)); k += 1
        data.append(zs)
    return data


def decode(u, N, literal):
    C, _, tau = R.ray(u, 2)
    alpha = ((C-7)//9) % 3
    has_head = int(u % 3 == 2)
    j = -1
    while R.shell_point(C, 4, j+1, 0) < 4*N:
        j += 1
    for i in (j-1, j):
        if i < 0:
            continue
        L, r = i//3+1, i % 3
        for mode in 'AB':
            if mode == 'A':
                he = sum(R.shell_point(C, 4, i, k) < 2*N for k in range(L))
                ho = sum(R.shell_point(2*C, 2, i, k) < 4*N for k in range(L))
                head_passed = None
                # k=-1 is integral only in the actual head case. Avoid a
                # negative exponent: C/32 is an integer in that case.
                if has_head:
                    head_passed = ((C//32-2)*4**(i+3)-3)//9 < 4*N
            else:
                he = sum(R.shell_point(C, 4, i+1, k) < 4*N for k in range(L))
                ho = sum(R.shell_point(2*C, 2, i, k) < 2*N for k in range(L))
                if has_head:
                    head_passed = ((C//32-2)*4**(i+3)-3)//9 < 2*N
            eps = he-ho; m = ho; n = L-he
            head = (1 if head_passed else 2) if has_head else 0
            if eps not in (0, 1) or n < 0 or head == 2 and (m or eps):
                continue
            args = mode, r, alpha, m, n, eps, head
            proposed = words(*args)
            if all(trim(a) == trim(b) for a, b in zip(literal, proposed)):
                return args
    # Before the first full paired shell, at most the odd head and even
    # source zero are present. Their Cauchy--Schwarz bound costs at most 3.
    need(sum(bool(Q(z)) for ray in literal for z in ray) <= 2,
         ('native word absent from grammar', u, N, literal, j))
    need(all(Q(z) <= 3 for ray in literal for z in ray), 'small-cutoff source bound')
    return None


def controls():
    count = 0; decoded = 0; modes = {}; peak = None; maxbits = 0
    targets = [u for u in range(1, 258) if R.primitive(u)]+[1009, 2**61-1]
    event_targets = {1, 3, 7, 9, 11, 17, 31, 127, 1009, 2**61-1}
    for u in targets:
        Ns = set(range(1, 513))
        if u in event_targets:
            for a in (1, 2):
                C, D, _ = R.ray(u, a)
                for i in list(range(24))+[35, 62, 95, 128, 191, 240]:
                    for k in sorted({0, i//6, i//3}):
                        p = R.shell_point(C, D-3, i, k)
                        for c in (2, 4):
                            Ns.update(n for n in (p//c, p//c+1) if n > 0)
        for N in sorted(Ns):
            actual = native_words(u, N)
            state = decode(u, N, actual)
            sc = score(*actual, int(u % 3 == 2))
            need(sc[0] <= 6, ('literal inequality control', u, N, sc))
            need(3*sc[1]-sc[2] <= 12, ('one-third literal inequality', u, N, sc))
            if state:
                need(score(*words(*state), state[-1]) == sc, 'score survives zero padding')
                decoded += 1; modes[str((state[0], state[-1]))] = modes.get(str((state[0], state[-1])), 0)+1
            count += 1; maxbits = max(maxbits, N.bit_length())
            if peak is None or sc[0] > peak['twice_excess']:
                peak = dict(target=u, N=str(N), twice_excess=sc[0], D_pair=sc[1], S=sc[2])
    # Independent source-first native enumeration: no inverse-ray formula.
    source_first = 0
    for N in (1, 2, 3, 7, 16, 31, 64, 128):
        buckets = {}
        for s in range(1, 4*N, 2):
            if not R.primitive(s):
                continue
            u, t, a = R.native(s)
            z = R.literal_charge(s, N)
            buckets.setdefault(u, {1: [], 2: []})[a].append((s, z if t % 2 == 0 else neg(z)))
        for u, rays in buckets.items():
            observed = [[z for _, z in sorted(rays[a])] for a in (1, 2)]
            need(observed == native_words(u, N), 'source-first inverse fibre completeness')
            source_first += 1
    # The two inherited hostile families: isolated-ray half bound fails;
    # paired debt is unbounded even though the new weighted bound holds.
    hostile = []
    for h in (0, 1, 5, 20, 100):
        for u, N, expected in (
            (1, (8*4**(9*h+10)-8186)//9, (4, 12*h+16)),
            (7, (22*4**(9*h+11)-65530)//9, (2*h+8, 12*h+14)),
            (7, (44*4**(9*h+2)-29)//9, (h, 3*h+3))):
            sc = score(*native_words(u, N), 0)
            need(sc[1:] == expected, 'inherited sharp-family formula')
            need(decode(u, N, native_words(u, N)) is not None, 'hostile belongs to grammar')
            hostile.append(dict(target=u, h=h, D_pair=sc[1], S=sc[2]))
    return dict(native_cutoffs=count, decoded_cutoffs=decoded, decoded_modes=modes,
                source_first_fibres=source_first, maximum_integer_bits=maxbits,
                sharp_control=peak, hostile_families=hostile)


def main():
    proof = certificate(); samples = controls()
    return dict(status='PROVED by event grammar and finite induction; controls FINITE-EXACT',
        theorem='D_pair <= min(S/2 + 3, S/3 + 4) for complete one-step native fibres',
                proof=proof, controls=samples, checks=CHECKS,
                scope='native unit-counting source charges and original terminal phase; no multi-depth or universal-entry conclusion')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
