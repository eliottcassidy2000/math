#!/usr/bin/env python3
"""FINITE-EXACT probe of a new one-step paired budget, not a universal proof."""
import json
import check_root_rays as r


def evaluate(u, N):
    data = []
    for a in (1, 2):
        C, D, tau = r.ray(u, a)
        zs = []; k = 0
        while r.source(C, D, k) < 4*N:
            zs.append(r.scale(r.literal_charge(r.source(C, D, k), N), (-1)**tau)); k += 1
        ds = []; tail = (0, 0)
        for z in reversed(zs):
            v = r.add(z, tail); ds.append(r.norm(v)-r.norm(z)-r.norm(tail)); tail = v
        data.append((zs, ds[::-1], tail))
    cross = r.norm(r.add(data[0][2], data[1][2]))-r.norm(data[0][2])-r.norm(data[1][2])
    S = sum(r.norm(z) for zs, _, _ in data for z in zs)
    off = int(u % 3 == 2)
    odd, even = data[0][1], data[1][1]
    debt = max(cross, 0)+(max(0, odd[0]) if off and odd else 0)
    debt += sum(max(0, (odd[k+off] if k+off < len(odd) else 0)+(even[k] if k < len(even) else 0))
                for k in range(max(len(odd)-off, len(even))))
    return debt, S


def main():
    peak = None; failures = []; checked = 0
    targets = [u for u in range(1, 258) if r.primitive(u)] + [1009]
    for u in targets:
        Ns = set(range(1, 2049))
        if u in (1, 3, 7, 9, 11, 17, 31, 127, 1009):
            for a in (1, 2):
                C, D, _ = r.ray(u, a)
                for i in range(37):
                    for k in range(i//3+1):
                        m = r.shell_point(C, D-3, i, k)
                        Ns.update(n for n in (m//2, m//2+1, m//4, m//4+1) if n > 0)
        for N in sorted(Ns):
            debt, S = evaluate(u, N); excess = 2*debt-S
            checked += 1
            if peak is None or excess > peak['twice_excess']:
                peak = dict(target=u, N=str(N), paired_debt=debt, source_norm=S, twice_excess=excess)
            if excess > 6 and len(failures) < 8:
                failures.append(dict(target=u, N=str(N), paired_debt=debt, source_norm=S))
    return dict(status='FINITE-EXACT; universal statement OPEN', candidate='D_pair <= S/2 + 3',
                targets=len(targets), dense_N=[1, 2048], event_shells=37,
                event_targets=[1, 3, 7, 9, 11, 17, 31, 127, 1009], tested_pairs=checked,
                peak=peak, recorded_failures=failures,
                first_arrival_scope='complete one-step native preimage fibres; ROOT source1 retains its stop tag')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
