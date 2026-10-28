#!/usr/bin/env python3
"""Hostile probe of iteration: original sources and clocks are retained."""
import json
import check_paired_budget as p
R = p.R


def evaluate(N):
    z = {}; mass = {}; ancestry = {}
    for s in range(1, 4*N, 2):
        if not R.primitive(s):
            continue
        v, t, _ = R.native(s)
        u, t2, _ = R.native(v)
        # A stopped source is carried with clock zero; its positive sign is
        # unchanged. This is an absorbing two-level transport, not a count
        # of distinct first-arrival depths.
        if s == 1:
            t = 0
        if v == 1:
            t2 = 0
        q = R.literal_charge(s, N)
        z[v] = p.add(z.get(v, (0, 0)), q if t % 2 == 0 else p.neg(q))
        mass[u] = mass.get(u, 0)+p.Q(q)
        if q != (0, 0):
            ancestry.setdefault(u, []).append(dict(source=s, child=v, target=u,
                first_clock=t, second_clock=t2, original_charge=q, original_energy=p.Q(q)))
    out = {}
    bound = max(z, default=0)
    for u in sorted({R.native(v)[0] for v, q in z.items() if q != (0, 0)}):
        rays = []; addresses = []
        for a in (1, 2):
            C, D, t = R.ray(u, a)
            k = 0; word = []; src = []
            while R.source(C, D, k) <= bound:
                v = R.source(C, D, k)
                q = z.get(v, (0, 0))
                word.append(q if t % 2 == 0 else p.neg(q)); src.append(v); k += 1
            rays.append(word); addresses.append(src)
        sc = p.score(*rays, int(u % 3 == 2))
        out[u] = dict(N=N, target=u, outer_debt=sc[1], child_energy=sc[2],
            original_energy=mass[u], outer_words=rays, child_addresses=addresses,
            endpoint_charge=p.add(*rays[0], *rays[1]),
            original_sources=ancestry[u])
    return out


def main():
    first = None; peak = None; cases = 0
    for N in range(1, 513):
        for u, row in evaluate(N).items():
            excess = 3*row['outer_debt']-row['child_energy']; cases += 1
            if excess > 12 and first is None:
                first = row
            if peak is None or excess > peak['thrice_excess']:
                peak = dict(N=N, target=u, thrice_excess=excess)
    control = evaluate(57)[1]
    p.need((first['N'], first['target'], first['outer_debt'], first['child_energy']) == (19, 1, 6, 3),
           'minimal finite two-level hostile')
    p.need((control['outer_debt'], control['child_energy'], control['original_energy']) == (27, 12, 16),
           'even original source energy does not repair naive iteration')
    p.need(p.Q(control['endpoint_charge']) == 39, 'actual two-level endpoint exceeds one-step bound')
    return dict(status='REFUTED naive iteration; FINITE-EXACT hostile search',
        universe='all positive N<=512, all nonzero targets of two absorbing native steps from sources<4N',
        cases=cases, first_one_third_failure=first, peak_thrice_excess=peak,
        original_weight_hostile=control,
        conclusion='One-step theorem remains valid. Child phase correlations and ancestry cannot be replaced by the child diagonal energies or the original total mass alone.')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
