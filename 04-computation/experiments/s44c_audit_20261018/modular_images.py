#!/usr/bin/env python3
"""Audit check (e)(v): projective images of the modular representations of M(2,5) and M(2,7) (numerical, exact up to
rounding in a finite group enumeration).  S and T from the standard minimal-model formulas
    S_(r,s),(r',s') = 2 sqrt(2/(p p')) (-1)^(1 + s r' + r s') sin(pi p r r'/p') sin(pi p' s s'/p),
    T = diag exp(2 pi i (h_(r,s) - c/24)),
with Kac-table representatives (r,s), 1 <= r <= p'-1, 1 <= s <= p-1, r p - s p' > 0 ... (one per pair).
The group <S, T> modulo scalars is enumerated by BFS with a rounded projective normal form.
"""
import cmath, math, itertools, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def minimal_model(p, pp):
    # Virasoro M(p, pp) with c = 1 - 6 (p - pp)^2 / (p pp); Kac labels (r, s), 1 <= r <= pp - 1, 1 <= s <= p - 1,
    # identified (r, s) ~ (pp - r, p - s)
    c = 1 - 6 * (p - pp) ** 2 / (p * pp)
    labels = []
    seen = set()
    for r in range(1, pp):
        for s in range(1, p):
            key = min((r, s), (pp - r, p - s))
            if key not in seen:
                seen.add(key)
                labels.append(key)
    h = {(r, s): ((r * p - s * pp) ** 2 - (p - pp) ** 2) / (4 * p * pp) for (r, s) in labels}
    n = len(labels)
    S = np.zeros((n, n), dtype=complex)
    for i, (r, s) in enumerate(labels):
        for j, (r2, s2) in enumerate(labels):
            S[i, j] = 2 * math.sqrt(2 / (p * pp)) * (-1) ** (1 + s * r2 + r * s2) * \
                math.sin(math.pi * p * r * r2 / pp) * math.sin(math.pi * pp * s * s2 / p)
    T = np.diag([cmath.exp(2j * math.pi * (h[l] - c / 24)) for l in labels])
    return labels, h, c, S, T


def proj_key(M, nd=6):
    flat = M.flatten()
    k = next(i for i, x in enumerate(flat) if abs(x) > 1e-7)
    ph = flat[k] / abs(flat[k])
    Mn = flat / ph
    return tuple((round(x.real, nd) + 0.0, round(x.imag, nd) + 0.0) for x in Mn)


def group_order(gens, limit=100000):
    n = gens[0].shape[0]
    I = np.eye(n, dtype=complex)
    seen = {proj_key(I): I}
    frontier = [I]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = s @ g
                k = proj_key(h)
                if k not in seen:
                    seen[k] = h
                    new.append(h)
                    if len(seen) > limit:
                        return None
        frontier = new
    return len(seen)


if __name__ == '__main__':
    out = {}
    for (p, pp) in ((2, 5), (2, 7), (2, 9), (3, 4)):
        labels, h, c, S, T = minimal_model(p, pp)
        unit = np.allclose(S @ S.conj().T, np.eye(len(labels)))
        S2 = np.allclose(S @ S, np.eye(len(labels)))
        ST3 = (S @ T) @ (S @ T) @ (S @ T)
        # projective relation (ST)^3 ~ scalar
        k = ST3[0, 0]
        st3_scalar = np.allclose(ST3, k * np.eye(len(labels)))
        tord = next(m for m in range(1, 1000) if np.allclose(np.linalg.matrix_power(T, m) / np.linalg.matrix_power(T, m)[0, 0],
                                                                np.eye(len(labels))))
        order = group_order([S, T])
        out[f'M({p},{pp})'] = {'c': c, 'labels': labels, 'h': {str(k): v for k, v in h.items()}, 'S_unitary': bool(unit),
                               'S2_identity': bool(S2), 'ST3_scalar': bool(st3_scalar), 'T_projective_order': tord,
                               'projective_image_order': order}
        print(f'M({p},{pp}): c = {c:.6f}, S unitary {unit}, S^2 = 1 {S2}, (ST)^3 scalar {st3_scalar}, '
              f'T projective order {tord}, |<S,T>/scalars| = {order}')
    json.dump(out, open(os.path.join(HERE, 'modular_images.json'), 'w'), indent=1, default=str)
