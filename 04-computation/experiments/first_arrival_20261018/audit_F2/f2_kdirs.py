#!/usr/bin/env python3
"""Audit F2, check (5) continued: the K-direction connector probe added by the session after the audit started
(kdirs.py, copy_connectors_kdirs.json, THM-4690 remark 'none whose first step lies in the 240 K-directions').

Independent code (f2_field.py; my own F_361 = F_19[s]/(s^2-2); my own places).  The four gamma's are rebuilt from the
(c0, c1) printed in copy_connectors_kdirs.out: gamma = j - c0 - c1 zeta.
  1. gamma vanishes at exactly one prime of each conjugate pair of K-primes above 19 (checked on all 16 homs of K);
  2. delta = conj(gamma)/gamma: unit vector; its valuations at the 24 primes of F above 19;
  3. the 240 directions zeta^l delta^(+-1), sorted as the session sorted them; recorded survivors' |t|^2 reproduced;
  4. every survivor, and every one of the 47515 + 33031 three-step tests, certified non-square at a degree-1 place
     where all residues involved are defined (places with p in a coordinate denominator are skipped).
"""
import sys, os, json, cmath, math, itertools, time
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_field import E, ONE, J, W, EXP, N, ZPOW, build, roots_of_unity, conj, norm_cc
from f2_connectors import places_all, _tab
import f2_valuations as V

HERE = os.path.dirname(os.path.abspath(__file__)); SESSION = os.path.dirname(HERE)
T0 = time.time()
def log(*a): print('[%5.0f s]' % (time.time() - T0), *a, flush=True)

GAM_C = [(9, -9), (-9, 9), (6, -2), (-6, 2)]          # (c0, c1) as printed in copy_connectors_kdirs.out


def red_or_none(x, p, emb):
    if x.d % p == 0:
        return None
    tab = _tab(p, emb)
    return sum(v * tab[k] for k, v in enumerate(x.n) if v) * pow(x.d, -1, p) % p


def main():
    O = build(); b, c, M, C = O['b'], O['c'], O['M'], O['C']
    roots = roots_of_unity()
    gammas = [J - E.scalar(c0) - ZPOW[1].scale(c1) for c0, c1 in GAM_C]
    deltas = [conj(g) / g for g in gammas]
    log('delta conj(delta) = 1 for all four:', all(norm_cc(d) == ONE for d in deltas))
    # 1. vanishing pattern on the 8 primes of K above 19 (homs of K: zeta -> z, j -> j0; Frobenius orbits)
    Khoms = [(z, j0) for z in V.Z15 for j0 in V.JR]
    Kprimes, seen = [], set()
    for hm in Khoms:
        if hm in seen: continue
        orb = (hm, (V.frob(hm[0]), V.frob(hm[1]))); seen.update(orb); Kprimes.append(hm)
    def kval(x, hm):                       # x in K (no w); evaluate with w -> anything (w-exponents absent)
        return V.hom_value(x, (hm[0], hm[1], (0, 0)))
    def conj_prime(hm):                    # the prime hit by x -> hm(conj x)
        z, j0 = hm; cz = V.fpow(z, 14); cj = V.fneg(j0)
        for q in Kprimes:
            if q == (cz, cj) or q == (V.frob(cz), V.frob(cj)): return q
    pairs = set(frozenset([q, conj_prime(q)]) for q in Kprimes)
    log('K: %d primes above 19 in %d conjugate pairs; no prime self-conjugate: %s'
        % (len(Kprimes), len(pairs), all(conj_prime(q) != q for q in Kprimes)))
    for (c0, c1), g in zip(GAM_C, gammas):
        zeros = [q for q in Kprimes if kval(g, q) == (0, 0)]
        log('gamma = j - (%d) - (%d) zeta: vanishes at %d K-prime(s); its conjugate prime not among them: %s'
            % (c0, c1, len(zeros), all(conj_prime(q) not in zeros for q in zeros)))
    covered = set()
    for g in gammas:
        for q in Kprimes:
            if kval(g, q) == (0, 0):
                covered.add(frozenset([q, conj_prime(q)]))
    log('conjugate pairs hit by the four gammas: %d of %d' % (len(covered), len(pairs)))
    # 2. valuations of delta at the 24 F-primes (my reduction; exact where decidable)
    homs = [(z, j0, w) for z in V.Z15 for j0 in V.JR for w in V.WR]
    Fprimes, seen = [], set()
    for hm in homs:
        if hm in seen: continue
        orb = (hm, tuple(V.frob(x) for x in hm)); seen.update(orb); Fprimes.append(hm)
    def typ(hm):
        z, j0, w = hm; hr = V.fadd(w, (3, 0))
        return '+j' if hr == j0 else ('-j' if hr == V.fneg(j0) else '7')
    for d in deltas:
        k = V.k19(d); dd = d.scale(19 ** k); di = d.inv(); ki = V.k19(di); ddi = di.scale(19 ** ki)
        rows = Counter()
        for hm in Fprimes:
            pole = V.hom_value(dd, hm) != (0, 0) and k == 1          # v(delta) = -1 exactly
            zero = V.hom_value(ddi, hm) != (0, 0) and ki == 1        # v(1/delta) = -1, i.e. v(delta) = +1
            rows[(typ(hm), 'v=-1' if pole else ('v=+1' if zero else 'v=0 or undecided'))] += 1
        log('delta: 19-power in denominator %d (inverse %d); F-primes by (type, valuation):' % (k, ki), sorted(rows.items()))
    # 3. the 240 directions, session order
    D240 = sorted({r * q for r in roots for q in deltas + [d.inv() for d in deltas]}, key=lambda e: e.fracs())
    log('|D240| = %d' % len(D240))
    wr = 1.5
    for _ in range(100): wr -= (wr ** 3 + 2 * wr ** 2 - 4 * wr - 4) / (3 * wr ** 2 + 4 * wr - 4)
    zc = cmath.exp(2j * cmath.pi * 2 / 15); jc = -1j * math.sqrt(11)
    pw = {(a, bb, cc_): zc ** a * jc ** bb * wr ** cc_ for (a, bb, cc_) in EXP.values()}
    def ev(x): return sum(v * pw[EXP[k]] for k, v in enumerate(x.n) if v) / x.d
    PL = places_all(40)
    Jn = json.load(open(os.path.join(SESSION, 'copy_connectors_kdirs.json')))
    allok = True; nsurv = 0
    for name, co in (('b', b), ('c', c)):
        for (i1, i2, a1, a2, qi, nt) in Jn[name]['three_step_survivors']:
            t = (C[a2] + co * M[i2]) - (C[a1] + co * M[i1]) - D240[qi]
            Nt = norm_cc(t); X = Nt * (Nt - 4)
            cert = None
            for p, embs in PL:
                if X.d % p == 0: continue
                for emb in embs:
                    r = red_or_none(X, p, emb)
                    if r and pow(r, (p - 1) // 2, p) == p - 1:
                        cert = (p, emb); break
                if cert: break
            nsurv += 1
            ok = cert is not None and abs(abs(ev(t)) ** 2 - nt) < 1e-9
            allok &= ok
            log('survivor %s (%d,%d) pts (%d,%d) q#%d: |t|^2 %.6f (recorded %.6f), certificate %s'
                % (name, i1, i2, a1, a2, qi, abs(ev(t)) ** 2, nt, cert[0] if cert else None))
    log('all %d recorded K-direction survivors reproduced and certified: %s' % (nsurv, allok))
    # 4. full recount with certification (residues; skip places where a residue is undefined)
    evD = [ev(q) for q in D240]
    P12 = PL[:12]
    for name, co in (('b', b), ('c', c)):
        pts = {(i, a): C[a] + co * M[i] for i in range(7) for a in range(len(C))}
        evp = {k: ev(v) for k, v in pts.items()}
        res = {}
        for p, embs in P12:
            for emb in embs:
                cemb = (pow(emb[0], 14, p), (-emb[1]) % p, emb[2])
                res[(p, emb)] = ({k: (red_or_none(v, p, emb), red_or_none(v, p, cemb)) for k, v in pts.items()},
                                 [(red_or_none(q, p, emb), red_or_none(q, p, cemb)) for q in D240])
        tested = 0; unc = []
        for (i1, i2) in itertools.combinations(range(7), 2):
            for a1 in range(len(C)):
                for a2 in range(len(C)):
                    t0 = evp[(i2, a2)] - evp[(i1, a1)]
                    if abs(t0) > 3 + 1e-9: continue
                    for qi, eq in enumerate(evD):
                        n1 = abs(t0 - eq) ** 2
                        if n1 > 4 + 1e-9 or n1 < 1e-12: continue
                        tested += 1; ok = False
                        for key, (rp, rd) in res.items():
                            p = key[0]
                            ry, rx, rq = rp[(i2, a2)], rp[(i1, a1)], rd[qi]
                            if None in ry or None in rx or None in rq: continue
                            nr = (ry[0] - rx[0] - rq[0]) * (ry[1] - rx[1] - rq[1]) % p
                            X = nr * (nr - 4) % p
                            if X and pow(X, (p - 1) // 2, p) == p - 1:
                                ok = True; break
                        if not ok: unc.append((i1, i2, a1, a2, qi))
        still = 0
        for (i1, i2, a1, a2, qi) in unc:
            t = pts[(i2, a2)] - pts[(i1, a1)] - D240[qi]; Nt = norm_cc(t); X = Nt * (Nt - 4)
            found = False
            for p, embs in PL:
                if X.d % p == 0: continue
                for emb in embs:
                    r = red_or_none(X, p, emb)
                    if r and pow(r, (p - 1) // 2, p) == p - 1: found = True; break
                if found: break
            if not found: still += 1
        log('K-directions, coefficient %s: three-step tests %d ; uncertified at 12 primes %d ; after 40 primes %d'
            % (name, tested, len(unc), still))


if __name__ == '__main__':
    main()
