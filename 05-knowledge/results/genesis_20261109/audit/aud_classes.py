"""AUDIT G1: THM-4985 (1) and (2) re-derived with own code (aud_core props + aud_iso).

(1) the 36 escaper extensions found by aud_core ext: classes, multiplicity, converse pairing, self-converse,
    |Aut|, H per class, lambda, twin distance, module defect.
(2) the 193 classes of escapers12_all_classes.txt: escaper status with S(10), S(11) from census T-lines,
    pairwise non-isomorphism, which are the 12 classes of (1), self-converse count, converse found,
    H statistics, (lambda, module defect, twin) triples, score sequences of the defect >= 2 classes.
"""
import os
import subprocess
from collections import Counter, defaultdict
from aud_iso import decode, encode, converse, refine, isomorphic, aut_count, classes

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'genesis_20261109'))


def S(n):
    return set(int(x) for x in open(os.path.join(HERE, 'S_%d.txt' % n)))


S10, S11 = S(10), S(11)


def props(codes, n=12):
    inp = ''.join('%d %s\n' % (n, c) for c in codes)
    r = subprocess.run([os.path.join(HERE, 'aud_core'), 'props'], input=inp, capture_output=True, text=True,
                       check=True).stdout.splitlines()
    res = []
    for line in r:
        a, b, c, d = line.split('|')
        t = a.split()
        H = int(t[1])
        kv = dict(x.split('=') for x in t[2:])
        dels = list(map(int, b.split()))
        sc = list(map(int, c.split()[1:]))
        res.append(dict(H=H, lam=int(kv['lam']), prime=int(kv['prime']), moddef=int(kv['moddef']),
                        twin=int(kv['twin']), dels=dels, sc=sc))
    return res


def escaper(p):
    return (p['H'] in S11 and all(d in S10 for d in p['dels']) and p['lam'] >= 3 and p['prime'] == 1)


def main():
    # ---------- (1)
    esc = [l.split()[1] for l in open(os.path.join(HERE, 'ext12_audit.out')) if l.startswith('E ')]
    print('(1) escaper extensions:', len(esc))
    outs = [decode(c, 12) for c in esc]
    P = props(esc)
    print('    all re-verified escapers:', all(escaper(p) for p in P))
    cl = classes(outs)
    print('    isomorphism classes:', len(cl), ' multiplicities:', sorted(Counter(len(c) for c in cl).items()))
    reps = [outs[c[0]] for c in cl]
    repP = [P[c[0]] for c in cl]
    rc = [refine(o) for o in reps]
    conv_of = {}
    for i, o in enumerate(reps):
        co = converse(o)
        cc = refine(co)
        js = [j for j in range(len(reps)) if isomorphic(co, reps[j], cc, rc[j])]
        conv_of[i] = js
    sc_cnt = sum(1 for i in conv_of if conv_of[i] == [i])
    print('    converse class found for every class:', all(len(v) == 1 for v in conv_of.values()),
          ' self-converse:', sc_cnt)
    pairs = set(tuple(sorted((i, conv_of[i][0]))) for i in conv_of)
    for a, b in sorted(pairs, key=lambda ab: repP[ab[0]]['H']):
        print('    pair H = %d / %d  |Aut| %d,%d  lambda %d,%d  twin %d,%d  moddef %d,%d  scores %s' % (
            repP[a]['H'], repP[b]['H'], aut_count(reps[a]), aut_count(reps[b]), repP[a]['lam'], repP[b]['lam'],
            repP[a]['twin'], repP[b]['twin'], repP[a]['moddef'], repP[b]['moddef'],
            ''.join(str(x) for x in sorted(repP[a]['sc']))))
    # ---------- (2)
    allc = [l.split()[1] for l in open(os.path.join(EXP, 'escapers12_all_classes.txt')) if l.startswith('12 ')]
    print()
    print('(2) classes listed:', len(allc), ' distinct codes:', len(set(allc)))
    Q = props(allc)
    ok = [escaper(p) for p in Q]
    print('    escaper status verified for all:', all(ok), ' failures:', ok.count(False))
    print('    lambda values:', Counter(p['lam'] for p in Q))
    o2 = [decode(c, 12) for c in allc]
    cl2 = classes(o2)
    print('    pairwise non-isomorphic (classes == listed):', len(cl2) == len(allc), len(cl2))
    c2 = [refine(o) for o in o2]
    # the 12 classes of (1) inside the list
    inside = 0
    for o, c in zip(reps, rc):
        if any(isomorphic(o, o2[j], c, c2[j]) for j in range(len(o2))):
            inside += 1
    print('    classes of (1) found in the list:', inside, 'of', len(reps))
    # converse
    selfc = 0
    conv_found = 0
    for i, o in enumerate(o2):
        co = converse(o)
        cc = refine(co)
        js = [j for j in range(len(o2)) if isomorphic(co, o2[j], cc, c2[j])]
        assert len(js) <= 1
        if js == [i]:
            selfc += 1
        if js:
            conv_found += 1
    print('    self-converse:', selfc, ' classes whose converse class is in the list:', conv_found)
    hv = defaultdict(list)
    for i, p in enumerate(Q):
        hv[p['H']].append(i)
    print('    H range', min(hv), max(hv), ' distinct H', len(hv), ' classes per H:', sorted(Counter(len(v) for v in hv.values()).items()))
    trip = Counter((p['lam'], p['moddef'], p['twin']) for p in Q)
    print('    (lambda, module defect, twin) ->', dict(trip))
    for i, p in enumerate(Q):
        if p['moddef'] >= 2:
            o = o2[i]
            js = [j for j in range(len(o2)) if isomorphic(converse(o), o2[j], refine(converse(o)), c2[j])]
            print('      moddef %d: H=%d scores %s self-converse %s |Aut| %d code %s' % (
                p['moddef'], p['H'], sorted(p['sc']), js == [i], aut_count(o), allc[i]))
    # sanity: dec = module defect <= lambda for every class
    print('    module defect <= lambda for all:', all(p['moddef'] <= p['lam'] for p in Q))
    with open(os.path.join(HERE, 'classes193_props.txt'), 'w') as g:
        for c, p in zip(allc, Q):
            g.write('%s H=%d lam=%d moddef=%d twin=%d sc=%s\n' % (c, p['H'], p['lam'], p['moddef'], p['twin'],
                                                               ''.join('%x' % x for x in sorted(p['sc']))))


if __name__ == '__main__':
    main()
