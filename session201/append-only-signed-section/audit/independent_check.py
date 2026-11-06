#!/usr/bin/env python3
"""Independent audit in nonnegative (a,f,e,p), not author's raw coordinates.
All checks survive python -O. Finite universes and conditional completeness
bounds are printed in JSON. The accompanying audit supplies the all-height proof.
"""
from collections import deque
from hashlib import sha256
from pathlib import Path
import argparse
import json

ROOT = (0, 0, 0, 1)
LETTERS = 'OZJ'

def need(ok, info):
    if not ok:
        raise RuntimeError(info)

def steps(s):
    a, f, e, p = s
    return ((a + e, a + 2*f + e + 2*p, 0, e + p),
            (2*a + f + e + p, f, e + p, 0),
            (f + p, 0, a + e + p, e))

def word(s, w):
    for letter in w:
        s = steps(s)[LETTERS.index(letter)]
    return s

def raw(s):
    a, f, e, p = s
    return (a + 2*e + p, e, a + f + 2*e + 2*p, e + p)

def positive(s):
    D, E, U, d = s
    return (D-E-d, U-D+E-d, E, d-E)

def raw_steps(s):
    D, E, U, d = s
    return ((D, 0, 2*U, d), (D-E+U, d, 2*U-d, d),
            (D-E+U, D-E, D+U, D))

def D(s):
    a, f, e, p = s
    return a + 2*e + p

def U(s):
    a, f, e, p = s
    return a + f + 2*e + 2*p

def deletion(s):
    return s[2] + s[3]

def section(s):
    _, f, e, p = s
    return f == 0 and e == 2*p

def target(s):
    return 6*U(s) - 7*D(s)

def bank(s, cap):
    need(U(s) <= cap, ('source outside cap', s, cap))
    seen = {s: ''}
    todo = deque([s])
    while todo:
        parent = todo.popleft()
        for letter, child in zip(LETTERS, steps(parent)):
            if U(child) <= cap and child not in seen:
                seen[child] = seen[parent] + letter
                todo.append(child)
    return seen

def short(s):
    return {target(word(s, 'J'*j)) for j in range(3)
            if section(word(s, 'J'*j))}

def bound(s, N):
    return max(34*N-42, 14*U(s)+6*N)

def formal_checks():
    for k in range(4):
        r = tuple(int(i == k) for i in range(4))
        s = positive(r)
        need(raw(s) == r, ('coordinate inverse', k))
        need(tuple(map(raw, steps(s))) == raw_steps(r), ('operation basis', k))
        Ds, Es, Us, ds = r
        w = raw(word(s, 'OJJJJ'))
        need(w == (5*Ds+12*Us, 2*Ds+4*Us, 6*Ds+14*Us, 3*Ds+6*Us), ('W basis', k))
        need(6*w[2]-7*w[0] == Ds, ('W target', k))
    # Independent full formal basis of the f=0 plane.
    for s in ((1,0,0,0), (0,0,1,0), (0,0,0,1)):
        r = word(s, 'Z')
        need(U(r) == D(r) and r[2] == deletion(s), ('Z balancing', s))
        q = word(s, 'ZJJJ')
        need(section(q) and target(q) == deletion(s), ('deletion realization', s))
        need(U(q) == 14*U(s)-8*deletion(s), ('deletion height', s))
    return {'raw_basis_vectors': 4, 'f_zero_basis_vectors': 3}

def geometry_checks(cap):
    b = bank(ROOT, cap)
    counts = {'profiles': len(b), 'f_zero': 0, 'f_positive': 0, 'section': 0,
              'short_hits_by_j': [0, 0, 0], 'edge_checks': 0}
    for s in b:
        need(min(s) >= 0 and deletion(s) >= 1, ('physical coordinates', s))
        if section(s):
            counts['section'] += 1
            need(target(s) > 0, ('nonpositive actual section', s))
        if s[1] == 0:
            counts['f_zero'] += 1
            a, f, e, p = s
            faces = (e == 2*p, e == a+p, a == 0)
            values = (p-a, a, p)
            for j in range(3):
                q = word(s, 'J'*j)
                need(section(q) == faces[j], ('short face', s, j))
                if faces[j]:
                    counts['short_hits_by_j'][j] += 1
                    need(target(q) == values[j] and values[j] > 0, ('short readout', s, j))
        else:
            counts['f_positive'] += 1
            need(short(s) <= {D(s)}, ('positive-f short anomaly', s))
        need(tuple(map(raw, steps(s))) == raw_steps(raw(s)), ('operation mismatch', s))
        for letter, v in zip(LETTERS, steps(s)):
            counts['edge_checks'] += 1
            need(U(v) > U(s), ('nonincreasing height', s, letter))
            need(D(v) >= D(s) and deletion(v) >= deletion(s), ('monotonicity', s, letter))
            if letter == 'Z':
                need(2*D(v)-U(v) >= 3*deletion(v) >= 3, ('Z core inequality', s))
            if letter == 'J':
                need(2*D(v)-U(v) >= 3, ('J core inequality', s))
    counts['height_cap'] = cap
    return counts

def classification_checks(seed_cap, N):
    seeds = bank(ROOT, seed_cap)
    total_forward = total_core = max_forward_cap = max_balanced_cap = 0
    positive_f = 0
    component_witnesses = {}
    larger_bound_checks = 0
    for index, (s, prefix) in enumerate(seeds.items()):
        C = short(s)
        m = D(s) if s[1] else min({deletion(s)} | C)
        forward_cap = bound(s, max(N,m))
        fwd = bank(s, forward_cap)
        total_forward += len(fwd)
        max_forward_cap = max(max_forward_cap, forward_cap)
        vals = {target(q) for q in fwd if section(q)}
        need(vals and min(vals) == m, ('minimum', prefix, m, sorted(vals)))
        observed = {n for n in vals if 1 <= n <= N}
        cores = bank(s, max(U(s), 2*N-3))
        total_core += len(cores)
        tails = {D(v) for v in cores if D(v) <= N}
        extra = {deletion(s)} if s[1] == 0 and deletion(s) <= N else set()
        expected = tails | extra | {n for n in C if n <= N}
        need(observed == expected, ('classification', prefix, sorted(observed ^ expected)))
        if s[1] > 0:
            positive_f += 1
            need(observed == tails, ('positive-f tail equality', prefix))
        # Realization checks give actual words for every RHS contribution.
        for v in cores:
            if D(v) <= N:
                q = word(v, 'OJJJJ')
                need(section(q) and target(q) == D(v) and U(q) <= bound(s,N), ('tail realization', prefix, v))
        if extra:
            q = word(s, 'ZJJJ')
            need(section(q) and target(q) == deletion(s) and U(q) <= bound(s,N), ('deletion witness', prefix))
        for j in range(3):
            q = word(s, 'J'*j)
            if section(q) and target(q) <= N:
                need(U(q) <= bound(s,N), ('short height', prefix, j))
        # Endpoint counts of balanced descendants, separate finite search.
        balanced_cap = max(5*N-6, N+2*U(s))
        max_balanced_cap = max(max_balanced_cap, balanced_cap)
        balanced = bank(s, balanced_cap)
        balanced_values = {v[2] for v in balanced if U(v) == D(v) and v[2] <= N}
        need(balanced_values == tails | extra, ('balanced scalar invariant', prefix))
        for label, isolated in [('optional_deletion', extra-tails-C), ('short_completions', C-tails-extra)]:
            isolated = {n for n in isolated if n <= N}
            if isolated and label not in component_witnesses:
                component_witnesses[label] = {'prefix_after_t':prefix, 'raw_source':raw(s), 'targets':sorted(isolated)}
        # A strictly larger forward bank is a noncircular stress test of B.
        if index < 12:
            larger = bank(s, 2*bound(s,N))
            larger_values = {target(q) for q in larger if section(q) and target(q) <= N}
            need(larger_values == expected, ('larger-bound stress', prefix))
            larger_bound_checks += 1
    return {'seed_height_cap':seed_cap, 'seed_count':len(seeds), 'target_window':[1,N],
            'source_target_decisions':len(seeds)*N, 'positive_f_sources':positive_f,
            'sum_forward_bank_sizes':total_forward, 'sum_core_bank_sizes':total_core,
            'largest_forward_height_cap':max_forward_cap, 'largest_balanced_height_cap':max_balanced_cap,
            'doubled_bound_sources':larger_bound_checks, 'necessary_component_controls':component_witnesses}

def hostile_checks():
    examples=[]
    for r in range(41):
        lam=3**r
        s=word(ROOT, 'Z'+'OJ'*r+'OZ')
        need(raw(s)==tuple(lam*x for x in (9,1,11,1)), ('hostile profile',r))
        need(s[1]==2*lam and target(s)==3*lam and target(s)<D(s), ('hostile obstruction',r))
        balanced=word(ROOT,'Z'+'OJ'*r)
        q=word(balanced,'OJJJJ')
        need(section(q) and target(q)==3*lam, ('already represented',r))
        if r<3:
            cap=bound(s,3*lam)
            b=bank(s,cap)
            need(not any(section(v) and target(v)==3*lam for v in b), ('finite hostile mismatch',r))
            examples.append({'r':r,'section_target':3*lam,'source_D':D(s),'height_cap':cap,'descendants':len(b)})
    return {'formula_range_r':[0,40], 'complete_bounded_obstruction_checks':examples}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed-cap',type=int,default=128)
    parser.add_argument('--target-cap',type=int,default=63)
    parser.add_argument('--geometry-cap',type=int,default=2048)
    args=parser.parse_args()
    result={'status':'PASS','representation':'independently derived (a,f,e,p)',
            'formal':formal_checks(), 'geometry':geometry_checks(args.geometry_cap),
            'classification':classification_checks(args.seed_cap,args.target_cap),
            'hostile':hostile_checks(),
            'scope':'Append-only classification for a supplied actual source. No all-target coverage claim.'}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
