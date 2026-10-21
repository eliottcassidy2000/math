#!/usr/bin/env python3
"""Independent exact audit; no imports from the producer and no assert statements."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import argparse
import ast
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PINS = {
    'README.md': 'cf49309a88e3113843a0b7cfc03ea3b15ec8fa789d98c85e314e8df3e9e7fdc3',
    'verify.py': '94e0c6b1a3bf5721c2cfffb24e74d5bad5895335c5eebdc26368c9337f5c35e0',
    'certificate.json': 'd02d1800e54f2dd2b150ca1ee4d658d1df5fe625d1f89a6b0f7ea4bb3bf520c7',
}
COUNTS = Counter()

def need(ok, label):
    if not ok:
        raise RuntimeError('AUDIT FAILURE: ' + label)

def equal(x, y, label):
    need(x == y, label + ': actual=' + repr(x) + '; expected=' + repr(y))

def bad(v, t, n):
    a, b = t.numerator, t.denominator
    z = (v * a) % b
    return n * min(z, b - z) < b

def margin(V, t):
    a, b = t.numerator, t.denominator
    return min(F(min((v*a)%b, (-v*a)%b), b) for v in V)

def unsafe(V, t):
    return [v for v in V if bad(v, t, len(V)+1)]

def isprime(q):
    return q >= 2 and all(q % p for p in range(2, isqrt(q)+1))

def inverse(v, q):
    return next(a for a in range(1, q) if (a*v)%q == 1)

def physical_rows(V):
    """Reconstruct from actual tooth intervals, not the producer's inequalities."""
    k = len(V); n = k+1; q = 2*k+1; delta = F(1, n)
    need(k >= 2 and len(set(V)) == k and min(V) > 0, 'valid distinct positive speeds')
    need(all(gcd(v,q) == 1 for v in V), 'unit speeds')
    grid = [a for a in range(1,q) if not unsafe(V,F(a,q))]
    classes = {min(v%q, q-v%q) for v in V}
    for a in range(1,q):
        equal(unsafe(V,F(a,q)), [v for v in V if (a*v)%q in (1,q-1)], 'grid occupancy')
    if len(classes) < k:
        need(grid, 'duplicate class must give grid witness')
        return {'q':q, 'grid':grid, 'repairs':[]}
    need(not grid and isprime(q), 'complete unit signed system forces prime')
    rows = []
    for v in V:
        a = inverse(v,q); clock = F(a,q); center = (v*clock).__floor__()
        equal(unsafe(V,clock), [v], 'unique clock owner')
        equal(v*clock-center, F(1,q), 'owner clock phase')
        for sign in (1,-1):
            t = (center + sign*delta)/v
            D = k if sign == 1 else 3*k+2
            equal(t-clock, sign*F(D,n*q*v), 'physical repair duration')
            actual = unsafe(V,t); corridor = []; first = []; times = []
            for w in V:
                if w == v:
                    continue
                j = (w*clock).__floor__()
                tooth = j+1 if sign == 1 else j
                low = (tooth-delta)/w; high = (tooth+delta)/w
                entry = low if sign == 1 else high
                need(sign*(entry-clock) > 0, 'first entry is in travel direction')
                r = (a*w)%q
                c = n*(q-r)-q if sign == 1 else n*r-q
                need(c >= 1 and 2 <= r <= q-2, 'positive first-entry numerator')
                equal(entry, clock + sign*F(c,n*q*w), 'first-entry formula')
                enters = sign*(t-entry) > 0
                equal(enters, D*w > c*v, 'corridor strict lower bound')
                in_first = low < t < high
                equal(in_first, c*v < D*w < (c+2*q)*v, 'first-tooth strict double bound')
                if enters:
                    corridor.append(w); times.append({'target':w,'first_entry':str(entry)})
                if in_first:
                    first.append(w)
                if w in actual and not in_first:
                    COUNTS['later_tooth_endpoint_blockers'] += 1
            need(set(first) <= set(actual) <= set(corridor), 'first/endpoint/corridor inclusions')
            need(v not in actual, 'owner is safe at endpoint')
            equal(unsafe(V, 1-t), actual, 'antipodal endpoint reflection')
            rows.append({'owner':v,'clock':str(clock),'sign':sign,'endpoint':str(t),
                         'actual_blockers':actual,'corridor_blockers':corridor,
                         'first_tooth_blockers':first,'timed_entries':times})
    forward = {r['owner']:r['corridor_blockers'] for r in rows if r['sign'] == 1}
    for v, targets in forward.items():
        for w in targets:
            need(v not in forward[w], 'no reciprocal forward corridor edges')
    COUNTS['physical_rows'] += 1
    COUNTS['physical_endpoints'] += len(rows)
    return {'q':q,'grid':grid,'repairs':rows}

def physical_tape(V):
    """Compute before/point/after occupancy directly in each open event cell."""
    n = len(V)+1
    times = sorted({F(n*j+s,n*v) for v in V for j in range(v) for s in (1,n-1)})
    out = []
    for p,t in enumerate(times):
        prev = times[p-1] if p else times[-1]-1
        nxt = times[p+1] if p+1 < len(times) else times[0]+1
        X = [v for v in V if (v*t)%1 == F(1,n)]
        I = [v for v in V if (v*t)%1 == F(n-1,n)]
        before = len(unsafe(V,(prev+t)/2)); at = len(unsafe(V,t)); after = len(unsafe(V,(t+nxt)/2))
        equal(at, before-len(X), 'simultaneous exit point equation')
        equal(at, after-len(I), 'simultaneous entry point equation')
        equal(after, before-len(X)+len(I), 'signed simultaneous transition')
        out.append({'time':str(t),'exit':X,'entry':I,'before':before,'point':at,'after':after})
    equal(out[0]['before'],len(V),'all initially unsafe')
    equal(out[-1]['after'],len(V),'all unsafe before period closes')
    COUNTS['event_tapes'] += 1; COUNTS['events'] += len(out)
    return out

def audit_certificate(cert):
    inverse_count = 0
    for q in range(5,102,2):
        n = (q+1)//2; k = n-1
        for a in range(2,q-1):
            if gcd(a,q) != 1:
                continue
            b = inverse(a,q); p = (n*a-q)*(n*b-q)
            need(p >= k*k, 'inverse product bound')
            equal(p == k*k, {a,b} == {2,n}, 'inverse product equality')
            inverse_count += 1
    equal(cert['inverse_pair_checks'], inverse_count, 'inverse pair count')
    expected_ap = []
    for k in range(2,14):
        tape = physical_tape(list(range(1,k+1)))
        need(all(r['after'] > 0 for r in tape), 'AP has no good open cell')
        row = next(r for r in tape if F(r['time']) == F(1,k+1))
        equal((row['exit'],row['entry'],row['before'],row['point'],row['after']),([1],[k],1,0,1),'AP tied isolated witness')
        expected_ap.append({'k':k,'tied_witness':row})
    equal(cert['tight_AP_controls'], expected_ap, 'all AP certificate fields')
    equal(cert['composite_unit_positive'],physical_rows([1,2,4,5]),'composite positive control')
    tri = physical_rows([11,19,50]); equal(cert['forward_triangle'],tri,'triangle fields')
    need(all(r['actual_blockers'] for r in tri['repairs'] if r['sign'] == 1),'triangle forward obstruction')
    need(not unsafe([11,19,50],F(27,200)), 'triangle backward witness')
    repaired = next(r for r in tri['repairs'] if r['owner'] == 11 and r['sign'] == -1)
    equal(repaired['endpoint'],'1/4','corridor-only example time')
    need(repaired['corridor_blockers'] and not repaired['actual_blockers'],'corridor is not endpoint graph')
    V = [4,5,9,12,14]; hostile = physical_rows(V)
    equal(cert['two_sided_hostile'],{'speeds':V,'witness':'1/24','certificate':hostile},'hostile complete certificate')
    need(all(r['actual_blockers'] for r in hostile['repairs']),'all ten actual endpoints fail')
    tape = physical_tape(V)
    equal(next(r['time'] for r in tape if r['point'] == 0),'1/24','hostile earliest good event')
    equal(margin(V,F(1,24)),F(1,6),'hostile actual witness margin')
    fwd = {r['owner']:r['actual_blockers'] for r in hostile['repairs'] if r['sign'] == 1}
    equal(fwd,{4:[14],14:[5],5:[12],12:[9],9:[4]},'actual five-cycle')
    edge = next(r for r in hostile['repairs'] if r['owner'] == 4 and r['sign'] == 1)
    arrival = edge['timed_entries'][0]['first_entry']
    equal(arrival,'23/84','physical arrival')
    need(F(arrival) != F(inverse(14,11),11),'arrival is not target own clock')
    expected_families = []
    for q in (11,19,43,59,67,83):
        k = (q-1)//2
        H = sorted({a*a%q for a in range(1,q)})
        need(isprime(q) and q%8 == 3,'QR family prime condition')
        equal(len(H),k,'QR count')
        equal({min(h,q-h) for h in H},set(range(1,k+1)),'QR complete signed classes')
        for M in (10,100,1000):
            V = [q*M+h for h in H]; L = min(V); U = max(V)
            need(F(U,L) < F(3*k+2,2*k+3),'QR strict spread')
            rows = physical_rows(V)['repairs']; by_res = {v%q:v for v in V}
            for row in rows:
                mult = -2 if row['sign'] == 1 else 4
                target = by_res[mult*row['owner']%q]
                need(target != row['owner'],'QR distinct target')
                need(target in row['first_tooth_blockers'],'QR first-tooth target blocks')
            t = F(1,2*U)
            equal(margin(V,t),F(L,2*U),'exact QR witness margin')
            need(margin(V,t) >= F(1,k+1),'QR actual lonely witness')
            expected_families.append({'q':q,'k':k,'M':M,'min':L,'max':U,
                'ordinary_witness':str(t),'margin':str(F(L,2*U)),'all_repair_endpoints_fail':True})
    equal(cert['quadratic_residue_families'],expected_families,'all QR certificate fields')
    equal(cert['scope'],'Exact certificate tests, not an LRC census or general proof.','scope label')
    equal(set(cert),{'inverse_pair_checks','tight_AP_controls','composite_unit_positive','forward_triangle','two_sided_hostile','quadratic_residue_families','scope'},'certificate schema')


def rational_cycle_lift(q, hs):
    """Exact rational construction: produce v_s(T)=q*T*W_s+h_s for all T>=T0."""
    k = (q-1)//2; n = k+1; ell = len(hs)
    cs = [n*(q-(hs[(s+1)%ell]*inverse(h,q))%q)-q for s,h in enumerate(hs)]
    need(all(c > 0 for c in cs),'positive cycle coefficients')
    need(prod(cs) < k**ell,'strict cycle product deficiency')
    lo = [F(c,k) for c in cs]; hi = [F(c+2*q,k) for c in cs]
    ratios = []; remaining = F(1)
    for s in range(ell-1):
        left = max(lo[s],remaining/prod(hi[s+1:])); right = min(hi[s],remaining/prod(lo[s+1:]))
        need(left < right,'nonempty feasible rational ratio interval')
        ratio = (left+right)/2
        ratios.append(ratio); remaining /= ratio
    ratios.append(remaining)
    need(all(L < r < U for L,r,U in zip(lo,ratios,hi)),'strict rational lower and upper ratios')
    equal(prod(ratios),F(1),'rational cycle closes')
    ws = [F(1)]
    for r in ratios[:-1]:
        ws.append(ws[-1]*r)
    den = lcm(*(w.denominator for w in ws))
    W = [int(w*den) for w in ws]
    T0 = 1
    for s,c in enumerate(cs):
        z = (s+1)%ell
        for A,B in ((q*(k*W[z]-c*W[s]),k*hs[z]-c*hs[s]),
                    (q*((c+2*q)*W[s]-k*W[z]),(c+2*q)*hs[s]-k*hs[z])):
            need(A > 0,'positive integer-lift slack slope')
            T0 = max(T0,(-B)//A+1)
    for T in (T0,T0+1,T0+37):
        vs = [q*T*w+h for w,h in zip(W,hs)]
        need(min(vs)>0 and len(set(vs)) == ell,'positive distinct lifts')
        for s,c in enumerate(cs):
            z = (s+1)%ell; v = vs[s]; w = vs[z]; a = inverse(v,q)
            endpoint = F(n*((a*v-1)//q)+1,n*v)
            center = (w*F(a,q)).__floor__()+1
            need(F(n*center-1,n*w)<endpoint<F(n*center+1,n*w),'constructed cycle actual first tooth')
            need(bad(w,endpoint,n),'constructed cycle direct unsafe phase')
            need(c*v < k*w < (c+2*q)*v,'constructed strict double inequality')
            COUNTS['constructed_cycle_edges'] += 1
    return {'q':q,'residues':list(hs),'coefficients':cs,'integer_weights':W,'T_min':T0}


def extra_controls():
    # Extended inverse lemma, including composite clocks and its exact equality case.
    for q in range(5,502,2):
        k = (q-1)//2; n = k+1
        for a in range(2,q-1):
            if gcd(a,q) != 1:
                continue
            b = pow(a,-1,q); ab = a*b
            need(ab >= q+1 and a+b <= F(ab,2)+2,'inverse lemma auxiliary inequalities')
            p = (n*a-q)*(n*b-q)
            need(p >= F(n*ab,2)-q >= n*n-q == k*k,'inverse product full chain')
            equal(p == k*k,{a,b} == {2,n},'extended equality classification')
            COUNTS['extended_inverse_pairs'] += 1
    # All signed systems and small height choices for q=5 and q=7.
    for q in (5,7):
        k = (q-1)//2
        for signs in product((1,-1),repeat=k):
            residues = [(s*h)%q for h,s in zip(range(1,k+1),signs)]
            for heights in product((0,1,3),repeat=k):
                physical_rows([h+q*height for h,height in zip(residues,heights)])
                COUNTS['exhaustive_signed_height_rows'] += 1
    # Every composite odd q through 51, with positive distinct unit speeds.
    for q in range(9,52,2):
        if isprime(q):
            continue
        k = (q-1)//2; units = [h for h in range(1,q) if gcd(h,q) == 1]
        V = [units[i%len(units)] + q*(i//len(units)) for i in range(k)]
        cert = physical_rows(V)
        need(cert['grid'] and not cert['repairs'],'composite all-unit grid discharge')
        COUNTS['composite_unit_rows'] += 1
    # Every directed residue cycle modulo rotation for q=5,7,11, with no antipodal duplicates.
    examples = []
    for q in (5,7,11):
        k = (q-1)//2; n = k+1
        for ell in range(2,k+1):
            for chosen in combinations(range(1,q),ell):
                if len({min(h,q-h) for h in chosen}) != ell:
                    continue
                for tail in permutations(chosen[1:]):
                    hs = (chosen[0],)+tail
                    rs = [(hs[(s+1)%ell]*inverse(h,q))%q for s,h in enumerate(hs)]
                    bs = [q-r for r in rs]; cs = [n*b-q for b in bs]
                    equal(prod(bs)%q,(-1)**ell%q,'residue cycle consistency')
                    COUNTS['residue_cycles'] += 1
                    if prod(cs) < k**ell:
                        ex = rational_cycle_lift(q,hs)
                        COUNTS['deficient_cycles_constructed'] += 1
                        if len(examples)<8:
                            examples.append(ex)
                    elif prod(cs) == k**ell:
                        COUNTS['cycle_equality_rejections'] += 1
                    else:
                        COUNTS['cycle_surplus_rejections'] += 1
    # Strict boundary controls: no corridor arrow at entry equality; no first-tooth
    # blocker at exit equality, even though the earlier corridor arrow exists.
    q=5; k=2; n=3; owner=2; a=inverse(owner,q)
    for target, label in ((1,'entry'),(11,'exit')):
        r=(a*target)%q; c=n*(q-r)-q
        endpoint=F(n*((a*owner-1)//q)+1,n*owner)
        if label=='entry':
            equal(k*target,c*owner,'entry equality fixture')
            need(not (k*target>c*owner),'entry equality is not an edge')
        else:
            equal(k*target,(c+2*q)*owner,'exit equality fixture')
            need(k*target>c*owner,'exit equality retains corridor edge')
        need(not bad(target,endpoint,n),'tooth boundary is safe')
    # Proper graph inclusion: the endpoint may be bad in a later tooth, not the first.
    later = next(r for r in physical_rows([1,18])['repairs']
                 if r['owner'] == 1 and r['sign'] == 1)
    equal(later['endpoint'],'1/3','later-tooth fixture endpoint')
    equal(later['actual_blockers'],[18],'later-tooth actual blocker')
    equal(later['first_tooth_blockers'],[],'later-tooth is not first tooth')
    equal(later['timed_entries'],[{'target':18,'first_entry':'11/54'}],
          'later-tooth first-entry sidecar')
    # q=11 threshold M=3 and larger, plus additional prime rows beyond the certificate.
    for q in range(11,252,8):
        if not isprime(q):
            continue
        k=(q-1)//2; H=sorted({a*a%q for a in range(1,q)})
        for M in (3,7):
            V=[q*M+h for h in H]; R=F(max(V),min(V))
            if R>=F(3*k+2,2*k+3):
                continue
            cert=physical_rows(V)
            need(all(r['actual_blockers'] for r in cert['repairs']),'extended QR endpoint failures')
            equal(margin(V,F(1,2*max(V))),F(min(V),2*max(V)),'extended QR exact margin')
            COUNTS['extended_QR_rows']+=1
    return examples


def rejection_controls(cert, source):
    # These corruptions leave JSON syntactically valid and must fail in normal and -O.
    for key in ('tie','endpoint','chronology','margin'):
        altered=deepcopy(cert)
        if key=='tie':
            altered['tight_AP_controls'][0]['tied_witness']['point']=1
        elif key=='endpoint':
            altered['two_sided_hostile']['certificate']['repairs'][0]['actual_blockers']=[]
        elif key=='chronology':
            altered['two_sided_hostile']['certificate']['repairs'][0]['timed_entries'][0]['first_entry']='4/11'
        else:
            altered['quadratic_residue_families'][0]['margin']='1/2'
        caught=False
        try:
            audit_certificate(altered)
        except RuntimeError:
            caught=True
        need(caught,'negative control rejected: '+key)
    # Ensure the revised producer retains its own explicit failures under optimization.
    command=("import runpy; m=runpy.run_path("+repr(str(source/'verify.py'))+"); "
             "m['critical']([1,1])")
    for flags in ([],['-O']):
        run=subprocess.run([sys.executable]+flags+['-c',command],capture_output=True,text=True)
        need(run.returncode!=0 and 'RuntimeError: verification failed:' in run.stderr,
             'producer rejects invalid duplicate speeds '+repr(flags))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-dir',type=Path,default=ROOT/'source')
    parser.add_argument('--certificate',type=Path)
    args=parser.parse_args(); source=args.source_dir.resolve()
    for name,digest in PINS.items():
        equal(sha256((source/name).read_bytes()).hexdigest(),digest,'frozen source pin '+name)
    for path in (source/'verify.py',Path(__file__)):
        need(not any(isinstance(n,ast.Assert) for n in ast.walk(ast.parse(path.read_text()))),
             'no optimization-disabled asserts: '+path.name)
    expected=(source/'certificate.json').read_bytes()
    for flags in ([],['-O']):
        run=subprocess.run([sys.executable]+flags+[str(source/'verify.py')],capture_output=True)
        equal(run.returncode,0,'producer replay exit '+repr(flags))
        equal(run.stdout,expected,'producer byte replay '+repr(flags))
    cert=json.loads((args.certificate or source/'certificate.json').read_text())
    audit_certificate(cert)
    examples=extra_controls()
    completed_counts=dict(COUNTS)
    rejection_controls(cert,source)
    print(json.dumps({'status':'PASS','source_sha256':PINS,'controls':completed_counts,
                      'rational_cycle_lift_examples':examples,
                      'negative_controls':['tie occupancy','endpoint blocker','arrival clock','witness margin',
                                           'producer duplicate speed normal','producer duplicate speed -O'],
                      'scope':'Exact certificate/obstruction audit; no general LRC proof or LRC(14) closure.'},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
