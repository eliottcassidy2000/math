#!/usr/bin/env python3
"""Semantic negative controls independent of byte-integrity rejection."""
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import independent as a

HERE = Path(__file__).resolve().parent

def reject(name, action, out):
    try:
        action()
    except (ValueError,KeyError,IndexError) as exc:
        out.append(dict(name=name,result='REJECTED',reason=str(exc)))
        return
    raise RuntimeError('hostile accepted: ' + name)

def relabel(rows, p):
    n = len(rows)
    result = [0] * n
    for i in range(n):
        for j in range(n):
            if (rows[i] >> j) & 1:
                result[p[i]] |= 1 << p[j]
    return result

def main():
    source = json.loads((HERE/'nonrigid'/'normal.json').read_text())
    out = []
    # These bypass integrity intentionally to prove semantic rejection.
    bad = copy.deepcopy(source)
    bad['q9_classes'].pop()
    reject('missing nonrigid class',lambda:a.audit_nonrigid(bad),out)
    bad = copy.deepcopy(source)
    record = copy.deepcopy(bad['q9_classes'][0])
    record['rows'] = relabel(record['rows'],[10]+list(range(10)))
    bad['q9_classes'].append(record)
    bad['nonrigid_classes'] += 1
    reject('nonliteral isomorphic duplicate',lambda:a.audit_nonrigid(bad),out)
    bad = copy.deepcopy(source)
    bad['q1_classes'][-1]['full_automorphism_order'] = 11
    reject('stabilizer order mistaken for full group',lambda:a.audit_nonrigid(bad),out)
    bad = copy.deepcopy(source)
    bad['q7_classes'][0]['vertex_orbit_sizes'] = [1,1,1,1,1,1,1,1,3]
    reject('wrong full vertex orbits',lambda:a.audit_nonrigid(bad),out)
    bad = copy.deepcopy(source)
    bad['q7_classes'][0]['error'] += 11
    reject('wrong nonrigid error',lambda:a.audit_nonrigid(bad),out)
    bad = copy.deepcopy(source)
    rows = bad['q9_classes'][0]['rows']
    rows[0] ^= 1 << 1
    reject('one-sided invalid edge',lambda:a.audit_nonrigid(bad),out)
    bad = copy.deepcopy(source)
    rows = bad['q9_classes'][0]['rows']
    rows[0] ^= 1 << 1
    rows[1] ^= 1
    reject('valid tournament made nonregular',lambda:a.audit_nonrigid(bad),out)
    good = a.audit_nonrigid(source)
    total = a.literal_count((5,)*11)
    reject('nonintegral rigid mass',lambda:a.derive_closure(total+1,good),out)
    reject('negative rigid mass',lambda:a.derive_closure(good['labelled_mass']-1,good),out)
    bad = copy.deepcopy(good)
    bad['pointed_classes'] += 1
    reject('pointed mass inconsistent with deficit',lambda:a.derive_closure(total,bad),out)
    # Positive relabelling control and orientation-preserving negative control.
    reps = a.nonrigid_records(source)
    first = reps[0]['rows']
    renamed = relabel(first,list(reversed(range(11))))
    a.require(len(a.all_isomorphisms(first,renamed)) == reps[0]['full_automorphism_order'], 'relabelled positive control')
    nonself = 0
    for item in reps:
        rows = item['rows']
        dual = [((1 << 11)-1) ^ (1 << i) ^ r for i,r in enumerate(rows)]
        if not a.all_isomorphisms(rows,dual,True):
            nonself += 1
    a.require(nonself > 0, 'ordinary isomorphism never folds all converses')
    # Executable logic mutations must be rejected in either interpreter mode.
    code = (HERE/'independent.py').read_text()
    mutations = [
        ('omit binomial labelled-choice multiplicity','weight * math.comb(m, x)','weight'),
        ('reverse residual-score update','r if j in out else r - 1','r - 1 if j in out else r'),
        ('quotient by repeated-score factorial','total += literal_count(residual)','total += literal_count(residual) // math.prod(math.factorial(z) for z in collections.Counter(residual).values())'),
        ('multiply by repeated-score factorial','total += literal_count(residual)','total += literal_count(residual) * math.prod(math.factorial(z) for z in collections.Counter(residual).values())'),
    ]
    mode = ['-O'] if not __debug__ else []
    with tempfile.TemporaryDirectory(prefix='s51-score-hostile-') as temp:
        packet = Path(temp)/'packet'
        packet.mkdir()
        shutil.copytree(HERE/'nonrigid',packet/'nonrigid',ignore=shutil.ignore_patterns('__pycache__'))
        unrelated = Path(temp)/'unrelated'
        unrelated.mkdir()
        for name,old,new in mutations:
            a.require(code.count(old) == 1,'unique mutation target')
            (packet/'independent.py').write_text(code.replace(old,new,1))
            run = subprocess.run([sys.executable]+mode+[str(packet/'independent.py')],cwd=unrelated,capture_output=True,text=True)
            a.require(run.returncode != 0,'executable hostile rejected: '+name)
            a.require('ValueError:' in run.stderr,'semantic rather than parse failure: '+name)
            reason = run.stderr.strip().splitlines()[-1].split('ValueError: ',1)[-1]
            out.append(dict(name=name,result='REJECTED',reason=reason))
    print(json.dumps(dict(status='PASS',semantic_rejections=len(out),controls=out,
                         positive_relabel_control=True,nonself_converse_representatives=nonself),sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
