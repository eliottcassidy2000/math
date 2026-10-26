#!/usr/bin/env python3
"""Independent, finite audit controls; standard library, explicit guards.

Does not decide a famous conjecture or validate an external paper's proof.
Producer code is locally authored session code inspected before this import.
"""
import ast
import hashlib
import itertools
import json
import math
from pathlib import Path
import runpy
import subprocess
from fractions import Fraction

ROOT = Path(__file__).resolve().parent
PRODUCER = ROOT.parent / 'session65-reductions'
P = runpy.run_path(str(PRODUCER / 'check.py'))

def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def own_encode(n, records):
    states = [('state', i) for i in range(n)]
    sends = [('send', i) for i in range(len(records))]
    receives = [('receive', i) for i in range(len(records))]
    verts = states + sends + receives
    arrows = set()
    for group in (states, sends, receives):
        arrows.update(itertools.combinations(group, 2))
    for e, (s, t, label) in enumerate(records):
        for v in range(n):
            arrows.add((states[v], sends[e]) if s == v else (sends[e], states[v]))
            arrows.add((receives[e], states[v]) if t == v else (states[v], receives[e]))
        for f in range(len(records)):
            arrows.add((sends[e], receives[f]) if e == f else (receives[f], sends[e]))
    index = {v: i for i, v in enumerate(verts)}
    rows = [sum(1 << index[b] for a, b in arrows if a == v) for v in verts]
    return {'n': n, 'm': len(records), 'rows': rows,
            'types': ['V'] * n + ['P'] * len(records) + ['Q'] * len(records),
            'labels': [None] * n + [x[2] for x in records] + [None] * len(records)}

def own_decode(c):
    n, m = c['n'], c['m']
    arcs = {(a, b) for a, row in enumerate(c['rows']) for b in range(n+2*m) if row & (1 << b)}
    result = []
    for i in range(m):
        p = n+i
        mates = [q for q in range(n+m,n+2*m) if (p,q) in arcs]
        need(len(mates) == 1, 'independent matching')
        q = mates[0]
        sources = [s for s in range(n) if (s,p) in arcs]
        targets = [t for t in range(n) if (q,t) in arcs]
        need(len(sources) == len(targets) == 1, 'independent incidence')
        result.append((sources[0],targets[0],c['labels'][p]))
    return result

def codecs():
    totals = []
    for n in (1,2):
        # Each read slot absent, or target/write/move fully specified.
        choices = [None] + list(itertools.product(range(n), (0,1), (-1,1)))
        count = 0
        for table in itertools.product(choices, repeat=2*n):
            records = [(slot//2, inst[0], [slot%2, inst[1], inst[2]])
                       for slot,inst in enumerate(table) if inst is not None]
            own = own_encode(n,records)
            their = P['codec'](n,records)
            need(own == their, 'independent encoding mismatch')
            need(P['decode'](own) == records, 'producer decode auditor encode')
            need(own_decode(their) == records, 'auditor decode producer encode')
            count += 1
        totals.append({'states':n,'full_read_write_move_target_tables':count})
    records = [(0,0,'a'),(0,0,'b'),(0,1,'c'),(1,0,'d')]
    c = own_encode(2,records)
    rejected = 0
    for a,b in itertools.combinations(range(len(c['rows'])),2):
        bad = {**c,'rows':c['rows'][:]}
        bad['rows'][a] ^= 1<<b
        bad['rows'][b] ^= 1<<a
        try:
            P['decode'](bad)
        except RuntimeError:
            rejected += 1
        else:
            raise RuntimeError('canonical single-arc mutation accepted')
    # A changed payload is a DIFFERENT valid graph, not an invalid code.
    changed = {**c,'labels':c['labels'][:]}
    changed['labels'][2] = 'changed-write-move-label'
    decoded = P['decode'](changed)
    need(decoded != records and decoded[0][2] == changed['labels'][2], 'payload not retained')
    return {'full_tables':totals,'canonical_arc_corruptions_rejected':rejected,
            'valid_payload_change_decodes_as_changed_source':True}

def adjacency(n,mask):
    result = set()
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        result.add((i,j) if mask & (1<<bit) else (j,i))
    return result

def covers(n,arcs):
    # Enumerate vertex permutations with cuts; divided by c! for unordered covers.
    ordered = [0]*(n+1)
    for order in itertools.permutations(range(n)):
        for cuts in range(1 << (n-1)):
            if all((cuts >> i)&1 or (order[i],order[i+1]) in arcs for i in range(n-1)):
                ordered[1+cuts.bit_count()] += 1
    return ordered[1:]

def hamiltonian(n,arcs):
    return sum(all((p[i],p[i+1]) in arcs for i in range(n-1)) for p in itertools.permutations(range(n)))

def continuation():
    result = []
    for mask in (69,76):
        arcs=adjacency(5,mask)
        profile=covers(5,arcs)
        # Five-vertex factor -> singleton 5 -> singleton 6 -> factor.
        context=arcs|{(v,5) for v in range(5)}|{(5,6)}|{(6,v) for v in range(5)}
        value=hamiltonian(7,context)
        need(value==3*profile[0]+profile[1], 'context response kernel')
        result.append({'mask':mask,'ordered_profile':profile,'two_vertex_context_H':value})
    need([r['ordered_profile'] for r in result]==[[15,78,198,240,120],[15,90,210,240,120]], 'profile hostile')
    need([r['two_vertex_context_H'] for r in result]==[123,135], 'context hostile')
    return result

def lrc():
    cases=0
    yes=0
    for k in range(1,4):
        for speeds in itertools.product(range(1,7),repeat=k):
            for eps in (Fraction(1,k+1),Fraction(1,2),Fraction(2,5)):
                grid=eps.denominator*math.lcm(*speeds)
                witnesses=[]
                for j in range(grid+1):
                    if all(min((v*j)%grid,grid-(v*j)%grid)*eps.denominator >= eps.numerator*grid for v in speeds):
                        witnesses.append(Fraction(j,grid))
                intervals=P['lrc_exact'](speeds,eps)
                need(bool(intervals)==bool(witnesses), 'integer-grid threshold disagreement')
                for a,b in intervals:
                    need(a<=b and all(min((v*a)%1,1-(v*a)%1)>=eps for v in speeds), 'endpoint not a witness')
                yes += bool(intervals)
                cases += 1
    return {'rational_threshold_tuple_cases':cases,'yes_cases':yes,'method':'exact common-denominator integer grid, k<=3 and speeds 1..6'}

def pin_check():
    frozen = {
        'REPORT.md': '434cba60cf4760eb9ddb0d9cf25500e54b1db321376eeefc6b4bdd5d946fc898',
        'check.py': 'a5bf1c7b398393416d660008b9f7804f77de2ca39dc9592c3beaa80e68626396',
        'normal.json': '153b7925cf69f895b634e94388e3f8f5a8708cad49f06b2f7c501e65c6cf16e3',
        'optimized.json': '153b7925cf69f895b634e94388e3f8f5a8708cad49f06b2f7c501e65c6cf16e3',
        'SOURCE_PINS.json': '1c203205287a884f5eb38091b9e0503b8bdb7ce21823c678950412087993e192',
    }
    for path, expected in frozen.items():
        need(hashlib.sha256((PRODUCER/path).read_bytes()).hexdigest() == expected,
             'audited producer freeze changed: '+path)
    manifest = json.loads((PRODUCER/'MANIFEST.json').read_text())
    for row in manifest['files']:
        b = (PRODUCER/row['path']).read_bytes()
        need(len(b) == row['bytes'] and hashlib.sha256(b).hexdigest() == row['sha256'],
             'producer manifest mismatch: '+row['path'])
    pins=json.loads((PRODUCER/'SOURCE_PINS.json').read_text())
    for row in pins['files']:
        content=(PRODUCER/row['snapshot']).read_bytes()
        need(hashlib.sha256(content).hexdigest()==row['sha256'], 'snapshot pin')
        actual=subprocess.check_output(['git','-C',str(ROOT.parent/'session65-repo'),'show',pins['commit']+':'+row['repo_path']])
        need(content==actual,'git provenance')
    # No Python assert: neither auditor nor producer gates can disappear under -O.
    for path in (Path(__file__),PRODUCER/'check.py'):
        need(not any(isinstance(n,ast.Assert) for n in ast.walk(ast.parse(path.read_text()))),'disabled assert gate')
    return {'git_commit':pins['commit'],'snapshot_and_git_blob_matches':len(pins['files']),
            'frozen_primary_files':len(frozen),'manifest_rows_verified':len(manifest['files']),
            'assert_statements':0}

print(json.dumps({'status':'PASS_INDEPENDENT_FINITE_AUDIT_ONLY','codec':codecs(),
                  'continuation_hostile':continuation(),'lrc_threshold_controls':lrc(),
                  'source_integrity':pin_check()},sort_keys=True,indent=2))
