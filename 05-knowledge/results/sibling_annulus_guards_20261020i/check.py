#!/usr/bin/env python3
"""Exact sibling interaction guards, not a universal annulus certificate."""
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('native', HERE.parent / 'session18_native_annulus/verify_annulus.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
checks = 0

def ck(ok, label):
    global checks
    checks += 1
    if not ok:
        raise ValueError(label)

def add(*vs):
    return tuple(sum(v[i] for v in vs) for i in range(2))

def band_charge(b, r):
    L, U = m.L, m.U
    return [add(L[r], L[(r+1)%3], L[(r+2)%3]),
            add(L[r], L[(r+1)%3], U[(r+2)%3]),
            add(L[r], L[(r+1)%3]),
            add(L[r], U[(r+1)%3]), L[r], U[r]][b]

def band(N, s):
    for b, hi in enumerate((8,16,32,64,128,256)):
        if s < hi*N:
            return b
    raise ValueError('outside annulus')

def source_row(N, s, routes):
    u,t,d = routes.get(s)
    a,b = m.charge(64*N,s)
    r15,r35 = s%15,s%35
    ck(r15%5 == r35%5 and (r15+15*((r35-r15)%7))%105 == s%105, 'joint register')
    return [s,r15,r35,u,t,d,a,b]

def run():
    table = []
    for b in range(5):
        row=[]
        for r in range(3):
            a,z = band_charge(b,r),band_charge(b+1,(2*r+1)%3)
            row.append([r,list(a),list(z),-m.beta(a,z)])
            ck(a != z, 'no same-charge near sibling')
        ck(sorted(v[3] for v in row) == ([0,0,0] if b==0 else [-3,0,3]), 'three-residue alphabet')
        table.append(row)

    # Independent direct B-chain and parent-group check, entire small universes.
    for N in range(1,65):
        parents=defaultdict(list)
        for s in range(4*N+1,256*N,2):
            if not m.source(s): continue
            ck(m.charge(64*N,s) == band_charge(band(N,s),s%3), 'six-band formula')
            u,t=m.edge(s);parents[u].append((s,t))
        for rows in parents.values():
            ck(len(rows)<=2,'at most two same-parent new sources')
            if len(rows)==2:
                (x,a),(y,b)=sorted(rows)
                near=(x%8==1 and y==2*x+1 and b-a==1)
                far=(x%8 in (3,7) and y==32*x+17 and b-a==5)
                ck(near or far, 'complete inverse classification')
                ck(m.charge(64*N,x)!=m.charge(64*N,y),'same-parent exact-zero obstruction')

    obj=json.loads((HERE.parent/'session18_native_annulus/frontier.json').read_text())
    V=set(obj['frontier']);routes=m.Routes(obj['frontier'],dict(obj['frontier_ranks']))
    N=4096; bounds=(4,8,16,32,64,128)
    counts=Counter();debit=[0]*257;retired=[0]*257;digest=hashlib.sha256();examples={}
    # Disjoint, consecutive triples within each exact band. Sources in V are
    # rejected before grouping: first arrival is otherwise at depth zero.
    for b,(lo,hi) in enumerate(zip(bounds,bounds[1:])):
        xs=[x for x in range(lo*N+1,hi*N,8) if x not in V and 2*x+1 not in V]
        # A deletion breaks a run; never bridge it and silently assume residues.
        runs=[];current=[]
        for x in xs:
            if current and x!=current[-1]+8:
                runs.append(current);current=[]
            current.append(x)
        if current:runs.append(current)
        for run in runs:
            counts['leftover_pairs']+=len(run)%3
            for j in range(0,len(run)-2,3):
                six=[];terms=[]
                for x in run[j:j+3]:
                    y=2*x+1;p=source_row(N,x,routes);q=source_row(N,y,routes)
                    ck(p[3]==q[3] and p[5]==q[5] and q[4]==p[4]+1,'actual sibling stop and clocks')
                    e=-m.beta(p[6:8],q[6:8]);terms.append((e,p[5]));six.extend((p,q))
                ck(len({v[0] for v in six})==6,'distinct source owners')
                counts['blocks']+=1
                neg=[d for e,d in terms if e==-3];pos=[d for e,d in terms if e==3]
                ck((not neg and not pos) or (len(neg)==len(pos)==1),'three-term pattern')
                paid=not pos or neg[0]<=pos[0]
                category='zero' if not pos else ('paid' if paid else 'unpaid')
                counts[category]+=1
                record={'band':b,'terms':terms,'sources':six,'paid':paid}
                digest.update(m.encode(record));examples.setdefault(category,record)
                for d in range(257):
                    value=sum(e for e,t in terms if t<=d)
                    debit[d]+=value
                    if paid:
                        ck(value<=0,'all-horizon selected-interaction payment')
                        retired[d]+=value
                if pos:
                    ck(paid == all(sum(e for e,t in terms if t<=d)<=0 for d in range(257)), 'necessary and sufficient finite guard')
    ck(counts['paid']>0 and counts['unpaid']>0,'positive and hostile native controls')
    mixed=[]
    for triple in ((1009,64625,4136049),(1015,64995,4159715)):
        rows=[source_row(16384,s,routes) for s in triple]
        ck(triple[0]<65536 and triple[1]<65536<=triple[2]<4194304,'two old and one new')
        ck(all(s not in V for s in triple),'mixed packet stopped-source guard')
        ck(len({(p[3],p[5]) for p in rows})==1,'synchronous actual arrival')
        vectors=[tuple((1-2*(p[4]%2))*a for a in p[6:8]) for p in rows]
        ck(add(*vectors)==(0,0),'literal mixed zero packet')
        v,w,z=vectors
        bill=m.beta(v,z)+m.beta(w,z)
        ck(bill==-2*m.Q(z) and m.Q(z)>0,'mixed selected old-new credit')
        ck(m.Q(add(v,w,z))-m.Q(add(v,w))==-m.Q(z),'old baseline retained')
        mixed.append({'scale':16384,'sources':rows,'selected_old_new_interaction':bill,
                      'packet_increment':-m.Q(z),'new_source_norm':m.Q(z)})
    receipt={'scope':'selected owned interactions only; all other pair terms retained',
             'band_table':table,'small_universe_N':[1,64], 'scale':N,
             'counts':dict(counts),'blocks_sha256':digest.hexdigest(),'examples':examples,
             'all_selected_debit_by_depth':debit,'retired_debit_by_depth':retired,
             'most_negative_retired':min(retired),'retired_at_final':retired[-1],
             'mixed_phase_controls':mixed,'routes':routes.audit(),'checks':checks}
    return receipt

if __name__=='__main__':
    result=run(); path=HERE/'receipt.json'
    data=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
    if '--write' in sys.argv:path.write_bytes(data)
    else:ck(path.read_bytes().replace(b'\r\n',b'\n')==data,'saved receipt agreement')
    print(json.dumps({'checks':checks,'sha256':hashlib.sha256(data).hexdigest(),
                      'counts':result['counts'],'most_negative_retired':result['most_negative_retired']}))
