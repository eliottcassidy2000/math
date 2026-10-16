#!/usr/bin/env python3
"""Finite exact annular certificates. Discovery never sees audit ROOT labels.

Run from the repo root; output is regenerated, not an input certificate.
All integers 1..8X receive the declared prefix observation budget.
"""
import itertools
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

from collatz_gram_slack_20261015b import (
    T, close_graph, digest, gram, independent_labels, need, observe,
    transport_slack,
)

ROOT = Path(__file__).resolve().parents[2]


def w(n):
    return 2 if n % 3 == 2 else -1


def add(row, n, coefficient, cache):
    e, s, _ = cache[n]
    row[e] += coefficient*s


def annular_rows(k, cache):
    """Literal source rows; fourfold cancellation checked, not assumed by code."""
    X, b = 2**k, 2**(k//2)
    h, d = defaultdict(int), defaultdict(int)
    for n in range(1, X+1):
        add(h, n, w(n), cache)
    for n in range(X+1, 4*X+1):
        add(d, n, w(n), cache)
    hs, ds = [], []
    for a in range(X, 2*X, b):
        hr, dr = defaultdict(int), defaultdict(int)
        for n in range(a, a+b):
            for e, v in h.items():
                hr[e] += v
            for e, v in d.items():
                dr[e] += v
            # Check the actual finite controller authenticates t(4m)=t(m).
            e, s, _ = cache[n+1]
            ee, ss, _ = cache[4*n+4]
            need((e, s) == (ee, ss), 'two-step native fourfold identity')
            add(h, n+1, w(n+1), cache)
            add(d, n+1, -w(n+1), cache)
            for j in range(1, 5):
                add(d, 4*n+j, w(4*n+j), cache)
        hs.append({e: v for e, v in hr.items() if v})
        ds.append({e: v for e, v in dr.items() if v})
    if k <= 12:
        # Independent coordinate formula counts N in [a,a+b) with N<n<=4N.
        for a,row in zip(range(X,2*X,b),ds):
            direct=defaultdict(int)
            for n in range(a+1,4*(a+b-1)+1):
                multiplicity=max(0,min(a+b-1,n-1)-max(a,(n+3)//4)+1)
                add(direct,n,multiplicity*w(n),cache)
            need({e:v for e,v in direct.items() if v}==row,'independent annular source multiplicities')
    return hs, ds


def push_rows(rows, controller):
    out = []
    for row in rows:
        new = defaultdict(int)
        for e, v in row.items():
            f, s = controller[e]
            new[f] += s*v
        out.append({e: v for e, v in new.items() if v})
    return out


def l1(G):
    return sum(abs(v) for row in G for v in row)


def measure(k, hs, ds):
    X, b = 2**k, 2**(k//2)
    cols = sorted({1} | {e for row in hs+ds for e in row})
    gs = {}
    for name, rows in (('H', hs), ('annular', ds),
                       ('sum', [{e: h.get(e, 0)+d.get(e, 0) for e in cols} for h, d in zip(hs, ds)])):
        gs[name] = gram(rows, cols)
    # Exact PSD Gram of the pair (H, Delta), including its unsymmetrized cross block.
    joint_rows = [{**{(0,e): v for e,v in h.items()}, **{(1,e): v for e,v in d.items()}}
                  for h,d in zip(hs, ds)]
    jcols = [(side,e) for side in (0,1) for e in cols]
    joint = gram(joint_rows, jcols)
    m = len(cols)
    for i in range(m):
        for j in range(m):
            need(gs['sum'][i][j] == gs['H'][i][j]+gs['annular'][i][j]+joint[i][m+j]+joint[m+i][j],
                 'joint cross blocks recover sum Gram exactly')
    # Annular Lipschitz envelope: 16 times variance of the integer ramp.
    fine = F(4,3)*X*(b*b-1)
    q = F(l1(gs['annular']), b)
    out = {'columns': len(cols), 'gram_hashes': {s:digest(g) for s,g in gs.items()},
           'coarse_H': [l1(gs['H']), b], 'coarse_annular': [l1(gs['annular']), b],
           'coarse_sum': [l1(gs['sum']), b], 'annular_fine': [fine.numerator, fine.denominator],
           'annular_upper_over_X2': float((q+fine)/(X*X)),
           'annular_coarse_over_X2': float(q/(X*X)),
           'log_Gram_target_passes': q <= F(k*k,16)*X*X,
           'marginal_Minkowski_budget_over_X2': float((F(l1(gs['H']),b)**F(1,2)+F(l1(gs['annular']),b)**F(1,2))**2/(X*X)),
           'joint_sum_coarse_over_X2': float(F(l1(gs['sum']), b*X*X)),
           'joint_Gram_hash':digest(joint)}
    return out, cols, gs


def run_case(k, depth):
    X, b, N = 2**k, 2**(k//2), 8*2**k+1
    edges = observe(N, depth)
    cache, cycles = close_graph(edges, range(1,N))
    need(not cycles, 'no nonROOT cycle in this declared census')
    hs, ds = annular_rows(k, cache)
    info, cols, gs = measure(k, hs, ds)
    info.update({'additional_steps':0, 'new_edges':0,
                 'unresolved_sources':sum(cache[n][0] != 1 for n in range(1,N)),
                 'compiled_max_native_clock':max(v[2] for v in cache.values())})
    records = [info]
    certificates = [(hs, ds, cache, info)]
    boundary = sorted({cache[n][0] for n in range(1,N)})
    schedule=(1,2,4,8,16) if (k,depth) in ((12,16),(14,16)) else (1,2,4,8)
    for q in schedule:
        enlarged = dict(edges)
        for e in boundary:
            x=e
            for _ in range(q):
                if x==1:
                    break
                enlarged[x]=T(x)
                x=T(x)
        bc, _ = close_graph(enlarged, boundary)
        controller = {e:bc[e][:2] for e in cols}
        ch, cd = push_rows(hs,controller), push_rows(ds,controller)
        ext, newcols, ng = measure(k,ch,cd)
        for name in gs:
            targetcols, transported, loss = transport_slack(gs[name],cols,controller)
            # Transport retains some zero columns; compare in a common index.
            index={e:i for i,e in enumerate(targetcols)}
            need(all(ng[name][i][j] == transported[index[e]][index[f]] for i,e in enumerate(newcols) for j,f in enumerate(newcols)),
                 'fresh Gram equals transported pair buckets')
            need(l1(transported)==l1(ng[name]), 'zero-column deletion exact')
            ext[name+'_paid_loss']=loss
        cc={n:(bc[e][0],s*bc[e][1],d+bc[e][2]) for n,(e,s,d) in cache.items()}
        ext.update({'additional_steps':q,'new_edges':len(enlarged)-len(edges),
                    'unresolved_sources':sum(cc[n][0]!=1 for n in range(1,N)),
                    'compiled_max_native_clock':max(v[2] for v in cc.values())})
        records.append(ext)
        certificates.append((ch,cd,cc,ext))
    frozen = digest(records)
    # Independent audit begins only after the complete controller schedule is frozen.
    labels=independent_labels(N)
    H=[0]*N
    for n in range(1,N):
        H[n]=H[n-1]+w(n)*labels[n]
    delta=[H[4*n]-H[n] for n in range(X,2*X)]
    need(all(abs(delta[j+1]-delta[j])<=4 for j in range(X-1)), 'literal increment bound')
    exact_D=sum(v*v for v in delta)
    exact_H=sum(H[n]*H[n] for n in range(X,2*X))
    exact_cross=sum(H[n]*delta[n-X] for n in range(X,2*X))
    for ch,cd,cc,rec in certificates:
        z={1:1}
        for n in range(1,N):
            e,s,d=cc[n]
            need(s==(-1 if d&1 else 1), 'native integer clock parity retained')
            value=s*labels[n]
            need(e not in z or z[e]==value, 'independent terminal label consistency')
            z[e]=value
        for j,a in enumerate(range(X,2*X,b)):
            need(sum(v*z[e] for e,v in ch[j].items())==sum(H[a:a+b]), 'actual H bin sum')
            need(sum(v*z[e] for e,v in cd[j].items())==sum(delta[a-X:a-X+b]), 'actual annular bin sum')
        upper=F(*rec['coarse_annular'])+F(*rec['annular_fine'])
        need(exact_D<=upper, 'certified actual consequence D bound')
    coarse_actual=sum(F(sum(delta[i:i+b])**2,b) for i in range(0,X,b))
    return {'k':k,'depth':depth,'source_count':N-1,'initial_observed_edges':len(edges),
            'frozen_certificate_hash':frozen,'controllers':records,
            'audit_only':{'D':exact_D,'E':exact_H,'cross':exact_cross,
                          'annular_coarse':[coarse_actual.numerator,coarse_actual.denominator],
                          'D_over_X2':float(F(exact_D,X*X))}}


def elementary_controls():
    count=0
    for L in range(1,9):
        for increments in itertools.product((-4,0,4),repeat=L-1):
            xs=[0]
            for inc in increments: xs.append(xs[-1]+inc)
            v=F(sum(x*x for x in xs))-F(sum(xs)**2,L)
            bound=F(4,3)*L*(L*L-1)
            need(v<=bound,'all small ramp hostile controls')
            for cut in range(1,L):
                a,b=xs[:cut],xs[cut:]
                refined=F(sum(a)**2,len(a))+F(sum(b)**2,len(b))+F(4,3)*(len(a)*(len(a)**2-1)+len(b)*(len(b)**2-1))
                need(refined<=F(sum(xs)**2,L)+bound,'nested upper envelope')
            count+=1
    old=json.loads((ROOT/'05-knowledge/results/collatz_discrepancy_20261013c.json').read_text())
    rows={r['k']:r for r in old['rows']}
    seeds=[]
    for k in (32,33):
        E=int(rows[k+1]['residues'][2][3])+rows[k]['residues'][2][0]**2-rows[k+1]['residues'][2][0]**2
        need(E<=F(k*k,16)*2**(2*k),'inherited exact energy seed')
        need(E<=(16*k-480)*2**(2*k),'linear-log propagation seed')
        seeds.append({'k':k,'E':E,'normalized':float(F(E,2**(2*k)))})
    return {'ramp_sequences':count,'inherited_energy_seeds':seeds}


def adaptive_probe(mode):
    """A declared heuristic, with lexicographic control; no optimality claim."""
    k,depth,q,width,cap=12,16,8,10,30
    X,b,N=2**k,2**(k//2),8*2**k+1
    base=observe(N,depth)
    edges=dict(base)
    cache,_=close_graph(edges,range(1,N))
    _,rows=annular_rows(k,cache)
    cols=sorted({1}|{e for row in rows for e in row})
    G=gram(rows,cols)
    initial_norm=l1(G)
    records=[]
    total_paid=0
    for step in range(cap):
        if F(l1(G),b)<=F(k*k,16)*X*X:
            break
        loads={e:sum(abs(v) for v in G[i]) for i,e in enumerate(cols) if e!=1}
        selected=sorted(loads,key=(lambda e:(-loads[e],e)) if mode=='weighted' else (lambda e:e))[:width]
        oldedges=len(edges)
        for e in selected:
            x=e
            for _ in range(q):
                if x==1: break
                edges[x]=T(x)
                x=T(x)
        bc,_=close_graph(edges,cols)
        controller={e:bc[e][:2] for e in cols}
        nc,ng,slack=transport_slack(G,cols,controller)
        rows=push_rows(rows,controller)
        need(gram(rows,nc)==ng,'adaptive direct Gram versus pair transport')
        cols,G=nc,ng
        total_paid+=slack['deleted_mass']+slack['cancelled_mass']
        records.append({'round':step+1,'selected':selected,'new_edges':len(edges)-oldedges,
                        'cumulative_new_edges':len(edges)-len(base),'columns':len(cols),
                        'coarse_over_X2':float(F(l1(G),b*X*X)),
                        'paid_loss':slack,'target_passes':F(l1(G),b)<=F(k*k,16)*X*X})
    need(initial_norm-l1(G)==total_paid,'adaptive cumulative opposite-sign invoice')
    frozen=digest(records)
    # Independent source labels only now; audit every retained source and final row.
    labels=independent_labels(N)
    final,_=close_graph(edges,range(1,N))
    z={1:1}
    for n in range(1,N):
        e,s,d=final[n]
        need(s==(-1 if d&1 else 1),'adaptive exact native clock')
        value=s*labels[n]
        need(e not in z or z[e]==value,'adaptive terminal consistency')
        z[e]=value
    H=[0]*N
    for n in range(1,N): H[n]=H[n-1]+w(n)*labels[n]
    for a,row in zip(range(X,2*X,b),rows):
        need(sum(v*z[e] for e,v in row.items())==sum(H[4*n]-H[n] for n in range(a,a+b)),
             'adaptive literal source mean')
    return {'mode':mode,'k':k,'depth':depth,'block_depth':q,'batch_width':width,'round_cap':cap,
            'frozen_hash':frozen,'rounds':records,
            'unresolved_sources':sum(final[n][0]!=1 for n in range(1,N)),
            'max_native_clock':max(v[2] for v in final.values())}


def main():
    controls=elementary_controls()
    cases=[]
    for k,d in ((10,8),(10,16),(12,16),(14,16),(14,24),(16,32)):
        case=run_case(k,d)
        cases.append(case)
        print(k,d,[(r['additional_steps'],round(r['annular_coarse_over_X2'],6),r['log_Gram_target_passes']) for r in case['controllers']],flush=True)
    adaptive=[adaptive_probe(mode) for mode in ('weighted','lexicographic')]
    for a in adaptive:
        print(a['mode'],len(a['rounds']),a['rounds'][-1]['cumulative_new_edges'],a['rounds'][-1]['coarse_over_X2'],flush=True)
    out={'status':'FINITE-EXACT; universal annular Gram budget OPEN','controls':controls,'cases':cases,'adaptive':adaptive}
    (ROOT/'05-knowledge/results/collatz_annular_overlap_20261015c.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__': main()
