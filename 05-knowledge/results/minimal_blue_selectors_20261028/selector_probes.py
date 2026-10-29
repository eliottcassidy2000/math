"""Endpoint decoder and measured selector hierarchy for the seven-vertex pair."""
from collections import Counter
from fractions import Fraction
from itertools import permutations,product
from pathlib import Path
import json
import blue_minimum as m
r=m.r;core=m.core;need=m.need
HERE=Path(__file__).resolve().parent

def sandwich(t):
    t=core.complement(t)
    for left in (True,True,False,False):t=r.adjoin(t,left,False)
    return core.complement(t)

def decode(t):
    N=len(t);left=[i for i,x in enumerate(t) if x.bit_count()==1]
    right=[i for i,x in enumerate(t) if N-1-x.bit_count()==1]
    middle=[i for i in range(N) if i not in left and i not in right]
    need(len(left)==len(right)==2 and not set(left)&set(right),'unique outer degree classes')
    starts={v for u in left for v in middle if core.edge(t,u,v)}
    ends={u for u in middle for v in right if core.edge(t,u,v)}
    need(len(starts)==len(ends)==1,'unique two rooted ports')
    rows=tuple(sum(core.edge(t,u,v)<<j for j,v in enumerate(middle)) for u in middle)
    return rows,middle.index(next(iter(starts))),middle.index(next(iter(ends)))

def decoder_audit():
    out=[]
    for n in range(1,8):
        N=1<<((n-1)*(n-2)//2)
        for x in range(N):
            t=core.presentation(n,x);s=sandwich(t)
            need(decode(s)==(t,0,n-1),'literal decoder')
            # A changing cyclic relabeling checks that recovery is structural.
            k=x%(n+4);p=tuple(range(k,n+4))+tuple(range(k));u=core.relabel(s,p)
            v,a,b=decode(u)
            labels=[p[j]-2 for j in range(n+4) if 2<=p[j]<n+2]
            need(v==core.relabel(t,labels) and labels[a]==0 and labels[b]==n-1,'equivariant decoder')
        out.append(dict(order=n,presentations=N,context_order=n+4))
    return out

def legal(t,p):return all(core.edge(t,p[i],p[i+1]) for i in range(len(t)-1))

def selector_audit():
    a,b=(core.presentation(7,x) for x in (1072,4176));pair=(a,b)
    maps=[[p for p in r.oracle.isomorphisms(t,u,True) if p[0]==0 and p[-1]==6]
          for t,u in (pair,tuple(map(core.complement,pair)))]
    marks=[]
    for i in range(7):marks.append(dict(position=i,separates_before=not any(p[i]==i for p in maps[0]),
                                      separates_after_tau=not any(p[i]==i for p in maps[1])))
    need(all(v['separates_before'] or v['separates_after_tau'] for v in marks[1:6]),'every interior mark separates')
    need(not any(v['separates_before'] or v['separates_after_tau'] for v in (marks[0],marks[-1])),'endpoint marks insufficient')
    guards=[]
    for p in permutations(range(7)):
        aa,bb=legal(a,p),legal(b,p)
        if aa!=bb:
            support=sum(i!=p[i] for i in range(7))
            moved=[i for i in range(7) if i!=p[i]]
            guards.append((support,max(moved)-min(moved)+1,p,aa,bb))
    guards.sort()
    fixed=[g for g in guards if g[2][0]==0 and g[2][-1]==6]
    q,blues=r.old.native_classes(7)
    reset=[]
    for before_tau in (False,True):
        banks=[]
        for t in pair:
            if before_tau:t=core.complement(t)
            paths=[p for p in core.paths(t) if p[0]==1]
            row=Counter()
            for p in paths:
                word=core.free_mask(core.relabel(t,p));row[q[word^32767],blues[word]]+=1
            banks.append((len(paths),row))
        normalized=[{k:Fraction(v,N) for k,v in row.items()} for N,row in banks]
        keys=set(normalized[0])|set(normalized[1]);different=sorted(k for k in keys if normalized[0].get(k,0)!=normalized[1].get(k,0))
        reset.append(dict(before_tau=before_tau,legal_path_counts=[x[0] for x in banks],
                          separating_partner_cells=[dict(class_id=k[0],blue=k[1],probabilities=[str(d.get(k,0)) for d in normalized]) for k in different],
                          full_rows=[[[k[0],k[1],v] for k,v in sorted(row.items())] for N,row in banks]))
    # A class-independent observation for the successor-filtered reset.
    H={c:r.oracle.hamilton(core.presentation(7,next(x for x,v in enumerate(q) if v==c))) for c in {v for row in reset for bank in row['full_rows'] for v,_,_ in bank}}
    for row in reset:
        hist=[]
        for bank in row['full_rows']:
            counts=Counter()
            for c,blue,mass in bank:counts[H[c]]+=mass
            hist.append(dict(sorted(counts.items())))
        row['partner_hamilton_histograms']=hist
    return dict(added_mark_tests=marks,all_guard_count=len(guards),minimum_guard=guards[0],
                minimum_endpoint_preserving_guard=fixed[0] if fixed else None,
                successor_start_reset=reset)


def carrier_census():
    old=json.loads((HERE.parent/'singleton_pair_reconstruction_20261028/rooted-seven.json').read_text(encoding='utf-8'))
    groups=old['nontrivial'];refined=old['cells'];remaining=[]
    for group in groups:
        need(len({m.state(core.presentation(7,x)) for x in group})==1,'B constant on old endpoint-pair/blue cell')
        cells=[]
        for x in group:
            t=core.presentation(7,x)
            for cell in cells:
                u=core.presentation(7,cell[0])
                if all(any(all(p[i]==i for i in (0,1,6)) for p in r.oracle.isomorphisms(v,w,True))
                       for v,w in ((t,u),(core.complement(t),core.complement(u)))):
                    cell.append(x);break
            else:cells.append([x])
        refined+=len(cells)-1;remaining.extend(c for c in cells if len(c)>1)
    return dict(order=7,minimal_cells=old['cells'],literal_states=32768,
                successor_port_cells=refined,remaining=remaining,
                scope='Both endpoint-rooted relations and B, then additionally mark current position 1 in both relations.')


def buffered_controls():
    maps=((0,1,4,2,3,5,6),(0,2,4,1,5,3,6));out=[]
    for k in range(6):
        a,b=(core.presentation(7,x) for x in (1072,4176))
        for side in (True,False):
            for _ in range(k):
                a=core.complement(r.adjoin(core.complement(a),side,False))
                b=core.complement(r.adjoin(core.complement(b),side,False))
        need(a!=b and all(core.is_strong(t) for t in (a,b)),'strong buffered controls')
        todo=[(a,b,k,0)];count=0
        for depth in range(4):
            nxt=[]
            for t,u,L,parity in todo:
                N=len(t);p=list(range(N));p[L:L+7]=[L+i for i in maps[parity]]
                need(r.is_map(t,u,p),'buffered actual isomorphism')
                need(all(p[i]==i for i in list(range(k+1))+list(range(N-k-1,N))),'all selected boundary ports fixed')
                need(m.state(t)==m.state(u)==('dead',),'blue-dead buffer')
                count+=1
                if depth<3:
                    for side,flip in product((False,True),repeat=2):
                        nxt.append((r.adjoin(t,side,flip),r.adjoin(u,side,flip),L+side,parity^flip))
            todo=nxt
        out.append(dict(boundary_width=k+1,order=7+2*k,contexts=count))
    return out

def main():
    result=dict(status='FINITE-EXACT',decoder=decoder_audit(),selectors=selector_audit(),
                carrier=carrier_census(),buffered=buffered_controls())
    (HERE/'selector-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(carrier=result['carrier'],buffered=result['buffered'],
                         minimum_guard=result['selectors']['minimum_endpoint_preserving_guard'],
                         reset_counts=[row['legal_path_counts'] for row in result['selectors']['successor_start_reset']])));return result

if __name__=='__main__':main()
