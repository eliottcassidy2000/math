"""Signed orbit spectra for completion cubes; no third-party dependencies."""
from collections import Counter
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,msg):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(msg)

def action(free,p):
    index={e:k for k,e in enumerate(free)}
    dest=[];flip=0
    for k,(u,v) in enumerate(free):
        a,b=p[u],p[v];j=index[tuple(sorted((a,b)))]
        dest.append(j)
        if a>b:flip|=1<<j
    return tuple(dest),flip

def linear(x,dest):
    out=0
    while x:
        b=x&-x;x-=b;out|=1<<dest[b.bit_length()-1]
    return out

def move(x,a):
    dest,flip=a
    return linear(x,dest)^flip

def signed_cycle_polynomial(a):
    dest,flip=a;seen=set();poly=[1]
    for i in range(len(dest)):
        if i in seen:continue
        j=i;length=0;parity=0
        while j not in seen:
            seen.add(j);length+=1;parity^=(flip>>j)&1;j=dest[j]
        new=poly+[0]*length
        sign=(-1)**parity
        for k,c in enumerate(poly):new[k+length]+=sign*c
        poly=new
    return poly

def spectrum(free,group):
    actions=[action(free,p) for p in group]
    total=[0]*(len(free)+1)
    for a in actions:
        poly=signed_cycle_polynomial(a)
        for k,c in enumerate(poly):total[k]+=c
    need(all(c%len(group)==0 for c in total),'integral invariant dimensions')
    answer=[c//len(group) for c in total]
    need(min(answer)>=0,'nonnegative invariant dimensions')
    return answer,actions

def literal_audit(actions,expected):
    """Independent support stabilizers, affine point orbits and quotient trace."""
    d=len(actions[0][0]);full=(1<<d)-1
    supports=set();by_degree=Counter()
    points=set();representatives=[];fixed_complement=0
    for x in range(1<<d):
        if x not in supports:
            orbit={linear(x,a[0]) for a in actions};supports.update(orbit)
            positive=all((x&a[1]).bit_count()%2==0 for a in actions if linear(x,a[0])==x)
            if positive:by_degree[x.bit_count()]+=1
        if x not in points:
            orbit={move(x,a) for a in actions};points.update(orbit);representatives.append(x)
            fixed_complement+=(x^full in orbit)
    observed=[by_degree[k] for k in range(d+1)]
    need(observed==expected,'literal signed support orbits equal polynomial spectrum')
    need(len(representatives)==sum(expected),'literal affine orbits equal spectrum at one')
    need(fixed_complement==sum((-1)**k*c for k,c in enumerate(expected)),
         'literal complement-fixed orbit trace')
    result={'literal_states':1<<d,'orbits':len(representatives),
            'complement_fixed_orbits':fixed_complement,'signed_support_spectrum':observed}
    if d<=6:
        ids={x:i for i,x in enumerate(representatives)}
        owner={move(x,a):ids[x] for x in representatives for a in actions}
        size=len(ids);matrix=[[0]*size for _ in range(size)]
        for x,i in ids.items():
            for j in range(d):matrix[i][owner[x^(1<<j)]]+=1
        power=[[int(i==j) for j in range(size)] for i in range(size)];traces=[]
        for degree in range(d+2):
            trace=sum(power[i][i] for i in range(size))
            need(trace==sum(c*(d-2*k)**degree for k,c in enumerate(expected)),
                 'literal weighted flip-matrix spectral moments')
            traces.append(trace)
            power=[[sum(power[i][h]*matrix[h][j] for h in range(size)) for j in range(size)] for i in range(size)]
        result['flip_matrix_trace_moments']=traces
    return result

def bigon(k):
    # Two directed paths u -> ... -> w, plus one isolated vertex.
    n=2*k+1;u=0;w=1;solo=n-1
    paths=[[u]+list(range(2,2+k-1))+[w],
           [u]+list(range(2+k-1,2+2*(k-1)))+[w]]
    arcs={(p[i],p[i+1]) for p in paths for i in range(k)}
    edges={tuple(sorted(e)) for e in arcs}
    swap=list(range(n))
    for a,b in zip(paths[0][1:-1],paths[1][1:-1]):swap[a],swap[b]=b,a
    group=[tuple(range(n)),tuple(swap)]
    need(all({(g[u],g[v]) for u,v in arcs}==arcs for g in group),'bigon automorphisms')
    free=[e for e in combinations(range(n),2) if e not in edges]
    return n,arcs,free,group

def theta(p,k):
    n=p*k+1;paths=[];nextv=2
    for j in range(p):
        paths.append([0]+list(range(nextv,nextv+k-1))+[1]);nextv+=k-1
    arcs={(path[i],path[i+1]) for path in paths for i in range(k)}
    g=list(range(n))
    for j in range(p):
        for a,b in zip(paths[j][1:-1],paths[(j+1)%p][1:-1]):g[a]=b
    need(len(arcs)==n-1,'theta excess minus one')
    need({(g[u],g[v]) for u,v in arcs}==arcs,'theta cyclic symmetry')
    free=[e for e in combinations(range(n),2) if e not in {tuple(sorted(a)) for a in arcs}]
    a=action(free,g)
    fixed=sum(signed_cycle_polynomial(a))
    exponent=k*(p*k-1)//2+(p*p-1)//2
    need(fixed==1<<exponent,'theta fixed-completion formula')
    possible=[]
    for f in range(0,p*k+1,p):
        candidate=(p*k*(p*k-1)+(p-1)*f*(f-1))//(2*p)
        need(candidate!=exponent,'all target C_p-set fixed counts excluded')
        possible.append([f,candidate])
    return {'p':p,'path_length':k,'vertices':n,'fixed_completion_exponent':exponent,
            'possible_target_fixed_vertices_and_exponents':possible}

def positive_carrier(n,arcs,name):
    """Full automorphism group, compared with its induced action on fixed arcs."""
    arcs=set(arcs);edge_list=sorted(arcs);edge_index={e:i for i,e in enumerate(edge_list)}
    group=[p for p in permutations(range(n)) if {(p[u],p[v]) for u,v in arcs}==arcs]
    edges={tuple(sorted(e)) for e in arcs}
    free=[e for e in combinations(range(n),2) if e not in edges]
    edge_group=[tuple(edge_index[(p[u],p[v])] for u,v in edge_list) for p in group]
    source,acts=spectrum(free,group)
    target,tacts=spectrum(list(combinations(range(n-1),2)),edge_group)
    need(source==target,'nonstar pair-level positive control')
    return {'name':name,'vertices':n,'arcs':edge_list,'automorphism_order':len(group),
            'spectrum':source,'source_literal':literal_audit(acts,source),
            'target_literal':literal_audit(tacts,target)}

def main():
    ordinary=[]
    for n in range(2,8):
        f=list(combinations(range(n),2));g=list(permutations(range(n)))
        coeff,acts=spectrum(f,g)
        row={'order':n,'sign_even_edge_spectrum':coeff,'tournament_classes':sum(coeff),
             'self_converse_classes':sum((-1)**k*c for k,c in enumerate(coeff))}
        if n<=5:row['literal']=literal_audit(acts,coeff)
        ordinary.append(row)
    need([r['tournament_classes'] for r in ordinary]==[1,2,4,12,56,456],'classical ordinary audit')
    need([r['self_converse_classes'] for r in ordinary]==[1,2,2,8,12,88],'classical converse audit')
    bigons=[]
    for k in range(2,9):
        n,arcs,free,group=bigon(k);coeff,acts=spectrum(free,group)
        d=k*(2*k-1)
        need(coeff[1]==k*(k-1)+2,'bigon degree-one invariant')
        need(sum(coeff)==1<<(d-1),'weak free-involution orbit count')
        need(sum((-1)**j*c for j,c in enumerate(coeff))==0,'bigon complement has no fixed orbit')
        targets=[]
        for r in range(k+1):
            p=list(range(2*k))
            for j in range(r):p[2*j],p[2*j+1]=2*j+1,2*j
            q,_=spectrum(list(combinations(range(2*k),2)),[tuple(range(2*k)),tuple(p)] if r else [tuple(p)])
            need(q[1]!=coeff[1],'degree-one spectrum excludes every target C2-set')
            if r:
                need(q[1]==k*(k-1)+(k-r)**2,'target degree-one formula')
                need(sum(q)==sum(coeff),'all faithful target actions have equal weak count')
            if r==k:
                need(sum((-1)**j*c for j,c in enumerate(q))==1<<(k*k-1),'edge-set complement discrepancy')
            targets.append({'transpositions':r,'degree_one':q[1],
                            'complement_fixed_orbits':sum((-1)**j*c for j,c in enumerate(q))})
        row={'k':k,'vertices':n,'free_bits':d,'sign_even_spectrum':coeff,
             'target_controls':targets}
        if k<=3:row['literal']=literal_audit(acts,coeff)
        bigons.append(row)
    # Positive linear trijection: an outward star versus its leaf tournament.
    positives=[]
    for leaves in range(2,6):
        group=[tuple(p)+(leaves,) for p in permutations(range(leaves))]
        free=list(combinations(range(leaves),2))
        a,_=spectrum(free,group)
        b,_=spectrum(free,list(permutations(range(leaves))))
        need(a==b,'star deletion preserves every spectral degree')
        positives.append({'leaves':leaves,'spectrum':a})
    nonstars=[positive_carrier(5,[(0,1),(0,2),(1,3),(2,4)],'two-arm outward tree')]
    nonstars += [positive_carrier(k+1,[(i,(i+1)%k) for i in range(k)],
                                 'directed cycle plus isolate') for k in (3,4,5)]
    out={'status':'FINITE-EXACT','ordinary':ordinary,'bigons':bigons,
         'theta_marks':[theta(p,k) for p in (3,5,7) for k in (2,3,4)],
         'positive_stars':positives,'positive_nonstars':nonstars,'checks':CHECKS}
    (HERE/'sign_spectrum.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','checks':CHECKS,'ordinary_counts':[r['tournament_classes'] for r in ordinary],
                      'bigon_k2':bigons[0],'bigon_k3_degree_one':bigons[1]['sign_even_spectrum'][1]}))

if __name__=='__main__':main()
