#!/usr/bin/env python3
"""Independent pointwise checks of saved partitions, without importing producer."""
from collections import Counter, defaultdict
from pathlib import Path
import json
ROOT = Path(__file__).resolve().parent

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def reflect(x, n):
    coords = [(i,j) for i in range(n) for j in range(i+2,n)]
    values = {p:(x>>k)&1 for k,p in enumerate(coords)}
    return sum(values[n-1-j,n-1-i]<<k for k,(i,j) in enumerate(coords))

def scores(x, n):
    out = [0]*n
    for i in range(n-1):out[i]+=1
    k=0
    for i in range(n):
        for j in range(i+2,n):
            out[i if (x>>k)&1 else j]+=1;k+=1
    return tuple(sorted(out))

def components(edges, active):
    adj = defaultdict(set)
    for (u,v),w in edges.items():
        if w:adj[u].add(v);adj[v].add(u)
    unseen=set(active);count=0
    while unseen:
        count+=1;stack=[unseen.pop()]
        while stack:
            for v in adj[stack.pop()] & unseen:
                unseen.remove(v);stack.append(v)
    return count

def inspect(g):
    n=g['n'];m=(n-1)*(n-2)//2;f=(n-1)//2
    total=1<<m;full=total-1;blue_mass=1<<((m+f)//2)
    p=g['merged_by_word'];c=[tuple(x) for x in g['ordinary_by_word']]
    need(len(p)==total==len(c),'word domain')
    fibers=Counter(p);ordinary=Counter(c)
    need(set(fibers)==set(range(len(g['degrees']))),'nonempty merged nodes')
    need(all(v%2 for v in ordinary.values()),'odd ordinary fibres')
    blue=Counter();black=Counter();bdeg=Counter();fix=0
    for x in range(total):
        rx=reflect(x,n);cx=x^full
        need(reflect(rx,n)==x,'R involution')
        need(reflect(cx,n)==(rx^full),'C,R commute')
        need(rx!=cx,'CR fixed-point-free')
        need(p[x]==p[rx],'merged R invariance')
        u=p[x]
        need(c[x][0]==u,'ordinary-to-merged map')
        if g['self_converse'][u]:
            need(c[rx]==c[x],'SC ordinary R invariance')
        else:
            need(c[rx]==(u,1-c[x][1]),'paired ordinary R action')
        if rx==x:fix+=1;bdeg[u]+=1
        if x<cx:
            pair=tuple(sorted((u,p[cx])))
            (blue if rx==x else black)[pair]+=1
    need(fix==blue_mass,'reflection fixed mass')
    need(all(w%2==0 for w in black.values()),'black twin multiplicities')
    need(blue==Counter({(u,v):w for u,v,w in g['blue_edges']}),'blue graph replay')
    need(black==Counter({(u,v):2*w for u,v,w in g['black_half_edges']}),'black graph replay')
    d=[fibers[u] for u in range(len(fibers))]
    need(d==g['degrees'],'degree equals fibre size')
    need(all((d[u]%2==1)==g['self_converse'][u]==bool(bdeg[u]) for u in range(len(d))),'SC characterization')
    need(all(d[u]%4==2 for u in range(len(d)) if not bdeg[u]),'NS mod4')
    need(all(b%2 for b in bdeg.values()),'blue degrees odd')
    V=len(d);S=len(bdeg);A=len(ordinary)
    need(A==2*V-S,'class count formula')
    need((S-4*V+3*sum(z*z for z in d))%8==0,'second moment congruence')
    for j in range(2,7):
        mod=1<<(j+2)
        need((sum(pow(z,1<<j,mod) for z in d)-S)%mod==0,'Euler moment congruence')
    # A concrete witness that this is not the tournament-isomorphism partition.
    representatives={};witness=None
    for x,key in enumerate(c):
        sig=scores(x,n)
        if key in representatives and representatives[key][1]!=sig:
            witness={'ordinary_fibre':list(key),'words':[representatives[key][0],x],
                     'score_multisets':[list(representatives[key][1]),list(sig)]}
            break
        representatives[key]=(x,sig)
    need(witness is not None,'explicit non-native witness')
    combined=blue+black
    result={'n':n,'N':total,'B':fix,'V':V,'S':S,'A_abstract':A,
            'M2':sum(z*z for z in d),'max_merged_degree':max(d),
            'blue_loops':sum(w for (u,v),w in blue.items() if u==v),
            'black_loops':sum(w for (u,v),w in black.items() if u==v),
            'blue_components_on_support':components(blue,set(bdeg)),
            'black_components_on_support':components(black,{v for e in black for v in e}),
            'full_components':components(combined,set(fibers)),
            'non_native_witness':witness}
    return result,tuple(sorted(combined)),Counter((d[u],bdeg[u]) for u in range(V))

def main():
    data=json.loads((ROOT/'certificate.json').read_text())
    star=[inspect(g) for g in data['models']]
    bounded=[inspect(g) for g in data['bounded_models']]
    need(star[0][1]==star[1][1],'same star support')
    for field in ['V','N','B','M2']:
        need(star[0][0][field]==star[1][0][field],'star aggregate collision')
        need(bounded[0][0][field]==bounded[1][0][field],'bounded aggregate collision')
    need([z[0]['S'] for z in star]==[2,10],'star S values')
    need([z[0]['S'] for z in bounded]==[12,20],'bounded S values')
    native=json.loads((ROOT/'native_reference.json').read_text())
    need([(r['n'],r['A'],r['V'],r['S']) for r in native]==[(3,2,2,2),(4,4,3,2),(5,12,10,8),(6,56,34,12)],'native small controls')
    native_hist=Counter((z['degree'],z['blue_degree']) for z in native[-1]['nodes'])
    need(bounded[0][2]==native_hist,'native degree-blue histogram match')
    need(set(bounded[1][2])<=set(native_hist),'bounded observed signature support')
    need({b for d,b in native_hist if b>0}=={1,3,5,7,9},'native blue palette')
    need({d for d,b in native_hist if b==1}=={1},'singleton-only blue degree one')
    ambiguity=[s for s in range(0,35,2) if 3*s-4<=64<=9*s-16 and (s-4*34+3*48724)%8==0]
    need(ambiguity==[12,20],'exact aggregate ambiguity set')
    for result,_,hist in bounded:
        need(hist[1,1]==2,'two native singleton pure-blue vertices')
        need(result['blue_loops']==native[-1]['blue_loops']==2 and result['black_loops']==native[-1]['black_loops']==24,'native total loop counts')
        need(result['blue_components_on_support']==result['black_components_on_support']==result['full_components']==1,'layer connectivity')
    # Validate all declared finite-moment signed relations, independently.
    for row in data['hierarchy']:
        degree=row['moment_degree'];rel={int(k):v for k,v in row['signed_relation'].items()}
        for q in range(degree+1):need(sum(c*d**q for d,c in rel.items())==0,'signed moment relation')
        need(row['degree_moments'][0]==row['degree_moments'][1],'hierarchy moments')
        need(abs(row['S'][0]-row['S'][1])==row['denominator'],'SC separation')
    archived=json.loads((ROOT/'archive'/'fibre_inverse.json').read_text())
    for local in native:
        old=next(x for x in archived['sizes'] if x['n']==local['n'])
        need(Counter(z['degree'] for z in local['nodes'])==Counter(z['predicted_tilings'] for z in old['node_records']),'historical degree histogram')
    print(json.dumps({'status':'PASS','star_countermodels':[z[0] for z in star],
                      'bounded_countermodels':[z[0] for z in bounded],
                      'finite_moment_orders_checked':list(range(1,9)),
                      'native_orders_checked':[3,4,5,6],
                      'bounded_exact_A_ambiguity_set':[48,56]},indent=2,sort_keys=True))

if __name__=='__main__':main()
