#!/usr/bin/env python3
"""Independent order-ten stable-cell classification audit. No producer imports."""
import argparse, itertools as it, json
from collections import Counter
from pathlib import Path

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def cyc(n): return tuple(sum(1<<((i+k)%n) for k in range(1,(n+1)//2)) for i in range(n))
def relabel(rows,p):
    out=[0]*len(rows)
    for i,row in enumerate(rows): out[p[i]]=sum(1<<p[j] for j in range(len(rows)) if row>>j&1)
    return tuple(out)
def stable(rows):
    c=(0,)*len(rows)
    while True:
        sig=tuple((c[i],)+tuple(sum(rows[i]>>j&1 for j in range(len(rows)) if c[j]==k) for k in sorted(set(c))) for i in range(len(rows)))
        names={s:k for k,s in enumerate(sorted(set(sig)))}; d=tuple(names[s] for s in sig)
        if len(set(c))==len(set(d)): return tuple(tuple(i for i in range(len(rows)) if d[i]==k) for k in sorted(set(d)))
        c=d
def strong(rows):
    n=len(rows)
    for start in range(n):
        seen=todo=1<<start
        while todo:
            bit=todo&-todo; todo^=bit; more=rows[bit.bit_length()-1]&~seen; todo|=more; seen|=more
        if seen!=(1<<n)-1: return False
    return True
def autos(rows):
    n=len(rows); cells=stable(rows); targets={v:c for c in cells for v in c}; order=sorted(range(n),key=lambda v:(len(targets[v]),v)); p=[-1]*n; out=[]
    def visit(k,used):
        if k==n: out.append(tuple(p)); return
        v=order[k]
        for w in targets[v]:
            if used>>w&1: continue
            if any((rows[v]>>u&1)!=(rows[w]>>p[u]&1) for u in order[:k]): continue
            p[v]=w; visit(k+1,used|1<<w); p[v]=-1
    visit(0,0); return tuple(out)
def error(rows,group):
    sizes={v:len(c) for c in stable(rows) for v in c}; seen=set(); total=0; orbits=[]
    for v in range(len(rows)):
        if v in seen: continue
        orbit={p[v] for p in group}; seen|=orbit; orbits.append(sorted(orbit)); total+=sizes[v]-len(orbit)
    return total,orbits
def graph(n,bits):
    rows=[0]*n
    for k,(i,j) in enumerate(it.combinations(range(n),2)):
        a,b=(i,j) if bits>>k&1 else (j,i); rows[a]|=1<<b
    return tuple(rows)
def partitions(n,lo=1):
    if not n: yield ()
    for a in range(lo,n+1,2):
        for tail in partitions(n-a,a): yield (a,)+tail

def five_rows(m):
    rows=list(cyc(5))+[r<<5 for r in cyc(5)]
    for i in range(5):
        for j in range(5):
            a,b=(i,5+j) if m[i]>>j&1 else (5+j,i); rows[a]|=1<<b
    return tuple(rows)
def move_matrix(m,a,b):
    out=[0]*5
    for i,row in enumerate(m): out[(i+a)%5]=sum(1<<((j+b)%5) for j in range(5) if row>>j&1)
    return tuple(out)
def five_sector():
    result={}; reps=[]
    for r in (1,2):
        choices=[sum(1<<j for j in c) for c in it.combinations(range(5),r)]
        matrices={x for x in it.product(choices,repeat=5) if all(sum(row>>j&1 for row in x)==r for j in range(5))}
        count=len(matrices); hist=Counter(); subtotal=0
        while matrices:
            m=min(matrices); orbit={move_matrix(m,a,b) for a in range(5) for b in range(5)}
            need(orbit<=matrices,'five orbit overlap'); matrices-=orbit; rows=five_rows(m)
            need(strong(rows),'five strong'); need(sorted(map(len,stable(rows)))==[5,5],'five stable')
            aut=autos(rows); e,orbits=error(rows,aut)
            need(len(aut)*len(orbit)==25,'five full aut versus rotations'); need(len(aut) in (1,5),'five group'); need(e==(0 if len(aut)==5 else 40),'five error')
            ct=tuple(sum(1<<j for j in range(5) if not(m[j]>>i&1)) for i in range(5))
            need(relabel(rows,tuple(range(5,10))+tuple(range(5)))==five_rows(ct),'five high degree swap')
            hist[(len(aut),e)]+=1; subtotal+=e; reps.append({'sector':'five','r':r,'rows':rows,'aut_order':len(aut),'error':e,'vertex_orbits':orbits})
        result[str(r)]={'presentations':count,'classes_by_aut_error':{str(k):v for k,v in sorted(hist.items())},'error':subtotal}
    need(result['1']['presentations']==120 and result['2']['presentations']==2040,'five universes')
    need(sum(v['error'] for v in result.values())==3360,'five total'); return result,reps

HOSTILE7=(74,84,97,22,37,11,56)
def substitute7(q):
    rows=list(HOSTILE7)+[0,0,0]; blocks=[tuple(range(7)),(7,),(8,),(9,)]
    for i,j in it.combinations(range(4),2):
        a,b=(i,j) if q[i]>>j&1 else (j,i)
        for x in blocks[a]:
            for y in blocks[b]: rows[x]|=1<<y
    return tuple(rows)
def seven_sector():
    labels={graph(4,b) for b in range(64)}; group=[(0,)+p for p in it.permutations((1,2,3))]; hist=Counter(); reps=[]
    while labels:
        q=min(labels); orbit={relabel(q,p) for p in group}; need(orbit<=labels,'seven quotient overlap'); labels-=orbit
        rows=substitute7(q); yes=strong(q); need(strong(rows)==yes,'seven strong equivalence')
        if not yes: hist['nonstrong']+=1; continue
        hist['strong']+=1; need(sorted(map(len,stable(rows)))==[1,1,1,7],'seven stable'); aut=autos(rows); e,orbits=error(rows,aut)
        need(len(aut)==3 and e==14,'seven inheritance'); reps.append({'sector':'seven','quotient_rows':q,'rows':rows,'aut_order':len(aut),'error':e,'vertex_orbits':orbits})
    need(hist['strong']==4,'seven rooted count'); return {'rooted_quotient_classes':dict(hist),'error':sum(x['error'] for x in reps)},reps

def triangle_rows(crossings,mask):
    rows=[r<<(3*i) for i in range(3) for r in cyc(3)]+[0]
    for (i,j),m in zip(it.combinations(range(3),2),crossings):
        for x in range(3):
            for y in range(3):
                a,b=(3*i+x,3*j+y) if m[x]>>y&1 else (3*j+y,3*i+x); rows[a]|=1<<b
    for i in range(3):
        for x in range(3):
            a,b=(3*i+x,9) if mask>>i&1 else (9,3*i+x); rows[a]|=1<<b
    return tuple(rows)
def triangle_group():
    for p in it.permutations(range(3)):
        for k in it.product(range(3),repeat=3): yield tuple(3*p[i]+(x+k[i])%3 for i in range(3) for x in range(3))+(9,)
def triangle_sector():
    pm=[tuple(1<<p[x] for x in range(3)) for p in it.permutations(range(3))]; choices=pm+[tuple(7^v for v in m) for m in pm]
    labels={triangle_rows(m,mask) for m in it.product(choices,repeat=3) for mask in range(1,7)}; count=len(labels); group=tuple(triangle_group()); reps=[]; hist=Counter(); skeletons=Counter()
    while labels:
        rows=min(labels); orbit={relabel(rows,p) for p in group}; need(orbit<=labels,'triangle orbit overlap'); labels-=orbit
        need(strong(rows),'triangle strong'); need(sorted(map(len,stable(rows)))==[1,3,3,3],'triangle stable'); aut=autos(rows); e,orbits=error(rows,aut)
        need(len(aut)*len(orbit)==162,'triangle hidden automorphism'); need(len(aut) in (1,3) and e==(18 if len(aut)==1 else 0),'triangle error')
        majority=tuple(sum(sum(rows[3*i]>>(3*j+y)&1 for y in range(3))==2 for j in range(3) if i!=j) for i in range(3)); kind='cyclic' if majority==(1,1,1) else 'transitive'
        hist[(kind,len(aut),e)]+=1; skeletons[(kind,rows[9].bit_count())]+=1; reps.append({'sector':'triangles','majority':kind,'rows':rows,'aut_order':len(aut),'error':e,'vertex_orbits':orbits})
    need(count==10368 and len(reps)==128,'triangle universes'); need(sum(x['error'] for x in reps)==576,'triangle total')
    # Exhaust every equitable three-triangle crossing and every uniform root
    # attachment, including omitted-support and regular-nine-overlap hostiles.
    all_choices=choices+[(0,0,0),(7,7,7)]
    coverage=Counter()
    for crosses in it.product(all_choices,repeat=3):
        for mask in range(8):
            rows=triangle_rows(crosses,mask)
            if not strong(rows):
                coverage['nonstrong_presentations']+=1
                continue
            coverage['strong_presentations']+=1
            need(sorted(map(len,stable(rows)))==[1,3,3,3],'strong three-cell completion unexpectedly coarsens')
            e,_=error(rows,autos(rows))
            if e:
                coverage['strong_positive_error_presentations']+=1
                need(mask in range(1,7) and all(c in choices for c in crosses),'omitted strong error outside complete support')
                need(e==18,'unexpected strong error in all crossing census')
    need(sum(coverage[k] for k in ('nonstrong_presentations','strong_presentations'))==14**3*8,'all crossing coverage')
    need(coverage['strong_positive_error_presentations']==5184,'all positive crossing coverage')
    controls={}; plus=pm[0]; minus=(1,4,2)
    for name,crosses,mask in [('unbalanced_mixed',(plus,plus,minus),1),('balanced_mixed',(plus,plus,plus),1),('missing_support_edge',((0,0,0),plus,minus),1),('unbalanced_sink_root',(plus,plus,minus),7),('unbalanced_source_root',(plus,plus,minus),0),('cyclic_unbalanced_mixed',(plus,tuple(7^v for v in plus),minus),1)]:
        rows=triangle_rows(crosses,mask); aut=autos(rows); e,_=error(rows,aut); controls[name]={'strong':strong(rows),'stable_sizes':sorted(map(len,stable(rows))),'aut_order':len(aut),'error':e,'deleted_root_stable_sizes':sorted(map(len,stable(tuple(v&511 for v in rows[:9]))))}
    need(controls['unbalanced_mixed']['error']==18 and controls['balanced_mixed']['error']==0 and controls['missing_support_edge']['error']==0,'signed controls')
    need(not controls['unbalanced_sink_root']['strong'] and not controls['unbalanced_source_root']['strong'],'strongness hostile')
    need(controls['cyclic_unbalanced_mixed']['strong'] and controls['cyclic_unbalanced_mixed']['error']==18 and controls['cyclic_unbalanced_mixed']['deleted_root_stable_sizes']==[9],'cyclic-majority overlap hostile')
    return {'all_equitable_crossings_coverage':dict(coverage),'presentations':count,'classes':len(reps),'classes_by_majority_aut_error':{str(k):v for k,v in sorted(hist.items())},'classes_by_majority_root_degree':{str(k):v for k,v in sorted(skeletons.items())},'error':sum(x['error'] for x in reps),'controls':controls},reps

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',default='receipt.json'); args=ap.parse_args()
    five,fr=five_sector(); seven,sr=seven_sector(); triangles,tr=triangle_sector(); se=sum(five[str(r)]['error'] for r in (1,2))+seven['error']+triangles['error']; a=[1,1,1,2]; square=[sum(a[i]*a[n-i] for i in range(n+1)) for n in range(4)]; full=se+square[1]*716+square[3]*14
    need(se==3992 and full==5508,'full correction'); result={'status':'independent finite-exact classification audit','odd_partitions_10':list(partitions(10)),'five_sector':five,'seven_sector':seven,'triangle_sector':triangles,'scc_transfer':{'A_coefficients_0_to_3':a,'A_squared_coefficients_0_to_3':square,'inherited_strong_error_7':14,'inherited_strong_error_9':716,'strong_error_10':se,'full_error_10':full},'representatives':fr+sr+tr}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n'); print(json.dumps({k:v for k,v in result.items() if k!='representatives'},indent=2,sort_keys=True))
if __name__=='__main__': main()
