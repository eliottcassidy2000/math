#!/usr/bin/env python3
"""Independent structural finite controls and full per-parent comparison.
No producer code is imported, and no known a9 is hard-coded.
"""
import csv, hashlib, itertools, json, math
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent
CHECKS=Counter()
def require(ok,msg):
    CHECKS[msg]+=1
    if not ok: raise RuntimeError(msg)
def write(name,obj):
    (HERE/name).write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
def decode(code,n):
    a=[0]*n; bit=n*(n-1)//2
    for i in range(n):
        for j in range(i+1,n):
            bit-=1
            if code>>bit&1:a[i]|=1<<j
            else:a[j]|=1<<i
    return tuple(a)
def relabel(a,p):
    return tuple(sum((1<<j) for j,w in enumerate(p) if a[v]>>w&1) for v in p)
def degree_key(a):
    cells={}
    for i,x in enumerate(a):cells.setdefault(x.bit_count(),[]).append(i)
    best=None
    for blocks in itertools.product(*(itertools.permutations(cells[d]) for d in sorted(cells))):
        p=sum(blocks,());z=0
        for i in range(len(p)):
            for j in range(i+1,len(p)):z=2*z+((a[p[i]]>>p[j])&1)
        if best is None or z<best:best=z
    return best

def stable(a):
    n=len(a);partition=[tuple(range(n))]
    while True:
        blocks=[]
        for block in partition:
            groups={}
            for v in block:
                incoming=tuple(sum((a[w]>>v)&1 for w in cell) for cell in partition)
                groups.setdefault(incoming,[]).append(v)
            blocks.extend(tuple(groups[k]) for k in sorted(groups))
        if len(blocks)==len(partition):return blocks
        partition=blocks

def exact_error(a,auts):
    cells=stable(a);sizes={v:len(c) for c in cells for v in c};remaining=set(range(len(a)));orbits=[]
    while remaining:
        v=min(remaining);o={g[v] for g in auts};require(o<=remaining,'aut orbits disjoint');remaining-=o;orbits.append(o)
    return sum(sizes[min(o)]-len(o) for o in orbits),sorted(map(len,orbits))

def regular7():
    n=7;out=set()
    def generate(v,a):
        if v==n:
            require(all(x.bit_count()==3 for x in a),'regular7 terminal degrees');out.add(tuple(a));return
        need=3-a[v].bit_count();future=list(range(v+1,n))
        if not 0<=need<=len(future):return
        for wins in itertools.combinations(future,need):
            b=a.copy();w=set(wins)
            for j in future:
                if j in w:b[v]|=1<<j
                else:b[j]|=1<<v
            if any(b[j].bit_count()>3 or b[j].bit_count()+n-v-2<3 for j in future):continue
            generate(v+1,b)
    generate(0,[0]*n)
    remaining=set(out);classes=[]
    while remaining:
        a=min(remaining);orbit=set();auts=[]
        for p in itertools.permutations(range(n)):
            b=relabel(a,p);orbit.add(b)
            if b==a:auts.append(p)
        require(orbit<=out,'regular7 permutation closure');require(orbit<=remaining,'regular7 class orbits disjoint');remaining-=orbit
        e,profile=exact_error(a,auts);require(len(orbit)*len(auts)==math.factorial(n),'regular7 orbit mass')
        classes.append({'rows':a,'aut_order':len(auts),'labelled_orbit':len(orbit),'error':e,'vertex_orbits':profile})
    require(sum(c['labelled_orbit'] for c in classes)==len(out),'regular7 complete mass')
    require(len(out)==2640 and len(classes)==3,'regular7 independent finite universe')
    require(sum(c['error'] for c in classes)==14,'regular7 error fourteen')
    return classes

def triangle_controls():
    perms=list(itertools.permutations(range(3)))
    choices=[('zero',None),('all',None)]+[(d,p) for d in (1,2) for p in perms]
    rot=[tuple(3*i+(x+k[i])%3 for i in range(3) for x in range(3)) for k in itertools.product(range(3),repeat=3)]
    exceptional={};patterns=Counter();regular=0
    for wiring in itertools.product(choices,repeat=3):
        a=[0]*9
        for i in range(3):
            for x in range(3):a[3*i+x]|=1<<(3*i+(x+1)%3)
        for ((i,j),(d,p)) in zip(itertools.combinations(range(3),2),wiring):
            for x in range(3):
                for y in range(3):
                    yes=d=='all' or d==1 and p[x]==y or d==2 and p[x]!=y
                    if yes:a[3*i+x]|=1<<(3*j+y)
                    else:a[3*j+y]|=1<<(3*i+x)
        a=tuple(a);partition=stable(a);shape=tuple(sorted(map(len,partition)));patterns[shape]+=1
        if shape==(9,):regular+=1;continue
        require(shape==(3,3,3),'triangle stable cells only three triples or one nine')
        # Stable cells are the original triangles individually, so the full
        # automorphism group is contained in these 27 rotations.
        aut=[p for p in rot if relabel(a,p)==a]
        e,profile=exact_error(a,aut)
        if not e:continue
        require(e==18 and len(aut)==1,'triangle positive error rigid eighteen')
        degrees=[a[3*i].bit_count() for i in range(3)]
        require(sorted(degrees)==[3,4,5],'positive error transitive majority')
        order=sorted(range(3),key=lambda i:degrees[i]);base=tuple(v for i in order for v in range(3*i,3*i+3));b=relabel(a,base)
        key=min(relabel(b,p) for p in rot)
        exceptional[key]=e
    require(sum(patterns.values())==14**3,'all equitable triangle crossings exhausted')
    require(len(exceptional)==4 and sum(exceptional.values())==72,'all nonregular triple exceptional classes')
    return {'labelled_wirings':sum(patterns.values()),'regular_presentations_excluded':regular,'partition_histogram':[[list(k),v] for k,v in sorted(patterns.items())],'exceptional_classes':[{'rows':k,'error':v} for k,v in sorted(exceptional.items())]}

def compare(producer):
    ours=json.loads((HERE/'compiled-normal.json').read_text());other=json.loads((producer/'normal.json').read_text());atlas=json.loads((producer/'parents-eight.json').read_text())
    require(ours['upper_trace']==other['upper_refinement_trace'],'independent full upper agreement')
    require(ours['pointed_states']==other['pointed_extension_states'],'independent pointed states agreement')
    require(ours['total_error']==other['total_error_correction'],'independent error agreement')
    require(ours['predicted_nine']==other['predicted_class_count'],'independent predicted count agreement')
    records={int(r['canonical_degree_code']):r for r in csv.DictReader((HERE/'compiled-normal-parents.tsv').open(),delimiter='\t')}
    for r,answer in zip(atlas['parents'],other['per_parent']):
        q=records[degree_key(r['rows'])]
        actual=[int(q[k]) for k in ['subset_orbits','upper','root1','root3','root5','root7','root9']]
        require(actual==answer,'all per-parent states upper and root cells agree')
        require(int(q['aut_order'])==len(r['actions']),'all parent automorphism orders agree')
        rows=tuple(r['rows']);groups={tuple(g) for g in r['actions']}
        require(len(groups)==len(r['actions']),'all supplied actions distinct')
        require(all(relabel(rows,g)==rows for g in groups),'all supplied actions preserve arcs')
    raw={}
    for q in csv.DictReader((HERE/'archive'/'merged_nodes.tsv.txt').open(),delimiter='\t'):
        for c in 'ab':raw[int(q['class_'+c])]=int(q['canon_'+c])
    require(all(0<=code<1<<28 for code in raw.values()),'all raw codes have twenty-eight bits')
    almost=[decode(code,8) for i,code in sorted(raw.items()) if all(x.bit_count() in (3,4) for x in decode(code,8))]
    ours_keys=[degree_key(r) for r in almost]
    actual_D={}
    for q in csv.DictReader((HERE/'compiled-normal-regular.tsv').open(),delimiter='\t'):
        actual_D[ours_keys[int(q['source'])],ours_keys[int(q['target'])]]=int(q['multiplicity'])
    reg=json.loads((producer/'regular-normal.json').read_text())
    producer_keys=[degree_key(q['parent']) for q in reg['states']]
    expected_D={(producer_keys[i],producer_keys[j]):weight for i,q in enumerate(reg['states']) for j,weight in q['reroot_row']}
    require(actual_D==expected_D,'all regular reroot matrix entries agree')
    require(sum(actual_D.get((k,k),0) for k in ours_keys)==ours['regular_trace'],'independent regular matrix trace')
    return {'compared_parent_classes':len(records),'producer_prediction':other['predicted_class_count'],'producer_upper_trace':other['upper_refinement_trace']}

def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--producer',type=Path,default=HERE/'archive'/'producer');p.add_argument('--output',default='normal.json');a=p.parse_args()
    r7=regular7();triangles=triangle_controls();comparison=compare(a.producer)
    receipt={'status':'PASS','regular_seven_classes':r7,'triangle_controls':triangles,'comparison':comparison,'checks':dict(CHECKS)}
    write(a.output,receipt);print(json.dumps({'status':'PASS','checks':sum(CHECKS.values()),'regular7_error':sum(c['error'] for c in r7),'triangle_error':sum(c['error'] for c in triangles['exceptional_classes']),'compared_parent_classes':comparison['compared_parent_classes']}))
if __name__=='__main__':main()
