"""Independent literal checks, inherited-basis comparison, marked core formulas."""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(test, message):
    if not test: raise RuntimeError(message)
def decode(n,c):
    a=[0]*n
    for b,(i,j) in enumerate(combinations(range(n),2)):
        if c>>b&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return tuple(a)
def encode(a,p):
    return sum(((a[p[i]]>>p[j])&1)<<b for b,(i,j) in enumerate(combinations(range(len(a)),2)))
def degree_perms(a):
    groups=defaultdict(list)
    for i,r in enumerate(a):groups[r.bit_count()].append(i)
    for parts in product(*(permutations(groups[d]) for d in sorted(groups))):
        yield tuple(x for part in parts for x in part)
def canon(a):return min(encode(a,p) for p in degree_perms(a))
def literal(a):
    return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)) for p in permutations(range(len(a))))
def hp(a):
    n=len(a)
    if not n:return 1
    counts={(1<<i,i):1 for i in range(n)}
    for size in range(1,n):
        nxt=defaultdict(int)
        for (s,i),c in counts.items():
            for j in range(n):
                if not(s>>j&1) and a[i]>>j&1:nxt[s|1<<j,j]+=c
        counts=nxt
    return sum(counts.values())
def induced(a,vertices):
    return tuple(sum(((a[i]>>j)&1)<<v for v,j in enumerate(vertices)) for i in vertices)
def deletes(a):return [hp(induced(a,[j for j in range(len(a)) if i!=j])) for i in range(len(a))]
def reach(a,i):
    seen={i};todo=[i]
    while todo:
        v=todo.pop()
        for j in range(len(a)):
            if a[v]>>j&1 and j not in seen:seen.add(j);todo.append(j)
    return seen
def scc(a):
    n=len(a);rs=[reach(a,i) for i in range(n)];unused=set(range(n));parts=[]
    while unused:
        i=min(unused);part={j for j in unused if j in rs[i] and i in rs[j]}
        unused-=part;parts.append(tuple(sorted(part)))
    return sorted(parts,key=lambda part:-sum(bool(a[part[0]]>>other[0]&1) for other in parts if part!=other))
def join(graphs):
    n=sum(map(len,graphs));a=[];offset=0
    for g in graphs:
        later=((1<<n)-1)^((1<<(offset+len(g)))-1)
        a += [(r<<offset)|later for r in g];offset+=len(g)
    return tuple(a)
def profile(a):
    # All degree-preserving vertex maps, independently of the C++ orbit code.
    degrees=[r.bit_count() for r in a];blocks=defaultdict(list)
    for i,d in enumerate(degrees):blocks[d].append(i)
    group=[]
    for images in product(*(permutations(vs) for vs in blocks.values())):
        mapping=[0]*len(a)
        for domain,image in zip(blocks.values(),images):
            for i,j in zip(domain,image):mapping[i]=j
        if all(((a[i]>>j)&1)==((a[mapping[i]]>>mapping[j])&1) for i in range(len(a)) for j in range(len(a))):group.append(mapping)
    one={tuple(sorted(p[i] for p in group)) for i in range(len(a))}
    two={tuple(sorted(p[i]*len(a)+p[j] for p in group)) for i in range(len(a)) for j in range(len(a)) if i!=j}
    return len(group),len(one),len(two)

def main():
    census=json.loads((HERE/'census.json').read_text(encoding='utf-8-sig'))
    records=[];atoms=[];checked=0
    for row in census['universes']:
        n=row['n'];need(row['frames']==2**((n-1)*(n-2)//2),'frame universe')
        mass=0;observed=[];induced_min=[]
        for item in row['classes']:
            a=decode(n,item['code']);h=literal(a);checked+=1
            need(h==hp(a)==item['H'],'literal paths')
            need(canon(a)==item['code'],'independent canonical representative')
            auto,r,q=profile(a);need((auto,r,q)==(item['aut'],item['roots'],item['ordered_distinct_ports']),'automorphism and marked counts')
            mass+=factorial(n)//auto;parts=scc(a)
            need((len(parts)==1)==item['strong'],'strongness')
            val=1
            for part in parts:val*=hp(induced(a,part))
            need(val==h,'SCC multiplicativity')
            ds=deletes(a);need(max(ds)==item['max_delete_H'],'deletion maxima')
            need(all(d<=h for d in ds),'insertion monotonicity')
            need((h>21 and max(ds)<=21)==item['minimal_above_21'],'minimality')
            rec={'n':n,'code':item['code'],'rows':list(a),'H':h,'aut':auto,'roots':r,'ordered_distinct_ports':q,'deletion_H':ds,'strong':len(parts)==1}
            records.append(rec)
            if item['strong'] and h>1 and h<=21:atoms.append(rec)
            if item['strong']:observed.append(h)
            if item['minimal_above_21']:induced_min.append(rec)
        need(mass==2**(n*(n-1)//2)==row['labelled_mass'],'orbit mass')
        need((min(observed) if observed else 0)==row['strong_floor'],'strong floor')
    need([r['class_count'] for r in census['universes']]==[1,1,2,4,12,56,456],'class counts')
    need([r['strong_floor'] for r in census['universes']][2:]==[3,5,9,15,25],'strong floors')

    # Complete minimal nonstrong possibilities follow from the written SCC proof.
    c3=next(r for r in atoms if r['n']==3);c4=next(r for r in atoms if r['n']==4)
    words=[(c3,c3,c3),(c4,c4)]
    for a in atoms:
        if a['n']==5:words += [(c3,a),(a,c3)]
    basis=[r for r in records if r['H']>21 and max(r['deletion_H'])<=21]
    for word in words:
        a=join([r['rows'] for r in word]);h=hp(a);ds=deletes(a)
        need(h>21 and max(ds)<=21,'nonstrong obstruction minimality')
        basis.append({'n':len(a),'code':canon(a),'rows':list(a),'H':h,'deletion_H':ds,'strong':False})
    need(len(basis)==56 and Counter(r['n'] for r in basis)=={6:29,7:13,8:13,9:1},'basis census')
    inherited=ROOT/'05-knowledge/results/session133_fractran_spectrum_refinement_20261101/spectrum/basis21.json'
    old=json.loads(inherited.read_text(encoding='utf8'));oldbasis=old['strong_obstructions']+old['nonstrong_obstructions']
    need({(r['n'],r['code']) for r in basis}=={(r['n'],canon(r['rows'])) for r in oldbasis},'complete inherited basis equality')
    # Explicit hostiles: adding a singleton is no longer minimal; triples show order9 is sharp.
    need(max(deletes(join([c3['rows']]*3)))==9,'triangle-chain deletion')
    padded=join([[0],basis[0]['rows']]);need(max(deletes(padded))>21,'padded obstruction hostile')

    cores=[]
    def visit(word,h):
        m=sum(r['n'] for r in word);k=len(word);r=sum(a['roots'] for a in word)
        q=sum(a['ordered_distinct_ports'] for a in word)+2*sum(a['roots']*b['roots'] for a,b in combinations(word,2))
        cores.append({'H':h,'m':m,'k':k,'root_orbits':r,'pair_orbits':q,'atoms':[(a['n'],a['code']) for a in word]})
        for a in atoms:
            if h*a['H']<=21:visit(word+[a],h*a['H'])
    visit([],1)
    # Weight each singleton allocation by its actual marked-orbit counts.
    formulas={str(h):[c for c in cores if c['H']==h] for h in range(1,22,2)}
    for n in range(1,8):
        for h in range(1,22,2):
            totals=[0,0,0]
            for c in formulas[str(h)]:
                j=n-c['m']
                if j<0:continue
                mult=comb(j+c['k'],c['k']);R=c['root_orbits'];Q=c['pair_orbits']
                weights=[1,R+j,Q+2*R*j+j*(j-1)]
                totals=[a+mult*b for a,b in zip(totals,weights)]
            actual=[sum(r[key] if key else 1 for r in records if r['n']==n and r['H']==h) for key in [None,'roots','ordered_distinct_ports']]
            need(totals==actual,'marked generating functions vs entire class census')
    need(profile(c3['rows'])==(3,1,2),'joint-port triangle hostile')
    for n in range(3,35):
        c=formulas['3'][0];j=n-c['m'];mult=comb(j+c['k'],c['k'])
        need(mult*(c['root_orbits']+j)==(n-2)**2,'rooted triangle family')
        need(mult*(c['pair_orbits']+2*j+j*(j-1))==(n-2)*((n-3)*(n-2)+2),'two-port triangle family')
    output={'status':'PASS','independent_literal_classes':checked,'basis_counts':dict(sorted(Counter(r['n'] for r in basis).items())),
            'inherited_basis':'exact graph-set equality after canonicalization','strong_atoms':atoms,'minimal_obstructions':basis,
            'core_formulas':formulas,'controls':['all vertex deletions','ordered triangle chain','nonminimal singleton padding','C3 one-port versus two-port orbit distinction']}
    (HERE/'audit.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf8')
    print(f'PASS: {checked} literal classes; 56 exact inherited obstructions; all-order marked core formulas')
if __name__=='__main__':main()
