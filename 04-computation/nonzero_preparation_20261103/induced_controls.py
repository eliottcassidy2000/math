#!/usr/bin/env python3
"""Independent literal birth controls and actual all-clock clone samples.

No C++ producer implementation or prior mathematical program is imported.
The all-order negative conclusion belongs to the complete C++ search plus
its analytic monotonicity and cited strong-order cutoff, not these samples.
"""
from itertools import permutations, combinations
from functools import cache
from pathlib import Path
import argparse, hashlib, json
GATES=0

def need(value, label):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(label)

def literal(rows):
    n=len(rows); P=[0]*(1<<n); P[0]=1
    E=[[0]*n for _ in P]; S=[[0]*n for _ in P]
    for length in range(1,n+1):
        for order in permutations(range(n),length):
            if all(rows[a]>>b&1 for a,b in zip(order,order[1:])):
                mask=sum(1<<v for v in order)
                P[mask]+=1; E[mask][order[-1]]+=1; S[mask][order[0]]+=1
    return P,E,S

def literal_h(rows):
    return sum(all(rows[a]>>b&1 for a,b in zip(order,order[1:]))
               for order in permutations(range(len(rows))))

def add(rows,incoming):
    n=len(rows)
    return tuple(row | (((incoming>>i)&1)<<n) for i,row in enumerate(rows))+(((1<<n)-1)^incoming,)

def strong(rows):
    for start in range(len(rows)):
        seen={start}; todo=[start]
        while todo:
            v=todo.pop()
            for w in range(len(rows)):
                if rows[v]>>w&1 and w not in seen:
                    seen.add(w); todo.append(w)
        if len(seen)!=len(rows):return False
    return True

def counted(rows):
    n=len(rows); full=(1<<n)-1
    @cache
    def end(mask,v):
        if not (mask>>v&1):return 0
        other=mask^(1<<v)
        if not other:return 1
        return sum(end(other,w) for w in range(n) if other>>w&1 and rows[w]>>v&1)
    @cache
    def start(mask,v):
        if not(mask>>v&1):return 0
        other=mask^(1<<v)
        if not other:return 1
        return sum(start(other,w) for w in range(n) if other>>w&1 and rows[v]>>w&1)
    P=[1]+[sum(end(mask,v) for v in range(n) if mask>>v&1) for mask in range(1,full+1)]
    F=[sum(end(mask,v)*P[full^mask] for mask in range(1,full+1) if mask>>v&1) for v in range(n)]
    starts=[sum(start(mask,v)*P[full^mask] for mask in range(1,full+1) if mask>>v&1) for v in range(n)]
    return P[full],F,starts

def core_signature(rows,port):
    P,E,S=literal(rows); full=len(P)-1
    return [P[full],sum(P[m]*P[full^m] for m in range(full+1)),
            sum(E[m][port]*P[full^m] for m in range(full+1)),
            sum(P[m]*S[full^m][port] for m in range(full+1)),
            sum(E[m][port]*S[full^m][port] for m in range(full+1)),
            E[full][port],S[full][port]]

def prepared(k):
    D=(2,4,57,3,11,27); root=k+7; u=k+6
    rows=[sum(1<<v for v in range(i+1,k)) | (127<<k) for i in range(k)]
    rows += [((d|64)<<k) | ((1<<root) if i==4 else 0) for i,d in enumerate(D)]
    rows += [1<<root,((1<<k)-1)|((127^80)<<k)]
    old=tuple(rows); x=len(old)
    first=add(old,((1<<len(old))-1)^(1<<u))
    final=add(first,(1<<len(old))-1)
    need(strong(old) and strong(first) and strong(final),'actual preparation strong')
    return final,root,x

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args()
    cases=0
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for bits in range(1<<len(pairs)):
            rows=[0]*n
            for q,(a,b) in enumerate(pairs):
                if bits>>q&1:rows[a]|=1<<b
                else:rows[b]|=1<<a
            P,E,S=literal(rows); full=(1<<n)-1
            K=[[sum(E[m][i]*S[full^m][j] for m in range(1,full)) for j in range(n)] for i in range(n)]
            for cut in range(1<<n):
                value=sum(E[full][i] if cut>>i&1 else S[full][i] for i in range(n))
                value+=sum(K[i][j] for i in range(n) for j in range(n) if cut>>i&1 and not(cut>>j&1))
                child=add(rows,cut)
                need(value==literal_h(child),'literal permutation arbitrary birth law')
                need(value>=P[full],'literal insertion monotonicity')
                cases+=1
    need(cases==1098,'complete small universe')
    source,root,x=prepared(1)
    census=json.loads(Path(__file__).with_name('induced_complete.json').read_text())
    need(list(source)==census['source'] and root==census['root'],'same actual source')
    for hit in census['hits']:
        rows=source
        for cut in hit['cuts']:rows=add(rows,cut)
        need(list(rows)==hit['rows'] and strong(rows),'literal strong endpoint and birth sequence')
        h,f,_=counted(rows)
        need((h,f[root])==(hit['H'],hit['F']),'independent endpoint deletion observers')
        need(tuple(r&2047 for r in rows[:11])==source,'old induced graph retained')
    need(len(census['hits'])==11,'all reported strong target controls')
    clocks=[]
    for k in (1,2,3):
        rows,r,x=prepared(k); old=rows; h,f,_=counted(rows)
        for _ in range(2):
            rows=add(rows,((1<<len(rows))-1)^(1<<x))
            need(strong(rows),'strong clone birth')
        hh,ff,ss=counted(rows)
        need(h==756*2**k-83 and f[r]==756*2**k+51,'actual prepared source family')
        need((hh,ff[r])==(3*h+4,3*f[r]),'uniform clone formula literal sample')
        need(ff[r]-hh==398,'cloned physical root defect')
        need(all(g-hh not in (1,2) for g in ss) and ss[x]==hh,'clone all-auxiliary C2 stopping control')
        need(tuple(v&((1<<len(old))-1) for v in rows[:len(old)])==old,'clone old arc preservation')
        clocks.append({'k':k,'H':hh,'F':ff[r],'start_defects':[g-hh for g in ss]})
    # An actual nontrivial-strong hypothesis hostile to the clone stopping lemma.
    one=(0,); aux=1; base=add(one,0); base=add(base,1); base=add(base,5)
    # This is a singleton B, x, followed by transitive Y2, with B->Y2->x->B.
    h,_,s=counted(base)
    need(strong(base) and h==5 and s[-1]-h==1,'singleton clone stopping hostile')
    sig=core_signature((2,4,57,3,75,27,47),4)
    need(sig==[41,712,47,128,0,2,9],'one-birth core hostile signature')
    need((2277*2-263,2277*2+123)==(4291,4677),'one-clock scalar coincidence')
    need(2277*4-263!=3*(756*4-83)+4,'next-clock hostile')
    result={'scope':'independent bounded controls; full census supplied separately',
            'gates':GATES,'literal_birth_cases':cases,'reported_endpoints_recounted':len(census['hits']),
            'clone_controls':clocks,'one_birth_core_signature':sig,
            'census_sha256':hashlib.sha256(Path(__file__).with_name('induced_complete.json').read_bytes()).hexdigest()}
    data=json.dumps(result,indent=2,sort_keys=True)+'\n';args.output.write_text(data)
    print(json.dumps({'gates':GATES,'bytes':len(data.encode()),'sha256':hashlib.sha256(data.encode()).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
