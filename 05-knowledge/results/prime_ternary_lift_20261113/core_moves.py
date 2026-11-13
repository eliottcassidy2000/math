"""Literal core replacements and universal two-run response certificates."""
from pathlib import Path
from itertools import permutations,combinations
import gzip,json
from frame_probe import clone,hp,two_run,strong,valid,need

HERE=Path(__file__).resolve().parent
SOURCE=((2,20,9,35,45,7),0)
PLUS2=((2,12,9,49,135,215,159,15),4)

def apex(word):
    n=len(word);rows=[]
    for i,b in enumerate(word):rows.append(sum(1<<j for j in range(i+1,n))|((b=='0')<<n))
    rows.append(sum((b=='1')<<i for i,b in enumerate(word)))
    return tuple(rows)

PLUS4=(apex('10010110'),4)

def skeleton(T,r,k):
    ext=[v for v in range(len(T)) if v!=r];count=0
    for order in permutations(ext):
        for slots in combinations(range(len(ext)+1),k):
            word=[];slots=set(slots)
            for j in range(len(ext)+1):
                if j in slots:word.append(r)
                if j<len(ext):word.append(order[j])
            if all(T[a]>>b&1 for a,b in zip(word,word[1:])):count+=1
    return count

def embedding(A,r,B,t):
    aa=[v for v in range(len(A)) if v!=r];bb=[v for v in range(len(B)) if v!=t]
    best=None
    for targets in permutations(bb,len(aa)):
        mapping=dict(zip(aa,targets));a=sum(((A[r]>>v)&1)!=((B[t]>>mapping[v])&1) for v in aa)
        b=sum(((A[v]>>w)&1)!=((B[mapping[v]]>>mapping[w])&1) for v,w in combinations(aa,2))
        key=(a,b,targets)
        if best is None or key<best[0]:best=(key,mapping)
    key,mapping=best
    return dict(per_module_vertex=key[0],fixed_core_flips=key[1],mapping=mapping,new_target_vertices=[v for v in bb if v not in mapping.values()])

def physical(A,r,B,t,m,E):
    """Keep every old label, retain the whole old TT_m, append new frame labels."""
    old=[v for v in range(len(A)) if v!=r];oldpos={v:i for i,v in enumerate(old)}
    module=list(range(len(old),len(old)+m));Nold=len(old)+m
    bytarget={E['mapping'][v]:[oldpos[v]] for v in old}
    bytarget[t]=module
    for j,v in enumerate(E['new_target_vertices']):bytarget[v]=[Nold+j]
    out=[0]*(len(B)-1+m)
    for v in range(len(B)):
        for w in range(len(B)):
            if v==w:
                if v==t:
                    for i,a in enumerate(module):
                        for b in module[i+1:]:out[a]|=1<<b
            elif B[v]>>w&1:
                for a in bytarget[v]:
                    for b in bytarget[w]:out[a]|=1<<b
    return tuple(out)

def run():
    profiles=[]
    for name,(T,r) in [('source',SOURCE),('plus2',PLUS2),('plus4',PLUS4)]:
        need(valid(T) and strong(T) and two_run(T,r),'valid strong two-run frame')
        a1=skeleton(T,r,1);a2=skeleton(T,r,2)
        need(hp(clone(T,r,1))==a1 and hp(clone(T,r,2))==a1+2*a2,'independent skeleton response')
        profiles.append(dict(name=name,rows=T,root=r,a1=a1,a2=a2))
    need([(r['a1'],r['a2']) for r in profiles]==[(29,15),(89,45),(91,45)],'exact universal coefficient table')
    moves=[]
    for digit,(T,r) in [(2,PLUS2),(4,PLUS4)]:
        E=embedding(*SOURCE,T,r);checks=[]
        for m in range(1,8):
            old=clone(*SOURCE,m);new=physical(*SOURCE,T,r,m,E)
            need(valid(old) and valid(new) and strong(new),'literal replacement stays a strong tournament')
            h=hp(old);H=hp(new)
            need(h==15*2**m-1 and H==3*h+digit,'exact target, not intermediate observer')
            flips=[(a,b) for a,b in combinations(range(len(old)),2) if ((old[a]>>b)&1)!=((new[a]>>b)&1)]
            need(len(flips)==E['per_module_vertex']*m+E['fixed_core_flips'],'literal old-arc change budget')
            need(all(((old[a]>>b)&1)==((new[a]>>b)&1) for a,b in combinations(range(5,5+m),2)),'old module retained')
            checks.append(dict(m=m,source_h=h,target_h=H,old_order=len(old),new_order=len(new),flips=flips,rows=new if m==1 else None))
        moves.append(dict(digit=digit,target_rows=T,target_root=r,embedding=E,checks=checks))
    result=dict(status='PROVED by universal two-run response, with independent finite certificates',scope='Every transitive module length m>=1; arbitrary module M gives scaled digit d*H(M). Literal checks m1..7, with all old vertices retained.',profiles=profiles,moves=moves)
    suffix='' if __debug__ else '.optimized'
    (HERE/('moves'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('UNIVERSAL',[(r['name'],r['a1'],r['a2']) for r in profiles],flush=True)
    print('REPLACEMENTS',[(r['digit'],r['embedding'],[(c['source_h'],c['target_h']) for c in r['checks'][:3]]) for r in moves],flush=True)

if __name__=='__main__':run()
