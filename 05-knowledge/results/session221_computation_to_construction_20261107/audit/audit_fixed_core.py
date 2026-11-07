#!/usr/bin/env python3
import sys,json,random
from pathlib import Path
from itertools import product
from functools import lru_cache
ROOT=Path(__file__).resolve().parent
producer=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT.parent/'automata'
if not producer.is_dir():producer=ROOT.parent/'session221-automata'
sys.path.insert(0,str(producer.resolve()))
import automata as a
def endpoint_deck(out,k):
    n=len(out)
    if n==0:return (1,(),(),())
    incoming=[sum(1<<v for v in range(n) if out[v]>>j&1) for j in range(n)]
    @lru_cache(None)
    def table(mask,end):
        rest=mask^(1<<end)
        if not rest:return tuple(int(v==end) for v in range(n))
        ans=[0]*n;prev=rest&incoming[end]
        while prev:
            b=prev&-prev;prev-=b
            z=table(rest,b.bit_length()-1)
            for i in range(n):ans[i]+=z[i]
        return tuple(ans)
    byend=[table((1<<n)-1,j) for j in range(n)]
    return (sum(map(sum,byend)),tuple(sum(v[i] for v in byend) for i in range(k)),tuple(sum(byend[j]) for j in range(k)),tuple(tuple(byend[j][i] for j in range(k)) for i in range(k)))

def main():
    cases=0;retained=0;rng=random.Random(22199);tests=[]
    for k in range(3):
        for corebits in range(1<<(k*(k-1)//2)):
            apex=[0]*k;q=0
            for i in range(k):
                for j in range(i+1,k):
                    if corebits>>q&1:apex[j]|=1<<i
                    else:apex[i]|=1<<j
                    q+=1
            for n in range(5):tests.extend((tuple(apex),w) for w in product(range(1<<k),repeat=n))
    for _ in range(60):
        k=rng.choice((3,4));apex=[0]*k
        for i in range(k):
            for j in range(i+1,k):
                if rng.randrange(2):apex[i]|=1<<j
                else:apex[j]|=1<<i
        tests.append((tuple(apex),tuple(rng.randrange(1<<k) for _ in range(rng.randrange(1,6)))))
    for apex,word in tests:
        out=a.tournament(apex,word);h,st,en,jt=endpoint_deck(out,len(apex));tf=a.Transfer(apex);state=tf.initial
        for mask in word:state=tf.step(state,mask)
        if (tf.output(state),*tf.endpoints(state),tf.joint_endpoints(state))!=(h,st,en,jt):raise ValueError(('transfer',apex,word))
        nw,removed=a.normalize(apex,word)
        if endpoint_deck(a.tournament(apex,nw),len(apex))[0]!=h or len(nw)>a.length_cap(h):raise ValueError('normal form')
        if True:
            ew,deleted=a.normalize_endpoints(apex,word)
            if endpoint_deck(a.tournament(apex,ew),len(apex))!=(h,st,en,jt) or len(ew)>a.length_cap(h)+2:raise ValueError(('endpoint normal form',apex,word))
            retained+=1
        cases+=1
    controls=0;profiles=set();apex=(2,0)
    for n in range(4):
        for word in product(range(4),repeat=n):
            h,st,en,jt=endpoint_deck(a.tournament(apex,word),2)
            if h<=15:profiles.add((h,st,en,jt))
    for h,st,en,jt in sorted(profiles):
        wanted={'start':list(st),'end':list(en),'joint':[list(row) for row in jt]}
        r=a.compile_target(apex,h,wanted)
        if r['status']!='FOUND' or endpoint_deck(r['outgoing'],2)!=(h,st,en,jt):raise ValueError(('constrained inverse',h,st,en,jt,r))
        controls+=1
    for value in (-1,2):
        if a.compile_target(apex,1,{'start':[value,0]})['status']!='ABSENT':raise ValueError('out of range')
    for value in (7,21):
        if a.compile_target(apex,value)['status']!='ABSENT':raise ValueError('scalar absence')
    result={'status':'PASS','fixed_core_literal_cases':cases,'endpoint_preserving_normal_forms':retained,'constrained_compiler_profiles':controls,'range_rejection_controls':2,'scalar_absence_controls':[7,21],'source_model':'Producer automata.py, tested against independently written backwards endpoint DP','dependency':'Discrete Busch floor as cited in THM4744 Setting; no independent external paper reproof.'}
    (ROOT/'fixed_core_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
