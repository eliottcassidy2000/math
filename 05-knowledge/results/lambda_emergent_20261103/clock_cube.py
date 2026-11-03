"""Commuting graph constructions, an idempotent repair, and a noncommuting peak clock."""
from itertools import permutations,product
from math import comb
from collections import Counter
import json
from effects import HERE,old,junction,need

def replay(w,u,word):
    T=old.source(w,u);m=len(w);base=list(range(m));a,b=m,m+1;ell=m-1
    native=helpers=0;peak=len(T);history=[]
    for command in word:
        if command=='P':
            for i in base[:-1]:
                if T[a]>>i&1:continue
                path=(i,a,b)+tuple(v for v in base if v!=i);n=len(T)
                need(all(T[x]>>y&1 for x,y in zip(path,path[1:])),'actual source guard')
                X=junction.ear(T,((1<<n)-1)^(1<<i));peak=max(peak,n+1)
                U=junction.native(junction.native(X,path+(n,)),path[:1]+(n,)+path[1:])
                T=tuple(row&((1<<n)-1) for row in U[:n]);native+=2;helpers+=1
        else:
            n=len(T);incoming=((1<<a)|(1<<b)) if command=='U' else (((1<<n)-1)^(1<<ell))
            T=junction.ear(T,incoming);peak=max(peak,len(T))
            if command=='U':base.insert(0,n)
            else:base.insert(len(base)-1,n)
        need(junction.strong(T),'all strong construction endpoints')
        history.append({'operation':command,'order':len(T),'native_so_far':native,'peak_so_far':peak})
    canonical=junction.sub(T,tuple(base+[a,b]))
    return T,a,b,canonical,{'word':''.join(word),'native':native,'helpers':helpers,'peak':peak,'history':history}

def main():
    examples=[]
    for w,u in (old.chart(3,4),((1,0,0,0,1,0),(1,1,1,1,1,0))):
        traces=[];target=None;z=w[:-1].count(0);n=len(w)+2
        for word in permutations('PUV'):
            T,a,b,C,receipt=replay(w,u,word)
            if target is None:target=C
            need(C==target,'same actual graph after all six construction orders')
            need(receipt['native']==2*z and receipt['helpers']==z,'additive clocks commute')
            need(receipt['peak']==n+2+int(word[-1]=='P' and z>0),'peak resource remembers evaluation order')
            need(junction.profile(T,a)[3]-junction.H(T)==1,'every ordering ends ready')
            traces.append(receipt)
        T,a,b,C,receipt=replay(w,u,'PPUV')
        need(C==target and receipt['native']==2*z,'idempotent normalization including total native cost')
        examples.append({'initial_words':[w,u],'initial_missing_arcs':z,'final_canonical_rows':target,'six_orders':traces})
    fibres=[]
    for m in range(2,9):
        counts=Counter(mid.count(0) for mid in product((0,1),repeat=m-2))
        need(all(counts[z]==comb(m-2,z) for z in range(m-1)),'exact labelled normalizer fibre clock polynomial')
        fibres.append({'base_size':m,'helper_clock_coefficients':[counts[z] for z in range(m-1)]})
    (HERE/'clock_cube.json').write_text(json.dumps({'status':'FINITE_EXACT_LITERAL_CONSTRUCTIONS','examples':examples,'normalizer_fibres':fibres},indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: idempotent repair commutes on actual graphs; native clocks agree and peak clock distinguishes order')
if __name__=='__main__':main()
