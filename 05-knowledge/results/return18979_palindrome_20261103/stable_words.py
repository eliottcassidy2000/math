"""Literal clone-stable sources, binary carry identities and converse transport."""
from pathlib import Path
import sys,json,subprocess
HERE=Path(__file__).resolve().parent
def need(b,s):
    if not b:raise ValueError(s)
def correction(v):
    need(v>=1,'nonempty strong word');out=2-v;q=v
    while q:
        if not q&1:out+=q>>1
        if q>1 and q&3==1:out+=1
        q//=2
    return out
def source(v):
    word=[0,0]+[(v>>i)&1 for i in range(v.bit_length())];n=len(word)+1;T=[0]*n
    # Apex=0; root=1; ordered transitive chain=1,...,n-1.
    for j,b in enumerate(word,1):
        if b:T[j]|=1
        else:T[0]|=1<<j
    for i in range(1,n):
        for j in range(i+1,n):T[i]|=1<<j
    return tuple(T)
def clone(T,m):
    # The first m chain vertices replace old root1; apex stays0.
    word=[0]*m+[int(T[i]&1) for i in range(2,len(T))];n=len(word)+1;U=[0]*n
    for i,b in enumerate(word,1):
        if b:U[i]|=1
        else:U[0]|=1<<i
        U[i]|=((1<<n)-1)^((1<<(i+1))-1)
    return tuple(U)
def dual(v):
    L=v.bit_length()+2;w=4*v;z=sum(((w>>i)&1)<<(L-i-1) for i in range(L));z=((1<<L)-1)^z
    return z//4 if z%4==0 else None
def birth(T,cut):
    n=len(T);return tuple(row|(((cut>>i)&1)<<n) for i,row in enumerate(T))+((((1<<n)-1)^cut),)
def induced(T,V):return tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(V)) for a in V)
def main():
    states=[];checks=0;dual_count=0;fixed={}
    for v in range(1,1024):
        T=source(v);states.append((v,1,T));c=correction(v)
        need(correction(2*v)==c,'even carry');
        if v>1:need(correction(2*v+1)==c-v-(v&1),'odd carry')
        for m in (2,3):
            U=clone(T,m);need(U==source(v<<(m-1)),'literal clone-to-word transport')
        w=dual(v)
        if w is not None:
            need(dual(w)==v,'word dual involution');U=source(w);p=[0]+list(range(len(T)-1,0,-1))
            need(all(((U[i]>>j)&1)==1-((T[p[i]]>>p[j])&1) for i in range(len(T)) for j in range(len(T)) if i!=j),'literal converse and reverse-chain map')
            need(4*v+c==4*w+correction(w),'dual exact scalar identity');dual_count+=1
            if w==v:fixed[len(T)]=fixed.get(len(T),0)+1
        checks+=1
    birth_receipts=[]
    for m in range(1,21):
        k=m-1;T=clone(source(7),m);n=len(T);c1=(1<<(k+4))-1;U=birth(T,c1)
        c2=((1<<len(U))-1)^(1<<(k+5));V=birth(U,c2)
        p=[0]+list(range(1,k+4))+[n,k+4,n+1,k+5]
        need(induced(V,p)==clone(source(21),m),'two literal births, all old names retained')
        if m<=6:
            states.extend([(7,m,T),(21,m,V)])
            birth_receipts.append(dict(m=m,source=T,cuts=[c1,c2],target=V,canonical_to_actual=p))
    for v in (7,21,149,4768):
        for m in (2,3):states.append((v,m,clone(source(v),m)))
    text=''.join(' '.join(map(str,[len(T),1,0,*T]))+'\n' for _,_,T in states)
    p=subprocess.run([sys.argv[1]],input=text,text=True,capture_output=True,check=True)
    for (v,m,T),line in zip(states,p.stdout.splitlines(),strict=True):
        z=list(map(int,line.split()));need(z[0]==4*v*2**(m-1)+correction(v) and z[5]-z[0]==1,'independent actual count and clone readiness')
    for n,c in fixed.items():need(n%2==1 and c==2**((n-5)//2),'apex-fixed self-converse census')
    hits={2:[],4:[]}
    for v in range(1,200001,2):
        delta=correction(3*v)-3*correction(v)
        if delta in hits:hits[delta].append(v)
    # Exact multiplication-by-three carry state; no inference of global nonexistence.
    for v in range(1,10001):
        A=correction(3*v)-3*correction(v);B=correction(3*v+1)-3*correction(v);C=correction(3*v+2)-3*correction(v);e=v&1
        def carry(x):return tuple(correction(3*x+j)-3*correction(x) for j in range(3))
        need(carry(2*v)==(A,A-3*v-e,B),'even triple carry')
        need(carry(2*v+1)==(B-2+4*e,C+3*v+3*e,C-2+2*e),'odd triple carry')
    out=dict(status='FINITE_EXACT_CHECKS_OF_ALL_SIZE_IDENTITIES',literal_sources=checks,dual_cases=dual_count,apex_fixed_census=fixed,
        stable_18977=dict(v=4768,coefficient=19072,correction=-95,triple_correction=correction(14304),uniform_offset=correction(14304)+285),
        finite_odd_core_search=dict(bound=200000,hits=hits,scope='not a proof of no offset2 or unique offset4 at all sizes'),
        all_size_plus4=dict(source_core=7,target_core=21,source='28*2^(m-1)-5',target='84*2^(m-1)-11',birth_receipts=birth_receipts,
            scope='Two actual births retaining every old name, apex and 1-ready root; uniform in m, but no iterated return from the output family'))
    (HERE/'stable_words.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS stable words, clone laws, converse and carry recurrences',flush=True)
if __name__=='__main__':main()
