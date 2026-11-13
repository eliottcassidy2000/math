"""Exact color-word compilation for an entry root and an output root in the context."""
from itertools import permutations,combinations
from math import factorial
import json,gzip
import numpy as np
from typed_glue import HERE,save,need,seven,glue,decode
from feedback_probe import all_profiles

# Feature order h,s,e,q,F,G,P of the arbitrary attached rooted tournament B.
UNIT=[tuple(int(i==j) for i in range(7)) for j in range(7)]


def feature(k,inc,out):
    z=[0]*7
    if k==1:z[1 if out else 2 if inc else 3]=1
    elif inc and out:z[0]=1;z[1]=z[2]=-1
    elif inc:z[4]=1;z[3]=z[2]=-1
    elif out:z[5]=1;z[3]=z[1]=-1
    else:z[6]=1;z[3]=-2
    return z


def coefficient(A,r,u,kind):
    """Count H, start(u), or first-bin-end(u); nonroot B vertices form <=2 runs."""
    n=len(A);z=[0]*7
    for k in (1,2):
        N=n+k
        for order in permutations(range(n)):
            for positions in combinations(range(N),k):
                slots=set(positions);it=iter(order);word=tuple(-1 if j in slots else next(it) for j in range(N))
                if kind=='s' and word[0]!=u:continue
                cut=word.index(u)+1 if kind=='F' else N
                inc=out=False;valid=True
                for j in range(N-1):
                    if j+1==cut:continue
                    a,b=word[j:j+2]
                    if a==-1:
                        if b!=r:valid=False;break
                        inc=True
                    elif b==-1:
                        if a==r:out=True
                    elif not(A[a]>>b&1):valid=False;break
                if valid:
                    v=feature(k,inc,out)
                    for j in range(7):z[j]+=v[j]
    return tuple(z)


def compiler(A,r,u):
    h=coefficient(A,r,u,'h');s=coefficient(A,r,u,'s');F=coefficient(A,r,u,'F')
    if u==r:
        z=seven(A,r);q=(0,0,0,z[3],0,0,0)
    else:
        order=[v for v in range(len(A)) if v!=u]
        B=tuple(sum(((A[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)
        q=coefficient(B,order.index(r),0,'h')
    rows=(h,s,q,F)
    need(all(row[j]==0 for row in rows for j in (2,5,6)),'four-coordinate closure')
    return tuple(tuple(row[j] for j in (0,1,3,4)) for row in rows)


def mul(M,z):return tuple(sum(a*b for a,b in zip(row,z)) for row in M)


def run():
    controls=0;small=[]
    for n in (2,3):
        for code in range(1<<(n*(n-1)//2)):
            A=decode(n,code)
            for r in range(n):
                for u in range(n):
                    M=compiler(A,r,u)
                    for m in (2,3,4):
                        for c in range(1<<(m*(m-1)//2)):
                            B=decode(m,c);Z=seven(B,0);z=(Z[0],Z[1],Z[3],Z[4])
                            T=glue(A,r,B,0);W=seven(T,u)
                            need(mul(M,z)==(W[0],W[1],W[3],W[4]),'literal arbitrary-output-root compiler')
                            controls+=1
                    small.append(dict(rows=A,entry=r,exit=u,matrix=M))
    print('SMALL PORT CONTROLS',controls,flush=True)
    d=json.loads(gzip.decompress((HERE/'reroot_probe.json.gz').read_bytes()))
    w=next(x for x in d['words'] if x['context']['output_slope']=='6/41')
    A=tuple(w['context']['rows']);r=w['context']['root'];u=6;M=compiler(A,r,u)
    v=(0,1,2,1);image=mul(M,v)
    print('RENEWAL MATRIX',M,'DIRECTION',image,flush=True)
    need(image[0]==0 and image[2]==2*image[1] and image[3]==image[1] and image[1]!=0,'universal renewable direction')
    source=w['reroot'];z0=[tuple(z[j] for j in (0,1,3,4)) for z in source['profiles']]
    iterates=[];zs=z0
    for k in range(13):
        need(zs[0][0]==zs[1][0] and zs[0][3]-zs[0][1]==zs[1][3]-zs[1][1],'all computed renewal observations')
        iterates.append(dict(k=k,order=11+7*k,profiles=zs,next_H=zs[0][0]+zs[0][3]-zs[0][1]))
        zs=[mul(M,z) for z in zs]
    need(iterates[1]['profiles']==[tuple(z[j] for j in (0,1,3,4)) for s in w['states'] if s['root']==u and s['mode']=='E' for z in s['profiles']],'independent literal order18 renewal')
    save('port_compiler',dict(status='PROVED-IN-NOTE; FINITE-EXACT certificate',scope='Color-word proof applies to every rooted tournament context with >=2 vertices and every output root in that context; numerical controls all labelled contexts2..3 and all labelled attachments2..4, all entry/exit vertices. Selected order8 matrix compiled from all abstract words, independent order18 literal check.',physical_controls=controls,small_contexts=small,renewal=dict(context_rows=A,entry=r,exit=u,matrix=M,input_direction=v,output_direction=image,initial_rows=next(x['rows'] for x in d['pairs'] if x['label']=='2/3>1/3'),initial_root=source['root'],initial_profiles=z0,iterates=iterates)))


if __name__=='__main__':run()
