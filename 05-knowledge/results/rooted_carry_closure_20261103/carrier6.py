"""Six signed polynomial coordinates for a change of the edge between two old ports."""
from collections import Counter
from itertools import product
from math import factorial
import json,sympy as s
from observer import HERE,old,glue,junction,observer,root_profile,need,records

FREE=((),((0,-1),),((1,-1),),((-1,0),),((-1,1),),((1,0),))
AB=((0,1),)

def difference(B,A):
    return Counter({key:B.get(key,0)-A.get(key,0) for key in B.keys()|A.keys() if B.get(key,0)!=A.get(key,0)})

def decode(V):
    out=Counter(V)
    for (sk,m),v in V.items():
        if sk==FREE[-1]:out[AB,m]=-v
    return out

def encode(V):
    need(all(sk in FREE or sk==AB for sk,m in V),'seven-skeleton support')
    degrees={m for sk,m in V}
    need(all(V.get((AB,m),0)==-V.get((FREE[-1],m),0) for m in degrees),'antisymmetric double-end coefficient')
    return Counter({key:v for key,v in V.items() if key[0]!=AB and v})

def apply_signed(V,ops):
    ans=Counter()
    for (state,m),value in V.items():
        for j,block in enumerate(ops[glue.STATES.index(state)]):
            if m<j:continue
            for (target,power),weight in block.items():
                ans[target,m-j+power]+=value*factorial(m)//factorial(m-j)*weight
    return Counter({key:v for key,v in ans.items() if v})

def birth_kernel(sk,m,mask,polarity):
    parts=list(sk)+[(-1,-1)]*m;n=len(parts)
    enters=lambda end:polarity if end==-1 else (mask>>end&1)
    leaves=lambda start:1-polarity if start==-1 else 1-(mask>>start&1)
    incoming=[None]+[i for i,p in enumerate(parts) if enters(p[1])]
    outgoing=[None]+[i for i,p in enumerate(parts) if leaves(p[0])]
    out=Counter()
    for a,b in product(incoming,outgoing):
        if a is not None and a==b:continue
        nxt={};prev={}
        if a is not None:nxt[a]=n;prev[n]=a
        if b is not None:nxt[n]=b;prev[b]=n
        ends=parts+[(-1,-1)];result=[]
        for i in range(n+1):
            if i in prev:continue
            j=i
            while j in nxt:j=nxt[j]
            result.append((ends[i][0],ends[j][1]))
        out[glue.split(result)]+=1
    return out

def birth_operator(mask,polarity):
    ans=[]
    for state in glue.STATES:
        A=birth_kernel(state,0,mask,polarity)
        B=difference(birth_kernel(state,1,mask,polarity),glue.shift(A,1))
        need(not difference(birth_kernel(state,2,mask,polarity),glue.shift(A,2)+glue.shift(B,1,2)),'first derivative suffices for a uniform birth')
        ans.append((A,B))
    return ans

def birth(A,ports,mask,polarity):
    names={p:i for i,p in enumerate(ports)}
    incoming=sum(1<<v for v in range(len(A)) if ((mask>>names[v]&1) if v in names else polarity))
    return junction.ear(A,incoming)

def reduced_operator(ops):
    blocks=[]
    for state in FREE:
        row=[]
        for j in range(len(ops[0])):
            raw=Counter(ops[glue.STATES.index(state)][j])
            if state==FREE[-1]:raw=difference(raw,ops[glue.STATES.index(AB)][j])
            reduced=encode(raw)
            row.append([{'target':FREE.index(sk),'power':power,'coefficient':v} for (sk,power),v in sorted(reduced.items())])
        blocks.append(row)
    return blocks

def converse(V):
    out=Counter()
    for (sk,m),v in V.items():out[tuple(sorted((b,a) for a,b in sk)),m]+=v
    return out

def main():
    G=glue.operator();births={(mask,p):birth_operator(mask,p) for mask,p in product(range(4),range(2))}
    controls=0
    for n in range(2,5):
        for A in junction.all_tournaments(n):
            if not(A[0]>>1&1):continue
            B=list(A);B[0]^=2;B[1]^=1;B=tuple(B)
            V,W=observer(A,(0,1)),observer(B,(0,1))
            need(V==glue.literal(A) and W==glue.literal(B),'independent permutation-and-cut observer')
            D=difference(W,V);need(decode(encode(D))==D,'lossless six-coordinate difference')
            controls+=1
    w,u=old.chart(3,4);A=old.source(w,u);B=old.source(w,u,False);ports=(7,8)
    V,W=observer(A,ports),observer(B,ports);D=difference(W,V)
    operations=[]
    for key,ops in births.items():
        mask,p=key;VA=observer(birth(A,ports,mask,p),ports);WB=observer(birth(B,ports,mask,p),ports)
        need(apply_signed(V,ops)==VA and apply_signed(W,ops)==WB,'literal uniform port birth')
        need(apply_signed(D,ops)==difference(WB,VA),'signed birth transport')
        encode(apply_signed(D,ops));operations.append({'mask':mask,'polarity':p,'operator':reduced_operator(ops)})
    GA,GB=observer(glue.G(A),ports),observer(glue.G(B),ports)
    need(apply_signed(V,G)==GA and apply_signed(W,G)==GB,'literal eighteen-polynomial gluing')
    need(apply_signed(D,G)==difference(GB,GA),'signed gluing transport')
    full=(1<<len(A))-1
    CA=tuple(full^row^(1<<i) for i,row in enumerate(A));CB=tuple(full^row^(1<<i) for i,row in enumerate(B))
    need(converse(D)==difference(observer(CB,ports),observer(CA,ports)),'literal converse intertwining')
    need(converse(apply_signed(D,G))==apply_signed(converse(D),G),'Frobenius-free physical converse commutation')
    progression=[]
    for depth in range(6):
        encode(D);progression.append({'depth':depth,'root_profiles_plus':[root_profile(V,i) for i in range(2)],
                                      'root_profiles_minus':[root_profile(W,i) for i in range(2)],'difference':records(encode(D))})
        V=apply_signed(V,G);W=apply_signed(W,G);D=apply_signed(D,G)
    # Six independent derivatives on strong sources of one common order7.
    raw=[(2,0,3,7,15),(6,0,2,7,15),(10,0,3,6,15),(18,0,3,7,14),(2,4,1,7,15),(2,8,3,5,15)]
    z=s.Symbol('z');matrix=[];witnesses=[]
    for T in raw:
        U=list(T);U[0]^=2;U[1]^=1;U=tuple(U);T=glue.G(T);U=glue.G(U)
        need(junction.strong(T) and junction.strong(U),'actual strong lower-bound sources')
        D=encode(difference(observer(U,(0,1)),observer(T,(0,1))))
        matrix.append([sum(v*z**m for (sk,m),v in D.items() if sk==state) for state in FREE])
        witnesses.append({'plus':T,'minus':U,'ports':[0,1]})
    det=s.factor(s.Matrix(matrix).det());need(det!=0,'sharp pointwise Q(z)-linear dimension six')
    out={'status':'PROVED_SCOPED','coordinate_order':FREE,'small_flip_controls':controls,
         'gluing_operator':reduced_operator(G),'birth_operators':operations,'H497_gluing_progression':progression,
         'strong_rank_witnesses':witnesses,'strong_rank_determinant':str(det),'strong_rank_matrix':[[str(x) for x in row] for row in matrix]}
    (HERE/'carrier6.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: six-polynomial edge-response carrier; eight births, G and converse; strong order7 rank certificate')
if __name__=='__main__':main()
