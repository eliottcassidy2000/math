"""Exact halving-fibre observers, positive parity laws and sharp attachment depth."""
from pathlib import Path
from itertools import permutations,combinations
import sys,json,subprocess,hashlib
from sympy import Matrix
HERE=Path(__file__).resolve().parent
from completion import base,need

def birth(T,mask):return tuple(x+((mask>>i&1)<<len(T)) for i,x in enumerate(T))+(((1<<len(T))-1)^mask,)
def relabel(T,p):return tuple(sum(1<<j for j in range(len(T)) if T[p[i]]>>p[j]&1) for i in range(len(T)))
def hp(T):
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    for v in range(n):D[1<<v][v]=1
    for m in range(1,1<<n):
      for v,c in enumerate(D[m]):
        if not c:continue
        left=T[v]&~m
        while left:
            bit=left&-left;left-=bit;D[m|bit][bit.bit_length()-1]+=c
    return sum(D[-1])
def ear_values(T,exe):
    p=subprocess.run([exe],input=' '.join(map(str,[len(T),*T])),capture_output=True,text=True)
    need(p.returncode==0,p.stderr);return json.loads(p.stdout)
def pair_flip(T,bits,N):return relabel(T,[2*b+(j^((bits>>b)&1)) for b in range(N) for j in range(2)]+list(range(2*N,len(T))))
def wrapped_source(N):
    # C3[TT_(2N),1,1]; the two old anchor ports are always retained.
    m=2*N;a=m;b=m+1
    return tuple([sum(1<<j for j in range(i+1,m))+(1<<a) for i in range(m)]+[1<<b,(1<<m)-1])
def wrapped_read(k,bits):
    # k ordinal four-vertex blocks, one fixed cut-5 reader in each; then C3 wrap.
    old=wrapped_source(2*k);T=list(pair_flip(old,bits,2*k));n=len(T)
    for g in range(k):
        mask=0
        for i in range(len(T)):
            if i<4*k:beat=i<4*g or (4*g<=i<4*g+4 and ((i-4*g)%2==0))
            elif i<4*k+2:beat=(i==4*k+1)
            else:beat=True # preceding readers belong to earlier ordinal blocks
            if beat:mask|=1<<i
        T=list(birth(tuple(T),mask))
    return tuple(T)
def main():
    exe=sys.argv[1];T=base(2);f=ear_values(T,exe);rows=[]
    for mask in range(1,63):
        row=[]
        for x in range(8):
            p=[2*b+(j^((x>>b)&1)) for b in range(3) for j in range(2)]
            row.append(f[sum(1<<p[i] for i in range(6) if mask>>i&1)])
        rows.append(row)
    parity=[(-1)**x.bit_count() for x in range(8)]
    need(Matrix(rows).rank()==7 and all(sum(a*b for a,b in zip(row,parity))==0 for row in rows),'all one-birth responses have exactly the parity blind line')
    basis=[];masks=[]
    for mask,row in enumerate(rows,1):
        if Matrix(basis+[row]).rank()>len(basis):basis.append(row);masks.append(mask)
    pair=[];pair_reversed=[]
    for x in range(8):
        U=pair_flip(T,x,3);V=birth(U,1)
        pair.append(hp(birth(V,20)));pair_reversed.append(hp(birth(V,84)))
    need(sum(a*b for a,b in zip(pair,parity))==-4,'two legal births reveal parity')
    need(len(set(pair))==8,'one depth-two context identifies every point state')
    need(sum(a*b for a,b in zip(pair_reversed,parity))==-4,'new-new edge does not alter top interaction')
    determinant=int(Matrix(basis+[pair]).det());need(determinant!=0,'full eight-response rank')
    laws=dict(first_cut=1,second_cut=20,responses=pair,even_sum=sum(pair[x] for x in range(8) if x.bit_count()%2==0),odd_sum=sum(pair[x] for x in range(8) if x.bit_count()%2),
        reversed_new_edge_responses=pair_reversed,one_birth_basis_masks=masks,completion_determinant=determinant,rank_filtration=[1,7,8],blind_dimensions=[7,1,0],holonomy_kernel_orders=[8,1,1])
    local=[]
    for x in range(4):
        p=[2*b+(j^((x>>b)&1)) for b in range(2) for j in range(2)]
        TT=tuple(sum(1<<j for j in range(i+1,4)) for i in range(4));U=birth(relabel(TT,p),5);h=q=0
        for w in permutations(range(5)):
            good=[bool(U[a]>>b&1) for a,b in zip(w,w[1:])];h+=all(good)
            q+=sum(all(good[:i]+good[i+1:]) for i in range(4))
        local.append([h,2*h+q])
    need(local==[[3,48],[5,64],[5,64],[9,84]],'literal local two-bit reader')
    sharp=[]
    for k in (1,2):
        counts=[]
        for bits in range(1<<(2*k)):
            h=f2=1
            for g in range(k):a,b=local[(bits>>(2*g))&3];h*=a;f2*=b
            literal=hp(wrapped_read(k,bits));need(literal==h+f2,'literal strong wrapper vs product formula');counts.append(literal)
        alternating=sum((-1)**x.bit_count()*v for x,v in enumerate(counts))
        need(alternating==2**k+4**k,'sharp parity depth for every-budget formula')
        sharp.append(dict(k=k,counts=counts,alternating_sum=alternating,source_order=4*k+2,target_order=5*k+2))
    sys.path.insert(0,str(HERE.parent/'overnight_rank_return_20261029'))
    from atom_precision import hostile
    host=[];ts=[]
    for flip in (False,True):
        S=tuple(hostile(flip));ts.append(S);v=ear_values(S,exe);need(hp(S)==462915,'actual hostile baseline');need(hp(birth(S,256))==v[256],'independent named hostile birth');host.append(v)
    differences=[m for m in range(1,4095) if host[0][m]!=host[1][m]]
    need(len(differences)==3584 and next(m for m in sorted(differences,key=lambda m:(m.bit_count(),m)))==256,'complete one-birth hostile scan')
    report=dict(status='FINITE_EXACT_CONTROLS_FOR_SCOPED_ALL_SIZE_PROOFS',parity= laws,local_reader=local,sharp_controls=sharp,
        hostile=dict(rows=ts,baseline=462915,contexts=4094,differing=len(differences),cut=256,responses=[x[256] for x in host],
            full_table_sha256=[hashlib.sha256(json.dumps(v).encode()).hexdigest() for v in host],automorphism_orders='81/27 inherited, not recomputed'))
    (HERE/'halving_response.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS halving response',laws,'sharp controls',len(sharp),flush=True)
if __name__=='__main__':main()
