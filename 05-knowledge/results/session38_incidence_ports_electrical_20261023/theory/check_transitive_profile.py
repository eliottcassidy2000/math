#!/usr/bin/env python3
"""Independent n15 exact-E/transitive-family obstruction and propagation.
No research producer imports. All subsets, root neighborhoods, and motifs
are counted exactly. Optimized Python retains every explicit check.
"""
import argparse,hashlib,itertools as it,json
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent
SETS=[[14,13,12,11,5,6,8],[14,13,12,4,10,9,8]]
CHECKS=0

def need(x,msg):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(msg)

def bits(x):
    while x:
        b=x&-x;yield b.bit_length()-1;x^=b

def circulant(n,S):return [sum(1<<((i+d)%n) for d in S) for i in range(n)]
def induce(A,vertices):return [sum(((A[v]>>w)&1)<<j for j,w in enumerate(vertices)) for v in vertices]

def F(A):
    n=len(A);ok=bytearray(1<<n);ok[0]=1;coeff=[0]*(n+1);coeff[0]=1
    for S in range(1,1<<n):
        k=S.bit_count()
        for v in bits(S):
            rest=S^(1<<v)
            if (A[v]&rest)==rest:
                ok[S]=ok[rest];break
        coeff[k]+=ok[S]
    while len(coeff)>1 and coeff[-1]==0:coeff.pop()
    return coeff

def rootF(A):return sorted(F(induce(A,list(bits(row)))) for row in A)

def ham_subsets(A):
    n=len(A);pred=[sum(((A[u]>>v)&1)<<u for u in range(n)) for v in range(n)];dp=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n)
    for v in range(n):dp[1<<v][v]=1;H[1<<v]=1
    for S in range(1,1<<n):
        if S&(S-1)==0:continue
        for v in bits(S):dp[S][v]=sum(dp[S^(1<<v)][u] for u in bits(pred[v]&S))
        H[S]=sum(dp[S])
    return H

def exactE(A):
    n=len(A);cv=[0]*n;ce=[[0]*n for _ in range(n)]
    for a,b,c in it.combinations(range(n),3):
        if ((A[a]>>b&1) and (A[b]>>c&1) and (A[c]>>a&1)) or ((A[b]>>a&1) and (A[c]>>b&1) and (A[a]>>c&1)):
            for v in (a,b,c):cv[v]+=1
            for u,v in [(a,b),(a,c),(b,c)]:ce[u][v]+=1;ce[v][u]+=1
    return sorted(sorted((cv[w],ce[v][w]) for w in bits(A[v])) for v in range(n))

def canon(A):
    n=len(A);deg=[a.bit_count() for a in A];blocks=[[v for v in range(n) if deg[v]==d] for d in sorted(set(deg))];pairs=list(it.combinations(range(n),2));best=None
    for choices in it.product(*(it.permutations(b) for b in blocks)):
        p=sum(choices,());code=sum(((A[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(pairs))
        if best is None or code<best:best=code
    return (tuple(sorted(deg)),best)

def motifs(A,k):return Counter(canon(induce(A,s)) for s in it.combinations(range(len(A)),k))
def motif_rows(c):return [{'scores':list(key[0]),'canonical_upper_bits':key[1],'count':v} for key,v in sorted(c.items())]

def polyadd(*xs):
    out=[0]*max(map(len,xs))
    for p in xs:
        for i,v in enumerate(p):out[i]+=v
    while len(out)>1 and out[-1]==0:out.pop()
    return out

def neg(x):return [-v for v in x]
def polymul(x,y):
    z=[0]*(len(x)+len(y)-1)
    for i,a in enumerate(x):
        for j,b in enumerate(y):z[i+j]+=a*b
    return z

def transitive(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def ordinal(A,B):
    n=len(A);m=len(B);return [a|(((1<<m)-1)<<n) for a in A]+[b<<n for b in B]
def source(A):return ordinal([0],A)
def lift(A,r):
    n=len(A);a=n;P=transitive(r);out=[v|(1<<a) for v in A]+[((1<<r)-1)<<(n+1)];out.extend((P[i]<<(n+1))|((1<<n)-1) for i in range(r));return out

def lift_poly(f,r):
    fa=[1,1];fp=[1]
    for _ in range(r):fp=polymul(fp,[1,1])
    return polyadd(polymul(f,fa),polymul(fa,fp),polymul(fp,f),neg(f),neg(fa),neg(fp),[1])

def lift_root_profiles(f,profiles,r):
    result=[polymul(p,[1,1]) for p in profiles];fp=[1]
    for _ in range(r):fp=polymul(fp,[1,1])
    result.append(fp);fp=[1]
    for j in range(r):result.append(polymul(fp,f));fp=polymul(fp,[1,1])
    return sorted(result)

def main(out):
    cases=[];rows=[circulant(15,S) for S in SETS]
    for number,A in enumerate(rows):
        need(all((A[i]>>i&1)==0 for i in range(15)),'zero diagonal');need(all((A[i]>>j&1)+(A[j]>>i&1)==1 for i,j in it.combinations(range(15),2)),'tournament')
        f=F(A);profiles=rootF(A);need(f==[1,15,105,315,375,150,15],'entire F15');need(profiles==[[1,7,21,25,10,1]]*15,'every entire outneighborhood F')
        for k in range(1,len(f)):need(f[k]==sum(p[k-1] if k-1<len(p) else 0 for p in profiles),'unique transitive source recurrence')
        E=exactE(A);need(E==[[(28,3)]*3+[(28,4)]*2+[(28,5),(28,6)]]*15,'exactE')
        neighborhoods=[]
        for v in range(15):
            U=induce(A,list(bits(A[v])));hs=ham_subsets(U);hist=Counter((U[a]&U[b]).bit_count() for a,b in it.combinations(range(7),2));B=sum(k*(k-1)//2*w for k,w in hist.items());m4=motifs(U,4);m5=motifs(U,5)
            need(hist==Counter({0:5,1:8,2:7,3:1}),'common-out histogram');need(B==10==F(U)[4],'local B/transitive4 count');need(hs[-1]==[101,89][number],'neighborhood H');regular5=sum(all((U[a]&sum(1<<b for b in S)).bit_count()==2 for a in S) for S in it.combinations(range(7),5));need(regular5==number,'direct regular-five neighborhood repair')
            neighborhoods.append({'root':v,'H':hs[-1],'common_out_histogram':sorted(hist.items()),'B':B,'regular_five_count':regular5,'motifs4':motif_rows(m4),'motifs5':motif_rows(m5)})
        need(all(x['motifs4']==neighborhoods[0]['motifs4'] and x['motifs5']==neighborhoods[0]['motifs5'] for x in neighborhoods),'vertex translation controls')
        hs=ham_subsets(A);F2=sum(hs[s]*hs[((1<<15)-1)^s] for s in range(1,(1<<15)-1));need(hs[-1]==[196264905,196277235][number],'global H');need(F2==[4857720210,4857619470][number],'global ordered2pathcovers')
        source_regular5=0
        for S in it.combinations(range(15),6):
            for v in S:
                rest=[w for w in S if w!=v];mask=sum(1<<w for w in rest)
                if A[v]&mask==mask and all((A[w]&mask).bit_count()==2 for w in rest):source_regular5+=1
        need(source_regular5==15*number,'literal source-plus-regular-five count')
        cases.append({'source_plus_regular_five_count':source_regular5,'steps':SETS[number],'rows':A,'F':f,'root_outneighbor_F':profiles,'E':E,'H':hs[-1],'ordered_two_path_covers':F2,'neighborhoods':neighborhoods,'global_motifs5':motif_rows(motifs(A,5))})
    need(cases[0]['neighborhoods'][0]['motifs4']==cases[1]['neighborhoods'][0]['motifs4'],'all four motifs agree locally')
    need(cases[0]['E']==cases[1]['E'] and cases[0]['F']==cases[1]['F'] and cases[0]['root_outneighbor_F']==cases[1]['root_outneighbor_F'],'simultaneous whole signature collision')
    neighborhood_differences=[];global_differences=[]
    for k,target in [(5,neighborhood_differences)]:
        a=motifs(induce(rows[0],list(bits(rows[0][0]))),k);b=motifs(induce(rows[1],list(bits(rows[1][0]))),k)
        for key in sorted(a.keys()|b.keys()):
            if a[key]!=b[key]:target.append({'scores':list(key[0]),'canonical_upper_bits':key[1],'counts':[a[key],b[key]]})
    a=motifs(rows[0],5);b=motifs(rows[1],5)
    for key in sorted(a.keys()|b.keys()):
        if a[key]!=b[key]:global_differences.append({'scores':list(key[0]),'canonical_upper_bits':key[1],'counts':[a[key],b[key]]})
    controls=[]
    for n in range(1,5):
        pairs=list(it.combinations(range(n),2))
        for code in range(1<<len(pairs)):
            A=[0]*n
            for k,(u,v) in enumerate(pairs):
                a,b=(u,v) if code>>k&1 else (v,u);A[a]|=1<<b
            f=F(A);rf=rootF(A);S=source(A);need(F(S)==polymul([1,1],f),'source F propagation');need(rootF(S)==sorted(rf+[f]),'source root profiles')
            for r in range(1,4):
                L=lift(A,r);need(F(L)==lift_poly(f,r),'cyclic F formula');need(rootF(L)==lift_root_profiles(f,rf,r),'cyclic root profiles')
        controls.append({'n':n,'all_labelled_tournaments':1<<len(pairs),'source_padding':True,'cyclic_r':[1,2,3]})
    seed_extension=[]
    for A in rows:
        f=F(A);rf=rootF(A);S=source(A);L=lift(A,1);need(F(S)==polymul([1,1],f),'literal n16 source control');need(rootF(S)==sorted(rf+[f]),'literal n16 root source control');need(F(L)==lift_poly(f,1),'literal n17 cyclic control');need(rootF(L)==lift_root_profiles(f,rf,1),'literal n17 rooted cyclic control');seed_extension.append({'source_n16_F':F(S),'cyclic_n17_F':F(L),'cyclic_n17_root_profiles':rootF(L)})
    result={'status':'FINITE-EXACT all-subset and all-neighborhood verification; propagation proved separately','incoming_commit':'f2cac6db1bb8ef780b9211d1cf3f7adc109c5fca','cases':cases,'local_nontransitive_five_motif_differences':neighborhood_differences,'global_nontransitive_five_motif_differences':global_differences,'small_exhaustive_propagation_controls':controls,'literal_seed_extension_controls':seed_extension,'checks':CHECKS,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'checks':CHECKS,'F':cases[0]['F'],'root_F':cases[0]['root_outneighbor_F'][0],'global_H':[x['H'] for x in cases],'neighborhood_H':[x['neighborhoods'][0]['H'] for x in cases],'local_five_differences':neighborhood_differences,'global_five_differences':global_differences},indent=2));print('ALL CHECKS PASSED')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=HERE/'transitive_profile_receipt.json');a=p.parse_args();main(a.out)
