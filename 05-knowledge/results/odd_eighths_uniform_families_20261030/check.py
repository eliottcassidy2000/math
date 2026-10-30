"""Exact odd-eighth, conditional-energy and finite-tail controls."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd
from pathlib import Path
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
GATES=Counter()

def need(ok,label):
    GATES[label]+=1
    if not ok: raise RuntimeError(label)

def compositions(A):
    for mask in range(1 << (A-1)):
        cuts=[0]+[j for j in range(1,A) if mask>>(j-1)&1]+[A]
        yield tuple(b-a for a,b in zip(cuts,cuts[1:]))

def address(w,q,b=1):
    P,Q,B=1,1,0
    for a in w:P,Q,B=q*P,(1<<a)*Q,q*B+b*Q
    return -B*pow(P,-1,Q)%Q

def parity(n,q,h):
    bits=[]
    for _ in range(h):
        bit=n&1;bits.append(bit)
        n=(q*n+1)//2 if bit else n//2
    return tuple(bits)

def addroot(out,exponent,weight,h):
    Q=1<<h;d=Q//2;exponent%=Q
    out[exponent%d]+=weight if exponent<d else -weight

def zero(h):return [F(0)]*(1<<(h-1))

def rootproduct(a,b,h,conjugate=False):
    out=zero(h)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if x and y:addroot(out,i-j if conjugate else i+j,x*y,h)
    return out

def short_vectors(q,c,h,M):
    out=[zero(h) for _ in range(M+1)]
    for A in range(1,h+1):
        for w in compositions(A):
            if len(w)<=M:addroot(out[len(w)],c*address(w,q)*(1<<(h-A)),F(1,1<<A),h)
    return out

def suffix_vectors(q,c,h):
    # Independent forward parity labels on all residues, including even starts.
    out=[zero(h) for _ in range(h+1)]
    for s in range(1<<h):addroot(out[sum(parity(s,q,h))],c*s,F(1),h)
    return out

def energy_gram(q,c,h,M,T=None):
    small=short_vectors(q,c,h,M)
    out=[[rootproduct(small[m],small[n],h,True) for n in range(M+1)] for m in range(M+1)]
    for j in range(1,M+1):
        weight=F(1,3**j) if T is None else sum((F(comb(L-1,j-1),4**L) for L in range(j,T-h+1)),F(0))
        weight/=4**h
        A=suffix_vectors(q,c*pow(q,-j,1<<h)%(1<<h),h)
        for a in range(min(h,M-j)+1):
            for b in range(min(h,M-j)+1):
                v=rootproduct(A[a],A[b],h,True)
                out[j+a][j+b]=[x+weight*y for x,y in zip(out[j+a][j+b],v)]
    return out

def direct_gram(q,c,h,M,T):
    Q=1<<T;div=1<<(T-h);terms=[]
    for A in range(1,T+1):
        for r in range(1,1<<A,2):
            m=sum(parity(r,q,A))
            if m<=M:terms.append((m,r*(1<<(T-A)),F(1,1<<A)))
    out=[[zero(h) for _ in range(M+1)] for _ in range(M+1)]
    for m,a,x in terms:
        for n,b,y in terms:
            diff=(a-b)%Q
            if diff%div==0:addroot(out[m][n],c*(diff//div),x*y,h)
    return out

def collatz_checks():
    records=[]
    for q in (3,5,7,9):
        for h in range(1,5):
            for c in range(1,1<<h,2):
                M,T=6,8
                actual=direct_gram(q,c,h,M,T)
                predicted=energy_gram(q,c,h,M,T)
                need(actual==predicted,'conditional finite energy identity')
                infinite=energy_gram(q,c,h,M)
                need(all(infinite[m][n]==zero(h) for m in range(1,M+1) for n in range(1,M+1) if abs(m-n)>h),'bandwidth')
                for j in range(1,M-h+1):
                    expected=zero(h)
                    r=next(a for a in range(1<<h) if all(parity(a,q,h)))
                    addroot(expected,-c*pow(q,-j,1<<h)*r,F(1,4**h*3**j),h)
                    need(infinite[j][j+h]==expected and expected!=zero(h),'sharp bandwidth and source decoder')
                for size in (h,h+1,h+3):
                    gram=energy_gram(q,c,h,size)
                    actual=[sum((gram[m][n][i] for m in range(1,size+1) for n in range(1,size+1)),F(0)) for i in range(1<<(h-1))]
                    actual[0]-=F(1,4)
                    expected=zero(h)
                    for b in range(h):
                        A=suffix_vectors(q,c*pow(q,-(size-b),1<<h),h)
                        B=[sum((A[l][i] for l in range(b+1)),F(0)) for i in range(1<<(h-1))]
                        norm=rootproduct(B,B,h,True)
                        expected=[x+F(3**b,4**h*3**size)*y for x,y in zip(expected,norm)]
                    need(actual==expected,'exact periodic conditional remainder')
                records.append({'q':q,'h':h,'c':c})
    for m in range(1,13):
        g=energy_gram(3,1,1,13)
        need(g[m][m]==[F(1,3**m)],'inherited odd Haar diagonal')
        need(g[m][m+1]==[-F(1,4*3**m)],'inherited odd Haar adjacent')
    polys=[]
    for c in (1,3,5,7):
        polys.append({'c':c,'A_ell':[list(map(str,a)) for a in suffix_vectors(3,c,3)]})
    remainder=[]
    for c in (1,3,5,7):
        for size in (3,4):
            g=energy_gram(3,c,3,size)
            coeff=[sum((g[m][n][i] for m in range(1,size+1) for n in range(1,size+1)),F(0)) for i in range(4)]
            coeff[0]-=F(1,4)
            sign=(-1)**size*(1 if c in (1,7) else -1)
            need([3**size*x for x in coeff]==[F(1,4),F(3*sign,64),F(0),F(-3*sign,64)],'explicit eighth-root remainder')
            remainder.append({'c':c,'M_parity':size%2,'coefficient_of_3_to_minus_M':list(map(str,[3**size*x for x in coeff]))})
    return {'finite_Gram_universes':records,'q3_h3_suffix_polynomials':polys,'q3_h3_remainder':remainder}

def transport_checks():
    rows=[]
    for q,p in product((3,5,7,9,11,13),repeat=2):
        for D in range(1,8):
            Q=1<<D;mapping={0:0}
            for A in range(1,D+1):
                for w in compositions(A):
                    mapping[address(w,q)*(1<<(D-A))]=address(w,p)*(1<<(D-A))
            need(set(mapping)==set(range(Q)) and len(set(mapping.values()))==Q,'full character permutation')
            u=q*pow(p,-1,Q)%Q
            linear=all(mapping[x]==u*x%Q for x in range(Q))
            need(linear==((p-q)%(1<<(D-1))==0),'all-depth multiplicativity classification')
            if D<=5:
                direct=all(mapping[(x+y)%Q]==(mapping[x]+mapping[y])%Q for x in range(Q) for y in range(Q))
                need(direct==linear,'independent additive-pair test')
        if p!=q:
            v=0;diff=abs(p-q)
            while diff%2==0:v+=1;diff//=2
            D=v+2;Q=1<<D;k=(1+2*pow(q,-1,Q))%Q
            w=(1,D-1)
            need(address(w,q)==k*address((D,),q)%Q,'source product relation')
            defect=(address(w,p)-k*address((D,),p))%Q
            need(defect==Q//2,'first failed product is antipodal')
            rows.append({'q':q,'p':p,'first_failed_depth':D,'defect':f'1/2'})
    for q,p,b,c,D in product((3,5,7),(3,5,7),(-1,1,3),(-1,1,5),range(1,6)):
        Q=1<<D;u=c*q*pow(b*p,-1,Q)%Q
        good=all(address(w,p,c)*(1<<(D-A))%Q==u*address(w,q,b)*(1<<(D-A))%Q
                 for A in range(1,D+1) for w in compositions(A))
        need(good==((p-q)%(1<<(D-1))==0),'odd offset transport classification')
    return {'first_failures':rows,'q3_cost3':[{'word':w,'address':str(F(address(w,3),8))} for w in compositions(3)]}

def load_lrc():
    path=ROOT/'05-knowledge/results/session114_continuation_symmetry_20261030/lrc/check_resonant_tail.py'
    spec=importlib.util.spec_from_file_location('old_lrc',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def lrc_checks():
    lrc=load_lrc();tables=[];banks=0;decisions=0
    for q in range(3,13):
        C=tuple(range(1,q))+(q+1,);s=next(a for a in range(2,q) if gcd(a,q)==1)
        for T in ((1,),(1,2),(2,3),(2,5),(3,4,7),(1,4,9),(2,4),(3,6,9)):
            V=C+tuple(q*m for m in T)
            for lam in (F(1,7),F(1,2),F(6,7),F(1)):
                r=F(1,q+2)+lam*(F(1,q)-F(1,q+2))
                tails=lrc.survivor(T,r)
                tau=tails[0][0] if tails else None
                want=tau is not None and q*r+s*tau<=1
                need(bool(lrc.survivor(V,r))==want,'all-tail first-arrival equivalence')
                decisions+=1
        for m in range(2,5):
            if 2*m<=s:continue
            D=q*m+s;a=pow(s,-1,q);K=a*m+(a*s-1)//q;r=F(m,D)
            base=C+(q*m,)
            need(gcd(K,D)==1,'primitive maximum clock')
            for v in range(1,2*D+1):
                good=m<=v*K%D<=D-m
                points=lrc.survivor(base+(v,),r)
                need(bool(points)==good,'complete arbitrary-tail equality classification')
                if good:
                    need(points==[(F(K,D),F(K,D)),(1-F(K,D),1-F(K,D))] or points==sorted([(F(K,D),F(K,D)),(1-F(K,D),1-F(K,D))]),'all maximizing points')
                decisions+=1
            good_v=[v for v in range(1,3*D+1) if m<=v*K%D<=D-m and v not in base]
            samples=(tuple(good_v[:4]),tuple(good_v[-4:]),(s+D,q*m+D,s+2*D,q*m+2*D))
            for Z in samples:
                V=base+Z;pmax=max([s]+[v for v in Z if v%D==s]);nmax=max([q*m]+[v for v in Z if v%D==q*m])
                B=2*(F(1,pmax)+F(1,nmax))
                eps=F(1,1000*max(V)**2)
                intervals=lrc.survivor(V,r-eps)
                measured=sum((b-a for a,b in intervals),F(0))/eps
                need(measured==B,'exact physical terminal mass slope')
                banks+=1
            if q in (8,12):tables.append({'q':q,'m':m,'M':str(r),'t':str(F(K,D)),'modulus':D,'forbidden_residues':[v for v in range(D) if not m<=v*K%D<=D-m]})
    # Direct global maximization is a separate path, at bounded heights.
    for q,m in ((3,2),(5,2),(8,2),(12,3)):
        C=tuple(range(1,q))+(q+1,);s=next(a for a in range(2,q) if gcd(a,q)==1);D=q*m+s
        for Z in ((s+D,),(q*m+D,),(q,)):
            actual,points,visits=lrc.maximum(C+(q*m,)+Z)
            K=pow(s,-1,q)*m+(pow(s,-1,q)*s-1)//q
            good=all(m<=v*K%D<=D-m for v in Z)
            need((actual==F(m,D))==good,'independent global maximum control')
    return {'total_tail_decisions':decisions,'multi_tail_slope_banks':banks,'examples':tables,'q8_core_owners':[[1,9,7],[3,3,5],[5,5,3],[7,7,9]]}

def inherited_controls():
    path=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec=importlib.util.spec_from_file_location('old_tournaments',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.hostile_checks()
    counts=Counter()
    edges=list(combinations(range(4),2))
    for mask in range(64):
        out=[set() for _ in range(4)]
        for i,(a,b) in enumerate(edges):
            if mask>>i&1:out[a].add(b)
            else:out[b].add(a)
        degrees=sorted(map(len,out))
        if degrees==[0,1,2,3]:label='transitive'
        elif degrees==[1,1,1,3]:label='source_plus_cycle'
        elif degrees==[0,2,2,2]:label='cycle_plus_sink'
        else:
            for start in range(4):
                reach={start}
                for _ in range(4):reach|={v for u in tuple(reach) for v in out[u]}
                need(len(reach)==4,'literal strong four-vertex control')
            label='strong'
        counts[label]+=1
    need(dict(counts)=={'transitive':24,'source_plus_cycle':8,'cycle_plus_sink':8,'strong':24},'inherited tournament eighth counts')
    return {'control_81_27':result,'four_vertex_counts':dict(counts),'inherited_gates':7}

def phase_bridge_checks():
    lrc=load_lrc();rows=[]
    for m in (2,3,6,10):
        V=tuple(range(1,7))+(8,)+tuple(7*m*j for j in range(1,8))
        points=lrc.survivor(V,F(1,8))
        need(bool(points) and all(a==b for a,b in points),'incoming harmonic literal maximizers')
        phases=[int(8*((7*m*a)%1)) for a,b in points]
        need(all(c in (1,3,5,7) for c in phases),'literal maximizing tail phase')
        for M in (3,4):
            K=zero(3)
            for c in phases:
                g=energy_gram(3,c,3,M)
                value=[sum((g[i][j][k] for i in range(1,M+1) for j in range(1,M+1)),F(0)) for k in range(4)]
                value[0]-=F(1,4)
                K=[x+3**M*y/len(phases) for x,y in zip(K,value)]
            delta=K[:];delta[0]-=F(1,4)
            correction=rootproduct(delta,[F(0),F(1),F(0),F(-1)],3)
            weight=[F(1024*(-1)**M,315)*x for x in correction];weight[0]+=F(88,105)
            B=[F(len(phases),7*m)*x for x in weight]
            eps=F(1,1000*max(V)**2)
            measured=sum((b-a for a,b in lrc.survivor(V,F(1,8)-eps)),F(0))/eps
            need(B==[measured,F(0),F(0),F(0)],'literal LRC mass from conditional Collatz energy')
            rows.append({'m':m,'M':M,'maximizer_count':len(phases),'mass_slope':str(measured),'scaled_average_energy':list(map(str,K))})
    return rows

def main():
    result={'collatz':collatz_checks(),'transport':transport_checks(),'lrc':lrc_checks()}
    result['inherited_controls']=inherited_controls()
    result['incoming_phase_bridge']=phase_bridge_checks()
    result['gates']=dict(GATES);result['total_gates']=sum(GATES.values())
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
