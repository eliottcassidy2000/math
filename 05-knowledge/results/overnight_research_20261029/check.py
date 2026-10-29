"""Exact producer checks; no third-party packages or disabled assertions."""
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
from functools import lru_cache
from math import gcd
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
GATES = 0
def require(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)

M = [[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]

def inverse(a):
    n=len(a)
    z=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for i in range(n):
        k=next(k for k in range(i,n) if z[k][i])
        z[i],z[k]=z[k],z[i]
        p=z[i][i]; z[i]=[x/p for x in z[i]]
        for j in range(n):
            if i!=j:
                p=z[j][i];z[j]=[x-p*y for x,y in zip(z[j],z[i])]
    return [row[n:] for row in z]

def pf(a):
    if not a:return 1
    return sum((-1)**(j+1)*a[0][j]*pf([[a[r][s] for s in range(1,len(a)) if s!=j]
                for r in range(1,len(a)) if r!=j]) for j in range(1,len(a)))

def sub(a, ids): return [[a[i][j] for j in ids] for i in ids]
def block(m,n,b):
    return [m[i]+b[i] for i in range(len(m))]+[
        [-b[j][i] for j in range(len(m))]+n[i] for i in range(len(n))]

def rank101(a):
    a=[[x%101 for x in row] for row in a];r=0
    for j in range(len(a[0])):
        k=next((k for k in range(r,len(a)) if a[k][j]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]
        u=pow(a[r][j],-1,101);a[r]=[(v*u)%101 for v in a[r]]
        for k in range(r+1,len(a)):
            u=a[k][j];a[k]=[(x-u*y)%101 for x,y in zip(a[k],a[r])]
        r+=1
    return r

def paired_invariants(a, total_pf=None):
    previous=1;out=[];n=len(a)
    for k in range(1,n//2+1):
        if k==n//2 and total_pf is not None:g=abs(total_pf)
        else:
            g=0
            for ids in combinations(range(n),2*k):
                g=gcd(g,pf(sub(a,ids)))
                if g==previous:break
        require(g%previous==0,"Pfaffian ideal divisibility")
        out.append(g//previous);previous=g
    require(all(out[i+1]%out[i]==0 for i in range(len(out)-1)),"invariant chain")
    return tuple(x for x in out if x!=1)

def tournament_census():
    ca=[[int(3*x) for x in row] for row in inverse(M)]
    rank_counts=Counter();dist=Counter();clock=Counter();reps={}
    perm=[0,2,3,1];pfs=set()
    for bits in range(1<<16):
        b=[[1 if (bits>>(4*i+j))&1 else -1 for j in range(4)] for i in range(4)]
        r=rank101(b);rank_counts[r]+=1
        if r>2:continue
        cb=[[sum(ca[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
        s=[[3*M[i][j]+sum(b[k][i]*cb[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
        p=pf(s)
        require(p%3==0,"integral Schur Pfaffian")
        p//=(-3);require(p%2!=0,"tournament odd Pfaffian")
        k=block(M,M,b); inv=paired_invariants(k,p)
        key=(r,inv);dist[str(key)]+=1;pfs.add(abs(p))
        if key not in reps:
            require(pf(k)==p,"direct Pfaffian control")
            reps[key]={"bits":bits,"B":b,"pfaffian":p,"paired_invariants":inv}
        if all(b[perm[i]][perm[j]]==b[i][j] for i in range(4) for j in range(4)):
            clock[str(key)]+=1
    require(sum(rank_counts.values())==65536,"entire labelled cut universe")
    return {"all_cut_ranks":dict(sorted(rank_counts.items())),"rank_le_two_groups":dict(sorted(dist.items())),
        "clock_preserving_groups":dict(sorted(clock.items())),"absolute_pfaffians":sorted(pfs),
        "representatives":{str(k):v for k,v in reps.items()}}

def integer_ports():
    count=0;special={}
    for d,e in product(range(1,10),repeat=2):
        for s in range(1,gcd(d,e)+1):
            if d%s or e%s:continue
            for n in range(1,61):
                b=[[s,0],[0,s*((d//s)*(e//s)-n)]]
                k=block([[0,d],[-d,0]],[[0,e],[-e,0]],b)
                require(pf(k)==s*s*n,"constructive two-port Pfaffian")
                require(paired_invariants(k)==tuple(x for x in (s,s*n) if x!=1),"two-port full group")
                count+=1
                if d==e==3 and s==1 and n in (25,27,57):special[n]={"B":b,"pfaffian":pf(k)}
    # Hostile: same order, different group, using the same starting blocks.
    require(paired_invariants(block([[0,3],[-3,0]],[[0,3],[-3,0]],[[0,0],[0,0]]))==(3,3),"H3 plus H3")
    require(paired_invariants(block([[0,3],[-3,0]],[[0,3],[-3,0]],[[1,0],[0,0]]))==(9,),"H9 hostile")
    return {"constructed_groups":count,"special":special}

def odd_step(n):
    z=3*n+1;a=(z&-z).bit_length()-1
    return z>>a,a

def parser(n,q0=1):
    source=n;q=q0;word=[];D=0;C=0;pow3=1;maxq=q;high=0
    while q>0:
        old=n;n,a=odd_step(n);oldq=q;q+=3*(a%3==0)-1
        if old>=9:
            if q-oldq==2: require(225*64*n<=224*25*old,"branch contraction")
            else: require(225*5*n<=224*8*old,"leaf contraction")
            high+=1
        C=3*C+(1<<D);D+=a;pow3*=3;word.append(a);maxq=max(maxq,q)
        require(len(word)<100000,"diagnostic guard, not the proof")
    require(pow3*source+C==n*(1<<D),"source endpoint affine identity")
    residue=((1<<D)-C)*pow(pow3,-1,1<<(D+1))%(1<<(D+1))
    require(source%(1<<(D+1))==residue,"exact valuation cylinder")
    return {"source":source,"endpoint":n,"steps":len(word),"clock":D,"carry":C,
            "queue_start":q0,"max_queue":maxq,"high_steps":high,"word":word}

def controller_census():
    maxima={"steps":(0,0),"max_queue":(0,0),"clock":(0,0)};growth=0;special={}
    for n in range(1,200000,2):
        r=parser(n)
        growth+=r["endpoint"]>=n
        for key in maxima:
            if r[key]>maxima[key][0]:maxima[key]=(r[key],n)
        if n in (1,3,5,7,9,11,25,27,57,343,1001):special[n]=r
    for n,q in product((1,3,5,7,9,11,25,27,57),range(1,41)):parser(n,q)
    require(special[3]["endpoint"]==5,"termination not descent hostile")
    # 3n-1: actual positive cycle 5 -> 7 -> 5, leaves both times.
    require((3*5-1)//2==7 and (3*7-1)//4==5,"minus map cycle")
    return {"positive_odd_sources":100000,"extra_source_queue_pairs":360,"nondecreasing_returns":growth,
            "maxima":maxima,"special":special}

def symplectic15():
    def b(x,y):return ((x&3)&(y>>2)).bit_count()%2 ^ (((x>>2)&(y&3)).bit_count()%2)
    lag={tuple(sorted((x,y,x^y))) for x in range(1,16) for y in range(x+1,16) if not b(x,y)}
    require(len(lag)==15,"15 Lagrangians")
    degrees=Counter(sum(bool(set(a)&set(c)) for c in lag if c!=a) for a in lag)
    require(degrees=={6:15},"doily incidence degrees")
    require(all(sum(x in a for a in lag)==3 for x in range(1,16)),"three lines through point")
    quadratic=[a for a in lag if all(((x&3)&(x>>2)).bit_count()%2==0 for x in a)]
    require(len(quadratic)==6,"quadratic refinement retains six")
    return {"lagrangians":[list(a) for a in sorted(lag)],"intersection_degrees":dict(degrees),
            "split_quadratic_lagrangians":len(quadratic),"symplectic_group_order":720}

def main():
    out={"integer_two_ports":integer_ports(),"tournament_two_ports":tournament_census(),
         "ternary_parser":controller_census(),"f2_geometry":symplectic15()}
    out["gates"]=GATES
    (ROOT/"check.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"gates":GATES,"integer_groups":out["integer_two_ports"]["constructed_groups"],
        "tournament":{k:v for k,v in out["tournament_two_ports"].items() if k!="representatives"},
        "parser":{k:v for k,v in out["ternary_parser"].items() if k!="special"}},indent=2))

if __name__=="__main__":main()
