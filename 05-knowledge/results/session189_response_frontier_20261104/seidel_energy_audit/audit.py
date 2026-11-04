#!/usr/bin/env python3
"""Independent exact rational/word checks. No imports from the producer packet."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from math import comb
import json, subprocess, hashlib, sys
P=Path(__file__).resolve().parent

def require(x,msg):
    if not x: raise RuntimeError(msg)
def solve(mat,rhs):
    n=len(mat); a=[[F(v) for v in row]+[F(rhs[i])] for i,row in enumerate(mat)]
    for col in range(n):
        piv=next((row for row in range(col,n) if a[row][col]),None)
        require(piv is not None,'singular matrix');a[col],a[piv]=a[piv],a[col]
        fac=a[col][col];a[col]=[v/fac for v in a[col]]
        for row in range(n):
            if row!=col:
                fac=a[row][col]
                if fac:a[row]=[u-fac*v for u,v in zip(a[row],a[col])]
    return [row[-1] for row in a]
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def matvec(A,x):return [dot(row,x) for row in A]
def trans(A):return list(map(list,zip(*A)))
def inv(A):return trans([solve(A,[int(i==j) for i in range(len(A))]) for j in range(len(A))])
def opposite(q):
    n=len(q);return tuple(((1<<n)-1)^q[i]^(1<<i) for i in range(n))
def graph(n,mask):
    q=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if mask>>k&1:q[i]|=1<<j
        else:q[j]|=1<<i
    return tuple(q)
def mod(x,M):
    x=F(x);require(x.denominator%2==1,'nonintegral 2-adic value');return x.numerator*pow(x.denominator,-1,M)%M
def val(x):
    x=F(x);require(x!=0,'valuation zero');return (abs(x.numerator)&-abs(x.numerator)).bit_length()-(x.denominator&-x.denominator).bit_length()
def basis(m):return [0]+[1<<i for i in range(m)]+[(1<<i)|(1<<j) for i,j in combinations(range(m),2)]
def flin(ids,r,sgn):
    p={0:F(1)}
    for i,x in zip(ids,r):p[1<<i]=sgn*4*x+4*x*x
    for j in range(len(ids)):
        for k in range(j):p[(1<<ids[j])|(1<<ids[k])]=8*r[j]*r[k]
    return p

def banks(q,check_inverse=True):
    m=len(q);bs=basis(m);even={(a,c):F(0) for a in bs for c in bs};odd=even.copy();energy=F(0);rs={};inv_checks=0
    for mask in range(1<<m):
        ids=[i for i in range(m) if mask>>i&1];n=len(ids)
        if not n:even[0,0]+=F(1,2**m);continue
        B=[[int(q[i]>>j&1) for j in ids] for i in ids]
        S=[[B[i][j]-B[j][i] for j in range(n)] for i in range(n)]
        M=[[int(i==j)+2*B[i][j] for j in range(n)] for i in range(n)]
        r=solve(M,[1]*n) if n%2 else solve(S,[1]*n)
        require(all(x.denominator%2 for x in r),'r odd denominators');require(sum(r)==(n%2),'r total');require(matvec(S,r)==([F(0)]*n if n%2 else [F(1)]*n),'S r')
        rs[mask]=(ids,r)
        if not n%2:energy+=sum(x*x for x in r)
        if check_inverse:
            for t in [F(0),F(1),F(-1,3)]:
                v=1+2*t;K=inv([[M[i][j]+2*t for j in range(n)] for i in range(n)]);left=matvec(K,[1]*n);right=matvec(trans(K),[1]*n)
                require(left==([x/v for x in r] if n%2 else r),'left inverse');require(right==([x/v for x in r] if n%2 else [-x for x in r]),'right inverse')
                if n%2:require(all(K[i][j]+K[j][i]==2*r[i]*r[j]/v for i in range(n) for j in range(n)),'odd sym inverse')
                else:
                    Sinv=inv(S);require(all(K[i][j]==Sinv[i][j]+v*r[i]*r[j] for i in range(n) for j in range(n)),'even inverse')
                # Direct segment formula with independently inverted K, including both port orders.
                for a,c in [(0,0),((1<<n)-1,1),(1,1<<(n-1))]:
                    av=[F(a>>i&1) for i in range(n)];cv=[F(c>>i&1) for i in range(n)]
                    La=1-2*dot(right,[x+t for x in av]);Lc=1-2*dot(right,[x+t for x in cv]);Ra=1-2*dot([1+t-x for x in av],left);Rc=1-2*dot([1+t-x for x in cv],left)
                    xy=1+t-2*dot([1+t-x for x in av],matvec(K,[x+t for x in cv]));yx=t-2*dot([1+t-x for x in cv],matvec(K,[x+t for x in av]))
                    got=La*xy*Rc+Lc*yx*Ra;sgn=-1 if n%2 else 1;expected=((1-2*dot(r,av))**2*(1-2*dot(r,cv))**2/(-v**3) if n%2 else v*(1+2*dot(r,av))**2*(1+2*dot(r,cv))**2)
                    require(got==expected,'walk resolvent constant');inv_checks+=1
        f=flin(ids,r,-1 if n%2 else 1);bank=odd if n%2 else even
        for a,u in f.items():
            for c,v in f.items():bank[a,c]+=u*v/F(2**m)
    energy*=F(64,2**m)
    return even,odd,energy,rs,inv_checks

def evalbank(bank,a,c):return sum(v for (u,w),v in bank.items() if u&a==u and w&c==w)
def encode(x):return str(x.numerator) if x.denominator==1 else str(x)
def main():
    # Every labeled quotient through order four, plus selected larger exact controls.
    qs=[]
    for m in range(2,5):qs.extend(graph(m,z) for z in range(1<<(m*(m-1)//2)))
    qs += [tuple(sum(1<<j for j in range(i+1,m)) for i in range(m)) for m in range(5,9)]
    qs += [(6,12,24,17,3),(14,28,56,112,97,67,7),(238,4,240,118,195,146,34,74)]
    qs=list(dict.fromkeys(qs));saved=[];requests=[];expected=[];inv_checks=0;coord=0;corners=0
    for ix,q in enumerate(qs):
        m=len(q);e,o,E,rs,ic=banks(q);ec,oc,Ec,_,_=banks(opposite(q),False);inv_checks+=ic;require(o==oc,'odd cancellation');require(E==Ec and E>0,'energy sign');bs=basis(m);d={k:e[k]+o[k]-ec[k]-oc[k] for k in e};coord+=len(d)
        direct=sum((1 if q[i]>>j&1 else -1)*d[(1<<i)|(1<<j),1<<i] for i in range(m) for j in range(m) if i!=j)
        require(E==direct,'energy factor 2^(6-m)');require(E.denominator%2==1 and val(E)>=0,'energy integrality')
        require(all(not z for (a,c),z in d.items() if a.bit_count()+c.bit_count()==4),'quartic cancellation')
        for i in range(m):
            for j in range(m):
                if i==j:continue
                formula=F(64,2**m)*sum(r[ids.index(i)]**2*r[ids.index(j)] for mask,(ids,r) in rs.items() if len(ids)%2==0 and i in ids and j in ids)
                require(d[(1<<i)|(1<<j),1<<i]==formula,'individual shared cubic constant')
        L=val(E)+1;M=1<<L
        require(mod(E,M)!=0,'finite modulus');chosen=next((i,j) for i in range(m) for j in range(m) if i!=j and mod(d[(1<<i)|(1<<j),1<<i],M))
        i,j=chosen;cube=[(a,c) for a in [0,1<<i,1<<j,(1<<i)|(1<<j)] for c in [0,1<<i]]
        require(any(mod(evalbank(d,a,c),M) for a,c in cube),'actual query cube');corners+=len(cube)
        # Full block-coordinate banks in both orientations, retaining all zeros.
        saved.append({'rows':q,'energy':encode(E),'L':L,'even':[encode(e[a,c]) for a in bs for c in bs],'odd':[encode(o[a,c]) for a in bs for c in bs],'converse_even':[encode(ec[a,c]) for a in bs for c in bs],'converse_odd':[encode(oc[a,c]) for a in bs for c in bs]})
        # All query vectors for small controls, and full finite-difference cube for larger controls.
        tests=[(a,c) for a in range(1<<m) for c in range(1<<m)] if (m<=3 or (m==4 and ix in [10,11,17,35,73])) else cube
        # Moduli through eight bits, exact complete profile coordinates up to the cutoff.
        for a,c in tests:
            for qq,ee,oo in [(q,e,o),(opposite(q),ec,oc)]:
                requests.append(' '.join(map(str,[m,L,a,c,*qq])))
                ev=evalbank(ee,a,c);od=evalbank(oo,a,c)
                expected.append([mod(ev*(1 if h==0 else 2 if h==1 else 0)+od*((-2)**h)*comb(h+2,2),M) for h in range(m+L+1)])
    (P/'run_input.txt').write_text('\n'.join(requests)+'\n')
    exe=P/'run_walk';subprocess.run(['g++','-O2','-std=c++17',str(P/'run_walk.cpp'),'-o',str(exe)],check=True)
    out=subprocess.run([str(exe)],input='\n'.join(requests)+'\n',text=True,capture_output=True,check=True).stdout
    actual=[list(map(int,line.split())) for line in out.splitlines()];require(actual==expected,'independent finite word profiles')
    (P/'run_output.txt').write_text(out);(P/'rational_banks.json').write_text(json.dumps(saved,separators=(',',':'))+'\n')
    # Direct Stirling recurrence checks even beyond physical order, at moduli beyond guaranteed bounds.
    stir_checks=0
    for L in range(1,9):
        M=1<<L
        for N in [M,2*M,3*M]:
            q=[0]*(N+L+2);q[0]=1
            for n in range(1,N+1):
                q=[0]+[(k*(q[k]+q[k-1]))%M for k in range(1,len(q))]
            for k in range(1,len(q)):
                require(q[k]==pow(-2,k-1,M),'uniform Stirling congruence');stir_checks+=1
    receipt={'status':'PASS','quotients':len(qs),'rational_named_coordinates_each_orientation':coord,'inverse_segment_identities':inv_checks,'finite_difference_corners':corners,'independent_complete_run_queries':len(requests),'run_coordinate_comparisons':sum(map(len,expected)),'stirling_congruences':stir_checks,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['audit.py','run_walk.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['rational_banks.json','run_input.txt','run_output.txt']}}
    (P/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
