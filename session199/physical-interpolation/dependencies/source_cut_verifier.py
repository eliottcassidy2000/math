#!/usr/bin/env python3
"""Exact source-cut and rearrangement checks, standard library only."""
import hashlib, itertools, json
from functools import lru_cache
from pathlib import Path

checks={}
def check(ok, label):
    if not ok: raise RuntimeError(label)
    checks[label.split(':')[0]]=checks.get(label.split(':')[0],0)+1

def step(s,x):
    D,E,U,d=s
    if x=='O':return D,0,2*U,d
    if x=='Z':return D-E+U,d,2*U-d,d
    if x=='J':return D-E+U,D-E,D+U,D
    raise ValueError(x)
def word(s,w):
    for x in w:s=step(s,x)
    return s

def bits(a,k):return tuple((a>>i)&1 for i in range(k))
def letters(a,k):return ''.join('Z' if e else 'O' for e in bits(a,k))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
@lru_cache(None)
def g(n):
    if n==1:return 1
    if n%2==0:return g(n//2)
    l=(n+1)//2
    return g(l)+l-l%2

def weights(k):return [(i,j,2**(j-i-1)+(j==i+1)) for i in range(k) for j in range(i+1,k)]
def cut(a,k):return sum(w for i,j,w in weights(k) if ((a>>i)^(a>>j))&1)
def block_formula(s,a,k):
    D,E,U,d=s;b=bits(a,k);R=sum(e*2**(k-1-i) for i,e in enumerate(b))
    C=sum(w*b[i]*b[j] for i,j,w in weights(k))
    return D+U*a-E*b[0]-d*C,d*b[-1],2**k*U-d*R,d

def cut_scalar(s,row,a,k):
    D,E,U,d=s;al,be,ga,et=row;b=bits(a,k)
    cs=[(i,j,al*d*w//2) for i,j,w in weights(k)]
    ls=[al*U*2**i-al*E*(i==0)+be*d*(i==k-1)-ga*d*2**(k-1-i) for i in range(k)]
    for i,j,c in cs:ls[i]-=c;ls[j]-=c
    return al*D+ga*2**k*U+et*d+dot(ls,b)+sum(c*abs(b[i]-b[j]) for i,j,c in cs)

P={
'N':[[6,4,10,9],[4,3,7,6],[10,7,15,14],[9,6,14,13]],
'M':[[11,8,20,18],[8,6,14,13],[20,14,34,31],[18,13,31,29]],
'A':[[14,10,21,20],[10,7,16,15],[21,16,33,31],[20,15,31,29]],
'B':[[14,10,26,23],[11,8,20,18],[20,15,35,32],[19,14,34,31]],
'C':[[23,18,44,40],[18,14,34,31],[44,34,80,73],[40,31,73,67]],
'J0':[[12,12,24,24],[12,12,21,21],[15,12,30,27],[15,12,30,27]],
'J1':[[12,12,21,21],[12,12,18,18],[18,15,32,30],[18,15,31,29]],
}
transpose=lambda p:[list(x) for x in zip(*p)]
P['T']=transpose(P['B']);P['J2']=transpose(P['J1']);P['J3']=transpose(P['J0'])

def y(s):
    D,E,U,d=s;return D-E-d,U-D+E-d,E,d-E

def readout(p,partner):
    r=[dot(row,y(partner)) for row in p]
    return r[0]-r[1],-r[0]+r[1]+r[2]-r[3],r[1],-r[0]-r[1]+r[3]

# Independent literal graph construction and subset path counting.
def birth(T,c,x):
    n=len(T);R=list(T)+[0]
    for i in range(n):
        if x!='O' and i==c:R[n]|=1<<i
        else:R[i]|=1<<n
    return R,n if x=='J' else c

def graphword(w):
    T,c=[2,0],0
    for x in w:T,c=birth(T,c,x)
    return T,c

def inventories(T):
    n=len(T);size=1<<n;end=[[0]*n for _ in range(size)];start=[[0]*n for _ in range(size)]
    for v in range(n):end[1<<v][v]=start[1<<v][v]=1
    for mask in range(1,size):
        for v in range(n):
            if mask>>v&1:
                old=mask^(1<<v)
                if old:
                    end[mask][v]=sum(end[old][u] for u in range(n) if old>>u&1 and T[u]>>v&1)
                    start[mask][v]=sum(start[old][u] for u in range(n) if old>>u&1 and T[v]>>u&1)
    H=[sum(row) for row in end];H[0]=1
    return H,end,start

def actual_profile(T,c):
    H,end,start=inventories(T);full=(1<<len(T))-1
    U=sum(H[full^q]*start[q][c] for q in range(full+1))
    check(start[full][c]==1,'literal_root_start')
    return H[full],end[full][c],U,H[full^(1<<c)]

def cycles(T):
    n=len(T);size=1<<n;dp=[[0]*n for _ in range(size)];dp[1][0]=1
    for mask in range(size):
        for u in range(n):
            if dp[mask][u]:
                for v in range(n):
                    if not(mask>>v&1) and T[u]>>v&1:dp[mask|1<<v][v]+=dp[mask][u]
    return sum(dp[-1][u] for u in range(n) if T[u]&1)

def main():
    tau=(1,0,2,1);F=(1,-1,1,0)
    certificates={}
    for name,p0 in P.items():
        for orient,p in [('left',p0),('right',transpose(p0))]:
            al=[p[0][j]-p[1][j] for j in range(4)]
            ga=p[1];ta=[p[0][j]+p[2][j]-p[3][j] for j in range(4)]
            extra=[p[2][j]-p[0][j] for j in range(4)]
            check(all(x>=0 for x in extra) and extra[2]>0 and extra[3]>0,'strong_cone_gate')
            check(all(x>=0 for x in al) and al[2]>0 and al[3]>0,'gate_positivity')
            check(all(x>0 for x in ga+ta),'gate_positivity')
            certificates[name+'_'+orient]={'alpha':al,'gamma':ga,'tau':ta,'r2_minus_r0':extra}
    # Complete physical source words through five extra vertices, not a class census.
    seeds={tau}
    for k in range(6):
        for letters0 in itertools.product('OZJ',repeat=k):
            w=''.join(letters0);s=word(tau,w);seeds.add(s)
            T,c=graphword(w)
            check(actual_profile(T,c)==s,'literal_profile')
            if k<=3:
                C,_=birth(T,c,'J');H,_,_=inventories(C)
                check(H[-1]==dot(F,s) and cycles(C)==1,'literal_c1')
    seedlist=sorted(seeds)
    # Exact source-independent identities on a vector-space basis, no physicality assertion here.
    basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
    for k in range(1,8):
        for a in range(1<<k):
            for s in basis:
                check(word(s,letters(a,k))==block_formula(s,a,k),'block_identity')
    # Exhaustive digit cuts through length fourteen, including equality classifications.
    cut_minima={}
    for k in range(1,15):
        vals=[];mins=[];mn=None
        for a in range(1<<k):
            K=cut(a,k);check(g(a+1)+g((1<<k)-a)-2==K,'cut_identity')
            if a not in (0,(1<<k)-1):
                if mn is None or K<mn:mn=K;mins=[a]
                elif K==mn:mins.append(a)
        if k>=2:
            expected=2 if k==2 else 2**((k-1)//2)+2**((k)//2)
            check(mn==expected,'cut_minimum')
            if k>=4:
                central={1<<((k-1)//2),1<<(k//2)};allbits=(1<<k)-1
                expected_args=central|{allbits^a for a in central}
                check(set(mins)==expected_args,'cut_equality')
        cut_minima[str(k)]={'minimum':mn,'minimizers':mins}
    rows=[F,(5,-2,16,-4),(3,-3,9,-6)]
    for partner in [tau,word(tau,'J'),word(tau,'ZO'),word(tau,'OZJ')]:
        for p in P.values():rows.extend([readout(p,partner),readout(transpose(p),partner)])
    rows=sorted(set(rows))
    def positive_row(row):
        al,be,ga,et=row
        return al+ga,ga,2*al+be+2*ga+et,al+2*ga+et
    def cone(r):return r[2]>r[0]>r[1]>0 and r[3]>0 and r[0]+r[2]>r[3]
    def transfer(r,x):
        a,b,c,d=r
        if x=='O':return a+b,2*b,a+b+d,2*b+d
        if x=='Z':return 2*a,a+b,a+c,a+c
        if x=='J':return c,a,c+d,a+c
        raise ValueError(x)
    for row in rows:
        check(cone(positive_row(row)),'strong_cone_initial')
        for k in range(6):
            for w0 in itertools.product('OZJ',repeat=k):
                w=''.join(w0);r=positive_row(row)
                for x in reversed(w):r=transfer(r,x)
                check(cone(r),'strong_cone_suffix')
                for s in [tau,word(tau,'OZJ')]:
                    actual=dot(row,word(s,w));check(actual==dot(r,y(s)),'strong_cone_transport')
                    lo=dot(row,word(s,'ZO'+w));hi=dot(row,word(s,'OZ'+w))
                    check(hi>lo,'arbitrary_suffix_exchange')
    # Exhaust every retained block through seven on selected genuinely different physical seeds.
    selected=[tau,word(tau,'J'),word(tau,'OZJ'),word(tau,'JOZ')]
    for s in selected:
        for k in range(1,8):
            states=[word(s,letters(a,k)) for a in range(1<<k)]
            for row in rows:
                values=[dot(row,v) for v in states]
                by_r={r:[a for a in range(1<<k) if a.bit_count()==r] for r in range(k+1)}
                for r,ids in by_r.items():
                    low=(1<<r)-1;high=((1<<r)-1)<<(k-r)
                    check(all(values[low]<=values[a]<=values[high] for a in ids),'rearrangement')
                    check(sum(values[a]==values[low] for a in ids)==1 and sum(values[a]==values[high] for a in ids)==1,'unique_extrema')
                for a in range(1<<k):check(values[a]==cut_scalar(s,row,a,k),'scalar_cut')
                # Independent closed interior formulas, with chronological word checks above.
                D,E,U,d=s;al,be,ga,et=row
                for r in range(1,k):
                    low=al*(D-E-U+2*d)+et*d+ga*2**k*(U-d)+al*2**r*(U-d)+ga*2**(k-r)*d
                    high=al*D+(al+ga)*2**k*U+(be+et+2*al+ga)*d-al*2**(k-r)*U-(al+ga)*2**r*d
                    check(low==dot(row,word(s,'Z'*r+'O'*(k-r))),'closed_extrema')
                    check(high==dot(row,word(s,'O'*(k-r)+'Z'*r)),'closed_extrema')
    mins={}
    for k in range(2,15):
        vs=[dot(F,word(tau,letters(a,k))) for a in range(1<<k)]
        actual=min(vs)
        max_expected=(2**(k//2+1)-1)*(2**((k+1)//2+1)-1)+2
        check(max(vs)==max_expected,'tau_maximum')
        expected=(2**(k//2)+1)*(2**((k+1)//2)+1)
        check(actual==expected,'tau_minimum');mins[str(k)]={'minimum':actual,'maximum':max(vs)}
    values=sorted({dot(F,word(tau,letters(a,4))) for a in range(16)})
    check(values==[25,27,31,33,37,39,45,47,51],'interior_gap_hostile')
    # Literal non-extremal witnesses ensure count optimization is not secretly only formal.
    literal=[]
    for k in range(2,9):
        w='Z'*(k//2)+'O'*((k+1)//2);T,c=graphword(w);C,_=birth(T,c,'J');H,_,_=inventories(C)
        expected=(2**(k//2)+1)*(2**((k+1)//2)+1)
        check(H[-1]==expected and cycles(C)==1,'literal_minimum')
        literal.append({'k':k,'word':w,'H':H[-1],'order':len(C)})
    out={'status':'PASS','checks':checks,'physical_seed_profiles':len(seeds),'readout_rows':len(rows),'gate_certificates':certificates,'cut_minima':cut_minima,'tau_minima':mins,'literal_minima':literal,'hostile_k4_image':values}
    data=json.dumps(out,sort_keys=True,indent=2)+'\n'
    Path(__file__).with_name('RESULT.json').write_text(data)
    print(data,end='')

if __name__=='__main__':main()
