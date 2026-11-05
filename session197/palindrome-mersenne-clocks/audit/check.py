"""Independent controls for the named palindrome clock snapshot, stdlib only."""
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x: raise RuntimeError(msg)

def f(n):
    return 1+sum((n>>(i+1))+((n>>(i+1))&1) for i in range(n.bit_length()) if not(n>>i&1))

def slot_H(n):
    # Independent literal Hamilton words in the transitive deletion. The new
    # first base vertex contributes the low zero of the incidence pattern.
    bits=n<<1; m=n.bit_length()+1; full=(1<<m)-1; h=0
    for before in range(1<<m):
        after=full^before
        if before and not(bits>>(before.bit_length()-1)&1): continue
        if after and bits>>((after&-after).bit_length()-1)&1: continue
        h+=1
    return h

def graph(n):
    # Apex label 0, unlike the submitted final-apex convention.
    m=n.bit_length()+1; rows=[0]*(m+1)
    for i in range(1,m+1):
        for j in range(i+1,m+1): rows[i]|=1<<j
        if i>=2 and n>>(i-2)&1: rows[i]|=1
        else: rows[0]|=1<<i
    return rows,list(range(1,m+1))

def cycle_count(rows):
    full=(1<<len(rows))-1
    def visit(mask,last):
        if mask==full: return (rows[last]&1)!=0
        out=rows[last]&(full^mask); total=0
        while out:
            bit=out&-out;out-=bit
            total+=visit(mask|bit,bit.bit_length()-1)
        return total
    return visit(1,0)

def born(rows,base,a,b,v):
    u=len(rows); full=(1<<u)-1
    S=rows+[full^(1<<a)^(1<<b)]
    S[a]|=1<<u; S[b]|=1<<u
    z=len(S)
    for i in range(z):
        if i!=v: S[i]|=1<<z
    S.append(1<<v)
    return S,[a,u]+base[1:-1]+[z,v]

def formula(c,p,q,t):
    reverse=int(bin(c)[2:][::-1],2)
    return (1<<(p+q+c.bit_length()+1+2*t))+((c<<(q+1))-(reverse<<p))*(1<<t)-c+f(c)+1

def main():
    counts=dict(formulas=0,port_returns=0,physical=0,offsets=0,mersennes=0,prime_divisors=0,minimality=0)
    digest=hashlib.sha256()
    for c in range(1,1024,2):
        rev=int(bin(c)[2:][::-1],2)
        for p in range(1,5):
            for q in range(4):
                B=(c<<(q+1))-(rev<<p)
                need((B==0)==(p==q+1 and rev==c),'iff pure cancellation')
                for t in range(13):
                    w='1'+'0'*(p+t)+bin(c)[2:]+'0'*(q+t)
                    n=int(w,2); h=n+f(n)+(n&1)
                    need(h==formula(c,p,q,t),'exact general two-block polynomial')
                    need(formula(c,p,q,t+1)-4*h==-2*B*(1<<t)+3*(c-f(c)-1),'return term including nonpure cases')
                    counts['formulas']+=1
                    digest.update(f'{c},{p},{q},{t},{h}\n'.encode())
                if c<64:
                    w='1'+'0'*p+bin(c)[2:]+'0'*q
                    rows,base=graph(int(w,2)); a=1;b=0;v=base[-1]
                    for t in range(1,6):
                        rows,base=born(rows,base,a,b,v)
                        w='1'+'0'*(p+t)+bin(c)[2:]+'0'*(q+t)
                        expected,_=graph(int(w,2)); mapping={b:0}
                        mapping.update({vertex:i+1 for i,vertex in enumerate(base)})
                        need((base[0],base[-1],b)==(a,v,0),'same physical ports')
                        need(all(((rows[x]>>y)&1)==((expected[mapping[x]]>>mapping[y])&1) for x in range(len(rows)) for y in range(len(rows))),'literal retained-port graph equality')
                        counts['port_returns']+=1
    # Every possible canonical strong feedback-one graph through order 12.
    maxima={}
    for n in range(1,1024):
        h=slot_H(n)
        need(h==n+f(n)+(n&1),'literal source count, no recurrence counter')
        maxima[n.bit_length()]=max(maxima.get(n.bit_length(),0),h)
        if n<128:
            rows,base=graph(n)
            need(cycle_count(rows)==1,'literal unique Hamilton cycle')
        counts['physical']+=1
    for k in range(4,12):
        need(max(v for width,v in maxima.items() if width<=k-1)<(1<<k)-1,'all smaller feedback-one classes')
        counts['minimality']+=1
    pals=[c for c in range(1,65536,2) if bin(c)[2:]==bin(c)[2:][::-1]]
    for delta in range(1,1000,2):
        bounded=[c for c in range(delta+2,2*delta+4,2) if bin(c)[2:]==bin(c)[2:][::-1] and c-f(c)-1==delta]
        need(bounded==[c for c in pals if c-f(c)-1==delta],'full bounded offset inverse')
        counts['offsets']+=1
    need([c for c in pals if c-f(c)-1==-1]==[1],'opposite unit offset')
    need([c for c in pals if c-f(c)-1==1]==[3,5],'Mersenne unit offset')
    for k in list(range(4,1025))+[2047,2048,4095,4096]:
        c=3 if k%2==0 else 5; t=(k-4)//2
        n=(1<<(k-1))+(c<<t)
        need(n.bit_length()==k and n+f(n)+(n&1)==(1<<k)-1,'arbitrarily large exact Mersenne formula')
        counts['mersennes']+=1
    witness=json.loads((HERE/'source'/'witness.json').read_text())
    for item in witness['witnesses']:
        n=item['L']-1; rows,base=graph(n); submitted=item['rows']; mapping={0:item['apex']}
        mapping.update({v:v-1 for v in base})
        need(all(((rows[x]>>y)&1)==((submitted[mapping[x]]>>mapping[y])&1) for x in range(len(rows)) for y in range(len(rows))),'literal published witness arcs')
        need(slot_H(n)==item['H'] and cycle_count(rows)==item['C']==1,'literal named witness counts')
        need(item['order']==item['exponent']+2,'sharp witness order')
        counts['physical']+=1
    for cert in witness['prime_certificates']:
        p=cert['prime']; expected=[[d,p%d] for d in range(2,math.isqrt(p)+1)]
        need(cert['remainders']==expected and all(r for _,r in expected),'complete literal prime certificate')
        counts['prime_divisors']+=len(expected)
    for k in range(4,1001):
        h=(1<<k)-1
        if h%7: continue
        j=k//3; a=(1<<j)-1; b=(1<<(2*j))+(1<<j)+1
        need(a*b==h and (a==7 or b==7)==(k==9),'all seven-prime factor controls')
    # The small non-feedback-one witness prevents an unrestricted minimum claim.
    R5=[6,12,24,17,3]
    hp=sum(all(R5[p[i]]>>p[i+1]&1 for i in range(4)) for p in itertools.permutations(range(5)))
    need(hp==15,'R5 order-minimality scope hostile')
    print(json.dumps(dict(status='PASS',counts=counts,digest=digest.hexdigest(),scope='Independent exact binary evaluation, literal slot counts and literal fixed-port arcs; no source implementation imported.'),sort_keys=True,indent=2))

if __name__=='__main__': main()
