#!/usr/bin/env python3
"""Independent cycle, path, and literal-family verification. No dependencies."""
import itertools, json
from pathlib import Path

def require(test, message):
    if not test: raise RuntimeError(message)

def comb(m, q, tau=0, reverse_apex=False):
    """Vertices x0,...,x(m+1), a=m+2, b=m+3; 0<=q<2**m."""
    require(m >= 1 and 0 <= q < (1 << m) and tau in (0,1), 'bad parameters')
    n=m+4; a=n-2; b=n-1
    A=[0]*n
    def edge(i,j): A[i]|=1<<j
    for i in range(m+2):
        for j in range(i+1,m+2): edge(i,j)
        edge(a,i) if i==1 else edge(i,a)
        loses = tau if i==1 else (q>>(i-2)&1 if i>=2 else 0)
        edge(i,b) if loses else edge(b,i)
    edge(b,a) if reverse_apex else edge(a,b)
    return A

def is_tournament(A):
    return all(not(A[i]>>i&1) for i in range(len(A))) and all(bool(A[i]>>j&1) != bool(A[j]>>i&1) for i in range(len(A)) for j in range(i))

def dp_cycles(A):
    n=len(A); full=(1<<n)-1
    if n<3:return 0
    dp={1:[1]+[0]*(n-1)}
    for mask in range(1,full+1,2):
        row=dp.get(mask)
        if row is None:continue
        for last,z in enumerate(row):
            if not z:continue
            avail=A[last]&~mask&full
            while avail:
                bit=avail&-avail; avail^=bit; j=bit.bit_length()-1
                new=dp.setdefault(mask|bit,[0]*n);new[j]+=z
    return sum(z for v,z in enumerate(dp.get(full,[0]*n)) if A[v]&1)

def dp_paths(A):
    n=len(A); full=(1<<n)-1; dp={1<<v:[int(i==v) for i in range(n)] for v in range(n)}
    for mask in range(1,full+1):
        row=dp.get(mask)
        if row is None:continue
        for last,z in enumerate(row):
            if not z:continue
            avail=A[last]&~mask&full
            while avail:
                bit=avail&-avail;avail^=bit;j=bit.bit_length()-1
                dp.setdefault(mask|bit,[0]*n)[j]+=z
    return sum(dp.get(full,[]))

def brute_cycles(A):
    n=len(A); out=[]
    for p in itertools.permutations(range(1,n)):
        cyc=(0,)+p
        if all(A[cyc[i]]>>cyc[(i+1)%n]&1 for i in range(n)):out.append(cyc)
    return out

def canonical(p):
    j=p.index(0); return tuple(p[j:]+p[:j])

def literal_cycles(m,q,tau=0,reverse_apex=False):
    a=m+2;b=m+3; chain=list(range(m+2)); result=[]
    if not reverse_apex:result.append(canonical([a,b]+chain))
    if tau:result.append(canonical([a,1,b]+[j for j in chain if j!=1]))
    for i in range(m):
        if not(q>>i&1):continue
        for mask in range(1<<i):
            chosen=[1]+[j+2 for j in range(i) if mask>>j&1]+[i+2]
            result.append(canonical([a]+chosen+[b]+[j for j in chain if j not in chosen]))
    return result

def strong(A):
    full=(1<<len(A))-1
    for v in range(len(A)):
        reached=1<<v
        while True:
            nxt=reached
            for w in range(len(A)):
                if reached>>w&1:nxt|=A[w]
            if nxt==reached:break
            reached=nxt
        if reached!=full:return False
    return True

def gates(A):
    n=len(A)
    return [(i,j) for i in range(n) for j in range(n) if A[i]==1<<j and sum(bool(A[v]>>j&1) for v in range(n))==1]

def target(k):
    require(k>=1,'positive target')
    q=max(0,k-2);tau=int(k>=2);m=max(1,q.bit_length())
    return comb(m,q,tau)

def all_tournaments(n):
    pairs=list(itertools.combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        A=[0]*n
        for i,(u,v) in enumerate(pairs):
            if mask>>i&1:A[u]|=1<<v
            else:A[v]|=1<<u
        yield A

def gate(A):
    n=len(A); B=A.copy()+[1<<(n+1),(1<<n)-1]
    for v in range(n):B[v]|=1<<n
    return B

def main():
    report={'status':'FINITE-EXACT verification of the all-order proof, not its replacement','convention':'Hamiltonian directed cycles modulo cyclic rotation; fixed starting vertex 0; n>=3','family_tests':[]}
    for m in range(1,9):
        checked=0; support=set()
        for q in range(1<<m):
            for tau in (0,1):
                for reverse in (False,True):
                    A=comb(m,q,tau,reverse); expect=q+tau+int(not reverse)
                    require(is_tournament(A),'not a tournament')
                    require(dp_cycles(A)==expect,('DP mismatch',m,q,tau,reverse))
                    require(strong(A)==bool(expect),'strongness mismatch')
                    if reverse:require(not gates(A),'reversed family has gate')
                    words=literal_cycles(m,q,tau,reverse)
                    require(len(words)==len(set(words))==expect,'literal count or duplicate')
                    require(all(len(p)==len(A) and len(set(p))==len(A) and all(A[p[i]]>>p[(i+1)%len(A)]&1 for i in range(len(A))) for p in words),'literal cycle invalid')
                    if m<=3:require(set(brute_cycles(A))==set(words),'permutation correspondence mismatch')
                    checked+=1
                    if not reverse:support.add(expect)
        require(support==set(range(1,(1<<m)+2)),'support mismatch')
        report['family_tests'].append({'m':m,'order':m+4,'parameter_cases':checked,'unreversed_support':[1,(1<<m)+1],'permutation_exhaustive':m<=3})
    witnesses={}
    for k in [1,2,3,7,21,64,65,127,128,129,255,256,257,1024,1025]:
        A=target(k); require(dp_cycles(A)==k,'target mismatch')
        witnesses[str(k)]={'order':len(A),'outrows':A,'cycles':k,'paths':dp_paths(A),'strong':strong(A),'gate_arcs':[[i,j] for i in range(len(A)) for j in range(len(A)) if A[i]==1<<j and sum(bool(A[v]>>j&1) for v in range(len(A)))==1]}
    report['witnesses']=witnesses
    growth_cases=0
    for m in range(1,6):
        for q in range(1<<m):
            for tau in (0,1):
                A=comb(m,q,tau,True)
                for bit in (0,1):
                    B=comb(m+1,2*q+bit,tau,True)
                    mapping=[0,1]+list(range(3,m+3))+[m+3,m+4]
                    require(all(bool(A[i]>>j&1)==bool(B[mapping[i]]>>mapping[j]&1) for i in range(len(A)) for j in range(len(A))), 'growth changes old arc')
                    require(dp_cycles(B)==2*dp_cycles(A)-tau+bit,'growth law mismatch')
                    growth_cases+=1
    report['growth_cases']=growth_cases
    controls=[]
    for n in range(1,5):
        count=0
        for A in all_tournaments(n):
            G=gate(A)
            require(dp_cycles(G)==dp_paths(A),'gate law mismatch')
            count+=1
        controls.append({'core_order':n,'labeled_tournaments':count})
    report['inherited_gate_controls']=controls
    report['zero_control']={'outrows':[6,4,0],'cycles':dp_cycles([6,4,0]),'strong':strong([6,4,0])}
    output=Path(__file__).with_name('RESULT.json');output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'family_cases':sum(r['parameter_cases'] for r in report['family_tests']),'witnesses':len(witnesses),'gate_controls':sum(r['labeled_tournaments'] for r in controls),'output':str(output)},indent=2))
if __name__=='__main__': main()
