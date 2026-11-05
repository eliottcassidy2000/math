#!/usr/bin/env python3
"""Standalone finite controls for universal zero slack and the grammar envelope."""
from itertools import combinations, permutations, product
from collections import deque
from pathlib import Path
import json, hashlib

HERE=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def bits(m):
    while m:
        b=m&-m;yield b.bit_length()-1;m-=b

def tournament(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return tuple(a)

def subgraph(a,labels):
    return tuple(sum(1<<j for j,v in enumerate(labels) if a[u]>>v&1) for u in labels)

def tables(a):
    n=len(a);E=[[0]*n for _ in range(1<<n)];S=[[0]*n for _ in E];H=[0]*len(E);H[0]=1
    for m in range(1,len(E)):
        for i in bits(m):
            prev=m^(1<<i)
            if not prev:E[m][i]=S[m][i]=1
            else:
                E[m][i]=sum(E[prev][j] for j in bits(prev) if a[j]>>i&1)
                S[m][i]=sum(S[prev][j] for j in bits(prev) if a[i]>>j&1)
        H[m]=sum(E[m]);need(H[m]==sum(S[m]),'start/end consistency')
    return H,E,S

def boundary(a,tab):
    H,E,S=tab;n=len(a);full=(1<<n)-1;B=[[0]*n for _ in a]
    for m in range(1,full):
        rest=full^m
        for y in bits(m):
            for x in bits(rest):B[y][x]+=E[m][y]*S[rest][x]
    return B

def profile(a,r,tab=None):
    H,E,S=tab or tables(a);full=len(H)-1
    R=sum(S[m][r]*H[full^m] for m in range(1,full+1) if m>>r&1)
    return H[full],E[full][r],R,H[full^(1<<r)]

def slack(q):
    h,e,R,d=q;return R+e-h-d

def path(a,w):return all(a[x]>>y&1 for x,y in zip(w,w[1:]))

def pair_injection(a,x,y,w):
    for i,z in enumerate(w):
        possible=[v for v in (x,y) if a[v]>>z&1]
        if possible:
            first=possible[0];last=y if first==x else x
            return w[:i]+(last,),(first,)+w[i:]
    return w+(x,),(y,)

def literal_slack(a,r):
    count=0
    for w in permutations(range(len(a))):
        cut=w.index(r)
        if 0<cut<len(a)-1 and a[r]>>w[cut-1]&1 and path(a,w[:cut]) and path(a,w[cut:]):count+=1
    return count

def unary(a,r,op):
    n=len(a)
    if op=='D':return tuple(x|(1<<n) for x in a)+(0,),r
    return tuple(x|((1<<n) if i!=r else 0) for i,x in enumerate(a))+(1<<r,),n if op=='A' else r

def change(q,op,bank):
    h,e,R,d=q
    if op=='K':return h+R-e,d,2*R-d,d
    if op=='D':return h,0,2*R,d
    g,s,L,E,Q,B,delta,dg=bank[op]
    return g*h-s*e+L*R,E*h-delta*e,Q*h+B*R,dg*h

def gadget(a,r):
    H,E,S=tables(a);full=len(H)-1;deleted=full^(1<<r)
    L=sum(E[m][r]*H[full^m] for m in range(1,full+1) if m>>r&1)
    R=sum(S[m][r]*H[full^m] for m in range(1,full+1) if m>>r&1)
    B=sum(H[m]*H[deleted^m] for m in range(full+1) if m&deleted==m)
    return H[-1],S[-1][r],L,E[-1][r],R,B,int(len(a)==1),H[deleted]

def check_normal_form(ops,states,bank,D,C):
    final=states[-1];need(slack(final)==0 and final[3]<=D,'normal-form domain')
    positions=[i for i,op in enumerate(ops) if isinstance(op,int)]
    if positions:
        last=positions[-1];z=states[last+1];prefix=states[last];g=bank[ops[last]]
        need(prefix[0]<=D and z[3]==g[7]*prefix[0],'bounded root-move input')
        k=0
        while last-k>0 and ops[last-k-1]=='D':k+=1
        base=states[last-k]
        need(base[2]<=C*D and prefix[2]==base[2]*(2**k),'bounded-count dyadic prefix')
        if k:need(prefix[1]==0,'sink prefix endpoint')
        suffix=ops[last+1:];suffix_states=states[last+1:]
    else:
        z=states[0];suffix=ops;suffix_states=states;k=0
    need(slack(z)==0,'positive slack cannot be repaired by K/D')
    d=z[3]
    if 'D' in suffix:
        i=suffix.index('D')
        need(suffix.count('D')==1 and all(x=='K' for x in suffix[:i]+suffix[i+1:]),'zero-slack suffix pattern')
        need(suffix_states[i]==(d,d,d,d),'only sink state can admit D at zero slack')
        j=len(suffix)-i-1
        expect=d if j==0 else d+d*(2**j)
    else:
        need(all(x=='K' for x in suffix),'suffix generator')
        j=len(suffix);expect=z[0] if j==0 else d+(2**j)*(z[0]-z[1])
        if positions and k:
            h,e,R,old_d=prefix;gg,s,L,E,Q,B,delta,dg=g;r0=base[2]
            A=(gg-dg)*h if j==0 else (gg-E)*h
            need(A>=0 and L*r0>0,'nonnegative sparse family')
            need(expect==d+(2**j)*(A+L*r0*(2**k)),'double-dyadic family')
    need(final[0]==expect,'normal form count')

def main():
    counts=dict(tournaments=0,rooted_checks=0,capacity_pairs=0,literal_pair_injections=0,literal_slack_checks=0,unit_deletion_zero_slack=0)
    cache={():([1],[[]],[[]],[])}
    for n in range(1,7):
        full=(1<<n)-1
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code);tab=tables(a);H,E,S=tab;counts['tournaments']+=1
            if n<=5:
                B=boundary(a,tab);cache[a]=(H,E,S,B)
                for x,y in combinations(range(n),2):
                    need(B[x][y]+B[y][x]>=H[full^(1<<x)^(1<<y)],'positive capacity lower bound')
                    counts['capacity_pairs']+=1
                    images=set()
                    others=[i for i in range(n) if i not in (x,y)]
                    for w in permutations(others):
                        if not path(a,w):continue
                        P,Q=pair_injection(a,x,y,w)
                        need(P and Q and path(a,P) and path(a,Q),'two-endpoint injected paths')
                        need({P[-1],Q[0]}=={x,y},'two-endpoint tags')
                        recover=tuple(z for z in P+Q if z not in (x,y))
                        need(recover==w,'two-endpoint recovery')
                        need((P,Q) not in images,'two-endpoint injection collision');images.add((P,Q))
                        counts['literal_pair_injections']+=1
                    need(len(images)==H[full^(1<<x)^(1<<y)],'injection source universe')
            for r in range(n):
                q=profile(a,r,tab);f=slack(q);labels=[i for i in range(n) if i!=r]
                small=subgraph(a,labels);Bs=cache[small][3];X=[i for i,u in enumerate(labels) if a[r]>>u&1]
                actual=sum(Bs[y][x] for x in X for y in X if x!=y)
                need(f==actual and f>=0,'exact slack boundary identity')
                need((f==0)==(a[r].bit_count()<=1),'zero slack iff outdegree')
                counts['rooted_checks']+=1
                if q[3]==1 and f==0:
                    h=q[0];need(h==1 or (h>2 and (h-1)&(h-2)==0),'unit-deletion sparse spectrum')
                    counts['unit_deletion_zero_slack']+=1
                if n<=4:
                    need(literal_slack(a,r)==f,'literal nonmergeable-cover slack');counts['literal_slack_checks']+=1
    positive=[]
    for m in range(1,9):
        n=m+1;chain=tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
        a,r=unary(chain,0,'A');q=profile(a,r)
        need(q==(2**m+1,1,2**m+1,1),'unit-deletion positive family');positive.append(dict(m=m,state=q))
    hostile=(10,12,9,16,7);need(tables(hostile)[0][-1]==15,'H15 normalization hostile')
    product_controls=[]
    for exponents in product(range(2),repeat=3):
        a=(0,);r=0;q=(1,1,1,1);expected=1
        for exponent,k in zip(exponents,(1,2,4)):
            for _ in range(exponent):
                old=q[0]
                for op in ['D']*k+['A']:a,r=unary(a,r,op)
                q=profile(a,r);expected*=1+2**k
                need(q==(expected,old,expected,old),'multiplicative actual graph macro')
        product_controls.append(dict(exponents=exponents,H=q[0],d=q[3],order=len(a)))
    bank=[gadget((0,),0),gadget((0,1),0),gadget((0,1,3),0),gadget((52,17,2,7,12,30),3),gadget((56,1,3,6,14,30),1)]
    bound=300;D=5;C=max([2]+[g[4]+g[5] for g in bank]);seed=(1,1,1,1)
    seen={seed:None};todo=deque([seed]);normal_checks=0;counts_flat=set()
    def history(q):
        qs=[q];ops=[]
        while seen[q] is not None:
            previous,op=seen[q];ops.append(op);qs.append(previous);q=previous
        return list(reversed(ops)),list(reversed(qs))
    while todo:
        q=todo.popleft()
        for op in ['K','D']+list(range(len(bank))):
            child=change(q,op,bank);h,e,R,d=child
            need(0<=e<=d<=h<=R and slack(child)>=0,'physical cone')
            if h<=bound and d<=D and slack(child)==0:
                ops,states=history(q);check_normal_form(ops+[op],states+[child],bank,D,C)
                normal_checks+=1;counts_flat.add(h)
            if R<=bound and child not in seen:seen[child]=(q,op);todo.append(child)
    receipt=dict(status='PASS',counts=counts,unit_deletion_positive_controls=positive,
                 product_controls=product_controls,H15_hostile_rows=hostile,
                 grammar_control=dict(physical_gadgets=len(bank),bound=bound,D=D,C=C,retained_states=len(seen),normal_form_checks=normal_checks,zero_slack_bounded_deletion_counts=sorted(counts_flat)),
                 scope='Universal slack identity and outdegree through order six; capacity/injection through order five; literal slack through order four. Finite grammar controls do not prove density statements.',
                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    rendered=json.dumps(receipt,indent=2,sort_keys=True)+'\n';(HERE/'receipt.json').write_text(rendered);print(rendered,end='')
if __name__=='__main__':main()
