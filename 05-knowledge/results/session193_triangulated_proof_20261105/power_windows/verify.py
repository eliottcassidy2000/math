#!/usr/bin/env python3
"""Standalone exact physical certificates, using standard-library integer DP."""
from itertools import product, permutations
from pathlib import Path
import hashlib,json
BASE=Path(__file__).resolve().parent

def need(test,msg):
    if not test:raise RuntimeError(msg)
def bits(m):
    while m:
        b=m&-m;yield b.bit_length()-1;m-=b

def validate(r):
    n=len(r);full=(1<<n)-1
    for i,a in enumerate(r):
        need(a>=0 and a&~full==0 and not(a>>i&1),'invalid outgoing row')
        for j in range(i):need((a>>j&1)+(r[j]>>i&1)==1,'not a tournament')
def tables(r):
    n=len(r);N=1<<n;E=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
    for i in range(n):E[1<<i][i]=1
    for M in range(1,N):
        for i in bits(M):
            a=E[M][i];H[M]+=a
            if a:
                for j in bits(r[i]&~M):E[M|(1<<j)][j]+=a
    return H,E

def hp(r):return tables(r)[0][-1]
def reversed_rows(r):return [((1<<len(r))-1)^a^(1<<i)for i,a in enumerate(r)]
def strong(r):
    for rr in (r,reversed_rows(r)):
        S=1
        while True:
            T=S
            for i in bits(S):T|=rr[i]
            if T==S:break
            S=T
        if S!=(1<<len(r))-1:return False
    return True

def build(r,mask,word=''):
    n=len(r);m=len(word);full=(1<<(n+m))-1
    out=[a|(full^((1<<n)-1)) for a in r]
    out += [full^((1<<(n+j+1))-1)for j in range(m)]
    S=mask+sum((int(b)<<(n+j))for j,b in enumerate(word))
    out=[a|((not(S>>i&1))<<(n+m))for i,a in enumerate(out)]+[S]
    validate(out);return out

def profile(r,mask):
    H,E=tables(r);_,starts=tables(reversed_rows(r));full=(1<<len(r))-1
    R=sum(starts[A][i]*H[full^A]for A in range(1,full+1)for i in bits(A&mask))
    e=sum(E[-1][i]for i in bits(mask))
    return H[-1],hp(build(r,mask)),R,e

def word_response(h,F,R,e,word):
    for b in word:
        if b=='0':F,R,e=F+R+e,2*R,0
        else:R,e=2*R+h,h
    return F,R,e

def literal(r):
    return sum(all(r[p[i]]>>p[i+1]&1 for i in range(len(r)-1))for p in permutations(range(len(r))))

def main():
    obj=json.loads((BASE/'seeds135.json').read_text());seeds=obj['seeds']
    need([z['offset']for z in seeds]==list(range(-115,154,2)),'exact 135 offset range')
    records=[];orders=[];physical=0
    for z in seeds:
        r=z['rows'];S=z['mask'];q=z['q'];b=z['offset'];validate(r)
        need(strong(r) and S>0,'strong core and nonempty attachment')
        h,F,R,e=profile(r,S)
        need(R==z['R']==1<<q,'power response')
        need(F-R+e==b,'intercept identity')
        need(R>=h and 0<=e<=h,'response bounds')
        need(q<=12,'k13 threshold')
        t=13-q;child=build(r,S,'0'*t);H=hp(child);physical+=1
        need(strong(child),'strong physical padded child')
        need(H==(1<<13)+b==F+(2**t-1)*R+e,'physical dyadic value')
        orders.append(len(child))
        records.append(dict(offset=b,h=h,F=F,R=R,e=e,core_order=len(r),child_order_at_k13=len(child),H_at_k13=H))
    # Independently evaluate the binary recurrence on explicit tournaments.
    wordchecks=0
    for b in (-115,-25,1,27,99,153):
        z=next(x for x in seeds if x['offset']==b);r=z['rows'];S=z['mask'];h,F,R,e=profile(r,S)
        for m in range(4):
            responses=set()
            for w in map(''.join,product('01',repeat=m)):
                FF,RR,ee=word_response(h,F,R,e,w)
                need(FF==hp(build(r,S,w)),'binary transfer physical H')
                need(RR==(1<<m)*R+h*(int(w,2)if w else 0),'integer carry response')
                need(RR not in responses,'same-length response injectivity');responses.add(RR)
                need((RR&(RR-1)==0)==(set(w)<={'0'}),'pure-power no-digit boundary')
                wordchecks+=1
    # Literal small controls do not share the subset recurrence.
    literalchecks=0
    for r,S in [([2,4,1],1),([20,1,2,7,14],10),([18,8,3,5,14],6)]:
        for w in ('','0','1'):
            rr=build(r,S,w);need(literal(rr)==hp(rr),'literal permutation control');literalchecks+=1
    report={'status':'PASS: finite exact physical seeds and controls; universal k follows from PROOF.md',
      'offsets':[-115,153,2],'number_of_offsets':len(seeds),'minimum_k':13,
      'core_order_range':[min(len(z['rows'])for z in seeds),max(len(z['rows'])for z in seeds)],
      'order_range_at_k13':[min(orders),max(orders)],'order_offsets_relative_to_k':sorted({n-13 for n in orders}),
      'physical_padded_children_checked':physical,'binary_word_checks':wordchecks,'literal_controls':literalchecks,
      'source_sha256':{p:hashlib.sha256((BASE/p).read_bytes()).hexdigest()for p in ['seeds135.json','verify.py']},'records':records}
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
