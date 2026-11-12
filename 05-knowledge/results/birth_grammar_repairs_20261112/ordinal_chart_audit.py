"""Independent slot-cover controls for the ordered-module extension of V2."""
from common import *
from itertools import product
from math import comb
from collections import defaultdict,Counter
from sparse_covers import kernel

def join(parts):
    total=sum(map(len,parts));rows=[];offset=0
    for part in parts:
        n=len(part);later=((1<<total)-1)^((1<<(offset+n))-1)
        rows.extend((r<<offset)|later for r in part);offset+=n
    return rows

def subset_H(rows):
    n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
    for mask in range(1,N):
        for v in range(n):
            if not mask>>v&1:continue
            rest=mask^(1<<v)
            dp[mask][v]=1 if not rest else sum(dp[rest][u] for u in range(n) if rows[u]>>v&1)
        h[mask]=sum(dp[mask])
    return h

def covers(rows,R=4):
    h=subset_H(rows);N=len(h);q=[0]*N;c=[0]*N;q[0]=c[0]=1;Q=[0];C=[0]
    for r in range(1,R+1):
        qq=[0]*N;cc=[0]*N
        for mask in range(N):
            sub=mask
            while True:
                qq[mask]+=q[mask^sub]*h[sub]
                if sub:cc[mask]+=c[mask^sub]*h[sub]
                if sub==0:break
                sub=(sub-1)&mask
        q,c=qq,cc;Q.append(q[-1]);C.append(c[-1])
        need(Q[r]==sum(comb(r,j)*C[j] for j in range(1,r+1)),'slot versus nonempty covers')
    return Q,C

def inflate(Q,modules):
    owner=[];offsets=[];rows=[]
    for i,A in enumerate(modules):offsets.append(len(owner));owner.extend([i]*len(A))
    for i,A in enumerate(modules):
        exterior=sum(1<<v for v,j in enumerate(owner) if Q[i]>>j&1)
        rows.extend((r<<offsets[i])|exterior for r in A)
    return rows

def run():
    atoms=[[0],[2,0],[2,4,1]];profiles=[covers(A) for A in atoms];stats=Counter();records=[]
    for length in range(1,4):
        for word in product(range(3),repeat=length):
            rows=join([atoms[i] for i in word]);Q,C=covers(rows)
            for r in range(1,5):
                expected=1
                for i in word:expected*=profiles[i][0][r]
                need(Q[r]==expected,'physical ordinal-product slot law');stats['slot_identities']+=1
                need(C[r]==sum((-1)**(r-j)*comb(r,j)*Q[j] for j in range(1,r+1)),'binomial inversion');stats['nonempty_identities']+=1
            stats['physical_words']+=1
    left=join([atoms[2],atoms[0]]);right=join([atoms[0],atoms[2]])
    need(covers(left)==covers(right),'same entire small cover profile for the hostile')
    need(sorted(r.bit_count() for r in left)!=sorted(r.bit_count() for r in right),'nonisomorphic order reversal hostile')
    outer=[2,4,1];D,_=kernel(outer,[0,1],[]);R=max(max(c) for c in D);aq,ac=covers(atoms[2],R);bq,bc=covers(atoms[0],R)
    exponential=defaultdict(int)
    for (r,s),coefficient in D.items():
        for i in range(1,r+1):
            for j in range(1,s+1):exponential[aq[i]*bq[j]]+=coefficient*(-1)**(r-i+s-j)*comb(r,i)*comb(s,j)
    exponential={b:c for b,c in exponential.items() if c}
    for a,b in product(range(1,4),range(1,4)):
        A=join([atoms[2]]*a);B=join([atoms[0]]*b);_,ac=covers(A,R);_,bc=covers(B,R)
        predicted=sum(c*ac[r]*bc[s] for (r,s),c in D.items());rows=inflate(outer,[A,B,[0]])
        actual=int(tables(np.array(rows,dtype=np.int64))[0]);need(predicted==actual,'whole graph substitution with nontransitive growing module')
        if a==b:need(actual==sum(c*base**a for base,c in exponential.items()),'whole exponential ordered-module diagonal')
        records.append({'a':a,'b':b,'order':len(rows),'H':actual})
    out={'scope':'All39 words of lengths1..3 in the fixed atoms singleton,TT2,C3; four slot/nonempty cover channels by independent subset convolution; nine literal C3[ordinal C3 powers,TT length,point] controls. These controls accompany the all-length binomial-inversion proof, not an extrapolation from finite data.','statistics':dict(stats),'atoms':atoms,'atom_slot_profiles':[r[0] for r in profiles],'order_loss_hostile':{'first':left,'second':right,'H':3},'diagonal_exponential_coefficients':exponential,'literal_controls':records}
    save('ordinal_chart_audit.json.gz',out);print(json.dumps(out),flush=True)

if __name__=='__main__':run()
