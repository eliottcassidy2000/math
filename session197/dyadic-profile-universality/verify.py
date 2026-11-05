import hashlib
import json
import random

def need(condition,message):
    if not condition:
        raise RuntimeError(message)

def constants(L):
    need(L>=1,'positive integer gadget')
    if L==1:return 1,1,1
    if L%2==0:
        g,Q,B=constants(L//2)
        return g,Q,2*B
    ell=(L+1)//2
    g,Q,B=constants(ell)
    return g+ell-ell%2,Q+B,2*B

def f(n):
    value=1
    while n:
        if n%2==0:value+=n//2+(n//2)%2
        n//=2
    return value

def evaluate(word):
    n=0;value=1
    for bit in word:
        if bit=='0':value+=n+n%2
        n=2*n+int(bit)
    return n,value

def revc(word):
    return int(''.join('1' if b=='0' else '0' for b in word[::-1]),2) if word else 0

def construct(k,u,z,q,extra=0):
    M=1<<k
    need(k>=1 and z%2 and q%2,'valid dyadic residues')
    target=format((q-1)%M,f'0{k}b')
    prefix=''.join('1' if b=='0' else '0' for b in target[::-1])
    need(prefix[0]=='1','minimal-width leading one')
    word=prefix+'0'*k+('1'+'0'*k)*extra
    _,initial=evaluate(word)
    r=(u-1)%M
    j=((initial+f(r)-1-z)//2)%(M//2)
    for i in range(k-1):
        if j>>i&1:
            if i:word+='1'*(i+1)+'0'*k
            word+='11'+'0'*k
    return word+format(r,f'0{k}b')

def graph(L):
    n=L-1;s=n.bit_length();rows=[0]*(s+1)
    for i in range(s):
        for j in range(i+1,s):rows[i]|=1<<j
        if n>>i&1:rows[i]|=1<<s
        else:rows[s]|=1<<i
    return rows,s

def seed_graph(L):
    G,b=graph(L);rows=[(((1<<len(G))-1)^(1<<b))<<1]
    rows.extend(row<<1 for row in G)
    rows[b+1]|=1
    return rows,b+1

def profile(rows,root):
    n=len(rows);size=1<<n;full=size-1
    ends=[[0]*n for _ in range(size)]
    starts=[[0]*n for _ in range(size)]
    for v in range(n):ends[1<<v][v]=starts[1<<v][v]=1
    for S in range(1,size):
        for v in range(n):
            if ends[S][v]:
                out=rows[v]&~S
                while out:
                    bit=out&-out;out-=bit;w=bit.bit_length()-1
                    ends[S|bit][w]+=ends[S][v]
            if starts[S][v]:
                unseen=full^S
                while unseen:
                    bit=unseen&-unseen;unseen-=bit;w=bit.bit_length()-1
                    if rows[w]>>v&1:starts[S|bit][w]+=starts[S][v]
    H=[sum(a) for a in ends];H[0]=1
    need(all(sum(starts[S])==H[S] for S in range(1,size)),'independent endpoint directions')
    R=sum(starts[S][root]*H[full^S] for S in range(size) if S>>root&1)
    L=sum(ends[S][root]*H[full^S] for S in range(size) if S>>root&1)
    return H[full],ends[full][root],starts[full][root],R,L,H[full^(1<<root)]

def sink(rows):
    n=len(rows)
    return [row|(1<<n) for row in rows]+[0]

def main():
    digest=hashlib.sha256();count=large=neutral=literal=padding=0;physical=set(range(2,130));examples=[]
    rng=random.Random(197)
    for k in range(1,10):
        M=1<<k
        triples=[(u,z,q) for u in range(M) for z in range(1,M,2) for q in range(1,M,2)] if k<=4 else [(rng.randrange(M),2*rng.randrange(M//2)+1,2*rng.randrange(M//2)+1) for _ in range(400)]
        for u,z,q in triples:
            word=construct(k,u,z,q);n,h=evaluate(word)
            need(n==int(word,2) and h==f(n),'two independent f evaluations')
            g,Q,B=constants(n+1)
            need((n+1)%M==u and g%M==z and Q%M==q and B%M==0,'joint physical profile residues')
            need(g==h and Q==1+revc(word) and B==1<<len(word),'physical response recurrences')
            need(len(word)<=4*k*k+3*k,'quadratic witness bound')
            other=construct(k,u,z,q,2);nn,hh=evaluate(other)
            need(nn>n and (nn+1)%M==u and hh%M==z and (1+revc(other))%M==q,'arbitrarily large neutral repetition')
            count+=1;neutral+=1
            digest.update(f'{k},{u},{z},{q},{n},{h},{Q}\n'.encode())
            if len(word)<=10:physical.add(n+1)
            if k==2:examples.append({'L':n+1,'g':g,'Q':Q,'B':B,'requested':[u,z,q],'word':word})
    # Very large exact integers; avoid recursion depth dependence in constants.
    for k in (16,32,64):
        M=1<<k
        for u,z,q in ((1,7,13),(M-1,M-1,M-3),(M//2,12345,54321)):
            word=construct(k,u,z,q);n,h=evaluate(word)
            need(h==f(n) and (n+1)%M==u%M and h%M==z%M and (1+revc(word))%M==q%M,'large exact profile residues')
            need(len(word)<=4*k*k+3*k,'large quadratic length bound')
            large+=1
    for L in sorted(physical):
        G,b=graph(L);T,r=seed_graph(L);g,Q,B=constants(L)
        need(profile(G,b)==(g,1,L%2,Q,L,1),'literal gadget response')
        h,e,s,R,_,d=profile(T,r)
        need((h,e,s,R,d)==(g+L-L%2,1,1,Q+B,1),'literal strong start-one response')
        order=list(range(len(T)))
        need(all(T[order[i]]>>order[(i+1)%len(T)]&1 for i in range(len(T))),'literal Hamilton cycle')
        need(all(T[i]>>j&1 for i in range(r) for j in range(i+1,r)),'literal transitive root deletion')
        if len(T)<=9:
            for t in (1,2):
                U=T
                for _ in range(t):U=sink(U)
                hh,ee,ss,RR,_,dd=profile(U,r)
                need((hh,ee,ss,RR,dd)==(h,0,1,(1<<t)*R,1),'literal retained-root sink transfer')
                padding+=1
        literal+=1
    print(json.dumps({'status':'PASS','joint_residue_triples':count,'neutral_extension_checks':neutral,'large_exact_triples':large,'literal_gadget_and_strong_pairs':literal,'literal_sink_controls':padding,'joint_witness_digest':digest.hexdigest(),'k2_examples':examples,'scope':'Simultaneous dyadic residues of actual physical responses; no exact-target or prime infinitude claim.'},sort_keys=True,indent=2))

if __name__=='__main__':main()
