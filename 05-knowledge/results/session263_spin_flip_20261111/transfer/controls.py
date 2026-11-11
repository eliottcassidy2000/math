#!/usr/bin/env python3
"""Exact tests of NEW word decoder / tomography; no inherited census audit."""
import itertools, json
from pathlib import Path

S0=(-1,3,1,1)
def step(s,a):
    r,M,h,E=s
    return (-r,a*(h+r)+(a*a+1)*M-E,r+(2*a-2)*M,a*M+r-h)
def guarded(s):
    r,M,h,E=s
    return r in (-1,1) and h>=1 and M>=2*h and 0<=E<=h and 0<=2*(h+r)<=M

def decode(s,initial=S0):
    word=[]
    while s!=initial:
        if not guarded(s) or s[1]<=initial[1]: raise ValueError('guard/initial bound')
        rp,Mp,hp,Ep=s; qp=hp+rp
        if 23*qp>=8*Mp and 19*qp<=8*Mp: a=3
        elif 57*qp>=16*Mp and 51*qp<=16*Mp: a=5
        else: raise ValueError('projective gap')
        if qp%(2*a-2): raise ValueError('nonintegral inverse')
        r=-rp; M=qp//(2*a-2); q=a*M+2*r-Ep
        E=a*q+(a*a+1)*M-Mp
        prev=(r,M,q-r,E)
        if not guarded(prev) or step(prev,a)!=s or 2*Mp<19*M:
            raise ValueError('inverse guard')
        word.append(a);s=prev
    return tuple(reversed(word))

def observe(s):
    return tuple(x[0]+2*x[1] for x in (s,step(s,3),step(s,5)))
def tomography(r,obs):
    A,B,C=obs
    assert (A-r)%2==0
    M=(A-r)//2
    assert (C-B-32*M)%4==0
    q=(C-B-32*M)//4
    assert (B+r)%2==0
    E=3*q+10*M-(B+r)//2
    return (r,M,q-r,E)

# Graphs use sets of outgoing vertices. All objects rebuilt literally.
def close(g,O,t):
    n=len(g); out=[set(x) for x in g]+[set() for _ in range(t+1)]
    v=n+t
    for x in range(n):
        out[x].update(range(n,n+t))
        if x in O: out[v].add(x)
        else: out[x].add(v)
    for x in range(n,n+t): out[x].update(range(x+1,v+1))
    return out

def physical(word):
    g=[{1},set()]; O={0,1}
    for a in word:
        n=len(g);t=2 if a==3 else 3
        g=close(g,O,t);O={n,n+t-1}
    return g,O

def all_ends(g):
    n=len(g);N=1<<n; incoming=[0]*n
    for u in range(n):
        for v in g[u]:incoming[v]|=1<<u
    dp=[[0]*n for _ in range(N)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,N):
        for v in range(n):
            if not mask>>v&1:continue
            old=mask^(1<<v)
            if not old:continue
            p=old&incoming[v];z=0
            while p:
                bit=p&-p;p-=bit;z+=dp[old][bit.bit_length()-1]
            dp[mask][v]=z
    return dp

def count(g):return sum(all_ends(g)[-1])
def recount(g,O):
    n=len(g);full=(1<<n)-1
    ed=all_ends(g)
    rev=[set() for _ in g]
    for u in range(n):
        for v in g[u]:rev[v].add(u)
    st=all_ends(rev)
    H=[sum(z) for z in ed];H[0]=1
    start=[sum(z[v] for v in O) for z in st]
    end=[sum(z[v] for v in O) for z in ed]
    M=sum(H[full^mask]*start[mask] for mask in range(1,full+1))
    W=sum(end[mask]*start[full^mask] for mask in range(1,full))
    return (H[full]-W,M,H[full],end[full])

def main():
    seen={}; obssets={};tested=0
    for k in range(11):
        obsset=set()
        for w in itertools.product((3,5),repeat=k):
            s=S0
            for a in w:
                old=s;s=step(s,a)
                assert guarded(s) and s[1]%2==s[3]%2==1
                assert 19*old[1]<=2*s[1]<=57*old[1]
            assert s not in seen
            seen[s]=w
            assert decode(s)==w
            obs=observe(s)
            assert tomography(s[0],obs)==s
            obsset.update(obs);tested+=1
        assert len(obsset)**3>=2**k
        obssets[k]=len(obsset)
    physical_receipts=[]
    words=[w for k in range(3) for w in itertools.product((3,5),repeat=k)]+[(3,3,3),(5,3,5)]
    for w in words:
        g,O=physical(w);s=S0
        for a in w:s=step(s,a)
        actual=recount(g,O);assert actual==s,(w,actual,s)
        # The initial and all length<=2 words get all three physical queries.
        values=[]
        if len(w)<=2:
            for suffix in ((),(3,),(5,)):
                G,P=physical(w+suffix)
                values.append(count(close(G,P,1)))
            assert tuple(values)==observe(s)
        physical_receipts.append({'word':w,'order':len(g),'state':s,'queries':values})
    rejected=0
    for s in [(1,29,11,8),(1,30,11,7),(-1,3,1,0),(1,11,1,1)]:
        try:decode(s)
        except ValueError:rejected+=1
        else:raise AssertionError(('forged state accepted',s))
    # Independent small determinant calculation, permutations only.
    C=[[1,2,0,0],[5,20,6,-2],[59,218,62,-20],[637,2368,674,-218]]
    det=0
    for p in itertools.permutations(range(4)):
        val=1
        for i in range(4):val*=C[i][p[i]]
        inv=sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
        det+=(-1)**inv*val
    assert det==16
    receipt={'status':'PASS','new_word_depth':10,'new_word_states':tested,
             'distinct_query_counts_by_depth':obssets,'physical_composites':physical_receipts,
             'hostile_rejections':rejected,'four_query_determinant':det,
             'scope':'New compiler/decoder/tomography controls only; not an inherited theorem or census re-audit.'}
    path=Path(__file__).with_name('receipt.json')
    path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','word_states':tested,'physical_composites':len(words),'hostile_rejections':rejected}))
if __name__=='__main__':main()
