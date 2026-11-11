"""All-length target saturation for the one-apex positive four-register grammar."""
from collections import deque
import json
from pathlib import Path

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def step(s,b,cap):
    h,o,e=s
    return (min(cap,h+2*o*(1-b)),min(cap,o+e+b),min(cap,e+o))

def compile_target(c):
    cap=c+1; start=(1,0,0); queue=deque([start]); pred={start:None}; edge={}
    witness=None
    while queue:
        s=queue.popleft()
        if s[0]==c and witness is None: witness=s
        for b in (0,1):
            t=step(s,b,cap)
            if t[0]>c: continue
            if t not in pred: pred[t]=s;edge[t]=b;queue.append(t)
    word=None
    if witness is not None:
        bits=[];s=witness
        while pred[s] is not None: bits.append(edge[s]);s=pred[s]
        word=bits[::-1]
    return dict(target=c,states=len(pred),word=word,complete=True)

def rows(word):
    # Apex 0, increasing chain vertices 1..m.
    n=len(word)+1; out=[0]*n
    for i in range(1,n):
        for j in range(i+1,n): out[i]|=1<<j
        if word[i-1]:out[0]|=1<<i
        else:out[i]|=1
    return out

def paths(out):
    n=len(out); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not dp[mask][v]:continue
            todo=out[v]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;dp[mask|bit][bit.bit_length()-1]+=dp[mask][v]
    return sum(dp[-1])

def zeckendorf(c):
    f=[1,2]
    while f[-1]<=c:f.append(sum(f[-2:]))
    digits=[0]*len(f);left=c
    for i in range(len(f)-1,-1,-1):
        if f[i]<=left:digits[i]=1;left-=f[i]
    while len(digits)>1 and digits[-1]==0:digits.pop()
    return digits

def digit_state(digits):
    # LSD first: previous digit, next position mod3, parity.
    p,r,z=0,0,0
    for b in digits:
        require(b in (0,1) and not(p and b),'Zeckendorf legality')
        z^=b*(r!=1);p=b;r=(r+1)%3
    return [p,r,int(z)]

def main():
    result=[]
    for c in range(1,64,2):
        row=compile_target(c);d=zeckendorf(c);row['zeckendorf_lsd']=d
        require(digit_state(d)[2]==c%2,'parity automaton')
        if row['word'] is not None:require(paths(rows(row['word']))==c,'independent witness count')
        result.append(row)
    controls=0
    for n in range(9):
        for mask in range(1<<n):
            w=[(mask>>i)&1 for i in range(n)];s=(1,0,0)
            for b in w:s=step(s,b,10**12)
            require(s[0]==paths(rows(w)),'all short words');controls+=1
    for c in (7,13,21):require(compile_target(c)['word'] is None,'hostile target')
    require(compile_target(23)['word'] is not None,'positive target')
    out={'status':'PASS','all_word_independent_checks':controls,'targets':result,
         'universe':'All binary words of arbitrary finite length for each target; finite saturation closure exhausted',
         'mathematical_scope':'One-apex target decision, not arbitrary tournament coverage'}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','checks':controls,'total_states':sum(r['states'] for r in result),'misses':[r['target'] for r in result if r['word'] is None]}))
if __name__=='__main__':main()
