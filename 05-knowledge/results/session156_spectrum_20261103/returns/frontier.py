"""Golden complement transfer: exact nonnegative path-component recurrence."""
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
# Most recent endpoint is 0; -1 is older than the retained interface.
E=()
STATES=(E, ((0,0),), ((-1,-1),), ((-1,0),), ((0,-1),),
        ((-1,-1),(0,-1)), ((-1,-1),(0,0)))
NAMES=("empty","single vertex","internal port","end port","start port",
       "two paths: port starts nontrivial","two paths: port isolated")
LETTERS={"s":(1,), "c":(1,1,3), "q0":(1,1,3,7), "q1":(1,1,3,5),
         "q2":(1,1,3,3), "q4":(1,1,1,7), "q5":(1,1,1,5)}
def need(x,msg):
    if not x:raise RuntimeError(msg)
def canonical(comps,window):
    def old(x):return -1 if x>=window else x
    result=tuple(sorted((old(a),old(b)) for a,b in comps))
    if sum(a==-1 for a,b in result)>1:return None
    return result
@lru_cache(None)
def step(state,backmask,window=3):
    # New vertex can take one incoming old tail and one outgoing old head.
    ans=Counter()
    incoming=[None]+[i for i,(a,b) in enumerate(state) if b==-1 or not (backmask>>b&1)]
    outgoing=[None]+[i for i,(a,b) in enumerate(state) if a>=0 and backmask>>a&1]
    for i,j in product(incoming,outgoing):
        if i is not None and i==j:continue  # closes a directed cycle
        advance=lambda x:-1 if x==-1 else x+1
        start=0 if i is None else advance(state[i][0])
        end=0 if j is None else advance(state[j][1])
        comps=[(advance(a),advance(b)) for k,(a,b) in enumerate(state) if k not in (i,j)]
        comps.append((start,end))
        result=canonical(comps,window)
        if result is not None:ans[result]+=1
    return ans
@lru_cache(None)
def block(state,letter):
    states=Counter({state:1})
    for mask in LETTERS[letter]:
        nxt=Counter()
        for s,w in states.items():
            for t,c in step(s,mask).items():nxt[t]+=w*c
        states=nxt
    ans=Counter()
    for s,w in states.items():
        t=canonical(s,1)
        if t is not None:ans[t]+=w
    need(all(t in STATES for t in ans),("boundary closure",letter,ans))
    return ans
def matrices():
    return {letter:[[block(s,letter).get(t,0) for t in STATES] for s in STATES]
            for letter in LETTERS}
M=matrices()
R=((1,0,0),(0,1,1),(1,0,0),(0,1,0),(1,0,1),(0,0,1),(0,0,1))
def product_matrix(a,b):
    return tuple(tuple(sum(row[k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for row in a)
SMALL={l:tuple(product_matrix(m,R)[i] for i in (2,3,5)) for l,m in M.items()}
for l in LETTERS:need(product_matrix(M[l],R)==product_matrix(R,SMALL[l]),('exact quotient',l))
def apply3(v,mat):
    a,b,c=v
    return tuple(a*mat[0][j]+b*mat[1][j]+c*mat[2][j] for j in range(3))
def vector3(word):
    v=(1,0,0)
    for letter in word:v=apply3(v,SMALL[letter])
    return v
def compile_descriptor(desc):
    matrix=[];window=max(1,len(desc)-1)
    for state in STATES:
        states=Counter({state:1})
        for mask in desc:
            nxt=Counter()
            for s,w in states.items():
                for t,c in step(s,mask,window).items():nxt[t]+=w*c
            states=nxt
        row=Counter()
        for s,w in states.items():
            t=canonical(s,1)
            if t is not None:row[t]+=w
        need(all(s in STATES for s in row),('general boundary',desc))
        matrix.append(tuple(row.get(s,0) for s in STATES))
    mr=product_matrix(matrix,R);small=tuple(mr[i] for i in (2,3,5))
    need(mr==product_matrix(R,small),('general block quotient',desc))
    return small
def multiply(v,letter):
    mat=M[letter]
    return tuple(sum(v[i]*mat[i][j] for i in range(len(v))) for j in range(len(v)))
def vector(word):
    v=(1,0,0,0,0,0,0)
    for letter in word:v=multiply(v,letter)
    return v
def h(word):
    return sum(vector(word)[:5])  # empty path convention included
def targets(word):
    # word is in PROCESSING order: reverse of source SCC order.
    masks=[m for letter in word for m in LETTERS[letter]]
    rows=[0]*len(masks)
    for v,mask in enumerate(masks):
        for u in range(v):
            back=(v-u<=3 and mask>>(v-u-1)&1)
            rows[v if back else u]|=1<<(u if back else v)
    return tuple(rows)
def subset_h(rows):
    n=len(rows)
    if not n:return 1
    dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,count in dp[mask].items():
            options=rows[v]&~mask
            while options:
                bit=options&-options;options-=bit;w=bit.bit_length()-1
                dest=dp[mask|bit];dest[w]=dest.get(w,0)+count
    return sum(dp[-1].values())
def words(n):
    if n==0:yield ()
    for letter,desc in LETTERS.items():
        if len(desc)<=n:
            for rest in words(n-len(desc)):yield rest+(letter,)
def main():
    checked=0
    for n in range(11):
        for word in words(n):
            need(h(word)==subset_h(targets(word)),("direct count",word,h(word),subset_h(targets(word))))
            checked+=1
    report={"states":NAMES,"matrices":M,"source_costs":{l:len(s) for l,s in LETTERS.items()},
            "minimal_three_matrices":SMALL,"quotient_matrix":R,
            "verified_through_order":10,"charts_checked":checked,
            "singleton_counts":[h(("s",)*n) for n in range(21)],
            "five_q4_responses":{l:{"isolated":h((l,)),"s_left":h(("s",l)),
                                "s_right":h((l,"s")),"both":h(("s",l,"s"))} for l in LETTERS if l[0]=="q"}}
    (HERE/'transfer.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','charts':checked,'three_matrices':SMALL,
                      'q4_responses':report['five_q4_responses']}))
if __name__=="__main__":main()

