import importlib.util
import itertools
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('inherited_compiler',SOURCE)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
C3=(2,4,1)

def require(p,msg):
    if not p:raise RuntimeError(msg)

def join(*parts):
    n=sum(map(len,parts));out=[];shift=0
    for part in parts:
        after=((1<<n)-1)^((1<<(shift+len(part)))-1)
        out.extend((x<<shift)|after for x in part);shift+=len(part)
    return out

def positive(sig):
    d,kp,km,ip,im,op,om,jp,jm=sig
    require(km==-d,'negative antipode')
    return (d+jp+jm,d,d+kp,ip,op,-im,-om)

def gap_state(r,M,k,l):
    Z=1<<r;R=Z-1-int(format(M,f'0{r}b')[::-1],2)
    alpha=M&1;beta=R&1;delta=lambda x:int(x in (1,2))
    h=9+6*(delta(k)+delta(l))+9*(c.g(M)-1)+9*k*R+3*k*beta+9*(3-l)*M+9*k*(3-l)*Z+3*alpha*(3-l)
    return (h,9,144*Z,36*(3-l)*Z+36*R+9*(3-k),36*k*Z+36*M+9*l,3*(3-k),3*l)

def count(s,A,B):
    h,d,L,i,o,start,end=s
    return h+d*(c.g(A)-1+c.g(B)-1)+L*A*B+i*A+o*B+start*(A&1)+end*(B&1)

def main():
    states=0;graphs=0;separators=0;first_quotients={}
    for r in range(1,5):
        middle=tuple(sum(1<<j for j in range(i+1,r)) for i in range(r))
        rows=join(C3,middle,C3);qs=set()
        for M in range(1<<r):
            middle_mask=sum(((M>>(r-1-j))&1)<<j for j in range(r))
            for k,l in itertools.product(range(4),repeat=2):
                O=(1<<k)-1 | (middle_mask<<3) | (((1<<l)-1)<<(3+r))
                sig=c.signature(rows,O);want=gap_state(r,M,k,l)
                require(positive(sig)==want,'middle-word observer')
                require(all(x%3==0 for x in want),'coefficient divisibility')
                qs.add(want[0]//3);states+=1
                tests=[(0,0)] if r<=2 else []
                if r==1:tests.extend([(0,1),(2,0),(1,2)])
                for A,B in tests:
                    got=c.hp_subsets(c.construct(rows,O,A,B))[-1]
                    require(got==count(want,A,B),'actual flank child')
                    graphs+=1
        first_quotients[str(r)]=sorted(qs)
    for n in range(1,4):
        for code in range(1<<(n*(n-1)//2)):
            D=c.tournament(n,code);rows=join(C3,D,C3)
            for M in range(1<<n):
                for k,l in itertools.product(range(4),repeat=2):
                    O=((1<<k)-1)|(M<<3)|(((1<<l)-1)<<(n+3))
                    state=positive(c.signature(rows,O))
                    require(all(x%3==0 for x in state),'arbitrary separator')
                    separators+=1
    rows=join(C3,C3);O=1|(1<<3)
    hostile=c.hp_subsets(c.construct(rows,O,0,0))[-1]
    require(hostile==41,'empty-gap hostile')
    return {'status':'PASS','middle_word_states':states,
            'literal_child_counts':graphs,'arbitrary_separator_states':separators,
            'empty_gap_hostile':hostile,'no_flank_quotients_by_gap':first_quotients,
            'scope':'r=1..4 all words and triangle cuts; labelled separator orders1..3 all cuts. General theorem uses inherited join/capacity laws.'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
