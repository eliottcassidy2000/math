"""Exact closed six-coordinate compiler for the inherited one-start grammar."""
from pathlib import Path
from functools import lru_cache
import argparse,json,sys,time
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'two_entry_funnels_20261103'))
import research as a
from verify import rooted_data

SINGLE=(1,1,1,2,1,1) # H,E,S(2),P(2),E(2),S(3)

@lru_cache(None)
def coefficients(word):
    one=a.tail_coefficients(word,1)
    two=a.tail_coefficients(word,2)
    three=a.tail_coefficients(word,3)
    return one,two,three

def advance(state,word):
    h,e,s2,p2,e2,s3=state;t=len(word)
    (l1,d1,b1),(l2,d2,b2),(l3,d3,_)=coefficients(word)
    return (l2*s2+b1*h-d1*e,
            l1*h-(0 if t else 1)*e,
            2**t*s2+d2*h,
            2*(l3*s3+b2*p2-d2*e2),
            l2*p2-e2,
            3**t*s3+d3*p2)

def decode(n,index):
    counts=a.class_counts(n)
    for m in [1]+list(range(3,n)):
        t=n-m-1;c=1 if t==0 else 1<<(t-1)
        core_count=1 if m==1 else counts[m]
        if index>=c*core_count:index-=c*core_count;continue
        word_index,rank=divmod(index,core_count)
        word=() if t==0 else tuple((word_index>>i)&1 for i in range(t-1))+(1,)
        return m,rank,word
    raise RuntimeError('decode exhausted')

def compile_rank(n,index):
    m,rank,word=decode(n,index)
    return advance(SINGLE if m==1 else compile_rank(m,rank),word)

def direct_state(T):
    p,e,s=a.bin_values(T,0,3)
    return (p[1],e[1],s[2],p[2],e[2],s[3])

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,default=15)
    parser.add_argument('--controls',type=int,default=7);args=parser.parse_args()
    targets={221:99,665:297,1997:891,5993:2673,5995:2673,667:297,1999:891,2003:891,2005:891,6013:2673}
    banks={1:[SINGLE]};hits={str(h):[] for h in targets};total=0;checks=0;t0=time.time()
    for n in range(3,args.order+1):
        states=[]
        for m in [1]+list(range(3,n)):
            t=n-m-1;words=[()] if t==0 else [tuple((j>>i)&1 for i in range(t-1))+(1,) for j in range(1<<(t-1))]
            for word in words:
                states.extend(advance(state,word) for state in banks[m])
        a.require(len(states)==a.class_counts(n)[n],'complete rank-ordered class count')
        for rank,state in enumerate(states):
            h=state[0]+state[4]-1
            if h in targets:hits[str(h)].append(dict(order=n,rank=rank,state=state))
            if n<=args.controls:
                a.require(state==direct_state(a.unrank(n,rank)),'independent complete six-state subset controls')
                a.tables.cache_clear();checks+=1
        banks[n]=states;total+=len(states)
        print('order',n,'classes',len(states),'hits',{h:len(v) for h,v in hits.items()},'elapsed',round(time.time()-t0,2),flush=True)
    data=dict(status='PASS',universe='Every strong unique-start rooted class through stated order',
              max_order=args.order,classes=total,independent_subset_controls=checks,hits=hits)
    (HERE/'funnel_six_state.json').write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':main()
