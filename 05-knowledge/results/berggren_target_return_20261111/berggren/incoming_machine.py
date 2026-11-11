"""Independent literal checks of the restricted two-link N/Z machine."""
from itertools import product
from pathlib import Path
import hashlib,json
from check import need,H,birth,root_profile,induced,strong

HERE=Path(__file__).resolve().parent
INCOMING=HERE.parent.parent/'odd_tables_signed_carries_20261111/history/NATIVE_ZERO.md'


def predecessor(T,r):
    incoming=[v for v in range(len(T)) if T[v]>>r&1]
    need(len(incoming)==1,'marked vertex has exactly one predecessor')
    return incoming[0]


def two_link(T,r):
    a=predecessor(T,r);p=predecessor(T,a)
    need(p not in(a,r),'three distinct interface vertices')
    return a,p


def triple(T,r):
    p=root_profile(T,r)
    return(p['H'],p['F']-p['H'],p['S'])


def step(T,r,letter):
    a,_=two_link(T,r)
    parent=r if letter=='N' else a
    U=birth(birth(T,parent),len(T));q=len(U)-1
    need(induced(U,list(range(len(T))))==T,'every older arc retained')
    two_link(U,q)
    if strong(T):need(strong(U),'strongness retained')
    return U,q


def expected(v,letter):
    h,b,s=v
    return(3*h+2*b,h-s,h+b) if letter=='N' else(3*h,2*s,h)


def tournament(n,mask):
    T=[0]*n;bit=0
    for u in range(n):
        for v in range(u+1,n):
            if mask>>bit&1:T[u]|=1<<v
            else:T[v]|=1<<u
            bit+=1
    return tuple(T)


def main():
    seeds=[tournament(n,mask) for n in range(1,4) for mask in range(1<<(n*(n-1)//2))]
    seeds.extend([(14,12,8,0),(2,4,9,3)])
    controls=[]
    for X in seeds:
        for r in range(len(X)):
            Y=birth(birth(X,r),len(X));root=len(Y)-1;v=triple(Y,root)
            need(v[0]-2*v[2]==H(X)>=1,'literal deleted-two-vertex identity')
            for length in(1,2):
                for letters in product('NZ',repeat=length):
                    T,q,state=Y,root,v
                    for letter in letters:
                        before=state;T,q=step(T,q,letter);state=triple(T,q)
                        need(state==expected(before,letter),'independent exact marked recurrence')
                        need(state[0]-2*state[2]==before[0],'strict two-link margin')
                        if before[2]>=2 and before[1]>=0:
                            need(state[1]>=3 and state[2]>=2,'forward-invariant excluded-ready region')
                    controls.append(dict(order=len(Y),source=v,word=''.join(letters),target=state))
    boundary=birth(birth((0,),0),1)
    need(triple(boundary,2)==(3,0,1),'necessary s>=2 boundary source')
    for letter in'NZ':
        U,r=step(boundary,2,letter)
        need(triple(U,r)==(9,2,3),'s=1 really can return b2')
    data=json.loads((HERE.parent/'reset/target_ready_return.json').read_text())
    actual=[]
    for rec in data['controls'][:2]:
        k=rec['k'];D=1<<k;Q=tuple(rec['ready665'])
        Y=birth(birth(Q,len(Q)-1),len(Q));r=len(Y)-1;v=triple(Y,r)
        need(v==(891*D+215,192*D+47,297*D+72),'actual all-k starting formula controls')
        for letter in'NZ':
            T,q=step(Y,r,letter);observed=triple(T,q)
            need(observed==expected(v,letter),'actual physical N/Z continuation')
            actual.append(dict(k=k,operation=letter,source=v,target=observed,rows=T,root=q))
    paths=[INCOMING,HERE.parent/'reset/TARGET_RETURN.md']
    out=dict(status='PASS',scope='All rooted labeled X of orders1..3 plus rooted TT4 and S4; both operations and all two-letter words; actual Q_k endpoint controls k1,2. No alternate-root or core-replacement claim.',
             literal_control_count=len(controls),controls=controls,actual=actual,
             boundary=dict(source=[3,0,1],N=[9,2,3],Z=[9,2,3]),
             source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (HERE/'incoming_machine.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('PASS: restricted two-link N/Z recurrence, invariant ready exclusion, actual Q1/Q2 controls, and s=1 hostile boundary')


if __name__=='__main__':main()
