"""Literal checks of the reflection exit on Session304's companion family."""
import json
from pathlib import Path
from mirror_module import inventory, opposite, need

def apex(w):
    n=len(w)
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]+[0]
    for i,b in enumerate(w):
        if b:rows[n]|=1<<i
        else:rows[i]|=1<<n
    return tuple(rows)

def reflect(w):
    return tuple(1-b for b in reversed(w))

def birth(T,r,retain):
    n=len(T);U=list(T);U[r]|=1<<n
    U.append(((1<<n)-1)^(1<<r))
    return tuple(U),r if retain else n

def state(T,r):
    H,E=inventory(T);_,S=inventory(opposite(T));mask=len(H)-1
    return (H[-1],S[-1][r],H[mask^(1<<r)],
            sum(E[m][r]*H[mask^m] for m in range(mask+1)))

def main():
    controls=[]
    for k in range(1,7):
        B=1<<k;w=(1,0,0)+(1,)*k+(0,0)
        T,U=apex(w),apex(reflect(w));h=15*B+3
        before,after=state(T,0),state(U,0)
        need(before==(h,3*B-1,3*B-1,h),'source full state')
        need(after==(h,5*B+4,5*B+5,h+1),'reflected full state')
        need(T[0]==U[0] and all((t&1)==(u&1) for t,u in zip(T,U)),
             'preserved complete root star')
        V,r=birth(U,0,False);V,r=birth(V,r,True)
        out=state(V,r);need(out==(3*h+2,h,h,3*h+2),'literal NR state')
        controls.append(dict(k=k,word=w,before=before,after=after,
                             output=out,rows=V,root=r))
    w=(1,0,1,0,0,1,0)
    hostile=[state(apex(x),0) for x in (w,reflect(w))]
    need(all(q-s==0 for h,s,q,f in hostile),'common fork is not mirror ready')
    receipt=dict(status='PASS',universe='k=1..6; both source mirrors and literal NR outputs; one common-fork guard hostile',
                 controls=controls,hostile=dict(word=w,states=hostile))
    Path(__file__).with_name('incoming_bridge_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS: six literal companion exits and the common-fork guard hostile')

if __name__=='__main__':main()
