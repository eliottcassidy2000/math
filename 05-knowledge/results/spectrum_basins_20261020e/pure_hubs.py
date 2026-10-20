"""Independent small replay of inherited pure-blue hubs and circulant witnesses."""
from itertools import permutations
from pathlib import Path
import json
from outer_probe import G

C3=(2,4,1)
C5=tuple(sum(1<<v for v in range(5) if (v-u)%5 in (1,2)) for u in range(5))
PATHS=[(0,1,2,3,4),(0,2,3,4,1),(0,2,4,1,3)]
BASE=[((26,20,9,18,4),(3,1,4,2,0),(0,1,2,3,4)),
      ((14,20,24,18,1),(1,2,3,4,0),(0,2,3,1,4)),
      ((22,12,24,17,2),(2,4,1,3,0),(1,3,0,2,4))]


def rigid_words(weight):
    if weight==0:yield ()
    for k in (1,3):
        if k<=weight:
            for w in rigid_words(weight-k):yield (k,)+w


def hp_count(a):
    return sum(G.is_path(a,p) for p in permutations(range(len(a))))


def pure_hubs():
    count=routes=0
    for weight in range(1,9):
        for word in rigid_words(weight):
            A=G.join([(0,) if k==1 else C3 for k in word]);r=len(A)
            U=G.join([A,C5,G.transpose(A)]);n=len(U)
            hub=list(U)
            for a,b in ((0,1),(0,n-1)):
                hub[a]^=1<<b;hub[b]^=1<<a
            hub=tuple(hub);profiles=[]
            for center in PATHS:
                p=tuple(range(r))+tuple(r+v for v in center)+tuple(range(r+5,n))
                G.need(G.is_path(U,p) and G.blue(U,p),'blue source path')
                S=G.complement(U,p);q=p[1:]+p[:1]
                G.need(G.is_path(S,q) and not G.blue(S,q),'black rotation')
                G.need(G.complement(S,q)==hub,'common labelled hub')
                ds=[x.bit_count() for x in S]
                G.need(set(ds[r:r+5])=={r+2} and max(ds[:r])<=r
                       and min(ds[r+5:])>=r+4,'intrinsic central score shell')
                profiles.append(hp_count(G.induced(S,tuple(range(r,r+5)))))
                routes+=1
            G.need(profiles==[11,15,13],'three distinct neighbor signatures')
            count+=1
    common=(26,20,25,18,0)
    for p,(expected,q,rename) in zip(PATHS,BASE):
        source=G.induced(C5,p)
        partner=G.complement(source,tuple(range(5)))
        G.need(partner==expected,'base partner matrix')
        G.need(G.is_path(partner,q) and not G.blue(partner,q),'base black path')
        marked=G.induced(partner,q)
        target=G.complement(marked,tuple(range(5)))
        G.need(G.induced(target,rename)==common,'base common hub')
    G.need(hp_count(common)==5,'base hub path count')
    return dict(outer_weights='1..8',outer_words=count,black_hub_routes=routes,
                empty_outer_routes=3,common_order_five_hub=common)


def circulant_witnesses():
    rows=[]
    for n in range(5,16,2):
        m=n//2;count=threshold=0
        for mask in range(1<<(m-1)):
            S={1}
            S.update(d if mask>>(d-2)&1 else n-d for d in range(2,m+1))
            B=tuple(sum(1<<v for v in range(n) if (v-u)%n in S) for u in range(n))
            pair=next(((t,k) for t in range(2,n-1) for k in range(t,n-1)
                       if t in S and k not in S),None)
            if pair:
                t,k=pair;u=k-t+2
                path=tuple(range(u,k+1))+(0,)+tuple(range(1,u))+tuple(range(k+1,n))
            else:
                threshold+=1
                G.need(S=={1}|set(range(m+1,2*m)),'threshold pattern')
                if m==2:path=(0,1,4,2,3)
                else:path=(0,1)+tuple(range(2,m))+(2*m,)+tuple(range(m,2*m))
            G.need(G.is_path(B,path),'constructed circulant word not directed')
            G.need(G.induced(B,path)!=B,'constructed word is automorphic')
            count+=1
        G.need(threshold==1,'threshold set count')
        rows.append(dict(order=n,marked_step_sets=count,threshold_sets=threshold))
    return rows


def main():
    result=dict(inherited_pure_hubs=pure_hubs(),
                elementary_rigidity_witnesses=circulant_witnesses())
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


if __name__=='__main__':main()
