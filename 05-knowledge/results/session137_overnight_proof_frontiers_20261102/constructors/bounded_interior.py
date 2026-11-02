#!/usr/bin/env python3
"""Complete small-source consumer for fixed-height/interior dyadic normal form."""
import json
from check import tournament,state,pos,HP,E,chain,need,induced,S0,multiply
from face_bases import sccs

def normalize(A,cut):
    """Delete neutral prefix singleton SCCs; move suffix singletons into S0 tail."""
    need(cut>0,'nonempty cut required')
    parts=sccs(A)
    hits=[i for i,P in enumerate(parts) if any(cut>>v&1 for v in P)]
    first,last=hits[0],hits[-1]
    kept=[];tail=0;r=0
    for i,P in enumerate(parts):
        if len(P)>1 or first<=i<=last:kept.extend(P)
        elif i>last:tail+=1
        if first<i<last and len(P)==1:r+=1
    kept.sort();B=induced(A,kept);newcut=sum(1<<j for j,v in enumerate(kept) if cut>>v&1)
    return B,newcut,tail,r,first==last

def main():
    tests=0;reencoded=0;maxL=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(n,code)
            for cut in range(1,1<<n):
                v=state(A,cut);C,L,g,t=pos(v);h=v[0]
                B,newcut,k,r,same=normalize(A,cut);base=state(B,newcut);w=base
                for j in range(k):w=multiply(S0,w)
                need(w==v,'SCC normalization changed actual state')
                need(pos(base)[1]==L and base[0]==h,'height/interior changed')
                bound=h+1 if L==0 else h+3+(L.bit_length()-1)
                need(len(B)<=bound,'base order bound')
                if L and not same and r:
                    need(2**(r-1)<=L,'even-path interior bound')
                if k:
                    need(v[1]==2**k*base[2]+h-2*L and v[3]==0,'dyadic output formula')
                    # Count the reencoded graph with k final unselected singleton blocks.
                    T=B;S=newcut
                    for j in range(k):
                        m=len(T);T=tuple(x|(1<<m) for x in T)+(0,)
                    need(state(T,S)==v,'literal suffix reencoding')
                    reencoded+=1
                tests+=1;maxL=max(maxL,L)
    # Every transitive base at h=1,L=1 up to the universal bound4.
    raw=[];terminal=[]
    for n in range(1,5):
        A=chain(n)
        for cut in range(1,1<<n):
            v=state(A,cut)
            if pos(v)[1]!=1:continue
            B,S,k,r,same=normalize(A,cut);base=state(B,S)
            raw.append((v[1],v[2],v[3]))
            terminal.append(base)
    actual=sorted(set(terminal));need(actual==[(1,1,3,1),(1,3,5,1)],'unit-height L1 base list')
    # Complete transitive cut control through n10 for the all-order scalar language.
    def belongs(F):
        if F in (1,3):return True
        for q in (3,5):
            z=F+1
            if z%q==0:
                u=z//q
                if u>=2 and not(u&(u-1)):return True
        return False
    L1=0
    for n in range(1,11):
        A=chain(n)
        for cut in range(1<<n):
            v=state(A,cut)
            if pos(v)[1]==1:
                need(belongs(v[1]),'fixed L1 output outside exact rays');L1+=1
    # Independently check the two empty-cut suffix operators commute on actual states.
    commuting=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            B=tournament(n,code);a=HP(B);e=E(B)
            M=((a,0,0,0),(0,a,e-a,a),(0,0,e,0),(0,0,0,0))
            for A,S in [(chain(1),0),(chain(1),1),(chain(2),1),(chain(2),3)]:
                v=state(A,S)
                need(multiply(M,multiply(S0,v))==multiply(S0,multiply(M,v)),'suffix commutation')
                commuting+=1
    print(json.dumps({'status':'PASS','complete_labelled_parent_orders':[1,2,3,4,5],
       'all_nonempty_cut_normalizations':tests,'literal_reencoded_suffix_controls':reencoded,
       'maximum_interior_L_in_universe':maxL,'suffix_operator_commutation_controls':commuting,
       'h1_L1_canonical_native_bases':actual,'h1_L1_transitive_cut_controls_through_order10':L1,
       'h1_L1_exact_all_order_H_language':'{1,3} union {3*2^k-1,5*2^k-1:k>=1}',
       'scope':'All-order statements use the proved fixed-(h,L) normal form; exact finite universes above are implementation controls.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
