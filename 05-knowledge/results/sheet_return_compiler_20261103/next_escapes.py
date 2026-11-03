"""Positive source-preserving continuation witnesses; not a closed bank."""
from bulk import target,need,HERE,tables,strong,native,birth
from escape_lifts import ear_values
from switch_search import hpaths
from return_program import field,opposite
import json


def profile(T,r=2):
    H,E=tables(T);_,S=tables(opposite(T));full=(1<<len(T))-1
    ends=E[-1][r];starts2=sum(S[mask][r]*H[full^mask] for mask in range(1,full+1))
    v=(H[-1],(starts2-ends)//2,(starts2+ends)//2)
    defects=[sum(E[mask][r]*H[full^mask] for mask in range(1,full+1))-H[-1] for r in range(len(T))]
    return dict(root=r,observer=v,phase=field(v),root_defects=defects)


def search(T):
    h=tables(T)[0][-1];wanted={3*h+2,3*h+4};found={};tested=0
    # A direct birth is a different and cheaper legal program.
    programs=[(None,T)]
    def choices():
        yield None,T
        for p in hpaths(T):yield p,native.native(T,p)
    for p,U in choices():
        tested+=1
        if tables(U)[0][-1]<=max(wanted):
            for mask,value in enumerate(ear_values(U)):
                if value in wanted-found.keys():
                    V=birth(U,mask)
                    if strong(V):
                        found[value]=dict(path=p,intermediate=U,intermediate_H=tables(U)[0][-1],mask=mask,target=V,profile=profile(V))
                        print('source',h,'target',value,'program',tested,'direct' if p is None else 'native',flush=True)
        tables.cache_clear()
        if set(found)==wanted:break
    return dict(source=T,H=h,profile=profile(T),programs_examined=tested,witnesses=found)


def run():
    rows=[]
    for kind in ('straight','revisit'):
        T,_=target(3,3,kind);row=search(T);row['kind']=kind;rows.append(row)
        (HERE/'next_escapes.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


if __name__=='__main__':run()
