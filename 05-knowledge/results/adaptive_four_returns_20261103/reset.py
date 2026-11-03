"""Seek an actual second complement for the sole missing one-move branch."""
from search import *


def run():
    T=bank()[1465];initial=profile(T);n=len(T);wanted=4397
    priority=tuple(sorted(range(n),key=lambda v:(initial['root_defects'][v],v)))
    candidates=[]
    for p in paths(T,priority):
        U=native.native(T,p);h=tables(U)[0][-1]
        candidates.append((h,p,U));tables.cache_clear()
    candidates.sort()
    print('first-complement census',len(candidates),'minimum',candidates[0][0],flush=True)
    tested=0
    for rank,(h,p,U) in enumerate(candidates):
        print('outer',rank,h,flush=True)
        up=profile(U);pr=tuple(sorted(range(n),key=lambda v:(up['root_defects'][v],v)))
        for q in paths(U,pr):
            V=native.native(U,q);hi=tables(V)[0][-1];tested+=1
            if hi<=wanted:
                for mask,value in enumerate(ear_values(V)):
                    if value==wanted:
                        W=birth(V,mask)
                        if strong(W):
                            result=dict(source=T,H=1465,first_path=p,first_target=U,first_H=h,
                                second_path=q,second_target=V,second_H=hi,mask=mask,target=W,
                                profile=profile(W),outer_rank=rank,second_programs_examined=tested,
                                first_census=len(candidates),first_minimum=candidates[0][0])
                            (HERE/'reset.json').write_text(json.dumps(result,indent=2)+'\n')
                            print('found',wanted,'via',h,hi,'mask',mask,flush=True)
                            return result
            tables.cache_clear()
    raise RuntimeError('two-complement language exhausted without witness')


if __name__=='__main__':run()
