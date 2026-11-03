"""Apply the inherited source-SCC recurrence to the classified quadratic graph."""
from functools import lru_cache
import json
from itertools import permutations


def count(F):
    n=len(F);E=list(enumerate(F))
    if any(a==b or F[b]==a for a,b in E):raise ValueError('not a partial orientation')
    @lru_cache(None)
    def all_count(U):
        k=U.bit_count()
        return 1<<(k*(k-1)//2-sum(bool(U>>a&1 and U>>b&1) for a,b in E))
    @lru_cache(None)
    def strong(U):
        if U==0:return 0
        result=all_count(U);A=(U-1)&U
        while A:
            B=U^A
            if all(not (B>>a&1 and A>>b&1) for a,b in E):result-=strong(A)*all_count(B)
            A=(A-1)&U
        return result
    return {'vertices':n,'all':all_count((1<<n)-1),'strong':strong((1<<n)-1)}


F=[6,3,1,0,0,1,3,6]
aut=[p for p in permutations(range(8)) if all(p[F[j]]==F[p[j]] for j in range(8))]
if len(aut)!=2:raise ValueError('functional automorphism count')
print(json.dumps({'full':count(F),'functional_automorphisms':aut,
                  'line_quotient':count([3,3,0,2])},indent=2,sort_keys=True))
