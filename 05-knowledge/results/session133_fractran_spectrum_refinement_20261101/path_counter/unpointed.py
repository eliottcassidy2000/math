"""Regenerate the unpointed fibre correction; independently test every apex <=9 vertices."""
import itertools
import json
from math import isqrt
from fibres import G,bounded_fibres,require
from check import rows_from_word,path_dp


def transitive_after_deletion(rows,x):
    vs=[v for v in range(len(rows)) if v!=x]
    for a,b,c in itertools.combinations(vs,3):
        if ((rows[a]>>b&1) and (rows[b]>>c&1) and (rows[c]>>a&1)) or \
           ((rows[a]>>c&1) and (rows[c]>>b&1) and (rows[b]>>a&1)):
            return False
    return True


def main():
    B=1000000;R=bounded_fibres(B)
    exceptions={2**k+1 for k in range(2,B.bit_length()) if 2**k+1<=B}
    fixed={G(q) for q in range(1,isqrt(B)+1) if G(q)<=B}
    for h in range(B+1):
        U=R[h]-int(h in exceptions)
        require(U>=0 and U%2==(int(h in fixed)+int(h in exceptions))%2,'unpointed parity correction')
    words=0
    for m in range(2,9):
        for middle in itertools.product((0,1),repeat=m-2):
            w=(1,)+middle+(0,);rows=rows_from_word(w)
            apices=[x for x in range(m+1) if transitive_after_deletion(rows,x)]
            expected=3 if m==2 else 2 if w in ((1,)+(0,)*(m-1),(1,)*(m-1)+(0,)) else 1
            require(len(apices)==expected,'deletion-apex classification');words+=1
        if m>=3:
            a=rows_from_word((1,)+(0,)*(m-1));b=rows_from_word((1,)*(m-1)+(0,))
            iso=(m,0)+tuple(range(1,m));anti=(1,0)+tuple(range(m,1,-1))
            require(all((a[u]>>v&1)==(b[iso[u]]>>iso[v]&1) for u in range(m+1) for v in range(m+1) if u!=v),'double-apex isomorphism')
            require(all((a[u]>>v&1)==(a[anti[v]]>>anti[u]&1) for u in range(m+1) for v in range(m+1) if u!=v),'double-apex anti-isomorphism')
            require(path_dp(a)==2**(m-1)+1,'exceptional H formula')
    return {'bound_H':B,'unpointed_strong_classes_total':sum(R)-len(exceptions),
            'unpointed_self_converse_classes_total':len(fixed)+len(exceptions),
            'unpointed_odd_fibres':len(fixed^exceptions),
            'double_apex_H_values':sorted(exceptions),
            'H_with_two_unpointed_self_converse_classes':sorted(fixed&exceptions),
            'finite_direct_apex_checks':{'reduced_words':words,'maximum_tournament_order':9},
            'pointed_versus_unpointed_hostile':{'H':5,'row_masks':[2,12,9,1],
                'pointed_words':['100','110'],'anti_isomorphism':[1,0,3,2],
                'R':2,'U':1,'pointed_self_converse':0,'unpointed_self_converse':1}}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
