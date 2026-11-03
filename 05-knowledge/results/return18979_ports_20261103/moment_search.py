"""Compile the inherited one-entry funnel grammar into requested moments only."""
from functools import lru_cache
import sys,random
from support import ROOT,HERE,need,save,junction
sys.path.insert(0,str(HERE.parent/'polarized_path_calculus_20261103'))
from recognize import class_counts,unrank
from explore import tail_coefficients

COUNTS=class_counts(20)
@lru_cache(None)
def decode(n,index):
    for m in [1]+list(range(3,n)):
        t=n-m-1;c=1 if not t else 1<<(t-1);u=1 if m==1 else COUNTS[m]
        if index>=c*u:index-=c*u;continue
        wordrank,rank=divmod(index,u)
        word=() if not t else tuple((wordrank>>i)&1 for i in range(t-1))+(1,)
        return m,rank,word
    raise RuntimeError('rank decode')
tail=lru_cache(None)(tail_coefficients)

@lru_cache(None)
def moment(n,index,m,kind):
    if n==1:return m if kind==0 else 1
    if kind==0 and m==0:return 0
    nn,j,word=decode(n,index);t=len(word);L,D,B=tail(word,m)
    if kind==0:return m*(tail(word,m+1)[0]*moment(nn,j,m+1,2)+B*moment(nn,j,m,0)-D*moment(nn,j,m,1))
    if kind==1:return L*moment(nn,j,m,0)-(m-1)**t*moment(nn,j,m,1)
    return m**t*moment(nn,j,m,2)+D*moment(nn,j,m-1,0)

def main():
    # Literal endpoint DP is independent of the funnel transfer calculation.
    checks=0
    for n in range(3,8):
        for j in range(COUNTS[n]):
            T=unrank(n,j);h,s,q,F=junction.profile(T,0)
            need((h,s,F)==(moment(n,j,1,0),moment(n,j,1,2),moment(n,j,2,1)),'independent parent moments');checks+=1
            junction.endpoints.cache_clear()
    lo=int(sys.argv[1]) if len(sys.argv)>1 else 14
    hi=int(sys.argv[2]) if len(sys.argv)>2 else 17
    searched=[];hit=None
    for n in range(lo,hi+1):
        ranks=list(range(COUNTS[n]));random.Random(18979+n).shuffle(ranks);examined=0
        for j in ranks:
            h=moment(n,j,1,0);F=moment(n,j,2,1);examined+=1
            if h+F-1==18979:
                T=unrank(n,j);need(junction.profile(T,0)[0::3]==(h,F),'actual positive parent profile')
                U=junction.ear(T,1);r=len(T);z=junction.profile(U,r)
                need(z[0]==18979 and z[3]-z[0]==1 and junction.strong(U),'actual ready18979')
                hit={'rows':U,'root':r,'profile':z,'parent_order':n,'budget':1,'rank':j};break
        searched.append({'parent_order':n,'class_count':COUNTS[n],'examined':examined})
        print('moment bank',n,examined,'hit',hit is not None,flush=True)
        if hit:break
    save('moment_search.json',{'status':'FINITE_EXACT_POSITIVE_WITNESS' if hit else 'FINITE_EXACT_ABSENCE_IN_SPECIFIED_ONE_ENTRY_BANKS','checks':checks,'banks':searched,'hit':hit,
      'source_formula':'polarized_path_calculus_20261103/explore.py:funnel_poly; selective integer evaluation, no polynomial interpolation',
      'scope':'positive finite search or listed one-entry banks only; no two-entry or all-order absence claim'})
    if hit:save('ready18979.json',{'status':'FINITE_EXACT_POSITIVE_WITNESS','hit':hit,'banks':searched,'scope':'first hit in declared shuffled one-entry banks; no minimality claim'})
    print('PASS moment search',None if hit is None else hit['profile'],flush=True)
if __name__=='__main__':main()
