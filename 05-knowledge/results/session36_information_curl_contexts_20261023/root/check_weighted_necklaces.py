#!/usr/bin/env python3
import itertools,json,math
from fractions import Fraction
from check_curl_decoder import matrix,key,encode,hp,check
from check_necklace_counts import permute,scc,phi

def logseries(a):
    out=[Fraction(0)]*len(a)
    for n in range(1,len(a)):
        out[n]=a[n]-sum(Fraction(j,n)*out[j]*a[n-j] for j in range(1,n))
    return out

def main():
    data=[([],[])];rows=[]
    for n in range(1,7):
        perms=list(itertools.permutations(range(n)));remain=set(range(1<<math.comb(n,2)));types=[];states={}
        while remain:
            b=min(remain);M=matrix(n,b)
            orbit=[permute(M,p) for p in perms]
            remain.difference_update(encode(X) for X in orbit)
            H=hp(M);types.append(H);K=min(key(X) for X in orbit)
            check(K not in states or states[K]==H,'H fibre clash');states[K]=H
        data.append((types,list(states.values())))
    errors=[]
    for s in (0,1,2):
        logs={d:logseries([1]+[sum(H**(d*s) for H in data[n][0]) for n in range(1,7)]) for d in range(1,7)}
        for n in range(1,7):
            expected=sum(Fraction(phi(d),d)*logs[d][n//d] for d in range(1,n+1) if n%d==0)
            actual=sum(H**s for H in data[n][1]);check(expected==actual,'weighted Adams formula')
            naive=sum(Fraction(phi(d),d)*logs[1][n//d] for d in range(1,n+1) if n%d==0)
            if naive!=actual:errors.append(dict(s=s,n=n,naive=str(naive),actual=actual))
            rows.append(dict(s=s,n=n,value=actual))
    check(any(x['s']==1 and x['n']==6 for x in errors),'missing repeated C3 hostile')
    print(json.dumps(dict(status='PASS',universe='all unlabelled tournament/curl types through n6, powers s0,1,2',rows=rows,omitted_Adams_hostiles=errors),sort_keys=True,indent=2))
if __name__=='__main__':main()
