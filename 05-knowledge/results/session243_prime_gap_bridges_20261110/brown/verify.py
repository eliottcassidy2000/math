#!/usr/bin/env python3
"""Finite exact Brown inverse controls; not an all-n conclusion."""
import argparse, json, math

def require(test, message):
    if not test: raise RuntimeError(message)

def factor(n):
    ans={}; d=2
    while d*d<=n:
        while n%d==0:
            ans[d]=ans.get(d,0)+1; n//=d
        d+=1
    if n>1: ans[n]=ans.get(n,0)+1
    return ans

def factorial_factors(n):
    out={}
    for k in range(2,n+1):
        for p,e in factor(k).items():out[p]=out.get(p,0)+e
    return out

def inverse(n):
    fs=factorial_factors(n); C=math.factorial(n)//4
    even=2**(fs[2]-2); blocks=[p**e for p,e in fs.items() if p!=2]
    products=[even]
    for b in blocks:products += [a*b for a in products]
    hits=[]
    for a in products:
        b=C//a
        require(a*b==C and math.gcd(a,b)==1,'invalid unitary partition')
        if abs(a-b)==1:hits.append({'root':a+b,'even':a,'odd':b})
    return len(products),hits

def signs(n,m):
    f=factorial_factors(n); e=f.pop(2)
    return (m%(2**(e-1)) in (1,2**(e-1)-1) and
            all(m%(p**v) in (1,p**v-1) for p,v in f.items()))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fail',action='store_true');args=ap.parse_args()
    require(not args.fail,'deliberate negative control')
    rows=[];checked=0; roots=[]
    for n in range(4,61):
        selections,hits=inverse(n);checked+=selections
        q=math.factorial(n)+1;m=math.isqrt(q)
        direct=[m] if m*m==q else []
        require([h['root'] for h in hits]==direct,'inverse/isqrt mismatch')
        for h in hits:
            require(signs(n,h['root']),'true Brown root failed sign test')
            roots.append([n,h['root']])
        rows.append({'n':n,'selections':selections,'hits':hits})
    require(roots==[[4,5],[5,11],[7,71]],'unexpected finite result')
    require(math.gcd(7,math.factorial(5))==1 and not signs(5,7),'wheel hostile failed')
    require(signs(4,23) and 23*23-1==22*math.factorial(4),'quotient hostile failed')
    congruences=0
    for n in range(4,15):
        for m in range(1,1001,2):
            require(signs(n,m)==((m*m-1)%math.factorial(n)==0),'sign iff failure')
            congruences+=1
    print(json.dumps({'status':'PASS','universe':'all n=4,...,60 and every unitary prime-power partition at each n','selection_count':checked,'sign_congruence_controls':congruences,'brown_roots':roots,'rows':rows,'hostiles':{'wheel_survivor':7,'sign_survivor_wrong_quotient':23,'nonphysical_square_root':[49,7]}},indent=2,sort_keys=True))
if __name__=='__main__':main()
