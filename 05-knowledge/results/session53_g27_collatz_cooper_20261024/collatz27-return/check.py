"""Exact fixed-prefix Collatz/character checks; no empirical asymptotic claims."""
import json
from collections import Counter

def require(ok, label):
    if not ok: raise RuntimeError(label)

def T(n): return (3*n+1)//2 if n%2 else n//2

def prefix(n,k):
    out=[n]; bits=[]; a=1;b=0
    for j in range(k):
        e=n%2;bits.append(e)
        if e:a*=3;b=3*b+2**j
        n=T(n);out.append(n)
    return out,bits,a,b

def orbit(n,ordinary=False):
    out=[n]
    while n!=1:
        if len(out)>10000:raise RuntimeError('bound')
        n=3*n+1 if ordinary and n%2 else T(n)
        out.append(n)
    return out

def oddorbit(n):
    out=[n]; vals=[]
    while n!=1:
        m=3*n+1;a=0
        while m%2==0:m//=2;a+=1
        vals.append(a);out.append(m);n=m
    return out,vals

def chi(x):
    r=x%7
    return 0 if r==0 else (1 if r in (1,2,4) else -1)

def run():
    c=orbit(27,True); t=orbit(27);u,vals=oddorbit(27)
    result={}
    for name,o in [('ordinary',c),('shortcut',t),('odd',u)]:
        result[name]={'total':len(o)-1,'first_descent':next(i for i in range(1,len(o)) if o[i]<27),'maximum':max(o),'first_peak_index':o.index(max(o))}
    require([len(c)-1,len(t)-1,len(u)-1]==[111,70,41],'clocks')
    require([result[x]['first_descent'] for x in ('ordinary','shortcut','odd')]==[96,59,37],'descents')
    require(sum(vals)==70,'division clock')
    result['odd_valuation_word']=vals
    result['odd_orbit']=u
    data={}
    for k in (59,70):
        o,w,a,b=prefix(27,k)
        require(o==t[:k+1],'prefix baseline')
        require(a*27+b==2**k*o[-1],'affine')
        mod_a=a*pow(2**k,-1,7)%7;mod_b=b*pow(2**k,-1,7)%7
        require((mod_a,mod_b)==(6,1 if k==59 else 0),'mod7 return')
        for h in list(range(128))+[2**100+9]:
            n=27+2**k*h;oo,ww,aa,bb=prefix(n,k)
            require(ww==w and oo[-1]==o[-1]+a*h,'lift')
            if k==59:require(next(i for i in range(1,len(oo)) if oo[i]<n)==59,'glide lift')
        data[str(k)]={'odd_count':sum(w),'multiplier_numerator':a,'carry':b,'endpoint':o[-1],'mod7_map':[mod_a,mod_b]}
    for j in range(1,59):
        o,w,a,b=prefix(27,j);require(a>2**j,'pre-descent coefficient')
    require(3**37<2**59,'descent coefficient')
    # Full shortcut cylinders are exactly all 2^k residue classes; exhaust all words through k=10.
    counts=0;edges=0
    for k in range(1,11):
        signatures={}
        for r in range(2**k):
            o,w,a,b=prefix(r,k);w=tuple(w)
            require(w not in signatures,'parity bijection');signatures[w]=r
            residues=[]
            for h in range(7):
                n=r+2**k*h;oo,ww,_,_=prefix(n,k)
                require(tuple(ww)==w,'CRT same word')
                residues.append(n%7)
                require(oo[-1]%7==(a*n+b)*pow(2**k,-1,7)%7,'endpoint residue')
                counts+=1
            require(sorted(residues)==list(range(7)),'CRT all residues')
            A=a*pow(2**k,-1,7)%7;B=b*pow(2**k,-1,7)%7
            require(chi(A)==(-1)**sum(w),'slope sign')
            for x in range(7):
                for y in range(7):
                    if x==y:continue
                    require(chi((A*y+B)-(A*x+B))==(-1)**sum(w)*chi(y-x),'Paley transport');edges+=1
    # Exact full odd word is one extra bit finer than the 70 shortcut steps.
    for h in range(32):
        n=27+2**71*h;cur=n;got=[]
        for v in vals:
            m=3*cur+1;a=0
            while m%2==0:m//=2;a+=1
            got.append(a);cur=m
        require(got==vals and cur==1+2*3**41*h,'exact odd cylinder')
    coarse=27+2**70
    oo,ww,a,b=prefix(coarse,70)
    require(oo[-1]==1+3**41 and oo[-1]%2==0,'coarse not exact odd endpoint')
    require(T(27)%7!=T(34)%7 and 27%7==34%7,'no autonomous mod7 map')
    # Completing the same word is not completing convergence at the same time.
    require(prefix(27+2**70,70)[0][-1]!=1,'no transferred total stopping claim')
    require(T(T(41))==31 and 27<31<41,'odd local descent not parent payment')
    result['return_words']=data
    result['bounded_verification']={'crt_integer_prefixes':counts,'directed_paley_pair_checks':edges,'all_words_max_length':10,'positive_lifts_per_return_word':129,'exact_odd_lifts':32}
    result['hostiles']={'coarse_odd_lift':coarse,'coarse_endpoint':oo[-1],'same_residue_sources':[27,34],'different_successors':[T(27),T(34)],'source_reset':[27,41,31]}
    return result
if __name__=='__main__': print(json.dumps(run(),indent=2,sort_keys=True))
