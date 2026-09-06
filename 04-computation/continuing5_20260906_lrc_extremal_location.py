"""A located-pair appendix to the correlated-overlap certificate.
Only the displayed16704 edge/word boundary is tested; no whole-clock sieve.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd,ceil,floor
from functools import reduce
from pathlib import Path
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
G=0

def need(ok,label):
    global G
    G+=1
    if not ok:raise ArithmeticError(label)

def atlas(s):
    if not 2<=s<=356:return False
    p=2
    while p*p<=s:
        k=0
        while s%p==0:s//=p;k+=1
        if k and (p%3!=2 or k>2):return False
        p+=1
    return s==1 or s%3==2

def intervals(p,q):
    # Literal rational arc intersections; unlike a length formula, retains position.
    def arcs(v):return [(max(F(0),F(14*j-1,14*v)),min(F(1),F(14*j+1,14*v))) for j in range(v+1)]
    a,b=arcs(p),arcs(q);i=j=0;out=[]
    while i<len(a) and j<len(b):
        lo=max(a[i][0],b[j][0]);hi=min(a[i][1],b[j][1])
        if lo<hi:out.append((lo,hi))
        if a[i][1]<b[j][1]:i+=1
        elif b[j][1]<a[i][1]:j+=1
        else:i+=1;j+=1
    need(out[0][0]==0 and out[-1][1]==1,'origin wrap split and rejoin')
    return [(out[-1][0]-1,out[0][1])]+out[1:-1]

def main():
    t,e,p,q=16704,4,3,308;n=t//e
    d=(12,16,72,58,64,9,9)
    need(t%e==0 and n==4176 and gcd(p,q)==1 and atlas(p+q),'literal clock divisor and strict actual ratio')
    need((e*gcd(n,p),e*gcd(n,q))==d[:2],'exact forced sheet margins')
    need(all(t%x==0 for x in d) and reduce(gcd,d)==1,'all word states divide clock and primitive word')
    name='overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    adjacent=Path(__file__).resolve().parent/name
    path=adjacent if adjacent.is_file() else Path('C:/w/s0905/04-computation')/name
    raw=path.read_bytes()
    need(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','inherited complete profile pin')
    profiles={int(k):{(c,tuple(w)) for c,w in v['profiles']} for k,v in json.loads(raw)['levels'].items()}
    checked=0
    for k in range(1,7):
        for I in combinations(range(7),k):
            c=reduce(gcd,(d[i] for i in I));word=tuple(sorted(gcd(c,d[j]) for j in range(7) if j not in I))
            need((c,word) in profiles[7-k],'every projected full complement word');checked+=1
    need(checked==126,'complete projected word universe, not complete actual profiles')
    excess=sum(x*ceil(F(t,7*x)) for x in d)-t
    need(excess==188,'literal attained full-word cost')
    I=intervals(p,q)
    lengths=[b-a for a,b in I]
    need(len(I)==45 and sum(lengths)==F(1,49),'complete literal clipped pair geometry')
    need(lengths.count(F(1,2156))==43 and lengths.count(F(1,4312))==2,'43 full and2 half components')
    need(e*sum(ceil(n*x)-1 for x in lengths)==172,'incoming separate component credit recovered')
    walls=sorted({(n*x)%1 for arc in I for x in arc})
    cells=[((a+walls[(i+1)%len(walls)]+(i==len(walls)-1))/2)%1 for i,a in enumerate(walls)]
    values=[]
    for alpha in walls+cells:
        count=sum(ceil(n*b-alpha)-floor(n*a-alpha)-1 for a,b in I)
        # Independent literal grid predicate; positions do not enter this path.
        den=n*alpha.denominator
        direct=0
        for j in range(n):
            num=alpha.numerator+j*alpha.denominator
            if all(14*min((v*num)%den,(-v*num)%den)<den for v in (p,q)):direct+=1
        need(count==direct,'every wall and open cell by literal strict danger predicates')
        values.append((count,alpha))
    need(len(walls)==90 and len(cells)==90,'complete discontinuity partition')
    need(min(values)==(84,F(53,539)) and max(values)==(87,F(74,77)),'exact located minimum and maximum with owners')
    need(e*84==336>excess and e*84-excess==148,'strict improved overlap credit and actual conditional consumer')
    # Open endpoints are necessary: contrast exact wall with an adjoining interior.
    need(any(sum(ceil(n*b-w)-floor(n*a-w)-1 for a,b in I)!=sum(ceil(n*b-c)-floor(n*a-c)-1 for a,b in I) for w,c in zip(walls,cells)),'boundary values genuinely affect the count')
    print('STATUS FINITE-EXACT selected edge boundary repair plus PROVED translated-grid transfer')
    print('INPUT t16704 e4 n4176 pair3,308 word',d,'forced12,16 profiles126 excess188')
    print('GEOMETRY intervals45 measure1/49 lengths43/2156+2/4312 walls90 cells90')
    print('LOCATED minimum84 at53/539 maximum87 at74/77; separate43; physical credits336 versus172')
    print('CONSUMER any actual selected edge with total excess<=188 gives at least148 weak-safe lifts')
    print('SCOPE no physical realization of the word and no elimination of whole clock16704; no other edge candidates tested')
    print('PASS',G,'always-active exact gates; raw LF')

if __name__=='__main__':main()
