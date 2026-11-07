"""Exact H27 leaf-population threshold; standard library, no imports from producer."""
import itertools,math,json,fractions
checks=0
def req(test,message):
 global checks
 checks+=1
 if not test:raise RuntimeError(message)
def N(n):
 return 2*(n-5)+3*max(n-6,0)+4*(math.comb(n-6,2) if n>=8 else 0)+(math.comb(n-6,3) if n>=9 else 0)
def S(n):return (n-5)//2 if n%2 else 0
def words(n):
 if n<9:return set()
 length=n-6
 return {tuple(3 if i in pos else 1 for i in range(length)) for pos in itertools.combinations(range(length),3)}
def A(n):return math.comb(n-6,3) if n>=9 else 0
def P(n):return (n-7)//2 if n%2 else 0
def L(n):return (A(n)-P(n))//2
def M(n):return (N(n)+S(n))//2
def E(n):return n**3-33*n*n+272*n-672
def O(n):return n**3-33*n*n+263*n-615
rows=[]
for n in range(6,81):
 req(N(n)==n*(n-5)*(n-4)//6,'literal SCC counts equal polynomial')
 req((N(n)+S(n))%2==0 and (A(n)-P(n))%2==0,'orbit integrality')
 if n<=40:
  ww=words(n);pal=sum(w==w[::-1] for w in ww)
  orbits={min(w,w[::-1]) for w in ww if w!=w[::-1]}
  req(len(ww)==A(n),'explicit word count')
  req(pal==P(n),'explicit palindrome count')
  req(len(orbits)==L(n),'explicit reversal orbit count')
 q=O(n) if n%2 else E(n)
 req(6*(2*A(n)-2*P(n)-N(n)-S(n))==q,'majority polynomial identity')
 req((2*L(n)>M(n))==(n>=23),'finite exact majority control')
 req(E(n+2)-E(n)==6*(n-10)**2-180,'even difference identity')
 req(O(n+2)-O(n)==6*(n-10)**2-198,'odd difference identity')
 rows.append({'n':n,'leaves':L(n),'merged_nodes':M(n),'fraction':str(fractions.Fraction(L(n),M(n))),'majority_sign':q})
req(all((O(n) if n%2 else E(n))<0 for n in range(6,23)),'entire initial range')
req(O(23)==144 and E(24)==672,'positive parity anchors')
req(6*(23-10)**2-198>0,'uniform tail difference lower bound')
print(json.dumps({'status':'PASS','active_checks':checks,'threshold':23,'scope':'Majority leaves in the fixed H27 fundamental paired-black stratum, not a global spectrum threshold. Counts inherit the complete low-H SCC inventory.','selected_rows':[r for r in rows if r['n'] in (14,21,22,23,24)],'all_initial_rows':[r for r in rows if r['n']<=24]},indent=2,sort_keys=True))
