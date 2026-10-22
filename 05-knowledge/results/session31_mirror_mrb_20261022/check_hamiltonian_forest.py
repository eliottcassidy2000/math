from pathlib import Path
import itertools, json, math
import sympy as s
x=s.symbols('x')
rows=[]
for n in range(2,6):
    edges=list(itertools.combinations(range(n),2)); sf=[0]*n; shf=[0]*n; sh=0
    for bits in range(1<<len(edges)):
        L=s.zeros(n); out=[0]*n
        for j,(a,b) in enumerate(edges):
            if bits>>j&1:a,b=b,a
            L[a,a]+=1; L[a,b]-=1;out[a]|=1<<b
        # coefficients det(I+x L) from characteristic polynomial, no interpolation
        cp=L.charpoly().all_coeffs()
        fs=[(-1)**j*cp[j] for j in range(n)]
        dp=[[0]*n for _ in range(1<<n)]
        for v in range(n):dp[1<<v][v]=1
        for mask in range(1,1<<n):
            for v in range(n):
                for w in range(n):
                    if not (mask>>w&1) and out[v]>>w&1:
                        dp[mask|1<<w][w]+=dp[mask][v]
        h=sum(dp[-1]);sh+=h
        for j in range(n):sf[j]+=fs[j];shf[j]+=h*fs[j]
    den=1<<len(edges);a=1+s.Rational(n,2)*x
    mean=sum(s.Rational(sf[j],den)*x**j for j in range(n))
    biased=sum(s.Rational(shf[j]*2**(n-1),den*math.factorial(n))*x**j for j in range(n))
    B=s.zeros(n)
    for j in range(n-1):
        B[j,j]+=s.Rational(1,2);B[j,j+1]-=s.Rational(1,2)
        B[j+1,j]+=s.Rational(1,2);B[j+1,j+1]-=s.Rational(1,2)
    assert s.expand(mean-a**(n-1))==0
    predicted=s.cancel((a*s.eye(n)+x*B).det()/a)
    assert s.expand(biased-predicted)==0
    assert s.expand(biased-mean).coeff(x,2)==s.Rational(n-2,4)
    assert sh==den*math.factorial(n)//2**(n-1)
    rows.append({'n':n,'tournaments':den,'mean':str(s.factor(mean)),'hamiltonian_biased':str(s.factor(biased)),'difference':str(s.expand(biased-mean))})
# Independent recurrence check including repeated zero roots in even orders.
z=s.symbols('z')
for n in range(2,13):
    Q=s.zeros(n)
    for j in range(n-1):Q[j,j]+=1;Q[j,j+1]-=1;Q[j+1,j]+=1;Q[j+1,j+1]-=1
    f0,f1=1,z
    for j in range(2,n):f0,f1=f1,s.expand(z*f1+f0)
    assert s.expand(Q.charpoly(z).as_expr()-z*f1)==0
result={'status':'exact rational arithmetic; all labelled tournaments n2..5; tridiagonal characteristic polynomials n2..12','rows':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
