import mpmath as m,json
from pathlib import Path
m.mp.dps=180
A=lambda x:2/m.pi*m.exp(-m.log(x)/x)*m.sin(m.pi/x)
C=m.nsum(lambda k:(-1)**k*m.expm1(m.log(k)/k),[1,m.inf])
r=C; rows=[]; sg=1; reconstruction=m.mpf(0)
for j in range(4):
 X=2/r;L=m.log(X);guess=X-L+(L-L**2/2-m.pi**2/6)/X
 # Initial moderate-scale root robust brackets; later residual normalized.
 x=m.findroot(lambda x:(A(x)-r)/r,(max(m.mpf(3),guess-1),guess+1))
 lo=int(m.floor(x));candidates=[lo,lo+1]
 k=min(candidates,key=lambda z:abs(r-A(z)))
 residual=r-A(k)
 reconstruction+=sg*A(k)
 sg*=1 if residual>=0 else -1
 nxt=abs(residual)
 assert nxt<=(A(lo)-A(lo+1))/2+m.mpf('1e-170')
 assert nxt<=1/m.mpf(lo)**2
 rows.append({'step':j+1,'index':str(k),'inverse':m.nstr(x,100),'residual':m.nstr(residual,85),'next_absolute_residual':m.nstr(nxt,85),'inverse_expansion_error':m.nstr(x-guess,20)})
 r=nxt
assert abs(C-reconstruction-sg*r)<m.mpf('1e-170')
out={'status':'180-digit numerical verification, not interval-certified integer decisions','rule':'nearest A-value with signed residual tracked','rows':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
