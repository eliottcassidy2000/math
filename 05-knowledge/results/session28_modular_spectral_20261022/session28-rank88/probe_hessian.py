import runpy,contextlib,io,mpmath as m,json
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):d=runpy.run_path(str(Path(__file__).with_name('check_hecke.py')))
m.mp.dps=60;N=48; f=d['f'];n=88
rows=[]
for tstr in ['0.5','0.75','0.99','0.999','1','1.001','1.01','1.25','1.5','2','3','4','5']:
 t=m.mpf(tstr);q=m.exp(-2*m.pi*t)
 Fp=sum(-2*m.pi*k*int(f[k])*q**k for k in range(1,N+1))
 Fpp=sum((2*m.pi*k)**2*int(f[k])*q**k for k in range(1,N+1))
 delta=q*m.fprod((1-q**k)**24 for k in range(1,200))
 a=2*t*t*Fpp/(n*(n+2))+t*Fp/n;b=(m.pi*t)**2*delta**4
 rows.append({'t':tstr,'a':m.nstr(a,30),'b':m.nstr(b,30),'a_div_b':m.nstr(a/b,30),'perfection_floor_bound':m.nstr(a-2723200*b,30)})
out={'status':'NUMERICAL only: truncated theta throughq48 and Delta product200; nearby values expose the self-dual local maximum, separately certified by scalar-followon. No uniform bound.','rows':rows}
print(json.dumps(out,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
