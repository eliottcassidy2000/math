from pathlib import Path
import subprocess,json
from find_intersections import g,rev,append_coeff
HERE=Path(__file__).resolve().parent
ns={'__name__':'compiler'}
exec(subprocess.check_output(['git','show','HEAD:05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'],text=True),ns)
cores=[('K1',(0,)),('C3',(2,4,1)),('S4',(2,4,9,3)),('S5a',(2,4,9,19,7)),('S5b',(2,4,25,3,11))]
sigs=[]
for name,core in cores:
 for cut in range(1<<len(core)):sigs.append((name,core,cut,ns['signature'](core,cut)))
hits={-33:[],-31:[]}
for p in range(78):
 lp=p.bit_length();rp=rev(p,lp)
 for l in range(8):
  for r in range(1<<l):
   aa,bb,cc=append_coeff(r,l)
   for digit in (0,1):
    if p==0 and digit==0:continue
    parity=(r&1) if l else digit
    if digit==0:
     X=p*(1<<l);Y=r;G=p*(1+aa);J=g(p)-p+(p&1)+cc
     revX=rev(r,l)*(1<<lp);revY=rp
    else:
     X=(p+1)*(1<<l);Y=r-(1<<l);G=aa*(p+1);J=g(p)-aa+bb+cc
     revX=(rev(r,l)+1)*(1<<lp);revY=rp-(1<<lp)
    for name,core,cut,s in sigs:
     d,kp,km,i,im,o,om,jp,jm=s;L=d+kp
     for shift in range(7):
      numer=(77<<shift)-d*G-o*X;denom=L*X
      if numer<0 or numer%denom:continue
      A=numer//denom;eff=L*A+o;constA=d*(g(A)-1)+i*A-im*(A&1)+jp+jm
      C=d*J+eff*Y-om*parity+constA
      if C in hits:
       Ux=(d+eff)*(1<<(lp+l))-d*revX;Uy=-d*revY
       rec={'shift':shift,'core':name,'rows':core,'cut':cut,'A':A,'p':p,'loop':digit,'suffix':r,'suffix_length':l,'B_coeff':[X,Y],'g_coeff':[G,J],'U_coeff':[Ux,Uy],'d':d}
       hits[C].append(rec)
# unique physical address coefficients percore/cut/A
for C in hits:
 hits[C]=list({(h['shift'],h['core'],h['cut'],h['A'],*h['B_coeff']):h for h in hits[C]}.values())
pairs=[]
for low in hits[-33]:
 for high in hits[-31]:
  if (low['d'],low['U_coeff'])==(high['d'],high['U_coeff']):pairs.append({'lower':low,'upper':high})
report={'scope':'All fixed five-core cuts; p digit^m suffix with0<=p<=77 (empty prefix only for one loops),suffixlength<=7; A solved from exponential coefficient77*2^shift, shift0..6; originalrayindex m=looplength+shift. Exact matches to77*2^m-33 / -31.', 'counts':{str(C):len(h) for C,h in hits.items()},'paired_matches':pairs,'hits':{str(C):h[:25] for C,h in hits.items()}}
(HERE/'ray_rebases.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
