#!/usr/bin/env python3
"""Exact rational certificate and compressed symbolic control; no full huge orbit."""
from fractions import Fraction as F
from math import gcd
import hashlib,json
checks=0

def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def logs(z,n=512):
    lo=2*sum((z**(2*i+1)/F(2*i+1) for i in range(n)),F(0))
    return lo,lo+2*z**(2*n+1)/(F(2*n+1)*(1-z*z))

l2,u2=logs(F(1,3));l3,u3=logs(F(1,2))
alo,ahi=l2/u3,u2/l3
transform=lambda a:(77-122*a)/(2*a-1)
tlo,thi=transform(ahi),transform(alo)
p,q,m=1003,9885,1116187
dlo,dhi=F(p,q)-thi,F(p,q)-tlo
need(F(1,2)<alo<ahi<F(77,122),'positive decreasing transform domain')
need(0<tlo<thi<1,'rotation interval')
need(gcd(p,q)==1 and m*p%q==1 and m>2*q,'residue geometry')
need(F(1,m*q)<dlo<dhi<F(1,(m-q)*q),'exact logarithmic bracket')
need(m*dhi<F(2,q),'all other residues cannot cross')
need(int(m*tlo)==int(m*thi)==113255,'last floor exact')
H=lambda n:(n*p)//q-int(n%q==0)
g=[H(n)-H(n-1) for n in range(1,q+1)]
need(g[0]==1 and sum(g)==p,'period boundary and count')
for n in range(1,q+1):
    need(g[n-1] in (0,1),'period is binary')
    need(H(n+q)==H(n)+p,'period height transport')
    need(H(n+q)-H(n+q-1)==g[n-1],'period increment transport')
for n in (1,q-1,q,q+1,m-q,m-1,m):
    lo,hi=int(n*tlo),int(n*thi)
    need(lo==hi,'selected exact floor')
    need(lo==H(n)-int(n==m),'selected floor formula')
need(H(m)-H(m-1)==1,'forced reset restores final model symbol')
need(int(m*tlo)-int((m-1)*tlo)==0,'true last conditional reset absent')
K='110110110110110101011011101010';V=K*4+'10'
rot=g[1:]+g[:1]
period=''.join(V+('10' if z else '') for z in rot)
P=122*q+2*p;h=77*q+p;a=122
need(len(period)==P and period.count('1')==h,'literal period coefficients')
need(3**h<2**P,'literal period strictly contracts')
L=122*m;N=L+2*(113255+1);r=(L+99)//100
need(r>a+P,'all long matches inherit continuation')
need(a+P+r<=N,'two full matching factors fit')
need(N<2**52,'concrete hostile excluded by primary short theorem')
cycles,rem=divmod(m-1,q)
tail_bits=cycles*P+sum(122+2*z for z in rot[:rem])
need(a+tail_bits==N,'compressed representation includes exact final reset')
small={'p':p,'q':q,'m':m,'N':N,'L':L,'r':r,'preperiod_bits':a,
 'period_bits':P,'period_odd_bits':h,'full_tail_periods':cycles,'remaining_iteration_images':rem,
 'period_sha256':hashlib.sha256(period.encode()).hexdigest(),'mismatch_margin':r-a-P}
p=8116910408926339217391896275524023351
q=79995672446647253540458867489349070728
m=15001586057656515185428419394326237781535
dlo,dhi=F(p,q)-thi,F(p,q)-tlo
need(0<p<q and gcd(p,q)==1 and m*p%q==1 and m>2*q,'large residue geometry')
need(F(1,m*q)<dlo<dhi<F(1,(m-q)*q),'large exact logarithmic bracket')
need(m*dhi<F(2,q),'large residue noncrossing margin')
last=int(m*tlo)
need(last==int(m*thi)==1522163965844627930660358762978517428452,'large exact final floor')
need(last==(m*p-1)//q-1,'large floor formula')
L=122*m;N=L+2*(last+1);r=(L+99)//100;P=122*q+2*p;h=77*q+p
need(N>2**52 and r>122+P,'large instance beyond short range with full mismatch margin')
need(122+P+r<=N,'large matching occurrences fit')
need(q*thi<p,'large period contracts by exact log sign')
cycles,rem=divmod(m-1,q)
# Sum g_2 through g_(rem+1) telescopes, no huge period materialization.
H=lambda n:(n*p)//q-int(n%q==0)
rem_ones=H(rem+1)-H(1)
need(122+cycles*P+122*rem+2*rem_ones==N,'large compressed word length')
large={'p':p,'q':q,'m':m,'N':N,'L':L,'r':r,'preperiod_bits':122,
 'period_bits':P,'period_odd_bits':h,'full_tail_periods':cycles,'remaining_iteration_images':rem,
 'mismatch_margin':r-122-P,'last_greedy_floor':last,'period_materialized':False}
# Bank-wide image types and lawful output controls.
def branch(w):
    aa,bb,mm=1,0,1
    for bit in w:
        if bit=='1':aa,bb=3*aa,3*bb+mm
        mm*=2
    return aa,bb,mm

def val3(z):
    if z==0:return 'infinity'
    z=abs(z);t=0
    while z%3==0:t+=1;z//=3
    return t

bankbits={'K':K,'V':V,'F':'10','H':'110101000'}
bank={z:branch(w) for z,w in bankbits.items()}
image_types={}
for z,t in [('K',3),('V',4),('H',1)]:
    aa,bb,mm=bank[z];hh=bankbits[z].count('1')
    need(val3(bb-mm)==t and hh>t,'exact bank image type from h greater than valuation')
    residue=-bb*pow(aa,-1,mm)%mm
    for j in range(4):
        x=residue+j*mm
        if x==0:x+=mm
        y=(aa*x+bb)//mm
        need(mm*y==aa*x+bb,'lawful integer output')
        need(val3(y-1)==t,'exact output type on lawful controls')
        need(y%(3**hh)==bb*pow(mm,-1,3**hh)%(3**hh),'full past-image guard')
    image_types[z]=t
for j in range(40):
    x=1+4*j;y=(3*x+1)//4
    old,new=val3(x-1),val3(y-1)
    need(new==('infinity' if old=='infinity' else old+1),'F finite increment and x=1 exception')
aa,bb,mm=bank['V']
rhoV=bb*pow(mm,-1,243)%243
rhoVF=(3*bb+mm)*pow(4*mm,-1,243)%243
need((rhoV,rhoVF)==(82,1),'V versus VF exact residue separation')
for rtest in range(1,65):
    residue=(rhoV-rhoVF)*pow(2**rtest,-1,243)%243
    need(min(residue,243-residue)==81,'CRT distance uniformly equals eighty-one')
# The universal convergent-family statement is proved by the matrix identity.
# These independent interval controls verify the first 39 upper instances,
# including the large displayed input, without constructing their words.
xl,xh=tlo,thi;cp0,cp1,cq0,cq1=0,1,1,0;previous_upper=None;cf_controls=[]
for index in range(79):
    digit=int(xl);need(digit==int(xh),'continued-fraction digit certified by interval')
    if previous_upper is not None:
        up,uq,prevq,ui=previous_upper
        um=(digit+1)*uq+prevq
        low,high=F(up,uq)-thi,F(up,uq)-tlo
        need(F(1,um*uq)<low<high<F(1,(um-uq)*uq),'upper-convergent family bracket')
        need(um>2*uq and um*up%uq==1,'upper family size and inverse residue')
        ur=(122*um+99)//100;uP=122*uq+2*up
        need(122+uP+ur<122*um,'matching future blocks fit for whole family control')
        cf_controls.append((ui,up,uq,um))
        previous_upper=None
    cp,cq=digit*cp1+cp0,digit*cq1+cq0
    if index%2==1:
        need(cp*cq1-cp1*cq==1,'upper determinant identity')
        need(F(cp,cq)>thi,'upper convergent orientation')
        previous_upper=(cp,cq,cq1,index)
    cp0,cp1,cq0,cq1=cp1,cp,cq1,cq
    xl,xh=1/(xh-digit),1/(xl-digit)
need(cf_controls[3][1:]==(small['p'],small['q'],small['m']),'small control is index-seven instance')
need(cf_controls[-1][1:]==(large['p'],large['q'],large['m']),'large control is index-seventy-seven instance')
image_types['F']='finite type plus one; infinity fixed'
print(json.dumps({'small_control':small,'large_control':large,'log_series_terms':512,'bank_image_types':image_types,
 'image_residues_mod_243':{'V':rhoV,'VF':rhoVF},'crt_distance_factor':81,
 'upper_convergent_controls':len(cf_controls),'status':'PASS image-guard descent controls and compressed structural hostiles',
 'checks':checks,'full_orbit_replay':False,'full_word_materialized':False},indent=2,sort_keys=True))
