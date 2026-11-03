#!/usr/bin/env python3
"""Exact arithmetic and F20 character checks; Baby Monster inputs remain CITED.
Run normally or with -O. --mutate bad-centralizer must fail.
"""
import argparse,json
parser=argparse.ArgumentParser();parser.add_argument('--mutate',choices=['bad-centralizer']);args=parser.parse_args()
checks=0
def require(value,message):
    global checks
    checks+=1
    if not value:raise SystemExit('FAIL: '+message)
M=1774080000;HS=M//40;ca=443520000;cb=6000000
if args.mutate:ca-=1
require(M%40==0 and HS==44352000,'ATLAS subgroup order and HS index-two-base order')
centralizing_subgroup=5*HS
require(centralizing_subgroup>cb,'normal C5 factor cannot belong to class 5B')
require(2*centralizing_subgroup==ca,'cited 5A centralizer order is twice the robust lower bound')
require(4*ca==M,'normalizer bound is attained by the listed subgroup')
H=[(a,b) for a in range(5) for b in range(4)]
def mul(x,y):return ((x[0]+pow(2,x[1],5)*y[0])%5,(x[1]+y[1])%4)
def inv(x):return ((-pow(pow(2,x[1],5),-1,5)*x[0])%5,(-x[1])%4)
def act(x,t):return (pow(2,x[1],5)*t+x[0])%5
require(len({tuple(act(x,t) for t in range(5)) for x in H})==20,'faithful affine F20 action')
require(all(mul(mul(x,y),z)==mul(x,mul(y,z)) for x in H for y in H for z in H),'affine multiplication associativity')
require(all(mul(x,inv(x))==(0,0) for x in H),'inverse formula')
require(all(act(mul(x,y),t)==act(x,act(y,t)) for x in H for y in H for t in range(5)),'affine action respects multiplication')
a=(1,0);b=(0,1)
require(mul(mul(b,a),inv(b))==(2,0),'b*a*b^-1=a^2')
# Gaussian integers are stored as integer pairs; no floating-point arithmetic.
def gaussmul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def conj(x):return (x[0],-x[1])
powers_i=[(1,0),(0,1),(-1,0),(0,-1)]
chars=[{x:powers_i[(t*x[1])%4] for x in H} for t in range(4)]
for ch in chars:
    require(all(ch[mul(x,y)]==gaussmul(ch[x],ch[y]) for x in H for y in H),'linear character multiplicativity')
W={x:(sum(act(x,k)==k for k in range(5))-1,0) for x in H};chars.append(W)
require(W[(0,0)]==(4,0),'deleted permutation representation degree4')
for x in H:
    require(W[x]==((4 if x==(0,0) else -1 if x[1]==0 else 0),0),'deleted permutation character formula')
for ch in chars:
    require(all(ch[mul(mul(y,x),inv(y))]==ch[x] for x in H for y in H),'character is a class function')
for r,ch in enumerate(chars):
    for s,psi in enumerate(chars):
        values=[gaussmul(ch[x],conj(psi[x])) for x in H]
        inner=(sum(z[0] for z in values),sum(z[1] for z in values))
        require(inner==((20 if r==s else 0),0),'five characters are mutually orthonormal')
degrees=[ch[(0,0)][0] for ch in chars]
require(degrees==[1,1,1,1,4] and sum(d*d for d in degrees)==20,'irreducible degrees exhaust F20')
d=4371;c=21
require((d+4*c)%5==0 and (d-c)%5==0,'integral C5 multiplicities')
m0=(d+4*c)//5;m=(d-c)//5
require((m0,m)==(891,870),'cited 5A restriction multiplicities')
for x in H:
    if x[1]==0:
        require(m0+m*W[x][0]==(d if x==(0,0) else c),'restriction reconstruction on normal C5')
require(4*m+m0==d,'F20 decomposition total dimension')
# Hostile control: same C5 data can have no trivial H constituent.
alternative={x:(m*W[x][0]+m0*chars[2][x][0],m0*chars[2][x][1]) for x in H}
require(alternative[(0,0)]==(4371,0),'alternative rational H-module has the same degree')
require(all(alternative[x]==((d if x==(0,0) else c),0) for x in H if x[1]==0),'alternative H-module has the same C5 character')
require((sum(alternative[x][0] for x in H),sum(alternative[x][1] for x in H))==(0,0),'alternative H-module has no trivial H constituent')
require(10 not in degrees,'no degree10 complex irreducible of F20')
print(json.dumps({'status':'PASS','checks_executed':checks,'normalizer_order':M,'centralizer_5A':ca,'centralizer_5B':cb,'HS_order':HS,'robust_centralizer_lower_bound':centralizing_subgroup,'F20_irreducible_degrees':degrees,'restriction':{'four_dimensional_multiplicity':m,'total_linear_multiplicity':m0,'rationality_constraint':'a1=a3; a0+a1+a2+a3=891','whole_permutation_modules_not_forced':True,'hostile_module_trivial_H_multiplicity':0},'boundary':'Abstract embedding conditional on cited ATLAS subgroup data; no explicit Baby Monster generators or intertwiner.'},indent=2))
