#!/usr/bin/env python3
"""Exact arithmetic and all 650 boundary inequalities for the E11 lift."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

def need(b,s):
    if not b:raise RuntimeError(s)

root=Path(__file__).resolve().parent
c=json.loads((root/'epoch2_9_audit.json').read_text())
need(c['status']=='PASS' and c['terminal_fee_two'],'certificate status')
eta=F(-c['terminal_numerator'],c['denominator']);gamma=eta-2
lam=F(5,8)*(1-F(1,5**11));beta=F(-c['maximum'],c['denominator'])
need(beta==2*lam and 0<gamma<beta,'global non-one margin')
points=list(product(range(5),repeat=2));phi=lambda x:F((-x[0])%5,5)
checks=0
for x in points:
 y=((x[0]+4)%5,(x[1]+4)%5)
 need(int(y==(0,0))+phi(y)-phi(x)<=F(1,5),'all-ones edge');checks+=1
 for y in points:
  need(int(y==(0,0))+phi(y)-phi(x)<=1,'arbitrary edge');checks+=1
T=2*5**10;reward=sum(T//(2*5**j) for j in range(11))
need(reward==T*lam,'sharp all-ones word')
nine=[(a,b) for b in range(1,4) for a in range(5-b)]
need(len(nine)==9,'nine count')
for a,b in nine:need(b*(eta-1)>=1-F(a,5),'native-price family')
seven=[(a,b) for b in range(1,5) for a in range(5) if F(a,5)+b*gamma<1]
need(seven==[(0,1),(1,1),(2,1),(3,1),(0,2),(1,2),(0,3)],'next family list')
out=dict(status='PASS',certificate='epoch2_9_audit.json',eta=str(eta),
         nonone_10_excursion_price=str(eta-1),lambda11=str(lam),
         nonone_S11_margin=str(gamma),global_F2_margin=str(beta),
         boundary_inequalities=checks,closed_E11_families=nine,
         sharp_word_clock=T,sharp_word_content=reward,
         remaining_E12_families=seven,all_layer_content='OPEN',universal_collatz='OPEN')
(root/'lift_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
