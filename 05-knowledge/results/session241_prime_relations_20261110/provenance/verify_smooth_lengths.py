"""Independent finite regression of the INHERITED smooth atom length law."""
import json
from pathlib import Path
B=12
atoms=((1,0),(0,2),(2,1),(0,3))
dp={(0,0):{0}}
for a in range(B+1):
    for b in range(B+1):
        if (a,b)==(0,0): continue
        dp[a,b]={k+1 for da,db in atoms if a>=da and b>=db for k in dp.get((a-da,b-db),set())}
        formula={a+v-w+t for v in range(b//2+1) for w in range(min(a//2,b)+1) for t in range(b//3+1) if 2*v+w+3*t==b}
        if dp[a,b]!=formula: raise RuntimeError((a,b,dp[a,b],formula))
        if bool(formula)!=((a,b) not in {(0,1),(1,1)}): raise RuntimeError('coverage')
        if formula and 2*max(formula)>5*min(formula): raise RuntimeError('elasticity')
if dp[4,2]!={2,5}: raise RuntimeError('sharp control')
print(json.dumps({'status':'PASS','universe':'0<=a,b<=12','nonunit_targets':168,'methods':['four-atom dynamic programming','three-counter exponent formula'],'sharp_control':{'value':3969,'lengths':sorted(dp[4,2])},'scope':'finite regression of inherited theorem; not a novelty claim'},indent=2))
