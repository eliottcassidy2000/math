"""Independent arithmetic controls for the structural addendum."""
from fractions import Fraction as F
import json

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

gamma=F(5878847,19531250)
D=8*5**10
Delta=F(1,2*5**11)
require(Delta==F(4,5*D),'clock correction')
require(1+2*gamma>F(9,8),'two-heavy exclusion')
expected=[0,0,1,2,2,3,4]
rows=[]
for j,maximum in zip(range(3,10),expected):
    require(gamma+F(maximum,5)<F(j,8),'retained count passes strict necessary test')
    require(gamma+F(maximum+1,5)>=F(j,8),'next integer count excluded')
    rows.append({'j':j,'maximum_allones_10_blocks':maximum})
# Algebraic conservation of the one terminal level-11 unit, for arbitrary
# synthetic sums. This does not assert actual native-word realization.
for A in range(1,101):
    for S10 in range(0,25):
        S11=S10+1
        l11=F(5,8)-F(1,D)
        l12=F(5,8)-F(1,5*D)
        require(l12*A-S11==(l11*A-S10)-1+Delta*A,'terminal conservation')
print(json.dumps({'status':'PASS','two_heavy_margin_over_9_over_8':str(1+2*gamma-F(9,8)),
 'native_cut_semantics':'proved in AUDIT.md under inherited assumptions',
 'bounds':rows,'additional_E12_families_closed':0},indent=2))
