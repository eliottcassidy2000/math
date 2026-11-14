# List collision pairs of F(u) = 908u + 9 g(u): u < u', F(u) = F(u'); show binary forms and d = u' - u.
def g(U):
    if U==0: return 1
    k=U.bit_length(); s=1
    for t in range(k):
        b=(U>>t)&1
        if not b: s+=U>>(t+1)
        if t+1<k and ((U>>(t+1))&1) and not b: s+=1
    return s
U=400000
vals={}
for u in range(1,U+1):
    F=908*u+9*g(u); vals.setdefault(F,[]).append(u)
pairs=[(v,us) for v,us in vals.items() if len(us)>1]
print("collision values up to u=%d: %d"%(U,len(pairs)))
from collections import Counter
dc=Counter()
for v,us in sorted(pairs)[:25]:
    a,b=us[0],us[1]
    print("F=%d: u=%d (%s, g=%d)  u'=%d (%s, g=%d)  d=%d=9*%d  g-drop %d"%(v,a,bin(a)[2:],g(a),b,bin(b)[2:],g(b),b-a,(b-a)//9,g(a)-g(b)))
for v,us in pairs: dc[us[1]-us[0]]+=1
print("most common d:",dc.most_common(12))
