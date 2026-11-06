"""Exact dyadic-class affine spanning controls; no external packages."""
from itertools import product
from functools import lru_cache
from pathlib import Path
import hashlib,json
from verify import need
@lru_cache(None)
def g(u):
    if u==1:return 1
    if u%2==0:return g(u//2)
    h=(u+1)//2;return g(h)+h-h%2

def block(word,h):
    u=h
    for letter in word:u=2*u if letter=='O' else 2*u-1
    return u

def affine(word):
    A,R,a,b,c=1,0,0,0,0;parity=None
    for letter in word:
        if letter=='Z':
            a+=A;c-=R
            if parity is None:b-=1
            else:c-=parity
            A,R=2*A,2*R+1;parity=1
        else:A,R=2*A,2*R;parity=0
    return A,R,a,b,c

def det(M):
    a,b,c=M[0];d,e,f=M[1];g,h,i=M[2]
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
words=span_controls=profile_controls=0
for k in range(9):
    residues=set()
    for letters in product('OZ',repeat=k):
        word=''.join(letters);A,R,a,b,c=affine(word)
        need(A==2**k and 0<=R<A and b in [0,-1],'affine block coefficients')
        residues.add((-R)%A);words+=1
        for h in range(2,60):
            u=block(word,h)
            need(u==A*h-R and g(u)==g(h)+a*h+b*(h%2)+c,'actual affine g formula')
            profile_controls+=1
        for m in range(1,8):
            hs=[2**m,2**(m+1),2**(m+2)-1]
            points=[(g(block(word,h)),block(word,h),1)for h in hs]
            base=[(g(h)+b*(h%2),h,1)for h in hs]
            need(det(points)==A*det(base) and det(points)!=0,'actual dyadic affine span')
            span_controls+=1
    need(residues==set(range(2**k)),'every dyadic residue class')
print(json.dumps(dict(status='PASS',max_block_length=8,binary_blocks=words,
    affine_profile_controls=profile_controls,exact_span_determinants=span_controls,
    power_parameter_range=[1,7],verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    scope='Every dyadic class has actual homogenized linear span3; excludes only independent-class even-scaling identities.'),indent=2,sort_keys=True))
