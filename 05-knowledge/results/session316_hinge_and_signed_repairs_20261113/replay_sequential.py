from itertools import product
from functools import cache
from replay_neutral_bank import R, rev, extend, profile, fuse

P=cache(profile)
def count(word, m):
    d=len(word)
    left=[P(R if i<m else extend(R,0),u) for i,u in enumerate(word)]
    right=[P(rev(R),31^u) for u in reversed(word)]
    return fuse(left+right)

for d in range(2,8):
    for m in range(d,-1,-1):
        fibers={}
        for word in product((3,5),repeat=d):
            key=(word.count(5),word[m:])
            h=count(word,m)
            if key in fibers and fibers[key]!=h:
                raise ValueError(('nonconstant',m,key))
            fibers[key]=h
        if m<d:
            for (weight, suffix),h in fibers.items():
                other=(weight,(5,)+suffix[1:])
                if suffix[0]==3 and other in fibers and not h<fibers[other]:
                    raise ValueError(('wrong sign',m,weight,suffix))
print('all neutral fibers and right-to-left signs through seven pairs passed')
