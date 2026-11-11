"""Finite controls on the proved all-word comparison; no coverage inference."""
from functools import lru_cache
from itertools import product
from pathlib import Path
import json

def move(s, letter):
    h,e,u,d=s
    return {'O':(h,0,2*u,d), 'Z':(h+u-e,d,2*u-d,d),
            'J':(h+u-e,h-e,h+u,h)}[letter]

def frame(s):
    h,e,u,d=s
    return h,u,u-d,h+u-e

def qmove(q, letter):
    h,u,v,t=q
    return {'O':(h,2*u,u+v,h+2*u), 'Z':(t,u+v,2*v,t+2*v),
            'J':(t,h+u,u,h+2*u)}[letter]

@lru_cache(None)
def g(n):
    if n<2:return 1
    a=n//2
    return g(a) if n&1 else g(a)+a+(a&1)

def apply(s, word):
    for letter in word:s=move(s,letter)
    return s

def word(n):
    return ''.join('Z' if n>>j&1 else 'O' for j in range(n.bit_length()))

def equal_ray(w):
    tail=w.lstrip('Z')
    return set(tail)<={'O'} or (tail.startswith('J') and set(tail[1:])<={'O'})

def main():
    linear=0
    for s in product((-2,0,3),repeat=4):
        for letter in 'OZJ':
            if frame(move(s,letter))!=qmove(frame(s),letter):raise RuntimeError(('frame',s,letter))
            linear+=1
    states=[(31,7,45,13),(211,13,347,13),(211,25,269,25),(0,8,8,8),(-3,9,-5,2)]
    tail_controls=0
    for s in states:
        h,e,u,d=s
        for n in range(4096):
            expected=h+d*(g(n)-1)+(u-d)*n+(d-e)*(n&1)
            if apply(s,word(n))[0]!=expected:raise RuntimeError(('tail',s,n))
            f=n+1-(n&1)-g(n)
            if min(g(n)-1,f)<0:raise RuntimeError(('positive weights',n))
            tail_controls+=1
    full_words=0
    for length in range(9):
        for w in product('OZJ',repeat=length):
            w=''.join(w)
            gap=apply((0,-12,78,-12),w)[0]
            if gap<0 or ((gap==0)!=set(w).issubset({'O'})):raise RuntimeError(('211 reset',w,gap))
            ray=apply((0,8,8,8),w)[0]
            if ray<0 or ((ray==0)!=equal_ray(w)):raise RuntimeError(('ray language',w,ray))
            full_words+=1
    # Each signed q difference violates only one coordinate; its witness
    # is a control of the linear comparison theorem, not a physical source.
    hostile=[]
    for q,w in [((-1,1,1,1),''),((1,-1,1,1),'OZ'),
                ((1,1,-1,1),'ZZ'),((1,1,1,-1),'Z')]:
        h,u,v,t=q;s=(h,h+u-t,u,u-v);out=apply(s,w)[0]
        if out>=0:raise RuntimeError(('necessity witness',q,w,out))
        hostile.append({'q_difference':q,'word':w,'output_difference':out})
    report={'status':'FINITE-EXACT controls; all-word theorem proved in PROOF.md; Chan-Pak OPEN',
            'linear_controls':linear,'tail_controls':tail_controls,'all_three_letter_words_through_length_8':full_words,
            'necessity_hostiles':hostile,
            'universes':'81 signed linear states x3 letters; 5 specified states x4096 addresses; all9841 O/Z/J words length<=8',
            'positive_controls':'actual211 reset state pair and H13-to-H5 ray difference',
            'negative_control_scope':'signed response differences only, not claimed physical tournaments'}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
