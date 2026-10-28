"""Independent least-prime-factor recurrence; imports no producer functions."""
from array import array
from collections import Counter
from pathlib import Path
from math import isqrt
import hashlib,json,sys

def need(b,s):
    if not b:raise RuntimeError(s)

def digest(a):
    need(a.itemsize==4,'uint32 encoding')
    a=array('I',a)
    if sys.byteorder!='little':a.byteswap()
    return hashlib.sha256(a.tobytes()).hexdigest()

def main():
    data=json.loads(Path(__file__).with_name('normal.json').read_text(encoding='utf-8-sig'))
    N=data['universe_max_source'];lp=array('I',[0])*(N+1)
    for p in range(2,N+1):
        if lp[p]:continue
        lp[p]=p
        for m in range(p*p,N+1,p):
            if not lp[m]:lp[m]=p
    sig=array('I',[0])*(N+1);tau=array('I',[0])*(N+1)
    power=array('I',[0])*(N+1);psum=array('I',[0])*(N+1);exp=array('B',[0])*(N+1)
    sig[1]=tau[1]=power[1]=psum[1]=1
    classes=Counter();unit=[];compressed=0;overlap=[];f343=[]
    for n in range(2,N+1):
        p=lp[n];m=n//p
        if m%p:
            power[n]=p;psum[n]=p+1;exp[n]=1;core=m
        else:
            power[n]=power[m]*p;psum[n]=psum[m]+power[n];exp[n]=exp[m]+1;core=n//power[n]
        sig[n]=sig[core]*psum[n];tau[n]=tau[core]*(exp[n]+1)
        v=sig[n]-n;d=(v-1)**2-4*n;c=d>=0 and isqrt(d)**2==d
        need(c==(tau[n]==4 or v==n+2),'compression iff')
        compressed+=c
        if v==n+2:unit.append(n)
        if c and tau[n]==6:overlap.append(n)
        if v==343:f343.append(n)
        classes['prime_power' if core==1 else ('factor_pair_compression' if c else 'remaining')]+=1
    need(digest(sig)==data['sigma_sha256_uint32le'],'every sigma value independently matched')
    need(digest(tau)==data['tau_sha256_uint32le'],'every divisor count independently matched')
    need(dict(classes)==data['class_counts'],'all classification counts')
    need(unit==data['unit_compression_sources'] and overlap==data['five_proper_divisor_overlap']==[20],'overlap lists')
    need(compressed==data['compressed_count'] and f343==data['complete_F343'],'compression count and inverse fibre')
    for t,values in data['finite_small_arity_overlaps'].items():
        need([n for n in unit if tau[n]-1==int(t)]==values,'independent small-arity overlaps')
    for v,values in data['complete_small_fibres'].items():
        v=int(v);need([n for n in range(1,(v-1)**2+1) if sig[n]-n==v]==values,'independent complete small fibre')
    print(json.dumps({'independent_range':N,'all_sigma_tau_values':'MATCH',
                      'class_counts':dict(classes),'unit_sources':unit,'five_term_overlap':overlap,
                      'F343_size':len(f343),'implementation':'least-prime-factor multiplicative recurrence'},sort_keys=True,indent=2))

if __name__=='__main__':main()
