#!/usr/bin/env python3
"""Monotone inverse of a marked palindrome energy; no pattern enumeration in inverse."""
import json,hashlib
from centered_fans import debt,patterns,fan,count

def pattern(d,index):
    if d<1 or not 0<=index<(1<<(d//2)):raise ValueError('pattern address')
    s={0,d};width=d//2
    for a in range(1,width+1):
        if (index>>(width-a))&1:s.update([a,d-a])
    return sorted(s)

def decode(d,target):
    if d<1:raise ValueError('positive marked diameter')
    lo=0;hi=(1<<(d//2))-1;steps=0
    while lo<=hi:
        mid=(lo+hi)//2;s=pattern(d,mid);v=debt(s);steps+=1
        if v==target:return mid,steps
        if v<target:lo=mid+1
        else:hi=mid-1
    return None,steps

def main():
    complete=roundtrips=0;maxsteps=0
    for d in range(1,21):
        oracle={debt(s):s for s in patterns(d)}
        for target,s in oracle.items():
            address,steps=decode(d,target)
            if address is None or pattern(d,address)!=s:raise RuntimeError('complete pattern inverse')
            maxsteps=max(maxsteps,steps);roundtrips+=1
        if d<=10:
            for target in range(1001):
                address,steps=decode(d,target)
                if (address is not None)!=(target in oracle):raise RuntimeError('complete negative query')
                complete+=1
    # Common binary source -> symmetric fan -> actual raw H -> original source.
    bridge=[]
    for bits in ['01010','1011010','001101100101']:
        d=2*len(bits)+1;K=d+4;address=int(bits,2);s=pattern(d,address)
        rows,_=fan(s,K);H=count(rows,[0]);D=(1<<(K-1))+1-H
        recovered,steps=decode(d,D)
        if recovered!=address or format(recovered,f'0{len(bits)}b')!=bits:
            raise RuntimeError('raw count bridge roundtrip')
        bridge.append({'bits':bits,'diameter':d,'vertices':len(rows),'raw_H':H,'decode_steps':steps})
    big=[]
    for width in [64,128,256]:
        bits=''.join(str(hashlib.sha256(f'{width}:{i}'.encode()).digest()[0]%2) for i in range(width))
        d=2*width+1;s=pattern(d,int(bits,2));value=debt(s);address,steps=decode(d,value)
        if address!=int(bits,2):raise RuntimeError('large exact energy inverse')
        big.append({'width':width,'diameter':d,'energy_bits':value.bit_length(),'decode_steps':steps})
    print(json.dumps({'status':'passed','complete_pattern_roundtrips':roundtrips,'complete_small_membership_queries':complete,
         'max_small_steps':maxsteps,'raw_count_bridge':bridge,'large_integer_roundtrips':big,
         'scope':'Inverse is conditional on a marked diameter; absent energies return None. Large tests certify arithmetic, not all tournament-spectrum coverage.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
