#!/usr/bin/env python3
"""Exact split/nonsplit clocks, with scalar fibres retained."""
import json
from pathlib import Path


def f32mul(a,b):
    c=0
    while b:
        if b&1:c^=a
        b>>=1;a<<=1
        if a&32:a^=0b100101 # z^5+z^2+1
    return c


def order(mul,x,one,limit):
    y=one
    for k in range(1,limit+1):
        y=mul(y,x)
        if y==one:return k
    raise RuntimeError('order bound')


def clock(q,add,mul,c,d,x):
    # u^2=c+d*u; product in the quadratic extension.
    def ext(a,b):
        ac=mul(a[0],b[0]);bd=mul(a[1],b[1])
        return add(ac,mul(c,bd)),add(add(mul(a[0],b[1]),mul(a[1],b[0])),mul(d,bd))
    one=(1,0);power=one;projective=None;scalar=None
    for k in range(1,q*q):
        power=ext(power,x)
        if power[1]==0 and projective is None:projective=k;scalar=power[0]
        if power==one:break
    else:raise RuntimeError('nonzero field element failed finite order')
    return {'q':q,'element':x,'linear_order':k,'projective_order':projective,
            'first_scalar':scalar,'scalar_fibre_order':k//projective},ext


if __name__=='__main__':
    add32=lambda a,b:a^b
    # Nonzero z has prime order 31, proving the degree-five polynomial irreducible.
    assert order(f32mul,2,1,31)==31
    assert all(f32mul(x,x)^x^1 for x in range(32))
    failed_generator,_=clock(32,add32,f32mul,1,1,(2,1))
    f32,_=clock(32,add32,f32mul,1,1,(6,1))
    add13=lambda a,b:(a+b)%13
    mul13=lambda a,b:a*b%13
    assert all((x*x-x-1)%13 for x in range(13))
    golden,_=clock(13,add13,mul13,1,1,(0,1))
    shifted,_=clock(13,add13,mul13,1,1,(2,1))
    assert golden['linear_order']==28 and golden['projective_order']==7
    assert shifted['linear_order']==56 and shifted['projective_order']==14
    assert order(mul13,2,1,12)==12
    assert f32['projective_order']==33
    assert failed_generator['projective_order']==11
    result={'split_32_order':31,'nonsplit_32':f32,'short_generator_32':failed_generator,'split_13_order':12,
            'golden_13':golden,'shifted_golden_13':shifted}
    Path(__file__).with_name('arithmetic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
