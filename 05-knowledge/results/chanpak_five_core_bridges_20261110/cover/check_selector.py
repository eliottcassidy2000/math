"""Independent direct-bit checks of mixed drift and the two-input identity."""
import itertools
import json
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent


def require(p, msg):
    if not p: raise RuntimeError(msg)


def g_direct(n):
    # Full binary inversion/adjacency sum, not the recursive compiler.
    bits=[int(x) for x in bin(n)[2:]]
    ans=1
    for i,x in enumerate(bits):
        if x:
            for j in range(i+1,len(bits)):
                if not bits[j]:
                    ans += 1 << (j-i-1)
                    if j==i+1: ans += 1
    return ans


def word_value(prefix,w,m,suffix):
    s=prefix+w*m+suffix
    return int(s or '0',2)


def drift(w):
    if not w: return Fraction(0)
    B=1<<len(w)
    r=int(w,2)
    x=-Fraction(r,B-1)
    p=int(w[-1]);D=Fraction(0)
    for b in w:
        if b=='0':D+=x+p
        x=2*x+int(b);p=int(b)
    return D


def apply_operator(values, roots):
    v=list(map(Fraction,values))
    for r in roots:
        v=[v[i+1]-r*v[i] for i in range(len(v)-1)]
    return v


def count(sig,a,b):
    d,kp,km,ip,im,op,om,jp,jm=sig
    return d*(g_direct(a)+g_direct(b)-1)+(d+kp)*a*b+(d+km)*(a&1)*(b&1)+ip*a-im*(a&1)+op*b-om*(b&1)+jp+jm


def main():
    words=[''.join(w) for k in range(1,7) for w in itertools.product('01',repeat=k)]
    one_checks=0
    for w in words:
        D=drift(w); B=1<<len(w)
        require((D>0)==('0' in w and '1' in w),'drift sign')
        for pre,suf in itertools.product(('','1','10'),('','0','11')):
            v=[g_direct(word_value(pre,w,m,suf)) for m in range(1,7)]
            got=apply_operator(v,[1,B])
            require(all(x==(1-B)*D for x in got),'one-input annihilator')
            one_checks+=1
    signatures=json.loads((HERE/'signatures.json').read_text())
    loops=['','0','1','01','10','001','101','0011']
    two_checks=0
    for w,z in itertools.product(loops,repeat=2):
        D=drift(w)+drift(z)
        roots=sorted({1,1<<len(w),1<<len(z),1<<(len(w)+len(z))})
        factor=1
        for r in roots:
            if r!=1:factor*=1-r
        for rec in signatures:
            vals=[count(rec['signature'],word_value('1',w,m,'10'),word_value('10',z,m,'1')) for m in range(1,8)]
            got=apply_operator(vals,roots)
            require(all(v==rec['signature'][0]*D*factor for v in got),'two-input annihilator')
            if D:require(any(got),'mixed drift lost')
            two_checks+=1
    return {'status':'PASS','single_loop_words':len(words),
            'single_loop_prefix_suffix_checks':one_checks,
            'two_input_signature_checks':two_checks,
            'scope':'Exact finite controls; general theorem is the written pumping proof.'}


if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
