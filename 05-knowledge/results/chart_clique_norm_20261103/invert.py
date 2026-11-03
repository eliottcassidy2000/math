"""Sharp ternary inverse for q^3+1-t^3+3qt when q=2 modulo3."""
import argparse,json

def value(q,t):return q**3+1-t**3+3*q*t

def recover(q,y,m):
    if q%3!=2 or m<1:raise ValueError('requires q=2 mod3 and at least one source digit')
    candidates=[r for r in range(3) if (value(q,r)-y)%9==0]
    if len(candidates)!=1:raise ValueError('residue is outside the realizable image modulo9')
    t=candidates[0]
    for j in range(1,m):
        candidates=[t+a*3**j for a in range(3) if (value(q,t+a*3**j)-y)%3**(j+2)==0]
        if len(candidates)!=1:raise RuntimeError('unit lifting failed')
        t=candidates[0]
    return t

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--q',type=int,default=5);p.add_argument('--residue',type=int,required=True)
    p.add_argument('--digits',type=int,default=2);a=p.parse_args()
    t=recover(a.q,a.residue,a.digits)
    print(json.dumps({'trace_residue':t,'source_modulus':3**a.digits,
                      'required_output_modulus':3**(a.digits+1),
                      'scope':'arithmetic inverse; elliptic-curve realization is separate'},indent=2))
