"""Independent polynomial-chart audit of the six port/Singer completions."""
from itertools import product
import json

def need(b,m):
    if not b:raise ArithmeticError(m)
def multiply(x,y,f):
    # Polynomial convolution followed by descending long division.
    r=0
    for i in range(5):
        if y>>i&1:r^=x<<i
    for i in range(8,4,-1):
        if r>>i&1:r^=f<<(i-5)
    return r

def main():
    polys=[41,55,47,37,61,59];tables=[]
    for f in polys:
        powers=[1]
        for _ in range(30):powers.append(multiply(powers[-1],2,f))
        need(len(set(powers))==31 and multiply(powers[-1],2,f)==1,'primitive orbit')
        coords=powers+[0];labels={x:i for i,x in enumerate(coords)}
        add=[[labels[coords[a]^coords[b]] for b in range(32)] for a in range(32)]
        mult=lambda a,b:31 if 31 in [a,b] else (a+b)%31
        for a,b in product(range(32),repeat=2):
            need(coords[mult(a,b)]==multiply(coords[a],coords[b],f),'field product')
            for c in range(32):need(mult(a,add[b][c])==add[mult(a,b)][mult(a,c)],'all distributivity')
        # All identity phases, not only one alternate phase.
        for t in range(31):
            transport=lambda a:31 if a==31 else (a+t)%31
            changed=lambda a,b:31 if 31 in [a,b] else (a+b-t)%31
            for a,b in product(range(32),repeat=2):
                need(changed(transport(a),transport(b))==transport(mult(a,b)),'phase transport isomorphism')
        tables.append(add)
    need(len({tuple(map(tuple,T)) for T in tables})==6,'six different additions')
    swap=lambda a:2 if a==1 else 1 if a==2 else a
    rot=lambda a:31 if a==31 else (a+1)%31
    bad=lambda a:swap(rot(swap(a)))
    orbit=[0]
    for _ in range(30):orbit.append(bad(orbit[-1]))
    logs={x:i for i,x in enumerate(orbit)}
    badmul=lambda a,b:31 if 31 in [a,b] else orbit[(logs[a]+logs[b])%31]
    A=tables[0];left=badmul(1,A[0][1]);right=A[badmul(1,0)][badmul(1,1)]
    need((left,right)==(16,6),'nonlinear clock hostile')
    print(json.dumps({'primitive_polynomials':polys,'identity_phases_per_field':31,
        'nonlinear_clock_witness':[1,0,1],'unequal_outputs':[left,right],
        'all_six_field_completions':True},sort_keys=True))
if __name__=='__main__':main()
