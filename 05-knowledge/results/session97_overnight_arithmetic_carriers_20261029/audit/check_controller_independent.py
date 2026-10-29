"""Independent integer-only audit; no imports from producer lanes."""
import json

def require(x, message):
    if not x:
        raise ArithmeticError(message)

def edge(n):
    z=3*n+1
    a=(z & -z).bit_length()-1
    return z>>a,a

def run(n):
    source=n; q=1; branches=steps=0
    while q:
        z,a=edge(n)
        branch=(a%3==0)
        if branch:
            require(n>=13 and a>=3,'branch source minimum')
            require(7*(z+1)<=3*(n+1),'branch ratio')
            branches+=1; q+=2
            # All branching outputs remain active. This certifies the bound
            # without logarithms and hence without rounding ambiguities.
            require(2*28**branches <= (source+1)*27**branches,
                    'universal branch bound')
        else:
            require(2*(z+1)<=3*(n+1),'leaf potential')
            require((a==1)==(2*(z+1)==3*(n+1)), 'leaf equality characterization')
            q-=1
        n=z;steps+=1
        require(steps <= 3*source+1, 'emergency finite test safety')
    require(steps==3*branches+1,'full ternary tree node count')
    return branches,steps,n

def main():
    max_b=(0,0);max_steps=(0,0)
    for n in range(1,100002,2):
        b,t,z=run(n)
        max_b=max(max_b,(b,n));max_steps=max(max_steps,(t,n))
    # sharp single-step branch contraction, and immediate growing leaf.
    require(edge(13)==(5,3),'sharp branch')
    require(run(3)==(0,1,5),'termination does not mean descent')
    require(run(1)==(0,1,1),'native fixed point is terminal leaf')
    print(json.dumps({'checked_positive_odd_seeds':50001,
       'max_branches_and_largest_seed':max_b,
       'max_steps_and_largest_seed':max_steps,
       'sharp_branch':[13,5,3], 'non_descent':[3,5]},sort_keys=True))
if __name__=='__main__':main()
