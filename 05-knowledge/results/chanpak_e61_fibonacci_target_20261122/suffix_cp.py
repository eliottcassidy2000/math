"""Exact unrestricted binary models of the two remaining order-123 suffix charts."""
import argparse
from ortools.sat.python import cp_model

from check import p as reading

CHARTS = {'singleton': (81,25,43,68,115), 'four': (1021,125,915,940,112)}
TARGET = (1 << 60)*((1 << 61)-1)-1


def build(chart, target, leading_zeros=0, trailing_ones=0):
    c,h,k0,k1,n = chart
    model = cp_model.CpModel()
    bit = [model.new_bool_var(f'b{i}') for i in range(n)]
    model.add(bit[-1]==1)
    for i in range(leading_zeros):
        model.add(bit[-2-i]==0)
    for i in range(trailing_ones):
        model.add(bit[i]==1)

    diagonals = [0]
    for d in range(1,n):
        pairs = []
        for j in range(n-d):
            x = model.new_bool_var('')
            model.add_bool_or([~bit[j+d],bit[j],x])
            model.add_implication(x,bit[j+d])
            model.add_implication(x,~bit[j])
            pairs.append(x)
        count = model.new_int_var(0,n-d,f'd{d}')
        model.add(count==sum(pairs))
        diagonals.append(count)
    # Work in base 128 so every linear equation fits the solver's int64 API.
    carry=0
    for s in range(0,(target-k0).bit_length()+8,7):
        pieces = []
        for shift in range(7):
            i=s+shift
            if i<n:
                pieces.append((1<<shift)*c*bit[i])
            if i+1<n:
                pieces.append((1<<shift)*h*diagonals[i+1])
        if s==0:
            pieces.extend([h*diagonals[1],(k1-k0)*bit[0]])
        next_carry=model.new_int_var(0,2*h*n+2*c+abs(k1-k0),f'carry{s+7}')
        model.add(carry+sum(pieces)==(((target-k0)>>s)&127)+128*next_carry)
        carry=next_carry
    model.add(carry==0)
    m=2*h
    lhs=sum(pow(2,i,m)*c*x for i,x in enumerate(bit))+(k1-k0)*bit[0]+k0+(-target)%m
    model.add_modulo_equality(0,lhs,m)
    return model,bit


def couple_ends(model, bit, chart, target, block_size, low_ones):
    c,h,k0,k1,n=chart
    width=block_size
    scale=1<<width
    low=model.new_int_var(0,scale-1,'low_word')
    high=model.new_int_var(scale//2,scale-1,'high_word')
    model.add(low==sum((1<<i)*bit[i] for i in range(width)))
    model.add(high==sum((1<<i)*bit[n-width+i] for i in range(width)))
    table=[]
    for r in range(scale):
        if (r+1)%(1<<low_ones):
            continue
        reverse=int(format(r,f'0{width}b')[::-1],2)
        den=c*scale+h*(scale-1-reverse)
        rem=target-(k1 if r%2 else k0)-c*r-h*reading(r)
        lo=max((rem+den+h-1)//(den+h),1<<(n-width-1))
        hi=min(rem//den,(1<<(n-width))-1)
        if lo<=hi:
            table.extend((r,t) for t in range(lo>>(n-2*width),(hi>>(n-2*width))+1))
    model.add_allowed_assignments([low,high],table)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('chart',choices=CHARTS)
    ap.add_argument('--seconds',type=float,help='optional solver time budget')
    ap.add_argument('--export',help='optional OR-Tools model filename (.txt for text)')
    ap.add_argument('--control',action='store_true',help='solve an 18-bit attainable test target')
    args=ap.parse_args()
    chart=CHARTS[args.chart]
    c,h,k0,k1,n=chart
    target=TARGET
    leading,tail=(8,5) if args.chart=='four' else (0,0)
    if args.control:
        n=18;chart=(c,h,k0,k1,n)
        z=(1<<(n-1))+(1<<(n//2))+31
        target=k1+c*z+h*reading(z)
        leading=tail=0
    model,bit=build(chart,target,leading,tail)
    if not args.control:
        couple_ends(model,bit,chart,target,15 if tail else 12,tail)
    if args.export:
        if not model.export_to_file(args.export):
            raise RuntimeError('model export failed')
    solver=cp_model.CpSolver()
    solver.parameters.num_search_workers=2
    solver.parameters.random_seed=61
    if args.seconds is not None:
        if args.seconds<=0:
            ap.error('--seconds must be positive')
        solver.parameters.max_time_in_seconds=args.seconds
    result=solver.solve(model)
    print(solver.status_name(result), 'conflicts',solver.num_conflicts,'branches',solver.num_branches)
    if result in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        y=sum(solver.value(x)<<i for i,x in enumerate(bit))
        assert (k1 if y%2 else k0)+c*y+h*reading(y)==target
        print('verified address',y,'lower',target)


if __name__=='__main__':
    main()
