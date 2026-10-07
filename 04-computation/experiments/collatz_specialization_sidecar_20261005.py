"""Fixed-degree spatial recovery from independently truthful localized moments.

No Collatz orbit, ROOT receipt, weight bank, or actual measurement oracle is
used. Synthetic laws test the exact source-indexed interval calculus only.
"""
from fractions import Fraction as F
from functools import lru_cache
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def natural(n):
    need(type(n) is int and n >= 0, "exact natural required")


def rational(x):
    need(type(x) in (int,F), "exact rational required")
    return F(x)


@lru_cache(None, typed=True)
def shell(d):
    natural(d)
    t = 1 << d
    return F(4*t,(1+t)**2)


def contraction(degree):
    natural(degree)
    need(degree >= 7, "degree at least seven required")
    return 2*(shell(1)**degree+shell(2)**degree+F(1,2**degree-1))


def tail_bound(degree, radius):
    contraction(degree); natural(radius)
    return 2*F(4**degree,2**(degree*(radius+1)))/(1-F(1,2**degree))


def parameters_to_error(degree, epsilon):
    """Finite measurement specification, conditional on truthful intervals.

    No target positivity or moment oracle is used in choosing these parameters.
    The returned address radius is a worst-case half-line stencil radius.
    """
    eta = contraction(degree)
    epsilon = rational(epsilon)
    need(0 < epsilon <= 1, "tolerance in (0,1] required")
    depth, power_bound = 1, eta**2
    while power_bound > epsilon/3:
        depth += 2
        power_bound *= eta**2
    radius = 0
    while 2*tail_bound(degree,radius)/(1-eta) > epsilon/3:
        radius += 1
    measurement_width = (1-eta)*epsilon/3
    bound = power_bound + (2*tail_bound(degree,radius)+measurement_width)/(1-eta)
    return {"degree":degree,"epsilon":epsilon,"depth":depth,"radius":radius,
            "address_radius":depth*radius,"measurement_width":measurement_width,
            "guaranteed_width":bound}


def row_step(row, degree, radius):
    out = {}
    for j, value in row.items():
        for d in range(1,radius+1):
            t = shell(d)**degree
            for k in (j-d,j+d):
                if k >= 0:
                    out[k] = out.get(k,F(0))+value*t
    return out


@lru_cache(None, typed=True)
def stencil(source, degree, radius, depth):
    natural(source); contraction(degree); natural(radius); natural(depth)
    row = {source:F(1)}
    coefficients = dict(row)
    for step in range(1,depth+1):
        row = row_step(row,degree,radius)
        for j,value in row.items():
            coefficients[j] = coefficients.get(j,F(0))+(-1)**step*value
    coefficients = tuple(sorted((j,c) for j,c in coefficients.items() if c))
    norm = sum((abs(c) for _,c in coefficients),F(0))
    remainder = sum(row_step(row,degree,radius).values(),F(0))
    return coefficients,norm,remainder


def certify(source, degree, radius, depth, intervals):
    """Conditional atom interval from a finite patch of H_degree(j) intervals.

    Passing checks does not authenticate their simultaneous provenance.
    """
    coeff,norm,remainder = stencil(source,degree,radius,depth)
    need(type(intervals) is dict, "indexed interval dictionary required")
    required = {j for j,_ in coeff}
    need(set(intervals) == required, "exactly the requested spatial addresses required")
    checked = {}
    for j,pair in intervals.items():
        natural(j)
        need(type(pair) is tuple and len(pair)==2,"interval pair required")
        lo,hi = map(rational,pair)
        need(0 <= lo <= hi <= 1,"ordered probability-moment interval required")
        checked[j] = (lo,hi)
    lo = sum((c*checked[j][0 if c>=0 else 1] for j,c in coeff),F(0))
    hi = sum((c*checked[j][1 if c>=0 else 0] for j,c in coeff),F(0))
    bill = norm*tail_bound(degree,radius)
    # S_depth(T_R)H = p-(-T_R)^(depth+1)p+S_depth(T_R)(T-T_R)p.
    lower = lo-bill-(remainder if depth%2==0 else 0)
    upper = hi+bill+(remainder if depth%2==1 else 0)
    lower,upper = max(F(0),lower),min(F(1),upper)
    need(lower <= upper,"packet contradicts a necessary law constraint")
    return {"source_index":source,"degree":degree,"radius":radius,"depth":depth,
            "lower":lower,"upper":upper,"raw_lower":lo,"raw_upper":hi,
            "coefficient_norm":norm,"tail_bill":bill,"power_remainder":remainder,
            "addresses":tuple(j for j,_ in coeff)}


def retained_bounds(results):
    need(type(results) is tuple and results,"nonempty result tuple required")
    source = results[0]["source_index"]
    lower,upper = F(0),F(1)
    history = []
    for result in results:
        need(result["source_index"]==source,"one fixed source required")
        lower,upper=max(lower,result["lower"]),min(upper,result["upper"])
        need(lower<=upper,"inconsistent retained atom intervals")
        history.append((lower,upper))
    return tuple(history)


def moments(law, degree, addresses):
    return {m:sum((mass*shell(abs(m-j))**degree for j,mass in law),F(0))
            for m in addresses}


def packet_for(law, source, degree, radius, depth, width=F(0)):
    coeff,_,_=stencil(source,degree,radius,depth)
    values=moments(law,degree,[j for j,_ in coeff])
    return {j:(max(F(0),v-width/2),min(F(1),v+width/2)) for j,v in values.items()}


def main():
    checks=0
    def check(ok,label):
        nonlocal checks
        need(ok,label);checks+=1
    for d in range(7,17):
        check(contraction(d)<1,"global contraction")
    check(2*(shell(1)**6+shell(2)**6)>1,"degree six fails uniform diagonal dominance")
    check(contraction(8)<=F(27,32),"practical degree eight bound")
    check(contraction(12)<F(1,2),"degree twelve faster inverse")
    plans = []
    for degree in (7,8,12):
        for epsilon in (F(1,2),F(1,16),F(1,256),F(1,65536)):
            plan=parameters_to_error(degree,epsilon)
            check(plan["guaranteed_width"]<=epsilon,"finite precision compiler")
            check(plan["depth"]%2==1,"lower-sided positive power remainder")
            check(plan["measurement_width"]>0,"positive requested measurement width")
            if epsilon==F(1,256): plans.append(plan)
    for m in range(7):
        for j in range(13):
            t=F(2)**(m-j)
            a,b=shell(abs(m-j)),shell(abs(m+1-j))
            check((8/b-4/a-2)/3==t,"two neighboring pointwise coordinates retain orientation")
    check(shell(1)==shell(abs(1-2)) and shell(2)!=1,"single-center reflection collision")
    laws=[((j,F(1)),) for j in range(10)]
    laws += [((a,F(1,4)),(b,F(3,4))) for a,b in ((0,2),(1,5),(2,7),(4,9))]
    cases=0
    for law in laws:
        for m in range(5):
            atom=sum((v for j,v in law if j==m),F(0))
            for degree in (7,8,12):
                for radius,depth in ((1,1),(2,2),(3,3),(3,5)):
                    result=certify(m,degree,radius,depth,packet_for(law,m,degree,radius,depth))
                    check(result["lower"]<=atom<=result["upper"],"synthetic atom bracket")
                    eta=contraction(degree)
                    check(result["coefficient_norm"]<=1/(1-eta),"uniform noise amplification")
                    check(result["power_remainder"]<=eta**(depth+1),"positive power remainder")
                    check(all(abs(j-m)<=radius*depth for j in result["addresses"]),"finite address radius")
                    cases+=1
    law=((0,F(1,16)),(2,F(15,16)))
    h=moments(law,8,(0,))[0]
    ghost_mass=(h-shell(3)**8)/(shell(1)**8-shell(3)**8)
    ghost=((1,ghost_mass),(3,1-ghost_mass))
    check(0<ghost_mass<1 and moments(ghost,8,(0,))[0]==h,
          "one target moment admits a zero-target completion")
    positive=certify(0,8,3,31,packet_for(law,0,8,3,31))
    check(positive["lower"]>0 and positive["lower"]<=F(1,16),
          "finite neighboring measurements give an independent synthetic floor")
    noisy=certify(0,8,3,31,packet_for(law,0,8,3,31,F(1,100000)))
    check(noisy["lower"]>0,"noise-robust positive synthetic floor")
    check(positive["raw_lower"]-noisy["raw_lower"]<=positive["coefficient_norm"]/100000,
          "noise bill")
    results=tuple(certify(0,12,r,k,packet_for(law,0,12,r,k))
                  for r,k in ((1,1),(2,3),(3,5),(3,9)))
    history=retained_bounds(results)
    check(all(history[i][0]<=history[i+1][0] and history[i][1]>=history[i+1][1]
              for i in range(len(history)-1)),"retained lower bounds survive refinement")
    delta=((0,F(1)),)
    check(all(v>0 for v in moments(delta,8,range(100)).values()),"smoothing is positive everywhere")
    absent=certify(3,12,3,15,packet_for(delta,3,12,3,15))
    check(absent["lower"]==0 and absent["upper"]<F(1,1000),"positive field does not imply positive atom")
    for bad in (lambda:contraction(6),lambda:contraction(True),
                lambda:stencil(True,8,3,3),lambda:stencil(0,8,-1,3),
                lambda:stencil(0,8,3,1.0),
                lambda:parameters_to_error(8,0),lambda:parameters_to_error(8,True),
                lambda:parameters_to_error(8,0.01),
                lambda:certify(0,8,3,3,{}),
                lambda:certify(0,8,0,0,{0:(0.0,1.0)}),
                lambda:retained_bounds((results[0],dict(results[0],source_index=1)))):
        try: bad()
        except (ValueError,TypeError):check(True,"exact API hostile")
        else:raise ValueError("hostile accepted")
    def compact(result):
        def display(x):
            if max(x.numerator.bit_length(),x.denominator.bit_length())<512:
                return str(x)
            scale=2**20
            floor=x.numerator*scale//x.denominator
            return {"dyadic_lower":str(F(floor,scale)),
                    "dyadic_upper":str(F(floor+1,scale)),
                    "numerator_bits":x.numerator.bit_length(),
                    "denominator_bits":x.denominator.bit_length()}
        return {k:display(result[k]) for k in ("lower","upper","coefficient_norm","tail_bill","power_remainder")}
    report={"status":"PROVED conditional spatial inversion; FINITE-EXACT synthetic laws only",
            "checks":checks,"synthetic_brackets":cases,
            "contraction_bounds":{str(d):str(contraction(d)) for d in (7,8,12)},
            "precision_plans_epsilon_1_over_256":[
                {key:(str(value) if type(value) is F else value)
                 for key,value in plan.items() if key!="guaranteed_width"} for plan in plans],
            "positive_control":compact(positive),"positive_control_addresses":len(positive["addresses"]),
            "noisy_positive_control":compact(noisy),"everywhere_positive_field_missing_atom":compact(absent),
            "scope":"Truthful simultaneous neighboring measurements remain a premise; no ROOT inputs or new Collatz coverage."}
    print(json.dumps(report,sort_keys=True,indent=2))


if __name__=="__main__":main()
