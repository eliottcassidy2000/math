"""Source floors, exact superlevels, and guarded deadline transport.

All threshold decisions terminate. No initial positive floor is assumed true
by a successful numerical approximation, and no universal support is claimed.
"""
from fractions import Fraction as F
from math import factorial, isqrt
import json

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def nat(n, positive=False):
    if type(n) is not int or n < int(positive):
        raise ValueError("exact natural integer required")


def odd(n):
    nat(n, True)
    if n % 2 != 1:
        raise ValueError("positive odd integer source required")


def floor_parameter(epsilon):
    if type(epsilon) not in (int, F) or not 0 < epsilon <= 1:
        raise ValueError("exact rational floor in (0,1] required")
    return F(epsilon)


def step(n):
    odd(n)
    if n == 1:
        raise ValueError("ROOT has no outgoing first-hit edge")
    v = 3*n+1
    a = (v & -v).bit_length()-1
    return v >> a, a


def weight(length, depth):
    nat(length)
    nat(depth)
    return F(2*factorial(depth)*factorial(length+1),
             factorial(length+depth+2))


def counter_ceiling(epsilon):
    epsilon = floor_parameter(epsilon)
    m = (2*epsilon.denominator)//epsilon.numerator
    return (isqrt(1+4*m)-3)//2


def largest_superlevel_source(epsilon):
    b = counter_ceiling(epsilon)
    return (4**(b+1)-1)//3


def word_weight(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("valuation tuple required")
    if not word:
        return F(1)
    return weight(len(word)-1, sum((a-1)//2 for a in word))


def threshold_receipt(source, epsilon):
    """Decide W(source)>=epsilon, including equality, by bounded actual replay.

    'below' does not mean unrooted. A complete word is retained whenever seen.
    The returned root word is strict, with no ROOT self-loop.
    """
    odd(source)
    epsilon = floor_parameter(epsilon)
    b = counter_ceiling(epsilon)
    n, word = source, []
    for used in range(b+1):
        if n == 1:
            exact = word_weight(tuple(word))
            return {"status": "met" if exact >= epsilon else "below",
                    "source":source, "word":tuple(word), "reached_root":True,
                    "weight":exact, "upper":exact, "counter_ceiling":b}
        if used == b:
            break
        n, a = step(n)
        word.append(a)
    return {"status":"below", "source":source, "word":tuple(word),
            "reached_root":False, "weight":None,
            "upper":F(2,(b+2)*(b+3)), "counter_ceiling":b}


def pointwise_interval(source, steps):
    """Enclose W(source) after a bounded actual forward prefix.

    Returning lower zero means no positive floor has yet been established.
    It does not identify the source as unrooted.
    """
    odd(source)
    nat(steps)
    current, word = source, []
    for used in range(steps+1):
        if current == 1:
            exact = word_weight(tuple(word))
            return {"lower":exact,"upper":exact,"reached_root":True,"word":tuple(word)}
        if used == steps:
            break
        current, a = step(current)
        word.append(a)
    return {"lower":F(0),"upper":F(2,(steps+2)*(steps+3)),
            "reached_root":False,"word":tuple(word)}


def rising(n, count):
    nat(n, True)
    nat(count)
    product = 1
    for j in range(count):
        product *= n+j
    return product


def transport_floor(epsilon, add_length, add_depth):
    """Conditional floor after a VERIFIED counter increment at a nonroot target."""
    epsilon = floor_parameter(epsilon)
    nat(add_length)
    nat(add_depth)
    if epsilon > F(1,3):
        raise ValueError("no nonroot source can have this floor")
    b = counter_ceiling(epsilon)
    return epsilon*F(factorial(add_depth+1)*factorial(add_length+1),
                     rising(b+3,add_length+add_depth))


def deadline_floor(epsilon, budget):
    """Uniform conditional floor when the whole increment r+k is <=budget."""
    epsilon = floor_parameter(epsilon)
    nat(budget)
    if epsilon > F(1,3):
        raise ValueError("nonroot anchor floor required")
    b = counter_ceiling(epsilon)
    r = budget//2
    k = budget-r
    return epsilon*F(factorial(r+1)*factorial(k+1),rising(b+3,budget))


def replay(source, word):
    odd(source)
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("valuation tuple required")
    n = source
    for requested in word:
        n, actual = step(n)
        if actual != requested:
            raise ValueError("source valuation guard failed")
    return n


def guarded_transport(source, word, epsilon):
    """Returns a conditional floor; epsilon must be justified at the actual endpoint.

    If the endpoint is ROOT, its actual source word gives an unconditional
    exact weight instead. Replaying the word verifies source identity/guards.
    """
    endpoint = replay(source,word)
    epsilon = floor_parameter(epsilon)
    if endpoint == 1:
        return {"endpoint":1,"floor":word_weight(word),"conditional":False,
                "r":len(word),"k":sum((a-1)//2 for a in word)}
    r, k = len(word), sum((a-1)//2 for a in word)
    return {"endpoint":endpoint,"floor":transport_floor(epsilon,r,k),
            "conditional":True,"r":r,"k":k}


def control_route(source, cap=512):
    n, word = source, []
    for _ in range(cap+1):
        if n == 1:
            return tuple(word)
        n, a = step(n)
        word.append(a)
    raise ValueError("finite control outside declared cap")


def main():
    report = {}
    floors = (F(1),F(1,2),F(1,3),F(1,6),F(1,12),F(1,60),F(1,420),F(1,1000))
    threshold_controls = 0
    interval_controls = 0
    literal_controls = 0
    for n in range(1,1024,2):
        word = control_route(n)
        exact = word_weight(word)
        if n > 1:
            length, depth = len(word)-1,sum((a-1)//2 for a in word)
            total = length+depth
            need(depth >= 1, "strict nonroot word has sibling fuel")
            need(exact <= F(2,(total+1)*(total+2)), "quadratic counter envelope")
            need(3*n+1 <= 4**(total+1), "source-height counter envelope")
            for cut in range(1,min(len(word),5)):
                prefix = word[:cut]
                endpoint = replay(n,prefix)
                endpoint_weight = word_weight(word[cut:])
                token = guarded_transport(n,prefix,endpoint_weight)
                need(token["conditional"] and token["endpoint"] == endpoint,
                     "actual target and nonroot guard retained")
                need(token["floor"] <= exact, "valid endpoint floor survives actual prefix")
                need(cut+counter_ceiling(endpoint_weight) >= len(word),
                     "structured finite ROOT deadline")
                need(token["r"]+token["k"] <= sum(prefix), "completed unary-bit deadline")
                literal_controls += 1
        for horizon in (0,1,2,4,8,16,32):
            enclosure=pointwise_interval(n,horizon)
            need(enclosure["lower"] <= exact <= enclosure["upper"],
                 "bounded pointwise enclosure")
            need(len(enclosure["word"]) <= horizon, "actual bounded work")
            if enclosure["reached_root"]:
                need(enclosure["lower"] == enclosure["upper"] == exact,
                     "exact enclosure after first ROOT")
            interval_controls += 1
        for epsilon in floors:
            result = threshold_receipt(n,epsilon)
            need((result["status"]=="met") == (exact >= epsilon), "total exact threshold decision")
            need(exact <= result["upper"], "valid upper endpoint including unresolved case")
            if result["status"]=="met":
                need(result["reached_root"] and replay(n,result["word"]) == 1,
                     "positive threshold returns actual first-hit ROOT certificate")
            if not result["reached_root"]:
                need(result["upper"] < epsilon, "strict threshold rejection at finite horizon")
            threshold_controls += 1
    report["finite_sources"] = {"positive_odds_below":1024,"count":512}
    report["threshold_controls"] = threshold_controls
    report["pointwise_interval_controls"] = interval_controls
    report["guarded_prefix_controls"] = literal_controls

    extrema = []
    for epsilon in floors:
        b = counter_ceiling(epsilon)
        need(epsilon*(b+1)*(b+2) <= 2 < epsilon*(b+2)*(b+3),
             "integer-only maximal counter ceiling")
        largest = largest_superlevel_source(epsilon)
        hit = threshold_receipt(largest,epsilon)
        need(hit["status"]=="met", "sharp actual largest-source witness")
        need(threshold_receipt(largest+2,epsilon)["status"]=="below",
             "next odd source lies below threshold")
        extrema.append({"epsilon":str(epsilon),"B":b,
                        "largest_source":largest,"witness_weight":str(hit["weight"])})
    report["exact_superlevel_extrema"] = extrema

    formal_controls = 0
    for epsilon in (F(1,3),F(1,12),F(1,60),F(1,420)):
        b = counter_ceiling(epsilon)
        eligible = [(l,k) for l in range(b+1) for k in range(1,b-l+1)
                    if weight(l,k) >= epsilon]
        need(bool(eligible),"nonempty formal envelope")
        for r in range(5):
            for kadd in range(5):
                transported = transport_floor(epsilon,r,kadd)
                sharp_formal = min(weight(l+r,k+kadd) for l,k in eligible)
                need(transported <= sharp_formal, "closed floor versus exact finite counter envelope")
                for l,k in eligible:
                    ratio=weight(l+r,k+kadd)/weight(l,k)
                    need(ratio == F(rising(k+1,kadd)*rising(l+2,r),
                                    rising(l+k+3,r+kadd)), "exact counter transport")
                formal_controls += 1
        for budget in range(9):
            common = deadline_floor(epsilon,budget)
            for r in range(budget+1):
                for kadd in range(budget-r+1):
                    need(common <= transport_floor(epsilon,r,kadd),
                         "whole-episode floor covers every allowed increment")
    report["formal_transport_envelopes"] = formal_controls

    # Existing source floors can stay unchanged as information is refined.
    # Changing source along inverse siblings is a different operation.
    escape = []
    for j in (0,1,2,4,8,16,32):
        n = 4**j*3+(4**j-1)//3
        word = control_route(n)
        need(word == (1+2*j,4), "actual sibling source with common endpoint5")
        exact = word_weight(word)
        need(exact == weight(1,j+1), "unbounded inverse refinement tends to zero")
        escape.append({"sibling_depth":j,"source":n,"weight":str(exact)})
    report["unbounded_changed_source_refinement"] = escape

    # AMM composition exchange does not preserve a fixed arithmetic endpoint.
    # Words12 and21 have equal counts but carry5 versus7 over common denominator8.
    need((8*4-5) % 9 == 0 and (8*4-7) % 9 != 0, "endpoint4 admits only inverseword12")
    need((8*2-7) % 9 == 0 and (8*2-5) % 9 != 0, "endpoint2 admits only inverseword21")
    need(replay(11,(1,2)) == 13 and replay(9,(2,1)) == 11,
         "actual positive ordered-word controls")
    report["composition_exchange_hostile"] = {
        "word12_endpoint_mod9":4,"word21_endpoint_mod9":2,
        "same_integer_endpoint_possible":False}
    # Scalar recompression is valid but can lose far more than retained bill.
    eps=F(1,3)
    combined=transport_floor(eps,2,0)
    sequential=transport_floor(transport_floor(eps,1,0),1,0)
    need(combined >= sequential, "retained aggregate bill beats this scalar reencoding")
    report["aggregation_example"] = {"anchor_floor":str(eps),
        "aggregate_two_steps":str(combined),"sequential_reencoding":str(sequential)}

    bad = [lambda: counter_ceiling(0.1),lambda: counter_ceiling(True),
           lambda: counter_ceiling(0),lambda: threshold_receipt(2,F(1,3)),
           lambda: threshold_receipt(True,F(1,3)),lambda: transport_floor(F(1,2),1,0),
           lambda: transport_floor(F(1,3),False,0),lambda: deadline_floor(F(1,3),-1),
           lambda: guarded_transport(11,(2,1),F(1,3)),
           lambda: guarded_transport(1,(2,),F(1,3)),
           lambda: pointwise_interval(3,True),lambda: pointwise_interval(3,-1)]
    for job in bad:
        try:
            job()
        except ValueError:
            need(True,"invalid type or guard rejected")
        else:
            raise ValueError("accepted hostile")
    report["malformed_controls"]=len(bad)
    report["exact_checks"]=CHECKS
    report["scope"]="Threshold decisions total; an initial floor for every source remains OPEN."
    print(json.dumps(report,indent=2,sort_keys=True))
    print("PASS")


if __name__=="__main__":
    main()
