"""Exact module-response, two-port obstruction, and induced ternary checks."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import runpy
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = ROOT/'05-knowledge/results/strategy_fibres_absorption_20261113/fibres/check.py'
base = runpy.run_path(str(BASE))
need, counts, profile = (base[k] for k in ('need','counts','profile'))
strong, permutation_count = (base[k] for k in ('strong','permutation_count'))
marked_parity = base['marked_parity']
A = (44,17,2,6,13,30,63)
B = (110,16,2,6,13,30,62)
SOURCE = (44,17,2,6,13,30)
PLUS2 = (110,92,24,112,97,70,4)
PLUS4 = (2,8,11,97,15,87,151,63)


def transitive(n):
    return tuple(((1 << n)-1) ^ ((1 << (i+1))-1) for i in range(n))


def covers(rows):
    return tuple(factorial(k)*d for k,d in enumerate(profile(rows)))


def acyclic(rows, ports):
    remaining = set(ports)
    while remaining:
        sources = [v for v in remaining if not any((rows[w] >> v) & 1 for w in remaining)]
        if not sources:
            return False
        remaining.remove(sources[0])
    return True


def response(rows, ports):
    """Run polynomial, exact when the substituted quotient is acyclic."""
    ports = tuple(ports)
    need(acyclic(rows,ports), 'cyclic port quotient has an unbounded run interface')
    exterior = [v for v in range(len(rows)) if v not in ports]
    full = (1 << len(exterior))-1
    zero = (0,)*len(ports)
    @lru_cache(None)
    def suffix(mask,last):
        out = defaultdict(int)
        if mask == full:
            out[zero] = 1
        choices = [(v,mask | (1 << i),zero) for i,v in enumerate(exterior)
                   if not (mask >> i) & 1]
        choices += [(v,mask,tuple(int(i==j) for j in range(len(ports))))
                    for i,v in enumerate(ports) if v != last]
        for vertex, nextmask, shift in choices:
            if last >= 0 and not (rows[last] >> vertex) & 1:
                continue
            for exponents, multiplicity in suffix(nextmask,vertex).items():
                out[tuple(a+b for a,b in zip(exponents,shift))] += multiplicity
        return dict(out)
    return {e:c for e,c in suffix(0,-1).items() if min(e) >= 1}


def evaluate(poly, modules):
    vectors = [covers(module) for module in modules]
    total = 0
    for exponents, coefficient in poly.items():
        term = coefficient
        for vector,k in zip(vectors,exponents):
            term *= vector[k] if k < len(vector) else 0
        total += term
    return total


def substitute(rows, replacements):
    blocks = [tuple(replacements.get(v,(0,))) for v in range(len(rows))]
    offsets = []
    total = 0
    for block in blocks:
        offsets.append(total)
        total += len(block)
    output = [0]*total
    for v,block in enumerate(blocks):
        for i,row in enumerate(block):
            output[offsets[v]+i] = row << offsets[v]
            for w,other in enumerate(blocks):
                if (rows[v] >> w) & 1:
                    output[offsets[v]+i] |= ((1 << len(other))-1) << offsets[w]
    return tuple(output), offsets


def marked_response(rows,port,root):
    """Coefficient vector for G=sum rooted-start path times complement H.

    Independent exterior-permutation / two-path-gap compiler. The first path
    starts at the fixed exterior root; the second path may be empty.
    """
    need(port != root, 'root must remain an exterior singleton')
    exterior = [v for v in range(len(rows)) if v != port]
    result = defaultdict(int)
    def gaps(path,root_first=False):
        if not path:
            return (0,1)
        forced = 0
        optional = int((rows[path[-1]] >> port) & 1)
        if not root_first:
            optional += int((rows[port] >> path[0]) & 1)
        for left,right in zip(path,path[1:]):
            direct = (rows[left] >> right) & 1
            through = ((rows[left] >> port) & 1) and ((rows[port] >> right) & 1)
            if not direct:
                if not through:
                    return None
                forced += 1
            elif through:
                optional += 1
        return forced,optional
    for path in permutations(exterior):
        if path[0] != root:
            continue
        for cut in range(1,len(path)+1):
            first,second = gaps(path[:cut],True),gaps(path[cut:])
            if first is None or second is None:
                continue
            forced,optional = first[0]+second[0],first[1]+second[1]
            for extra in range(optional+1):
                result[forced+extra] += comb(optional,extra)
    return {(k,):v for k,v in result.items() if k >= 1}


def serialize(poly):
    return [{'runs':list(e),'coefficient':c} for e,c in sorted(poly.items())]


def shifted_second_cover(poly,coordinate):
    # Replace K2=2*K1+L, with L>=0 for a nonsingleton tournament.
    shifted = defaultdict(int)
    for exponents,c in poly.items():
        shifted[exponents] += c
        if exponents[coordinate] == 2:
            lowered = list(exponents)
            lowered[coordinate] = 1
            shifted[tuple(lowered)] += 2*c
    return dict(shifted)


def two_port_audit():
    parity_pairs = {(0,4),(0,6),(1,4),(1,6)}
    positive_when_second_large = {(2,4),(2,6),(3,4),(3,6),(5,6)}
    modules = ((0,),transitive(2),(2,4,1))
    result = []
    literal_checks = 0
    for pair in combinations(range(7),2):
        left,right = response(A,pair),response(B,pair)
        delta = {e:right.get(e,0)-left.get(e,0) for e in left.keys() | right.keys()
                 if right.get(e,0) != left.get(e,0)}
        if pair in parity_pairs:
            category = 'nontrivial first module: difference2mod4'
            for hx,k3x,hy,k2y,k3y in product((1,3),(0,2),(1,3),(0,2),(0,2)):
                x,y = {1:hx,2:2,3:k3x},{1:hy,2:k2y,3:k3y}
                residue = sum(c*x.get(i,0)*y.get(j,0) for (i,j),c in delta.items()) % 4
                need(residue == 2,'two-adic obstruction')
            need(all(c > 0 for (i,j),c in delta.items() if i==1), 'singleton first endpoint')
        elif pair in positive_when_second_large or pair == (4,5):
            coordinate = 0 if pair == (4,5) else 1
            shifted = shifted_second_cover(delta,coordinate)
            need(min(shifted.values()) >= 0 and shifted[(1,1)] > 0,'physical cone positivity')
            category = 'positive after K2>=2H in coordinate'+str(coordinate)
        elif pair == (4,6):
            need(min(delta.values()) > 0,'positive pair')
            category = 'all coefficients positive'
        else:
            need(max(delta.values()) < 0 and delta[(1,2)] < 0 and delta[(2,1)] < 0,
                 'negative pair')
            category = 'all coefficients negative'
        for mx,my in product(modules,repeat=2):
            for rows,poly in ((A,left),(B,right)):
                graph,_ = substitute(rows,dict(zip(pair,(mx,my))))
                need(counts(graph)[-1] == evaluate(poly,(mx,my)), 'bivariate compiler versus literal DP')
                literal_checks += 1
        result.append({'ports':pair,'delta':serialize(delta),'certificate':category})
    return {'pairs':result,'literal_DP_checks':literal_checks}


def embeddings():
    maps = [({2:6,0:4,1:2,3:5,4:3,5:0},PLUS2,6),
            ({2:2,0:6,1:3,3:4,4:5,5:7},PLUS4,2)]
    result = []
    for mapping,frame,port in maps:
        need(all(((SOURCE[u] >> v) & 1) == ((frame[mapping[u]] >> mapping[v]) & 1)
                 for u in range(6) for v in range(6)), 'induced frame embedding')
        newborns = sorted(set(range(len(frame)))-set(mapping.values()))
        result.append({'source_to_target':mapping,'newborns':newborns,'old_arc_changes':0,'module_port':port})
    return result


def ternary_audit():
    source_poly = response(SOURCE,(2,))
    plus2_poly = response(PLUS2,(6,))
    plus4_poly = response(PLUS4,(2,))
    marked_poly = marked_response(SOURCE,2,1)
    need(source_poly == {(1,):25,(2,):12},'source run vector')
    need(plus2_poly == {(1,):77,(2,):36},'plus2 run vector')
    need(plus4_poly == {(1,):79,(2,):36},'plus4 run vector')
    need(marked_poly == {(1,):27,(2,):12},'marked preparation')
    for rows,port,poly in ((SOURCE,2,source_poly),(PLUS2,6,plus2_poly),(PLUS4,2,plus4_poly)):
        gap = base['module_coefficients'](rows,port)
        need({(k,):c for k,c in enumerate(gap) if k >= 1 and c} == poly,
             'independent exterior-gap compiler')
    for rows in (SOURCE,PLUS2,PLUS4):
        need(strong(rows),'strong base frame')
        need(counts(rows)[-1] == permutation_count(rows),'independent base permutation count')
    literal = []
    for m in range(1,9):
        module = transitive(m)
        source,offsets = substitute(SOURCE,{2:module})
        left,_ = substitute(PLUS2,{6:module})
        right,_ = substitute(PLUS4,{2:module})
        h,l,r = (counts(rows)[-1] for rows in (source,left,right))
        need(h == 12*2**m+1 and l == 3*h+2 and r == 3*h+4,'exact ternary family')
        need(all(strong(rows) for rows in (source,left,right)),'inflated strongness')
        h0,e,o = marked_parity(source,offsets[1])
        need(h0 == h and e+o-h == 2,'literal marked source')
        literal.append({'m':m,'source':h,'plus2':l,'plus4':r,'orders':[len(source),len(left),len(right)]})
    cyclic = (2,4,1)
    ch = [counts(substitute(rows,{port:cyclic})[0])[-1]
          for rows,port in ((SOURCE,2),(PLUS2,6),(PLUS4,2))]
    need(ch == [147,447,453],'nontransitive digit scaling hostile')
    base_capabilities = {}
    for name,rows in (('plus2',PLUS2),('plus4',PLUS4)):
        defects = []
        for root in range(len(rows)):
            h,e,o = marked_parity(rows,root)
            defects.append(e+o-h)
        need(2 not in defects,'selected child must fail automatic marked return')
        base_capabilities[name] = defects
    return {'source_rows':SOURCE,'source_port':2,'plus2_rows':PLUS2,'plus2_port':6,
            'plus4_rows':PLUS4,'plus4_port':2,'run_vectors':[serialize(p) for p in (source_poly,plus2_poly,plus4_poly)],
            'marked_source_root':1,'marked_source_run_vector':serialize(marked_poly),
            'embeddings':embeddings(),'literal_family':literal,'cyclic_module_counts':ch,
            'base_child_G_minus_H_by_root':base_capabilities}


def triple_audit():
    ports = (0,1,6)
    left,right = response(A,ports),response(B,ports)
    hx,px = counts(SOURCE)[-1],base['p2'](SOURCE)
    need((hx,px)==(25,304),'old six-core response')
    replacements = {0:(0,),1:transitive(2),6:SOURCE}
    ga,_ = substitute(A,replacements)
    gb,_ = substitute(B,replacements)
    ha,hb = counts(ga)[-1],counts(gb)[-1]
    need((ha,hb)==(2725,8179) and hb == 3*ha+4,'three-port plus4 witness')
    need(strong(gb) and not strong(ga),'three-port SCC boundary')
    need(all(8179 % d for d in range(2,91)),'8179 primality')
    examples = []
    for b in range(1,6):
        for module in ((0,),(2,4,1),SOURCE):
            h,p = counts(module)[-1],base['p2'](module)
            source = (14*3**b-17)*h
            target = 2*(4*2**b-3)*p + (3*2**b-1)*h
            modules = ((0,),transitive(b),module)
            need(evaluate(left,modules)==source and evaluate(right,modules)==target,'closed triple family')
            examples.append({'b':b,'module_H':h,'module_P2':p,'source':source,'target':target})
    return {'source_coefficients':serialize(left),'target_coefficients':serialize(right),
            'formula_examples':examples,'witness':{'source_rows':ga,'target_rows':gb,'H':[ha,hb],
                                                  'orders':[len(ga),len(gb)],'target_prime':True}}


def source_transfer():
    five = (2,12,1,21,7)
    child = (2,12,33,53,7,19)
    incoming = (2,20,9,35,45,7)
    iso = (0,1,2,4,5,3)
    need(all(((child[i] >> j)&1) == ((incoming[iso[i]] >> iso[j])&1)
             for i in range(6) for j in range(6)), 'marked incoming source isomorphism')
    need(all((child[i]&31) == five[i] for i in range(5)), 'five-core induced birth')
    need(response(five,(0,)) == {(1,):9,(2,):5}, 'inherited five-core chart')
    need(marked_response(five,0,1) == {(1,):11,(2,):5}, 'five-core marked capability')
    need(response(child,(0,)) == {(1,):29,(2,):15}, 'five-to-fifteen transfer')
    family = []
    for k in range(1,6):
        source,_ = substitute(five,{0:transitive(k)})
        target,_ = substitute(child,{0:transitive(k)})
        next_source,_ = substitute(incoming,{0:transitive(k)})
        mapping = tuple(range(k))+tuple(k+iso[i]-1 for i in range(1,6))
        need(all(((target[i]>>j)&1) == ((next_source[mapping[i]]>>mapping[j])&1)
                 for i in range(len(target)) for j in range(len(target))), 'inflated marked isomorphism')
        h,ht = counts(source)[-1],counts(target)[-1]
        need((h,ht)==(5*2**k-1,15*2**k-1), 'five-to-fifteen literal family')
        need(strong(source) and strong(target),'strong source transfer')
        family.append({'k':k,'source':h,'target':ht})
    ready,_ = substitute(five,{0:transitive(4)})
    edges = [(i,j) for j in range(8) for i in range(j)]
    minimum = 29
    minimizers = []
    for mapping in permutations(range(8)):
        flips = []
        for i,j in edges:
            if ((PLUS4[i] >> j)&1) != ((ready[mapping[i]] >> mapping[j])&1):
                flips.append((i,j))
                if len(flips) > minimum:
                    break
        if len(flips) < minimum:
            minimum,minimizers = len(flips),[(mapping,flips)]
        elif len(flips) == minimum:
            minimizers.append((mapping,flips))
    need(minimum==2 and len(minimizers)==1,'minimum chosen-target core repair')
    mapping,flips = minimizers[0]
    need(mapping==(2,3,1,4,0,5,6,7) and flips==[(0,3),(5,6)],'exact neutral repair')
    repaired = tuple(sum(1 << j for j in range(8) if (ready[mapping[i]] >> mapping[j])&1)
                     for i in range(8))
    need(repaired==(10,8,11,96,15,23,183,63),'repaired literal rows')
    need(counts(repaired)[-1]==79 and strong(repaired),'repaired exact target')
    before_poly,after_poly = response(PLUS4,(2,)),response(repaired,(2,))
    need(after_poly=={(1,):79,(2,):40} and before_poly=={(1,):79,(2,):36},
         'universal four-K2 repair boundary')
    root = 3
    need(marked_parity(repaired,root)[1:]==[40,41],'literal ready root')
    one_step_counts = []
    for i,j in flips:
        single = list(PLUS4)
        single[i] ^= 1 << j
        single[j] ^= 1 << i
        one_step_counts.append(counts(single)[-1])
    next_cut = 159
    h239 = tuple(row | ((1 << 8) if not (next_cut >> i)&1 else 0)
                 for i,row in enumerate(repaired)) + (next_cut,)
    need(counts(h239)[-1]==239 and strong(h239),'physical repaired-source birth')
    incoming239,_ = substitute(incoming,{0:transitive(4)})
    to_incoming = (2,3,1,4,0,5,7,8,6)
    need(all(((h239[i]>>j)&1) == ((incoming239[to_incoming[i]]>>to_incoming[j])&1)
             for i in range(9) for j in range(9)), 'physical H239 incoming interface')
    word = '10010110'
    frame = []
    for i,color in enumerate(word):
        row = ((1 << 8)-1) ^ ((1 << (i+1))-1)
        if color=='0':
            row |= 1 << 8
        frame.append(row)
    frame.append(sum(1 << i for i,color in enumerate(word) if color=='1'))
    canonical721,_ = substitute(tuple(frame),{4:transitive(4)})
    physical_to_canonical = (6,7,5,11,4,3,0,2,1,8,9,10)
    h721 = tuple(sum(1 << j for j in range(12)
                    if (canonical721[physical_to_canonical[i]] >> physical_to_canonical[j])&1)
                 for i in range(12))
    old_changes = [(i,j) for j in range(9) for i in range(j)
                   if ((h239[i]>>j)&1) != ((h721[i]>>j)&1)]
    need(old_changes==[(5,8)],'incoming plus4 old-arc certificate')
    need(counts(h721)[-1]==721 and strong(h721),'literal composed H721 output')
    return {'inherited_five_core':five,'induced_child':child,'incoming_source':incoming,
            'marked_isomorphism':iso,'family':family,
            'neutral_H79_repair':{'source_rows':PLUS4,'target_rows':repaired,
                                 'minimum_distance_to_A5_TT4_class':minimum,'minimizing_bijections':1,
                                 'bijections_checked':40320,'old_to_A5_TT4':mapping,'flipped_arcs':flips,
                                 'single_flip_counts':one_step_counts,'module_order':[4,2,0,1],'ready_root':root,
                                 'before_response':serialize(before_poly),'after_response':serialize(after_poly)},
            'literal_composition':{'counts':[25,79,79,239,721],'orders':[6,8,8,9,12],
                                   'H239_rows':h239,'H239_newborn_outmask':next_cut,
                                   'H239_to_incoming_S15_TT4':to_incoming,'H721_rows':h721,
                                   'H239_to_H721_old_changes':old_changes,'H721_newborns':[9,10,11]}}


def main():
    result = {'status':'FINITE-EXACT','two_port_obstruction':two_port_audit(),
              'induced_ternary_family':ternary_audit(),'three_port_selector':triple_audit(),
              'actual_source_transfer':source_transfer(),
              'inherited_helper_pin':{'path':str(BASE.relative_to(ROOT)).replace('\\','/'),
                                     'sha256_lf':hashlib.sha256(BASE.read_bytes().replace(b'\r\n',b'\n')).hexdigest()}}
    suffix = '.optimized' if sys.flags.optimize else ''
    (HERE/('receipt'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'two_port_pairs':21,
                      'bivariate_literal_checks':result['two_port_obstruction']['literal_DP_checks'],
                      'ternary_family':result['induced_ternary_family']['literal_family'],
                      'three_port_witness':[2725,8179]}))


if __name__ == '__main__':
    main()
