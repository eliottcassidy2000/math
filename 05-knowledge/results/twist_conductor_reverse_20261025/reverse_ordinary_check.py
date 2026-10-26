"""Normalized ordinary A000568 residues on the exact structural transfer DAG.

Uses the audited transition helpers in modular.py. Graph construction never
uses a truncated residue. The graph retains stage and every exact defect key.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

_spec=spec_from_file_location('ordinary_modular_transition',Path(__file__).resolve().parent.parent/'port_residual_integral_20261025/counting/modular.py')
_bank=module_from_spec(_spec);_spec.loader.exec_module(_bank)


def build_graph(n):
    if type(n) is not int or n<0:raise ValueError('n must be a nonnegative integer')
    top=n if n%2 else n-1
    initial=(top,n,(0,)*((top+1)//2))
    terminal=(-1,0,())
    if n==0:initial=terminal
    order=[];edges={};layer={initial};widths=[]
    for ell in range(top,0,-2):
        widths.append((ell,len(layer)));following=set()
        for node in sorted(layer):
            stage,s,defect=node
            if stage!=ell:raise RuntimeError('stage lost')
            order.append(node);out=[]
            b=_bank.bound(defect)
            selected=defect[(ell-1)//2] if ell<=s else 0
            for k in range(s//ell+1):
                rest=s-ell*k
                if ell==1 and rest:continue
                dn=_bank.next_defect(ell,rest,defect,k)
                bp=_bank.bound(dn)
                J=_bank.intrinsic(ell,k)-(k-k.bit_count())+k*(ell*b-selected)+rest*(k+b-bp)
                if J<0:raise RuntimeError('negative original edge cost')
                destination=(ell-2,rest,dn)
                out.append((destination,J,ell,k,n-s));following.add(destination)
            edges[node]=tuple(out)
        layer=following
    if n and layer!={terminal}:raise RuntimeError('wrong terminal layer')
    order.append(terminal);edges[terminal]=()
    V={terminal:0}
    for node in reversed(order[:-1]):V[node]=min(J+V[w] for w,J,ell,k,u in edges[node])
    B=V[initial]
    prefix={initial:0}
    for node in order:
        if node not in prefix:raise RuntimeError('unreachable structural state')
        if prefix[node]+V[node]<B:raise RuntimeError('invalid prefix lower bound')
        for w,J,ell,k,u in edges[node]:
            if J+V[w]-V[node]<0:raise RuntimeError('negative reduced cost')
            value=prefix[node]+J
            prefix[w]=min(prefix.get(w,value),value)
    return {'n':n,'initial':initial,'terminal':terminal,'order':tuple(order),
            'edges':edges,'suffix':V,'prefix':prefix,'baseline':B,'layer_widths':tuple(widths)}


def normalized_target(n,K,graph=None):
    """Return a(n)/2**B mod 2**K, where B is the termwise minimum valuation."""
    if type(n) is not int or n<0 or type(K) is not int or K<1:
        raise ValueError('n>=0 and K>=1 must be integers')
    graph=build_graph(n) if graph is None else graph
    if graph['n']!=n:raise ValueError('graph belongs to a different source order')
    modulus=1<<K;odd=_bank.odd_factorials(n,modulus)
    V=graph['suffix'];prefix=graph['prefix'];B=graph['baseline']
    values={graph['initial']:1};active_nodes=active_edges=zero_nodes=0
    coefficient_cache={}
    for node in graph['order']:
        if prefix[node]+V[node]-B>=K:continue
        active_nodes+=1
        value=values.get(node,0)
        zero_nodes+=int(value==0)
        for w,J,ell,k,u in graph['edges'][node]:
            reduced=J+V[w]-V[node]
            # This is a structural cost cut, not a residue-based support cut.
            if prefix[node]+V[node]-B+reduced>=K:continue
            active_edges+=1
            key=(ell,k,u)
            if key not in coefficient_cache:
                coefficient_cache[key]=_bank.integral_coeff(odd,u,ell,k,0,modulus)
            T=coefficient_cache[key]
            values[w]=(values.get(w,0)+(T*value<<reduced))%modulus
    result=values.get(graph['terminal'],0)*pow(odd[n],-1,modulus)%modulus
    return result,{'n':n,'normalized_bits':K,'baseline':B,
                   'structural_nodes':len(graph['order']),
                   'structural_edges':sum(map(len,graph['edges'].values())),
                   'active_nodes':active_nodes,'active_edges':active_edges,
                   'retained_zero_residue_nodes':zero_nodes,
                   'odd_coefficient_evaluations':len(coefficient_cache)}
