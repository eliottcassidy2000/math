#!/usr/bin/env python3
"""Exact terminal closure of the degree-22 full-support BCDEW chart.

The fixed quintic admits two possible proper factor types: a line and its
quartic complement over the degree-five root field, or a quadratic and its
cubic complement over the degree-ten unordered-root-pair field.  For each
field this script lifts the factorization through orders eleven and twelve.

After the torus substitution A=E/C and S=C^2, four independent order-eleven
rows are linear in W and S.  The six Cramer pivots give quadratic Delta
polynomials, two compatibility equations in K[A,D], and a degree-eleven
resultant H(A).  Exact saturation by Delta leaves degree seven.  Each of the
five order-twelve rows supplies two degree-twenty-one projections.  Direct
number-field gcds show that the saturated H and all ten projections generate
the unit ideal on every pivot chart, while the six Delta polynomials also have
gcd one.  Thus the order-eleven/order-twelve terminal system has no torus
point over an algebraic closure of either factor field.

For independent replay, named coefficient-convolution minors are reduced at
simple good residues F_113 (root generator 24) and F_103 (pair generator 61).
Every exact coefficient denominator is checked before reduction, exact and
independently recomputed residue polynomials are compared up to units, and
every selected square determinant is nonzero.  Saturation is always performed
in characteristic zero before reduction: this explicitly avoids the rejected
F_109 root-field artifact where reduction creates an extra Delta gcd.

This is only the full-support terminal calculation in the inherited genuine
nonsplit degree-twenty-two branch; it makes no claim about split/even edges,
integral raising, other degree branches, JC(2), or DC(2).
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_bcd_weighted_hensel_kummer_thm2636.py"
BASE_SHA256 = "0866a29f665aedc6d2c226f35943852e56907ff821e705a0dbca2651e71fa15c"
require(
    hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
    "audited THM-2636 dependency changed",
)
spec = importlib.util.spec_from_file_location("base", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load base")
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)

t, v, zeta = base.t, base.v, base.zeta
c, d, e, w = base.c, base.d, base.e, base.w


class GenericHensel(base.HenselSystem):
    """The audited two-parameter engine with arbitrary exponent-tuple length."""

    def __init__(self, name, field, q0, s0, params):
        self.params = tuple(params)
        self.zero_exp = (0,) * len(self.params)
        super().__init__(name, field, q0, s0)

    def pc(self, terms=None):
        result = {}
        for exponent, coefficient in (terms or {}).items():
            require(len(exponent) == len(self.params), "bad exponent length")
            reduced = self.field.reduce(coefficient)
            if not self.field.zero(reduced):
                result[tuple(exponent)] = reduced
        return result

    def pc_const(self, coefficient):
        return self.pc({self.zero_exp: s.sympify(coefficient)})

    def pc_mul(self, left, right):
        result = {}
        for exponent1, coefficient1 in left.items():
            for exponent2, coefficient2 in right.items():
                exponent = tuple(a + b for a, b in zip(exponent1, exponent2))
                result[exponent] = self.field.reduce(
                    result.get(exponent, 0) + coefficient1 * coefficient2
                )
        return self.pc(result)

    def v_divmod_monic(self, dividend, divisor):
        divisor = self.vtrim(divisor)
        require(divisor[-1] == {self.zero_exp: s.Integer(1)}, "divisor not monic")
        remainder = self.vtrim(dividend[:])
        quotient = [{} for _ in range(max(1, len(remainder) - len(divisor) + 1))]
        while len(remainder) >= len(divisor) and not (
            len(remainder) == 1 and not remainder[0]
        ):
            shift = len(remainder) - len(divisor)
            leading = remainder[-1]
            quotient[shift] = self.pc_add(quotient[shift], leading)
            subtraction = [{} for _ in range(shift)] + [
                self.pc_mul(coefficient, leading) for coefficient in divisor
            ]
            remainder = self.v_add(remainder, self.v_neg(subtraction))
        return self.vtrim(quotient), self.vtrim(remainder)

    def inverse_s0_mod_q0(self, q0_field, s0_field):
        q0 = self.v_from_field(q0_field)
        s0 = self.v_from_field(s0_field)
        columns = []
        for exponent in range(self.k):
            monomial = [{} for _ in range(exponent)] + [self.pc_const(1)]
            remainder = self.v_remainder(self.v_mul(s0, monomial), q0)
            entries = []
            for index in range(self.k):
                coefficient = remainder[index] if index < len(remainder) else {}
                require(set(coefficient).issubset({self.zero_exp}), "fixed inverse has params")
                entries.append(coefficient.get(self.zero_exp, s.Integer(0)))
            columns.append(entries)
        if self.k == 1:
            inverse_coefficients = [self.field.inverse(columns[0][0])]
        else:
            a00, a10 = columns[0]
            a01, a11 = columns[1]
            determinant = self.field.reduce(a00 * a11 - a01 * a10)
            inverse_coefficients = [
                self.field.reduce(a11 * self.field.inverse(determinant)),
                self.field.reduce(-a10 * self.field.inverse(determinant)),
            ]
        inverse = self.v_from_field(inverse_coefficients)
        require(
            self.v_equal(self.v_remainder(self.v_mul(s0, inverse), q0), self.v_from_field([1])),
            "fixed inverse failed",
        )
        return inverse

    def expression_to_vpoly(self, expression):
        poly = s.Poly(expression, v, *self.params, domain=s.QQ)
        degree_v = poly.degree(v)
        result = [{} for _ in range(max(0, degree_v) + 1)]
        for monomial, coefficient in poly.terms():
            iv, *parameter_exponents = monomial
            result[iv] = self.pc_add(
                result[iv], self.pc({tuple(parameter_exponents): coefficient})
            )
        return self.vtrim(result)

    def lift(self, r_polys, through):
        require(self.v_equal(self.v_mul(self.q0, self.s0), r_polys[0]), "fixed product failed")
        qs, ss = [self.q0], [self.s0]
        zero = self.v_from_field([0])
        for n in range(1, through + 1):
            rn = r_polys[n] if n < len(r_polys) else zero
            convolution = zero
            for index in range(1, n):
                convolution = self.v_add(convolution, self.v_mul(qs[index], ss[n-index]))
            residual = self.v_add(rn, self.v_neg(convolution))
            qn = self.v_remainder(self.v_mul(residual, self.inv_s0_mod_q0), self.q0)
            numerator = self.v_add(residual, self.v_neg(self.v_mul(qn, self.s0)))
            sn, remainder = self.v_divmod_monic(numerator, self.q0)
            require(all(not coefficient for coefficient in remainder), f"division failed {n}")
            qs.append(qn)
            ss.append(sn)
            reconstructed = zero
            for index in range(n + 1):
                reconstructed = self.v_add(reconstructed, self.v_mul(qs[index], ss[n-index]))
            require(self.v_equal(reconstructed, rn), f"product failed {n}")
        return qs, ss

    def field_matrix_rank(self, matrix):
        a = [[self.field.reduce(x) for x in row] for row in matrix]
        rows = len(a)
        cols = len(a[0]) if rows else 0
        pivots = []
        r = 0
        for col in range(cols):
            pivot = next((i for i in range(r, rows) if not self.field.zero(a[i][col])), None)
            if pivot is None:
                continue
            a[r], a[pivot] = a[pivot], a[r]
            inv = self.field.inverse(a[r][col])
            a[r] = [self.field.reduce(x * inv) for x in a[r]]
            for i in range(rows):
                if i == r or self.field.zero(a[i][col]):
                    continue
                scale = a[i][col]
                a[i] = [self.field.reduce(x - scale*y) for x, y in zip(a[i], a[r])]
            pivots.append(col)
            r += 1
            if r == rows:
                break
        return r, pivots, a

charts = [
    ("BCDEW", base.R_universal, (c,d,e,w)),
]


for name, raw, params in charts:
    leading = s.Poly(raw, v).coeff_monomial(v**5)
    require(leading != 0 and not leading.free_symbols, f"{name}: bad leading")
    P = s.cancel(raw / leading)
    require(s.Poly(P, v).LC() == 1, f"{name}: not monic")
    D = s.Poly(P, t).degree()
    require(D == 10 and len(s.Poly(raw,t,v,*params).terms()) == 60,
            "full-support eliminant shape changed")
    coefficients = [s.Poly(P, t).coeff_monomial(t**n) for n in range(D+1)]
    require(s.expand(coefficients[0] - base.P5_expr) == 0, f"{name}: fixed section changed")
    print(f"CHART {name}: terms={len(s.Poly(raw,t,v,*params).terms())} D={D} N={D+1} params={params}", flush=True)
    for field_name, field, q0, s0 in (
        ("root", base.root_field, base.root_q0, base.root_s0),
        ("pair", base.pair_field, base.pair_q0, base.pair_s0),
    ):
        system = GenericHensel(name+"_"+field_name, field, q0, s0, params)
        r_polys = [system.expression_to_vpoly(x) for x in coefficients]
        qs, ss = system.lift(r_polys, D+2)
        order_equations = {}
        for n in (D+1,D+2):
            equations = qs[n] + ss[n]
            order_equations[n] = equations
            require(len(equations) == 5, f"{name}/{field_name}/{n}: equation count")
            support = sorted(set().union(*(set(x) for x in equations)))
            matrix = [[equation.get(m,0) for m in support] for equation in equations]
            rank, pivots, rref = system.field_matrix_rank(matrix)
            print(
                f" {field_name} order={n} support_count={len(support)} rank={rank} "
                f"pivots={pivots} support={support}", flush=True,
            )
        if name == "BCDEW":
            support11 = sorted(set().union(*(set(x) for x in order_equations[D+1])))
            support12 = sorted(set().union(*(set(x) for x in order_equations[D+2])))
            require(
                support11 == [
                    (0,0,1,0),(0,0,1,1),(0,1,1,0),
                    (1,0,0,0),(1,0,0,1),(1,1,0,0),
                    (1,2,0,0),(2,0,1,0),(3,0,0,0),
                ],
                "order-eleven monomial support changed",
            )
            require(
                support12 == [
                    (0,0,0,0),(0,0,0,1),(0,0,0,2),(0,0,2,0),
                    (0,1,0,0),(0,1,0,1),(0,2,0,0),(0,3,0,0),
                    (1,0,1,0),(1,1,1,0),(2,0,0,0),(2,0,0,1),
                    (2,1,0,0),(4,0,0,0),
                ],
                "order-twelve monomial support changed",
            )
            matrix11 = [[equation.get(m,0) for m in support11]
                        for equation in order_equations[D+1]]
            selected_rows = None
            for candidate_rows in combinations(range(5),4):
                candidate = [matrix11[i] for i in candidate_rows]
                candidate_rank,_,_ = system.field_matrix_rank(candidate)
                if candidate_rank == 4:
                    selected_rows = candidate_rows
                    break
            require(selected_rows is not None,"full support terminal rank below four")
            require(selected_rows == (0,1,2,3),
                    "canonical four-row terminal chart changed")
            selected_raw = [matrix11[i] for i in selected_rows]

            class FastField:
                def __init__(self, slow_field):
                    self.generator=slow_field.generator
                    self.n=slow_field.modulus.degree()
                    self.mod=tuple(self._fraction(slow_field.modulus.nth(i))
                                   for i in range(self.n))
                    self.zero=(Fraction(0),)*self.n
                    self.one=(Fraction(1),)+(Fraction(0),)*(self.n-1)
                    self.inverse_cache={}

                @staticmethod
                def _fraction(value):
                    value=s.Rational(value)
                    return Fraction(int(value.p),int(value.q))

                def from_sympy(self,value):
                    poly=s.Poly(s.expand(value),self.generator,domain=s.QQ)
                    return tuple(self._fraction(poly.nth(i)) for i in range(self.n))

                def iszero(self,value):
                    return value==self.zero

                def add(self,left,right):
                    return tuple(a+b for a,b in zip(left,right))

                def neg(self,value):
                    return tuple(-a for a in value)

                def scale(self,value,scalar):
                    scalar=Fraction(scalar)
                    return tuple(scalar*a for a in value)

                def sub(self,left,right):
                    return tuple(a-b for a,b in zip(left,right))

                def mul(self,left,right):
                    work=[Fraction(0)]*(2*self.n-1)
                    for i,a in enumerate(left):
                        if a:
                            for j,b in enumerate(right):
                                work[i+j]+=a*b
                    for degree in range(2*self.n-2,self.n-1,-1):
                        coefficient=work[degree]
                        if coefficient:
                            shift=degree-self.n
                            for i,m in enumerate(self.mod):
                                work[shift+i]-=coefficient*m
                    return tuple(work[:self.n])

                def inv(self,value):
                    require(not self.iszero(value),"fast inverse of zero")
                    if value in self.inverse_cache:
                        return self.inverse_cache[value]
                    # Multiplication matrix columns value*x^j, solve M b=1.
                    matrix=[]
                    columns=[]
                    for j in range(self.n):
                        basis=tuple(Fraction(1) if i==j else Fraction(0)
                                    for i in range(self.n))
                        columns.append(self.mul(value,basis))
                    for i in range(self.n):
                        matrix.append([columns[j][i] for j in range(self.n)]
                                      +[Fraction(1) if i==0 else Fraction(0)])
                    for column in range(self.n):
                        pivot=next(i for i in range(column,self.n)
                                   if matrix[i][column])
                        matrix[column],matrix[pivot]=matrix[pivot],matrix[column]
                        scale=matrix[column][column]
                        matrix[column]=[x/scale for x in matrix[column]]
                        for row in range(self.n):
                            if row==column or not matrix[row][column]:
                                continue
                            scale=matrix[row][column]
                            matrix[row]=[x-scale*y for x,y in zip(matrix[row],matrix[column])]
                    result=tuple(matrix[i][-1] for i in range(self.n))
                    require(self.mul(value,result)==self.one,"fast inverse failed")
                    self.inverse_cache[value]=result
                    return result

            fast=FastField(field)
            selected=[[fast.from_sympy(x) for x in row] for row in selected_raw]

            def kclean(poly):
                return {m:x for m,x in poly.items() if not fast.iszero(x)}

            def kadd(left,right):
                result = dict(left)
                for m,x in right.items():
                    result[m] = fast.add(result.get(m,fast.zero),x)
                return kclean(result)

            def kneg(poly):
                return kclean({m:fast.neg(x) for m,x in poly.items()})

            def kmul(left,right):
                result = {}
                for (a,d0),x in left.items():
                    for (b,d1),y in right.items():
                        m=(a+b,d0+d1)
                        result[m]=fast.add(result.get(m,fast.zero),fast.mul(x,y))
                return kclean(result)

            def kdet(matrix0):
                size=len(matrix0)
                if size==0:
                    return {(0,0):fast.one}
                total={}
                for col in range(size):
                    minor=[row[:col]+row[col+1:] for row in matrix0[1:]]
                    term=kmul(matrix0[0][col],kdet(minor))
                    total=kadd(total,term if col%2==0 else kneg(term))
                return total

            def linA(constant,linear):
                return kclean({(0,0):constant,(1,0):linear})

            # support11 = E,EW,DE,C,CW,CD,CD^2,C^2E,C^3.
            components=[]
            for row in selected:
                wpoly=linA(row[4],row[1])
                spoly=linA(row[8],row[7])
                bpoly=kclean({
                    (0,0):row[3], (1,0):row[0], (0,1):row[5],
                    (1,1):row[2], (0,2):row[6],
                })
                components.append((wpoly,spoly,bpoly))
            delta=kadd(kmul(components[0][0],components[1][1]),
                       kneg(kmul(components[1][0],components[0][1])))
            compat=[]
            for index in (2,3):
                compat.append(kdet([
                    list(components[0]),list(components[1]),list(components[index])
                ]))
            print(
                f"  terminal_selected_rows={selected_rows} delta_support={sorted(delta)} "
                f"compat_supports={[sorted(poly) for poly in compat]}",flush=True,
            )

            def coefficient_in_D(poly, exponent):
                return {(a,0):x for (a,d0),x in poly.items() if d0==exponent}

            def resultant_in_D(left,right):
                dl=max(d0 for a,d0 in left)
                dr=max(d0 for a,d0 in right)
                lc=[coefficient_in_D(left,k) for k in range(dl,-1,-1)]
                rc=[coefficient_in_D(right,k) for k in range(dr,-1,-1)]
                size=dl+dr
                syl=[]
                for shift in range(dr):
                    syl.append([{}]*shift+lc+[{}]*(size-shift-len(lc)))
                for shift in range(dl):
                    syl.append([{}]*shift+rc+[{}]*(size-shift-len(rc)))
                return kdet(syl)

            terminal_resultant=resultant_in_D(compat[0],compat[1])
            require(all(d0==0 for a,d0 in terminal_resultant),"resultant retained D")
            print(
                f"  terminal_resultant_A_degree={max(a for a,d0 in terminal_resultant)} "
                f"terminal_resultant_terms={len(terminal_resultant)} "
                f"terminal_resultant_support={sorted(terminal_resultant)}",flush=True,
            )

            def udegree(poly):
                return max((a for (a,d0) in poly),default=-1)

            def udivmod(dividend,divisor):
                require(divisor,"zero univariate divisor")
                quotient={}
                remainder=dict(dividend)
                divisor_degree=udegree(divisor)
                divisor_lc=divisor[(divisor_degree,0)]
                divisor_lc_inv=fast.inv(divisor_lc)
                while remainder and udegree(remainder)>=divisor_degree:
                    shift=udegree(remainder)-divisor_degree
                    scalar=fast.mul(remainder[(udegree(remainder),0)],divisor_lc_inv)
                    quotient[(shift,0)]=scalar
                    subtraction={(a+shift,0):fast.mul(x,scalar)
                                 for (a,d0),x in divisor.items()}
                    remainder=kadd(remainder,kneg(subtraction))
                return kclean(quotient),kclean(remainder)

            def ugcd(left,right):
                left,right=kclean(left),kclean(right)
                while right:
                    _,remainder=udivmod(left,right)
                    left,right=right,remainder
                if not left:
                    return {}
                degree=udegree(left)
                inverse=fast.inv(left[(degree,0)])
                return kclean({m:fast.mul(x,inverse) for m,x in left.items()})

            def d_slice(poly,exponent):
                return {(a,0):x for (a,d0),x in poly.items() if d0==exponent}

            w_numerator=kadd(kmul(components[0][1],components[1][2]),
                             kneg(kmul(components[1][1],components[0][2])))
            s_numerator=kadd(kmul(components[1][0],components[0][2]),
                             kneg(kmul(components[0][0],components[1][2])))
            projected_bad=[]
            if field_name == "root":
                infinity_gcd=ugcd(d_slice(compat[0],2),d_slice(compat[1],2))
                d_zero_gcd=ugcd(d_slice(compat[0],0),d_slice(compat[1],0))
                delta_gcd=ugcd(terminal_resultant,delta)
                terminal_derivative={(a-1,0):fast.scale(x,a)
                                     for (a,d0),x in terminal_resultant.items() if a}
                repeated_gcd=ugcd(terminal_resultant,terminal_derivative)
                print(
                    f"  bad_gcd_degrees=infinity:{udegree(infinity_gcd)},D0:{udegree(d_zero_gcd)},"
                    f"delta:{udegree(delta_gcd)},repeated:{udegree(repeated_gcd)} "
                    f"H_constant_nonzero={not fast.iszero(terminal_resultant[(0,0)])}",
                    flush=True,
                )
                for label,badpoly in (("W0",w_numerator),("S0",s_numerator)):
                    left_bad=resultant_in_D(compat[0],badpoly)
                    right_bad=resultant_in_D(compat[1],badpoly)
                    common_bad=ugcd(left_bad,right_bad)
                    projected_bad.append(common_bad)
                    print(
                        f"  {label}_resultant_gcd_degree={udegree(common_bad)} "
                        f"left_degree={udegree(left_bad)} right_degree={udegree(right_bad)}",flush=True,
                    )

            def kpow(poly,exponent):
                result={(0,0):fast.one}
                factor=poly
                power=exponent
                while power:
                    if power&1:
                        result=kmul(result,factor)
                    factor=kmul(factor,factor)
                    power//=2
                return result

            n12_numerators=[]
            for equation in order_equations[D+2]:
                transformed={}
                for (ic,id0,ie,iw),coefficient in equation.items():
                    require((ic+ie)%2==0,"N12 odd C power after E=A*C")
                    is_power=(ic+ie)//2
                    denominator_power=iw+is_power
                    require(denominator_power<=2,"N12 denominator exceeds Delta^2")
                    term={(ie,id0):fast.from_sympy(coefficient)}
                    term=kmul(term,kpow(w_numerator,iw))
                    term=kmul(term,kpow(s_numerator,is_power))
                    term=kmul(term,kpow(delta,2-denominator_power))
                    transformed=kadd(transformed,term)
                n12_numerators.append(transformed)

            # Exact optimized K[A,D] arithmetic using SymPy's ANP domain.
            # This is the proof path; the residue calculations below are only
            # independent finite-exact controls.
            algebraic_field=QQ.alg_field_from_poly(field.modulus)
            polynomial_ring,D_ring,A_ring_generator=ring("D,A",algebraic_field)

            def to_anp(value):
                coefficients=[QQ.convert(x) for x in reversed(value)]
                while len(coefficients)>1 and not coefficients[0]:
                    coefficients.pop(0)
                return algebraic_field.new(coefficients)

            def to_ring2(poly):
                return polynomial_ring.from_dict({(d0,a):to_anp(x)
                                                  for (a,d0),x in poly.items()})

            def to_ringA(poly,target_ring):
                return target_ring.from_dict({(a,):to_anp(x)
                                              for (a,d0),x in poly.items()})

            def exact_saturate(poly,bad):
                result=poly
                while True:
                    common=result.gcd(bad)
                    if common.degree()<=0:
                        return result
                    result=result.exquo(common)

            exact_delta_polys=[]
            exact_pivot_results=[]
            common_A_ring=None
            for pivot_pair in combinations(range(4),2):
                other_rows=tuple(index for index in range(4) if index not in pivot_pair)
                first,second=pivot_pair
                delta_pair=kadd(
                    kmul(components[first][0],components[second][1]),
                    kneg(kmul(components[second][0],components[first][1])),
                )
                compat_pair=[]
                for index in other_rows:
                    compat_pair.append(kdet([
                        list(components[first]),list(components[second]),
                        list(components[index]),
                    ]))
                compat_ring=[to_ring2(poly) for poly in compat_pair]
                h_pair=compat_ring[0].resultant(compat_ring[1])
                if common_A_ring is None:
                    common_A_ring=h_pair.ring
                delta_A=to_ringA(delta_pair,h_pair.ring)
                exact_delta_polys.append(to_ringA(delta_pair,common_A_ring))
                physical_pair=exact_saturate(h_pair,delta_A)
                require(h_pair.degree()==11 and physical_pair.degree()==7,
                        f"exact terminal degree changed at {pivot_pair}")
                nw_pair=kadd(
                    kmul(components[first][1],components[second][2]),
                    kneg(kmul(components[second][1],components[first][2])),
                )
                ns_pair=kadd(
                    kmul(components[second][0],components[first][2]),
                    kneg(kmul(components[first][0],components[second][2])),
                )
                row_degrees=[]
                exact_projection_pairs=[]
                for equation in order_equations[D+2]:
                    numerator_pair={}
                    for (ic,id0,ie,iw),coefficient in equation.items():
                        require((ic+ie)%2==0,"exact N12 odd C power")
                        is_power=(ic+ie)//2
                        denominator_power=iw+is_power
                        term={(ie,id0):fast.from_sympy(coefficient)}
                        term=kmul(term,kpow(nw_pair,iw))
                        term=kmul(term,kpow(ns_pair,is_power))
                        term=kmul(term,kpow(delta_pair,2-denominator_power))
                        numerator_pair=kadd(numerator_pair,term)
                    numerator_ring=to_ring2(numerator_pair)
                    left_projection=compat_ring[0].resultant(numerator_ring)
                    right_projection=compat_ring[1].resultant(numerator_ring)
                    exact_projection_pairs.append((left_projection,right_projection))
                    row_degrees.append((left_projection.degree(),right_projection.degree()))
                require(all(pair==(21,21) for pair in row_degrees),
                        f"exact order-twelve projection degree changed at {pivot_pair}")
                exact_unit_gcd=physical_pair
                for exact_pair in exact_projection_pairs:
                    for exact_projection in exact_pair:
                        exact_unit_gcd=exact_unit_gcd.gcd(exact_projection)
                exact_pivot_results.append({
                    "pivot":pivot_pair,
                    "H":h_pair,
                    "physical_H":physical_pair,
                    "delta":delta_A,
                    "projection_pairs":tuple(exact_projection_pairs),
                    "projection_degrees":tuple(row_degrees),
                    "unit_gcd":exact_unit_gcd,
                })
                require(exact_unit_gcd.degree()==0,
                        f"exact terminal/N12 ideal not yet the unit ideal at {pivot_pair}")
                print(
                    f"  EXACT pivot_pair={pivot_pair} H_degree={h_pair.degree()} "
                    f"saturated_H_degree={physical_pair.degree()} "
                    f"N12_projection_degree_pairs={row_degrees} "
                    f"all_projection_gcd_degree={exact_unit_gcd.degree()}",flush=True,
                )
            exact_all_delta_gcd=exact_delta_polys[0]
            for delta_A in exact_delta_polys[1:]:
                exact_all_delta_gcd=exact_all_delta_gcd.gcd(delta_A)
            require(exact_all_delta_gcd.degree()==0,
                    "exact Cramer pivots have a common root")
            print(
                f"  EXACT all_pivot_deltas_gcd_degree={exact_all_delta_gcd.degree()} "
                f"exact_resultants_constructed=True",flush=True,
            )
            # The exact gcd above is the primary characteristic-zero proof.
            # Independently, reduce named exact coefficient matrices at fixed
            # simple residues.  A nonzero residue determinant certifies that
            # the corresponding exact determinant is nonzero.
            def elem_mod(value,prime,residue):
                total=0
                power=1
                for coefficient in value:
                    require(coefficient.denominator%prime,"bad residue denominator")
                    total=(total+coefficient.numerator*pow(coefficient.denominator,-1,prime)*power)%prime
                    power=power*residue%prime
                return total

            all_elements=[]
            for poly in compat+[delta,w_numerator,s_numerator]+n12_numerators:
                all_elements.extend(poly.values())
            mod,residue_root={"root":(113,24),"pair":(103,61)}[field_name]
            require(
                all(
                    coefficient.denominator % mod
                    for value in all_elements + [fast.mod]
                    for coefficient in value
                ),
                "selected residue has a zero coefficient denominator",
            )
            modulus_coefficients=[
                coefficient.numerator
                * pow(coefficient.denominator,-1,mod) % mod
                for coefficient in fast.mod
            ]
            modulus_value=(pow(residue_root,fast.n,mod)+sum(
                modulus_coefficients[index]*pow(residue_root,index,mod)
                for index in range(fast.n)
            ))%mod
            modulus_derivative=(fast.n*pow(residue_root,fast.n-1,mod)+sum(
                index*modulus_coefficients[index]
                * pow(residue_root,index-1,mod)
                for index in range(1,fast.n)
            ))%mod
            require(modulus_value==0 and modulus_derivative!=0,
                    "selected residue is not a simple modulus root")

            def mclean(poly):
                return {m:x%mod for m,x in poly.items() if x%mod}

            def mconvert(poly):
                return mclean({m:elem_mod(x,mod,residue_root) for m,x in poly.items()})

            def madd(left,right):
                result=dict(left)
                for monomial,value in right.items():
                    result[monomial]=(result.get(monomial,0)+value)%mod
                return mclean(result)

            def mneg(poly):
                return mclean({m:-x for m,x in poly.items()})

            def mmul(left,right):
                result={}
                for (a,d0),x in left.items():
                    for (b,d1),y in right.items():
                        monomial=(a+b,d0+d1)
                        result[monomial]=(result.get(monomial,0)+x*y)%mod
                return mclean(result)

            def mdet(matrix0):
                size=len(matrix0)
                if size==0:
                    return {(0,0):1}
                total={}
                for column in range(size):
                    minor=[row[:column]+row[column+1:] for row in matrix0[1:]]
                    term=mmul(matrix0[0][column],mdet(minor))
                    total=madd(total,term if column%2==0 else mneg(term))
                return total

            def mcoeffD(poly,exponent):
                return {(a,0):x for (a,d0),x in poly.items() if d0==exponent}

            def mresultantD(left,right):
                dl=max(d0 for a,d0 in left)
                dr=max(d0 for a,d0 in right)
                lc=[mcoeffD(left,k) for k in range(dl,-1,-1)]
                rc=[mcoeffD(right,k) for k in range(dr,-1,-1)]
                size=dl+dr
                syl=[]
                for shift in range(dr):
                    syl.append([{}]*shift+lc+[{}]*(size-shift-len(lc)))
                for shift in range(dl):
                    syl.append([{}]*shift+rc+[{}]*(size-shift-len(rc)))
                return mdet(syl)

            def mudegree(poly):
                return max((a for a,d0 in poly),default=-1)

            def mudivmod(left,right):
                quotient={}
                remainder=dict(left)
                dr=mudegree(right)
                inverse=pow(right[(dr,0)],-1,mod)
                while remainder and mudegree(remainder)>=dr:
                    shift=mudegree(remainder)-dr
                    scalar=remainder[(mudegree(remainder),0)]*inverse%mod
                    quotient[(shift,0)]=scalar
                    subtraction={(a+shift,0):x*scalar%mod
                                 for (a,d0),x in right.items()}
                    remainder=madd(remainder,mneg(subtraction))
                return mclean(quotient),mclean(remainder)

            def mugcd(left,right):
                while right:
                    _,remainder=mudivmod(left,right)
                    left,right=right,remainder
                if not left:
                    return {}
                degree=mudegree(left)
                inverse=pow(left[(degree,0)],-1,mod)
                return mclean({m:x*inverse%mod for m,x in left.items()})

            compat_mod=[mconvert(poly) for poly in compat]
            print(f"  residue_certificate=F_{mod}:generator={residue_root}:simple=True",flush=True)
            terminal_mod=mconvert(terminal_resultant)
            delta_mod=mconvert(delta)

            def msaturate(poly,bad):
                result=dict(poly)
                while True:
                    common=mugcd(result,bad)
                    if mudegree(common)<=0:
                        return result
                    quotient,remainder=mudivmod(result,common)
                    require(not remainder,"modular saturation division failed")
                    result=quotient

            physical_terminal_mod=msaturate(terminal_mod,delta_mod)
            if projected_bad:
                bad_union_mod=mmul(mconvert(projected_bad[0]),mconvert(projected_bad[1]))
                torus_bad_cover=mugcd(physical_terminal_mod,bad_union_mod)
                bad_cover_text=str(mudegree(torus_bad_cover))
            else:
                bad_cover_text="not_computed"
            print(
                f"  mod_terminal_degree={mudegree(terminal_mod)} "
                f"mod_delta_saturated_degree={mudegree(physical_terminal_mod)} "
                f"mod_W_or_S_bad_cover_degree={bad_cover_text}",flush=True,
            )
            aggregate12=dict(physical_terminal_mod)
            for index,numerator12 in enumerate(n12_numerators):
                numerator_mod=mconvert(numerator12)
                left12=mresultantD(compat_mod[0],numerator_mod)
                right12=mresultantD(compat_mod[1],numerator_mod)
                common12=mugcd(left12,right12)
                aggregate12=mugcd(aggregate12,common12)
                print(
                    f"  N12_row={index} numerator_support={len(numerator12)} "
                    f"D_degree={max(d0 for a,d0 in numerator12)} "
                    f"mod_common_gcd_degree={mudegree(common12)} "
                    f"mod_left_A_degree={mudegree(left12)} mod_right_A_degree={mudegree(right12)}",flush=True,
                )
            print(f"  all_N12_rows_terminal_projection_gcd_degree={mudegree(aggregate12)}",flush=True)

            def mpow(poly,exponent):
                result={(0,0):1}
                factor=poly
                power=exponent
                while power:
                    if power&1:
                        result=mmul(result,factor)
                    factor=mmul(factor,factor)
                    power//=2
                return result

            def anp_coefficient_mod(value):
                total=0
                for coefficient in value.rep:
                    numerator=int(coefficient.numerator)
                    denominator=int(coefficient.denominator)
                    require(denominator%mod,"exact ANP certificate has bad denominator")
                    total=(total*residue_root
                           +numerator*pow(denominator,-1,mod))%mod
                return total

            def exact_univariate_mod(poly):
                return mclean({(monomial[0],0):anp_coefficient_mod(coefficient)
                               for monomial,coefficient in poly.items()})

            def mmonic(poly):
                degree=mudegree(poly)
                require(degree>=0,"cannot monic-normalize zero polynomial")
                inverse=pow(poly[(degree,0)],-1,mod)
                return mclean({m:value*inverse%mod for m,value in poly.items()})

            def independent_row_labels(rows,labels,column_count):
                basis={}
                selected=[]
                for original,label in zip(rows,labels):
                    row=[value%mod for value in original]
                    for pivot in sorted(basis):
                        basis_row=basis[pivot]
                        if row[pivot]:
                            scale=row[pivot]
                            row=[(x-scale*y)%mod for x,y in zip(row,basis_row)]
                    pivot=next((index for index,value in enumerate(row) if value),None)
                    if pivot is None:
                        continue
                    inverse=pow(row[pivot],-1,mod)
                    row=[value*inverse%mod for value in row]
                    for old_pivot,basis_row in list(basis.items()):
                        if basis_row[pivot]:
                            scale=basis_row[pivot]
                            basis[old_pivot]=[(x-scale*y)%mod
                                              for x,y in zip(basis_row,row)]
                    basis[pivot]=row
                    selected.append(label)
                    if len(basis)==column_count:
                        break
                return selected

            def determinant_mod(matrix):
                work=[[value%mod for value in row] for row in matrix]
                size=len(work)
                determinant=1
                for column in range(size):
                    pivot=next((row for row in range(column,size)
                                if work[row][column]),None)
                    require(pivot is not None,"selected modular minor is singular")
                    if pivot!=column:
                        work[column],work[pivot]=work[pivot],work[column]
                        determinant=(-determinant)%mod
                    pivot_value=work[column][column]
                    determinant=determinant*pivot_value%mod
                    inverse=pow(pivot_value,-1,mod)
                    work[column]=[value*inverse%mod for value in work[column]]
                    for row in range(column+1,size):
                        if not work[row][column]:
                            continue
                        scale=work[row][column]
                        work[row]=[(x-scale*y)%mod
                                   for x,y in zip(work[row],work[column])]
                require(determinant,"selected modular determinant vanished")
                return determinant

            def convolution_rank_certificate(polys,poly_labels):
                max_degree=max(mudegree(poly) for poly in polys)
                for bound in range(max_degree,max_degree+41):
                    rows=[]
                    labels=[]
                    for poly,label in zip(polys,poly_labels):
                        degree=mudegree(poly)
                        for shift in range(bound-degree+1):
                            row=[0]*(bound+1)
                            for (exponent,d0),value in poly.items():
                                row[exponent+shift]=value
                            rows.append(row)
                            labels.append((label,shift))
                    selected=independent_row_labels(rows,labels,bound+1)
                    if len(selected)==bound+1:
                        row_lookup={label:row for label,row in zip(labels,rows)}
                        square=[row_lookup[label] for label in selected]
                        return bound,selected,determinant_mod(square)
                raise RuntimeError("univariate Macaulay rank did not fill")

            components_mod=[tuple(mconvert(poly) for poly in component)
                            for component in components]
            exact_by_pivot={record["pivot"]:record for record in exact_pivot_results}
            pivot_deltas=[]
            for pivot_pair in combinations(range(4),2):
                other_rows=tuple(index for index in range(4) if index not in pivot_pair)
                first,second=pivot_pair
                delta_pair=madd(
                    mmul(components_mod[first][0],components_mod[second][1]),
                    mneg(mmul(components_mod[second][0],components_mod[first][1])),
                )
                pivot_deltas.append(delta_pair)
                compat_pair=[]
                for index in other_rows:
                    compat_pair.append(mdet([
                        list(components_mod[first]),list(components_mod[second]),
                        list(components_mod[index]),
                    ]))
                h_pair=mresultantD(compat_pair[0],compat_pair[1])
                physical_pair=msaturate(h_pair,delta_pair)
                nw_pair=madd(
                    mmul(components_mod[first][1],components_mod[second][2]),
                    mneg(mmul(components_mod[second][1],components_mod[first][2])),
                )
                ns_pair=madd(
                    mmul(components_mod[second][0],components_mod[first][2]),
                    mneg(mmul(components_mod[first][0],components_mod[second][2])),
                )
                aggregate_pair=dict(physical_pair)
                aggregate_full=dict(h_pair)
                row_degrees=[]
                projection_pairs_mod=[]
                for equation in order_equations[D+2]:
                    numerator_pair={}
                    for (ic,id0,ie,iw),coefficient in equation.items():
                        require((ic+ie)%2==0,"mod N12 odd C power")
                        is_power=(ic+ie)//2
                        denominator_power=iw+is_power
                        coefficient_mod=elem_mod(fast.from_sympy(coefficient),mod,residue_root)
                        term={(ie,id0):coefficient_mod}
                        term=mmul(term,mpow(nw_pair,iw))
                        term=mmul(term,mpow(ns_pair,is_power))
                        term=mmul(term,mpow(delta_pair,2-denominator_power))
                        numerator_pair=madd(numerator_pair,term)
                    left_pair=mresultantD(compat_pair[0],numerator_pair)
                    right_pair=mresultantD(compat_pair[1],numerator_pair)
                    projection_pairs_mod.append((left_pair,right_pair))
                    row_common=mugcd(left_pair,right_pair)
                    row_degrees.append(mudegree(row_common))
                    aggregate_pair=mugcd(aggregate_pair,row_common)
                    aggregate_full=mugcd(aggregate_full,row_common)
                print(
                    f"  pivot_pair={pivot_pair} delta_degree={mudegree(delta_pair)} "
                    f"H_degree={mudegree(h_pair)} saturated_H_degree={mudegree(physical_pair)} "
                    f"N12_row_gcd_degrees={row_degrees} "
                    f"full_aggregate_degree={mudegree(aggregate_full)} "
                    f"saturated_aggregate_degree={mudegree(aggregate_pair)}",
                    flush=True,
                )
                exact_record=exact_by_pivot[pivot_pair]
                exact_physical_mod=exact_univariate_mod(exact_record["physical_H"])
                exact_h_mod=exact_univariate_mod(exact_record["H"])
                require(
                    mmonic(exact_h_mod)==mmonic(h_pair),
                    f"exact/modular unsaturated H mismatch at pivot {pivot_pair}",
                )
                # Reduction can create an extra gcd with Delta.  The proof row
                # is therefore the direct reduction of the exact saturated H,
                # not a saturation performed after reduction.
                if mmonic(exact_physical_mod)!=mmonic(physical_pair):
                    print(
                        f"  residue_extra_delta_gcd pivot_pair={pivot_pair} "
                        f"exact_Hsat_mod_degree={mudegree(exact_physical_mod)} "
                        f"post_reduction_Hsat_degree={mudegree(physical_pair)}",
                        flush=True,
                    )
                for exact_pair,mod_pair in zip(
                    exact_record["projection_pairs"],projection_pairs_mod
                ):
                    for exact_projection,mod_projection in zip(exact_pair,mod_pair):
                        require(
                            mmonic(exact_univariate_mod(exact_projection))
                            ==mmonic(mod_projection),
                            f"exact/modular N12 projection mismatch at {pivot_pair}",
                        )
                bezout_polys=[exact_physical_mod]
                bezout_labels=["Hsat"]
                for row_index,(left_projection,right_projection) in enumerate(projection_pairs_mod):
                    bezout_polys.extend((left_projection,right_projection))
                    bezout_labels.extend((f"N12_{row_index}_L",f"N12_{row_index}_R"))
                bezout_bound,bezout_rows,bezout_determinant=convolution_rank_certificate(
                    bezout_polys,bezout_labels
                )
                require(bezout_bound==21,
                        f"lifted unit-ideal certificate bound changed at {pivot_pair}")
                print(
                    f"  EXACT_BEZOUT pivot_pair={pivot_pair} bound={bezout_bound} "
                    f"selected_rows={bezout_rows} det_mod_{mod}={bezout_determinant}",
                    flush=True,
                )
            all_delta_gcd=pivot_deltas[0]
            for delta_pair in pivot_deltas[1:]:
                all_delta_gcd=mugcd(all_delta_gcd,delta_pair)
            print(f"  all_pivot_deltas_gcd_degree={mudegree(all_delta_gcd)}",flush=True)
            exact_delta_mods=[exact_univariate_mod(poly)
                              for poly in exact_delta_polys]
            for exact_delta_mod,delta_mod_control in zip(exact_delta_mods,pivot_deltas):
                require(mmonic(exact_delta_mod)==mmonic(delta_mod_control),
                        "exact/modular pivot delta mismatch")
            delta_bound,delta_rows,delta_determinant=convolution_rank_certificate(
                exact_delta_mods,[f"Delta_{pair}" for pair in combinations(range(4),2)]
            )
            require(delta_bound==2,"lifted Delta-cover certificate bound changed")
            print(
                f"  EXACT_BEZOUT all_pivot_deltas bound={delta_bound} "
                f"selected_rows={delta_rows} det_mod_{mod}={delta_determinant}",flush=True,
            )

print("DONE")
