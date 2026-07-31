#!/usr/bin/env python3
"""Exact conditional second-edge circuit certificate for THM-2997.

The three reduced resultant log-jet numerators below are a frozen exact output
of the classical 36-by-36 Macaulay chart minus its 10-by-10 extraneous minor.
An independent raw-chart rebuild is recorded in the theorem audit.  This
companion performs every wall, residue, rational-box, finite-prefix, and
redundant cleared-circuit check from those exact coefficient tables.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx


M, U, V = sp.symbols("M U V")
X, Y, T = sp.symbols("X Y T")
CTX = fmpq_mpoly_ctx.get(["X", "Y", "T"])


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


# Each row is (M exponent, U exponent, V exponent, numerator, denominator).
# Canonical compact JSON digest: cfb36557e1d54a0328a309375a948ace99c78e0688a54a014aef0906c1b90513.
JET_DATA_TEXT = r'''[[[2,4,0,23,1],[2,3,0,93,1],[2,2,1,-87,1],[2,2,0,170,1],[2,1,1,-288,1],[2,1,0,-120,1],[2,0,2,207,1],[2,0,1,114,1],[2,0,0,23,1],[1,4,0,-12,1],[1,3,0,225,1],[1,2,1,-327,1],[1,2,0,-174,1],[1,1,1,288,1],[1,1,0,144,1],[1,0,2,-99,1],[1,0,1,-42,1],[1,0,0,-3,1],[0,4,0,66,1],[0,3,0,-900,1],[0,2,1,576,1],[0,2,0,1794,1],[0,1,1,-1152,1],[0,1,0,-1152,1],[0,0,2,90,1],[0,0,1,396,1],[0,0,0,282,1]],[[4,7,0,225,1],[4,6,1,-289,1],[4,6,0,-1647,1],[4,5,1,36,1],[4,5,0,2133,1],[4,4,2,-2244,1],[4,4,1,-2643,1],[4,4,0,-765,1],[4,3,2,1035,1],[4,3,1,-4122,1],[4,3,0,-81,1],[4,2,3,-2601,1],[4,2,2,-12441,1],[4,2,1,1437,1],[4,2,0,-27,1],[4,1,3,-5292,1],[4,1,2,180,1],[4,1,1,-612,1],[4,1,0,-36,1],[4,0,3,576,1],[4,0,2,-192,1],[4,0,1,64,1],[3,8,0,-368,3],[3,7,0,-4046,1],[3,6,1,16714,3],[3,6,0,52142,3],[3,5,1,-21408,1],[3,5,0,-40530,1],[3,4,2,30672,1],[3,4,1,99230,1],[3,4,0,21242,1],[3,3,2,-72414,1],[3,3,1,-17628,1],[3,3,0,1146,1],[3,2,3,49794,1],[3,2,2,-27630,1],[3,2,1,-37802,1],[3,2,0,-14206,3],[3,1,3,22536,1],[3,1,2,-1944,1],[3,1,1,9048,1],[3,1,0,824,1],[3,0,4,-9936,1],[3,0,3,-12480,1],[3,0,2,-5472,1],[3,0,1,-3520,3],[3,0,0,-368,3],[2,8,0,96,1],[2,7,0,17055,1],[2,6,1,-26639,1],[2,6,0,-225123,1],[2,5,1,298104,1],[2,5,0,352719,1],[2,4,2,-220014,1],[2,4,1,-689217,1],[2,4,0,-208197,1],[2,3,2,765549,1],[2,3,1,404106,1],[2,3,0,70281,1],[2,2,3,-198567,1],[2,2,2,-447039,1],[2,2,1,-1917,1],[2,2,0,-11517,1],[2,1,3,145260,1],[2,1,2,123372,1],[2,1,1,-60636,1],[2,1,0,-3708,1],[2,0,4,7128,1],[2,0,3,-48672,1],[2,0,2,-47664,1],[2,0,1,5216,1],[2,0,0,24,1],[1,8,0,-94,3],[1,7,0,-30634,1],[1,6,1,297710,3],[1,6,0,2318848,3],[1,5,1,-1264092,1],[1,5,0,-2509746,1],[1,4,2,496986,1],[1,4,1,4726486,1],[1,4,0,2189080,1],[1,3,2,-3041514,1],[1,3,1,-4672788,1],[1,3,0,-660162,1],[1,2,3,520398,1],[1,2,2,3236550,1],[1,2,1,2097098,1],[1,2,0,281446,3],[1,1,3,-583848,1],[1,1,2,-1436328,1],[1,1,1,-527928,1],[1,1,0,-23672,1],[1,0,4,-2664,1],[1,0,3,239904,1],[1,0,2,230928,1],[1,0,1,150112,3],[1,0,0,-184,3],[0,8,0,-11172,1],[0,7,0,105552,1],[0,6,1,-90432,1],[0,6,0,-1332888,1],[0,5,1,1611648,1],[0,5,0,5774736,1],[0,4,2,-443376,1],[0,4,1,-7810848,1],[0,4,0,-10621332,1],[0,3,2,3401568,1],[0,3,1,14634432,1],[0,3,0,9000096,1],[0,2,3,-443136,1],[0,2,2,-7456752,1],[0,2,1,-10644192,1],[0,2,0,-3728688,1],[0,1,3,1537920,1],[0,1,2,4527360,1],[0,1,1,3067776,1],[0,1,0,764928,1],[0,0,4,-119448,1],[0,0,3,-616992,1],[0,0,2,-745488,1],[0,0,1,-290592,1],[0,0,0,-70680,1]],[[6,11,0,375,1],[6,10,1,-4913,9],[6,10,0,-12429,1],[6,9,1,28616,3],[6,9,0,75606,1],[6,8,2,-61258,3],[6,8,1,-785750,9],[6,8,0,-128550,1],[6,7,2,70454,1],[6,7,1,-4592,3],[6,7,0,70737,1],[6,6,3,-127650,1],[6,6,2,-2567558,3],[6,6,1,291961,9],[6,6,0,-26625,1],[6,5,3,-89532,1],[6,5,2,253758,1],[6,5,1,-77504,1],[6,5,0,7950,1],[6,4,4,-183774,1],[6,4,3,-583542,1],[6,4,2,-768482,1],[6,4,1,127646,9],[6,4,0,-3000,1],[6,3,4,-171873,1],[6,3,3,-1411452,1],[6,3,2,190122,1],[6,3,1,-55892,3],[6,3,0,-177,1],[6,2,5,-44217,1],[6,2,4,-736023,1],[6,2,3,221958,1],[6,2,2,-270506,3],[6,2,1,22811,9],[6,2,0,-171,1],[6,1,5,-74088,1],[6,1,4,62280,1],[6,1,3,-34128,1],[6,1,2,5200,1],[6,1,1,-1624,3],[6,1,0,-24,1],[6,0,5,4608,1],[6,0,4,-7680,1],[6,0,3,3072,1],[6,0,2,-2560,3],[6,0,1,512,9],[5,11,0,-5445,1],[5,10,1,93755,9],[5,10,0,244071,1],[5,9,1,-1246928,3],[5,9,0,-1575936,1],[5,8,2,1415198,3],[5,8,1,34165484,9],[5,8,0,2891988,1],[5,7,2,-3057304,1],[5,7,1,-23059924,3],[5,7,0,-2247831,1],[5,6,3,3036744,1],[5,6,2,41278888,3],[5,6,1,38183321,9],[5,6,0,1056411,1],[5,5,3,-4122468,1],[5,5,2,-9625362,1],[5,5,1,-1165288,1],[5,5,0,-349506,1],[5,4,4,4234134,1],[5,4,3,7933626,1],[5,4,2,3267094,1],[5,4,1,302398,9],[5,4,0,80460,1],[5,3,4,272745,1],[5,3,3,-3953988,1],[5,3,2,-4175850,1],[5,3,1,-448396,3],[5,3,0,-20151,1],[5,2,5,841653,1],[5,2,4,-1406877,1],[5,2,3,-4456350,1],[5,2,2,1207522,3],[5,2,1,-362591,9],[5,2,0,-585,1],[5,1,5,-563976,1],[5,1,4,-1731960,1],[5,1,3,311472,1],[5,1,2,-191984,1],[5,1,1,-19832,3],[5,1,0,-792,1],[5,0,5,32256,1],[5,0,4,30720,1],[5,0,2,-1024,3],[5,0,1,5632,9],[4,12,0,-1472,3],[4,11,0,27219,1],[4,10,1,-263147,3],[4,10,0,-3042605,1],[4,9,1,4850424,1],[4,9,0,22829718,1],[4,8,2,-4769686,1],[4,8,1,-158064634,3],[4,8,0,-44671558,1],[4,7,2,47802102,1],[4,7,1,123686000,1],[4,7,0,36142797,1],[4,6,3,-29952474,1],[4,6,2,-168577474,1],[4,6,1,-104980551,1],[4,6,0,-53001251,3],[4,5,3,101167380,1],[4,5,2,162545814,1],[4,5,1,58156176,1],[4,5,0,6795510,1],[4,4,4,-40357062,1],[4,4,3,-59236038,1],[4,4,2,-93696066,1],[4,4,1,-46141414,3],[4,4,0,-1738576,1],[4,3,4,49304187,1],[4,3,3,21873780,1],[4,3,2,42869970,1],[4,3,1,1554900,1],[4,3,0,312171,1],[4,2,5,-6030153,1],[4,2,4,-10116855,1],[4,2,3,-58554,1],[4,2,2,-19773038,1],[4,2,1,-2771975,3],[4,2,0,-102635,1],[4,1,5,7014600,1],[4,1,4,-1969704,1],[4,1,3,-9426672,1],[4,1,2,1550064,1],[4,1,1,-22616,1],[4,1,0,1464,1],[4,0,6,-357696,1],[4,0,5,-1876608,1],[4,0,4,-1374912,1],[4,0,3,699648,1],[4,0,2,-300224,1],[4,0,1,-3712,1],[4,0,0,-1472,3],[3,12,0,512,1],[3,11,0,-51867,1],[3,10,1,5792461,9],[3,10,0,15596421,1],[3,9,1,-28425560,1],[3,9,0,-180723864,1],[3,8,2,26638018,1],[3,8,1,3624690268,9],[3,8,0,446825604,1],[3,7,2,-352610312,1],[3,7,1,-1188743332,1],[3,7,0,-448929417,1],[3,6,3,157811832,1],[3,6,2,1364633144,1],[3,6,1,11218306543,9],[3,6,0,253013453,1],[3,5,3,-890897004,1],[3,5,2,-1564814766,1],[3,5,1,-773850888,1],[3,5,0,-92387262,1],[3,4,4,192907230,1],[3,4,3,1158029238,1],[3,4,2,871166850,1],[3,4,1,3116080370,9],[3,4,0,25556232,1],[3,3,4,-446830425,1],[3,3,3,-555644412,1],[3,3,2,-258715446,1],[3,3,1,-87793852,1],[3,3,0,-4235193,1],[3,2,5,34039683,1],[3,2,4,276250245,1],[3,2,3,78367566,1],[3,2,2,51272666,1],[3,2,1,93339575,9],[3,2,0,3345,1],[3,1,5,-25866648,1],[3,1,4,-44072424,1],[3,1,3,-6198768,1],[3,1,2,-16794832,1],[3,1,1,-2553080,1],[3,1,0,-65736,1],[3,0,6,342144,1],[3,0,5,6863616,1],[3,0,4,5516160,1],[3,0,3,-1024512,1],[3,0,2,-230528,1],[3,0,1,1425152,9],[3,0,0,128,1],[2,12,0,-752,3],[2,11,0,512022,1],[2,10,1,-29864254,9],[2,10,0,-53660786,1],[2,9,1,305567896,3],[2,9,0,723097884,1],[2,8,2,-276254576,3],[2,8,1,-13995936676,9],[2,8,0,-2709868228,1],[2,7,2,1415091604,1],[2,7,1,19976673056,3],[2,7,0,3658229634,1],[2,6,3,-435372420,1],[2,6,2,-21515671060,3],[2,6,1,-85760208010,9],[2,6,0,-7662880970,3],[2,5,3,3514009464,1],[2,5,2,11890010220,1],[2,5,1,6876725984,1],[2,5,0,1056891516,1],[2,4,4,-527467584,1],[2,4,3,-7275455316,1],[2,4,2,-9431297716,1],[2,4,1,-24667077308,9],[2,4,0,-265063708,1],[2,3,4,1792917990,1],[2,3,3,5951812680,1],[2,3,2,3845086788,1],[2,3,1,1753935896,3],[2,3,0,43521414,1],[2,2,5,-119821950,1],[2,2,4,-1520283090,1],[2,2,3,-2160702540,1],[2,2,2,-2793969868,3],[2,2,1,-333598598,9],[2,2,0,-4959818,1],[2,1,5,73494432,1],[2,1,4,526118112,1],[2,1,3,419229504,1],[2,1,2,172025024,1],[2,1,1,-32209568,3],[2,1,0,-313248,1],[2,0,6,-191808,1],[2,0,5,-29621376,1],[2,0,4,-71479488,1],[2,0,3,-49420032,1],[2,0,2,-73173056,3],[2,0,1,6481792,9],[2,0,0,-1472,3],[1,11,0,596448,1],[1,10,1,12782472,1],[1,10,0,119827224,1],[1,9,1,-899815840,3],[1,9,0,-1561807968,1],[1,8,2,485587456,3],[1,8,1,3554787696,1],[1,8,0,7885330128,1],[1,7,2,-2698602528,1],[1,7,1,-57510031904,3],[1,7,0,-16954580256,1],[1,6,3,610179360,1],[1,6,2,52075706624,3],[1,6,1,43224852648,1],[1,6,0,17559960408,1],[1,5,3,-6454584576,1],[1,5,2,-43048711008,1],[1,5,1,-45712431872,1],[1,5,0,-9726094752,1],[1,4,4,761102256,1],[1,4,3,19710550752,1],[1,4,2,45992593760,1],[1,4,1,27233453856,1],[1,4,0,2916344304,1],[1,3,4,-3886406928,1],[1,3,3,-20971704960,1],[1,3,2,-27162241056,1],[1,3,1,-30573877760,3],[1,3,0,-421530192,1],[1,2,5,271124424,1],[1,2,4,3928431864,1],[1,2,3,11576018256,1],[1,2,2,28376674640,3],[1,2,1,2530051176,1],[1,2,0,19792152,1],[1,1,5,-180877536,1],[1,1,4,-1970106912,1],[1,1,3,-3168771264,1],[1,1,2,-1854878784,1],[1,1,1,-1175852320,3],[1,1,0,-1370016,1],[1,0,5,73285632,1],[1,0,4,349317120,1],[1,0,3,333637632,1],[1,0,2,471574528,3],[1,0,1,27359232,1],[0,12,0,-1123568,1],[0,11,0,9483552,1],[0,10,1,-10809216,1],[0,10,0,-226998672,1],[0,9,1,360112896,1],[0,9,0,2024434752,1],[0,8,2,-128782368,1],[0,8,1,-3471050304,1],[0,8,0,-10731247536,1],[0,7,2,1869717888,1],[0,7,1,20833668864,1],[0,7,0,31005021600,1],[0,6,3,-299109888,1],[0,6,2,-15057157824,1],[0,6,1,-62945770752,1],[0,6,0,-50641475120,1],[0,5,3,4684545792,1],[0,5,2,50309179776,1],[0,5,1,100232522496,1],[0,5,0,49081798656,1],[0,4,4,-517353408,1],[0,4,3,-18944008704,1],[0,4,2,-80036800416,1],[0,4,1,-88346214720,1],[0,4,0,-29248583520,1],[0,3,4,3177523584,1],[0,3,3,31451348736,1],[0,3,2,63705864960,1],[0,3,1,44703680256,1],[0,3,0,10786772352,1],[0,2,5,-187414272,1],[0,2,4,-6259559616,1],[0,2,3,-21850223616,1],[0,2,2,-26524333440,1],[0,2,1,-12832051968,1],[0,2,0,-2388115392,1],[0,1,5,708113664,1],[0,1,4,3389536512,1],[0,1,3,6827217408,1],[0,1,2,5484704256,1],[0,1,1,1916156160,1],[0,1,0,290886912,1],[0,0,6,-54793152,1],[0,0,5,-191940480,1],[0,0,4,-612044352,1],[0,0,3,-773906688,1],[0,0,2,-444290112,1],[0,0,1,-111484800,1],[0,0,0,-15647168,1]]]'''
EXPECTED_JET_DIGEST = "cfb36557e1d54a0328a309375a948ace99c78e0688a54a014aef0906c1b90513"


def load_jets():
    data = json.loads(JET_DATA_TEXT)
    canonical = json.dumps(data, separators=(",", ":"))
    require(hashlib.sha256(canonical.encode("ascii")).hexdigest() == EXPECTED_JET_DIGEST, "jet table digest")
    expressions = []
    for rows in data:
        terms = {
            (m_power, u_power, v_power): sp.Rational(numerator, denominator)
            for m_power, u_power, v_power, numerator, denominator in rows
        }
        expressions.append(sp.Poly.from_dict(terms, (M, U, V), domain=sp.QQ).as_expr())
    require([len(sp.Poly(value, M, U, V).terms()) for value in expressions] == [27, 122, 333], "jet term counts")
    return tuple(expressions)


def sum_power(upper, power):
    if power == 0:
        return upper
    if power == 1:
        return upper * (upper + 1) * sp.Rational(1, 2)
    if power == 2:
        return upper * (upper + 1) * (2 * upper + 1) * sp.Rational(1, 6)
    if power == 3:
        return (upper * (upper + 1) * sp.Rational(1, 2)) ** 2
    raise ValueError(power)


def interval_power(low, high, power):
    if high == low - 1:
        return sp.Integer(0)
    return sum_power(high, power) - sum_power(low - 1, power)


def wall_stat(residue, power):
    local_t = sp.symbols("T", integer=True, nonnegative=True)
    width = 30 * local_t + residue
    a = 10 * local_t + residue // 3
    b = 15 * local_t + residue // 2
    c = 20 * local_t + (2 * residue) // 3
    value = 6 * 1**power + 21 * 2**power
    value += 26 * interval_power(3, a, power)
    value += 24 * interval_power(a + 1, b, power)
    value += 20 * interval_power(b + 1, c, power)
    value += 19 * interval_power(c + 1, width, power)
    value += 3 * interval_power(2, b, power)
    value += 4 * interval_power(b + 1, width - 1, power)
    value += 2 * width**power
    value -= interval_power(3, b, power)
    value -= width**power
    if residue % 10 == 1:
        value += ((4 * width + 1) * sp.Rational(1, 5)) ** power
    value += 6 if power == 0 else sp.Rational(6, 2**power)
    return local_t, sp.expand(value)


def to_flint(expression):
    poly = sp.Poly(expression, X, Y, T, domain=sp.QQ)
    return CTX.from_dict(
        {
            tuple(map(int, exponent)): fmpq(int(sp.numer(coefficient)), int(sp.denom(coefficient)))
            for exponent, coefficient in poly.terms()
        }
    )


def as_fraction(value):
    return Fraction(int(value.numer()), int(value.denom()))


def coefficient_list(poly):
    degree = max(poly)
    return [poly.get(power, Fraction(0)) for power in range(degree, -1, -1)]


def residue_polynomials(jets, residue):
    P1, P2, P3 = jets
    D = U**2 + 3 * U - 3 * V - 1
    local_t, wall_degree = wall_stat(residue, 0)
    _, wall1 = wall_stat(residue, 1)
    _, wall2 = wall_stat(residue, 2)
    _, wall3 = wall_stat(residue, 3)
    replacements = {local_t: T}
    wall_degree, wall1, wall2, wall3 = [sp.expand(value.subs(replacements)) for value in (wall_degree, wall1, wall2, wall3)]
    width = 30 * T + residue
    substitution = {M: width, U: 2**residue * X, V: 3**residue * Y}
    Ds = D.subs(substitution)
    p1 = to_flint(P1.subs(substitution))
    p2 = to_flint(P2.subs(substitution))
    p3 = to_flint(P3.subs(substitution))
    df = to_flint(46 * width - 26 - wall_degree)
    Df = to_flint(Ds)
    w1, w2, w3 = to_flint(wall1), to_flint(wall2), to_flint(wall3)
    A = p1 - w1 * Df**2
    B = p2 / 8 + w2 * Df**4
    C = -3 * p3 / 128 - w3 * Df**6
    return {
        "D": Df,
        "P1": p1,
        "P2": p2,
        "P3": p3,
        "A": A,
        "B": B,
        "C": C,
        "d": df,
        "w1": w1,
        "w2": w2,
        "w3": w3,
        "degree_text": str(sp.expand(46 * width - 26 - wall_degree)),
    }


def certify_positive(poly, minimum_t, name):
    grouped = {}
    for (xp0, yp0, tp0), coefficient in poly.to_dict().items():
        xp, yp, tp = int(xp0), int(yp0), int(tp0)
        grouped.setdefault((xp, yp), {})[tp] = as_fraction(coefficient)
    bases = sorted(grouped, key=lambda pair: 2**pair[0] * 3**pair[1], reverse=True)
    require(len(bases) >= 2, (name, "needs two exponential bases"))
    dominant, runner = bases[:2]
    dominant_base = 2**dominant[0] * 3**dominant[1]
    runner_base = 2**runner[0] * 3**runner[1]
    dominant_coeffs = grouped[dominant]
    dominant_degree = max(dominant_coeffs)
    leading = dominant_coeffs[dominant_degree]
    require(leading > 0, (name, "dominant leading coefficient", leading))
    negative_terms = tuple(
        (power, -coefficient)
        for power, coefficient in dominant_coeffs.items()
        if power < dominant_degree and coefficient < 0
    )
    lower_polys = [grouped[base] for base in bases[1:]]
    lower_degree = max(max(one) for one in lower_polys)
    degree_delta = max(0, lower_degree - dominant_degree)
    absolute_tariff = sum(
        (sum((abs(coefficient) for coefficient in one.values()), Fraction(0)) for one in lower_polys),
        Fraction(0),
    )
    ratio = Fraction(runner_base, dominant_base)
    tail_start = minimum_t
    while True:
        floor = leading - sum(
            (coefficient / tail_start ** (dominant_degree - power) for power, coefficient in negative_terms),
            Fraction(0),
        )
        margin = floor - absolute_tariff * tail_start**degree_delta * ratio ** (30 * tail_start)
        if floor > 0 and margin > 0:
            break
        tail_start += 1
    monotone_step = 2**degree_delta * ratio**30
    require(monotone_step < 1, (name, "tail majorant not monotone", monotone_step))
    prefix = []
    for value in range(minimum_t, tail_start):
        exact = as_fraction(poly(2 ** (30 * value), 3 ** (30 * value), value))
        require(exact > 0, (name, value, "finite prefix not positive", exact))
        prefix.append((value, exact))
    grouped_text = "\n".join(
        f"{xp},{yp}:" + ",".join(str(value) for value in coefficient_list(grouped[xp, yp]))
        for xp, yp in bases
    )
    prefix_text = "\n".join(f"{value}:{exact}" for value, exact in prefix)
    return {
        "name": name,
        "terms": len(poly),
        "bases": len(grouped),
        "grouped_sha256": hashlib.sha256(grouped_text.encode("ascii")).hexdigest(),
        "dominant": dominant,
        "runner": runner,
        "base_ratio": str(ratio),
        "dominant_degree": dominant_degree,
        "lower_degree": lower_degree,
        "tail_start": tail_start,
        "tail_margin": str(margin),
        "monotone_step": str(monotone_step),
        "finite_prefix_count": len(prefix),
        "finite_prefix_sha256": hashlib.sha256(prefix_text.encode("ascii")).hexdigest(),
    }


def circuit_polynomial(row):
    A, B, C, d = row["A"], row["B"], row["C"], row["d"]
    return d * (d - 2) * (A**2 + B) ** 3 - A**3 * (d - 1) ** 2 * (A**3 + 3 * A * B + 2 * C)


def box_record(row, residue):
    A, B, C, d = row["A"], row["B"], row["C"], row["d"]
    minimum_t = max(1, (43 - residue + 29) // 30)
    polynomials = (
        ("A_positive", A),
        ("x_ge_129_over_100", -100 * d * B - 129 * A**2),
        ("x_le_2", 2 * A**2 + d * B),
        ("z_nonnegative", C),
        ("z_le_39_over_20", 39 * A**3 - 20 * d**2 * C),
    )
    return {
        "residue": residue,
        "minimum_T": minimum_t,
        "minimum_M": 30 * minimum_t + residue,
        "certificates": [certify_positive(poly, minimum_t, name) for name, poly in polynomials],
    }


def h_record(row, residue):
    minimum_t = 2 if residue < 4 else 1
    return {
        "residue": residue,
        "degree": row["degree_text"],
        "certificate": certify_positive(circuit_polynomial(row), minimum_t, "cleared_H"),
    }


def unreduced_r4_record(row):
    D, P1, P2, P3, d = row["D"], row["P1"], row["P2"], row["P3"], row["d"]
    delta = -2 * D
    N1 = 16 * D**2 * P1
    N2 = 16 * D**4 * P2
    N3 = -32 * D**6 * P3
    Au = N1 - row["w1"] * delta**4
    Bu = 2 * N2 + row["w2"] * delta**8
    Cu = 3 * N3 - row["w3"] * delta**12
    unreduced = d * (d - 2) * (Au**2 + Bu) ** 3 - Au**3 * (d - 1) ** 2 * (Au**3 + 3 * Au * Bu + 2 * Cu)
    return certify_positive(unreduced, 1, "unreduced_r4_H")


def direct_boundary_records(rows):
    records = []
    for width in range(33, 43):
        residue = width % 30
        value = (width - residue) // 30
        row = rows[residue]
        xv, yv = 2 ** (30 * value), 3 ** (30 * value)
        A = as_fraction(row["A"](xv, yv, value))
        B = as_fraction(row["B"](xv, yv, value))
        C = as_fraction(row["C"](xv, yv, value))
        D = as_fraction(row["D"](xv, yv, value))
        d = as_fraction(row["d"](xv, yv, value))
        H = as_fraction(circuit_polynomial(row)(xv, yv, value))
        require(A > 0 and D != 0 and d > 2, (width, "direct normalization"))
        xnorm = -d * B / A**2
        znorm = d**2 * C / A**3
        quadratic = (
            d**2 * (3 * xnorm**2 - 2 * znorm - 1)
            + d * (-xnorm**3 - 6 * xnorm**2 + 3 * xnorm + 4 * znorm)
            + 2 * xnorm**3
            - 2 * znorm
        )
        require(quadratic == d**2 * H / A**6, (width, "quadratic/H identity"))
        require((quadratic < 0) if width == 33 else (quadratic > 0), (width, "sharp boundary sign"))
        u = A / D**2
        ell2 = B / D**4
        ell3 = C / D**6
        h1 = u / d
        h2 = (u**2 + ell2) / (d * (d - 1))
        h3 = (u**3 + 3 * u * ell2 + 2 * ell3) / (d * (d - 1) * (d - 2))
        R1 = h1**2 / h2
        R2 = h2**2 / (h1 * h3)
        require((R2 < R1) if width == 33 else (R2 > R1), (width, "ratio boundary sign"))
        records.append(
            {
                "M": width,
                "degree": int(d),
                "x": str(xnorm),
                "z": str(znorm),
                "quadratic": str(quadratic),
                "R2_minus_R1": str(R2 - R1),
            }
        )
    text = "\n".join(json.dumps(record, sort_keys=True, separators=(",", ":")) for record in records)
    return {
        "M33_sign": "NEGATIVE",
        "M34_through_M42_sign": "STRICTLY_POSITIVE",
        "M34_degree": records[1]["degree"],
        "M34_R2_minus_R1_float": format(float(Fraction(records[1]["R2_minus_R1"])), ".17g"),
        "records_sha256": hashlib.sha256(text.encode("ascii")).hexdigest(),
    }


def asymptotic_record(jets):
    P1, P2, P3 = jets
    D = U**2 + 3 * U - 3 * V - 1
    limits = (
        sp.cancel(sp.LC(sp.Poly(P1.subs(V, 0), U)) / sp.LC(sp.Poly(D**2, U))),
        sp.cancel(sp.LC(sp.Poly(P2.subs(V, 0), U)) / (16 * sp.LC(sp.Poly(D**4, U)))),
        sp.cancel(-sp.LC(sp.Poly(P3.subs(V, 0), U)) / (128 * sp.LC(sp.Poly(D**6, U)))),
    )
    expected = (
        23 * M**2 - 12 * M + 66,
        -(184 * M**3 - 144 * M**2 + 47 * M + 16758) / 24,
        (92 * M**4 - 96 * M**3 + 47 * M**2 + 210669) / 24,
    )
    require(all(sp.expand(left - right) == 0 for left, right in zip(limits, expected)), "dominant-U log jets")
    wall_limits = []
    for power in range(4):
        one = []
        for residue in range(30):
            local_t, expression = wall_stat(residue, power)
            degree = power + 1 if power > 0 else 1
            one.append(sp.cancel(sp.Poly(expression, local_t).LC() / 30**degree))
        require(len(set(one)) == 1, (power, "wall leading residue mismatch"))
        wall_limits.append(one[0])
    require(wall_limits == [sp.Rational(76, 3), sp.Rational(145, 12), sp.Rational(2551, 324), sp.Rational(1681, 288)], "wall leading moments")
    dlead = 46 - wall_limits[0]
    ulead = 23 - wall_limits[1]
    ell2lead = 2 * sp.LC(sp.Poly(expected[1], M)) + wall_limits[2]
    ell3lead = 3 * sp.LC(sp.Poly(expected[2], M)) - wall_limits[3]
    require((dlead, ulead, ell2lead, ell3lead) == (sp.Rational(62, 3), sp.Rational(131, 12), -sp.Rational(2417, 324), sp.Rational(1631, 288)), "core leading jets")
    xlimit = sp.cancel(-dlead * ell2lead / ulead**2)
    zlimit = sp.cancel(dlead**2 * ell3lead / ulead**3)
    curvature = sp.cancel(3 * xlimit**2 - 2 * zlimit - 1)
    require(xlimit == sp.Rational(599416, 463347), "x limit")
    require(zlimit == sp.Rational(12539128, 6744273), "z limit")
    require(curvature == sp.Rational(21630685837, 71563480803) > 0, "curvature limit")
    return {
        "resultant_L1": str(expected[0]),
        "resultant_L2": str(expected[1]),
        "resultant_L3": str(expected[2]),
        "d_over_M": str(dlead),
        "u_over_M2": str(ulead),
        "ell2_over_M3": str(ell2lead),
        "ell3_over_M4": str(ell3lead),
        "x_limit": str(xlimit),
        "z_limit": str(zlimit),
        "curvature_limit": str(curvature),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    jets = load_jets()
    asymptotic = asymptotic_record(jets)
    rows = [residue_polynomials(jets, residue) for residue in range(30)]
    degree_minima = []
    for residue, row in enumerate(rows):
        minimum_t = 2 if residue < 4 else 1
        degree_here = as_fraction(row["d"](1, 1, minimum_t))
        degree_next = as_fraction(row["d"](1, 1, minimum_t + 1))
        require(degree_next - degree_here == 620, (residue, "degree slope"))
        degree_minima.append(degree_here)
    require(min(degree_minima) == 701, "minimum degree in M>=34 scope")
    require(2 * Fraction(923, 10000) * 701 - 26 > 0, "quadratic floor increases after d=701")
    box = [box_record(rows[residue], residue) for residue in range(30)]
    h_records = [h_record(rows[residue], residue) for residue in range(30)]
    unreduced = unreduced_r4_record(rows[4])
    direct = direct_boundary_records(rows)

    box_lines = [json.dumps(record, sort_keys=True, separators=(",", ":")) for record in box]
    h_lines = [json.dumps(record, sort_keys=True, separators=(",", ":")) for record in h_records]
    box_prefixes = sum(certificate["finite_prefix_count"] for record in box for certificate in record["certificates"])
    h_prefixes = sum(record["certificate"]["finite_prefix_count"] for record in h_records)
    require(box_prefixes == 17 and h_prefixes == 32, "prefix census")
    require(all(record["certificate"]["terms"] == 2529 and record["certificate"]["bases"] == 169 for record in h_records), "cleared H census")
    require(unreduced["terms"] == 9369 and unreduced["bases"] == 625, "unreduced r4 census")
    transcript = "\n".join(
        [
            "FIRST-GAP WALL-STRIPPED ALL-WIDTH SECOND-EDGE CIRCUIT CERTIFICATE",
            f"jet_table_sha256={EXPECTED_JET_DIGEST};terms=27,122,333",
            "response_D=U^2+3U-3V-1;log_jets=P1/D^2,P2/(16D^4),-P3/(128D^6)",
            "discrete_reframing=R_k>=R_(k-1)_iff_Delta^3(log_h)_(k-2)<=0",
            "normalized_circuit=F=d^2(3x^2-2z-1)+d(-x^3-6x^2+3x+4z)+2x^3-2z",
            "box=129/100<=x<=2;0<=z<=39/20;scope=M>=43",
            "box_implies=F>=923*d^2/10000-26*d-39/10",
            "minimum_degree_M34_plus=701;floor_at_701=271264123/10000",
            "asymptotic=" + json.dumps(asymptotic, sort_keys=True, separators=(",", ":")),
            "direct_boundary=" + json.dumps(direct, sort_keys=True, separators=(",", ":")),
            *["box_residue=" + line for line in box_lines],
            "box_record_digest=" + hashlib.sha256("\n".join(box_lines).encode("ascii")).hexdigest(),
            f"box_finite_prefix_count={box_prefixes}",
            *["cleared_H_residue=" + line for line in h_lines],
            "cleared_H_record_digest=" + hashlib.sha256("\n".join(h_lines).encode("ascii")).hexdigest(),
            f"cleared_H_finite_prefix_count={h_prefixes}",
            "unreduced_r4=" + json.dumps(unreduced, sort_keys=True, separators=(",", ":")),
            "encoded_continuation_scope=M_GE_34",
            "encoded_second_edge_sign=STRICTLY_POSITIVE",
            "actual_core_scope=M34_PROVED;M35_PLUS_REQUIRES_STATED_CONTINUATION",
            "full_no_return_full_ULC_arbitrary_radial_GMC2=NOT_ASSERTED",
            "all_exact_checks=PASS",
        ]
    ) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
