"""Exact finite full-record separator compiler; no assertions as checks.

Input allowed tables are the complete conjunction of factors assigned to a bag.
Records can be any hashable values. Charges use tuples with modulus zero for Z.
"""
from itertools import product


class Group:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        if any(m < 0 for m in self.moduli):
            raise ValueError('negative modulus')
        self.zero = (0,) * len(self.moduli)

    def add(self, a, b):
        return tuple((x+y) % m if m else x+y
                     for x, y, m in zip(a, b, self.moduli))


def merge(a, b):
    if any(k in a and a[k] != v for k, v in b.items()):
        raise ValueError('inconsistent actual separator records')
    return a | b


class TreeCarrier:
    def __init__(self, bags, parents, domains, charges, costs, group):
        self.bags = [tuple(b) for b in bags]
        self.parents = list(parents)
        self.domains, self.charges, self.costs, self.group = domains, charges, costs, group
        n = len(bags)
        if not n or len(parents) != n or parents[0] != -1:
            raise ValueError('one nonempty rooted bag tree required')
        self.children = [[] for _ in bags]
        for v in range(1, n):
            if not 0 <= parents[v] < v:
                raise ValueError('parents must precede children')
            self.children[parents[v]].append(v)
        if any(len(b) != len(set(b)) for b in self.bags):
            raise ValueError('duplicate slot in a bag')
        slots = set().union(*map(set, self.bags))
        if slots != set(domains):
            raise ValueError('domains must name exactly the covered slots')
        self.slots = sorted(slots)
        self.separators = [()]
        for v in range(1, n):
            self.separators.append(tuple(i for i in self.bags[v]
                                         if i in self.bags[parents[v]]))
        self.owned = [tuple(i for i in b if i not in self.separators[v])
                      for v, b in enumerate(self.bags)]
        for i in slots:
            holders = {v for v, b in enumerate(self.bags) if i in b}
            first = min(holders)
            if any(v != first and parents[v] not in holders for v in holders):
                raise ValueError('running intersection fails')
            if sum(i in o for o in self.owned) != 1:
                raise ValueError('slot charge/cost has no unique owner bag')
            for r in domains[i]:
                if len(charges[i][r]) != len(group.moduli):
                    raise ValueError('charge type mismatch')
            if any(costs[i][a, b] < 0 for a, b in product(domains[i], repeat=2)):
                raise ValueError('negative replacement cost')

    def local_charge(self, v, record):
        q = self.group.zero
        for i in self.owned[v]:
            q = self.group.add(q, self.charges[i][record[i]])
        return q

    def validate_tables(self, tables):
        if len(tables) != len(self.bags):
            raise ValueError('one complete local relation per bag required')
        for b, allowed in zip(self.bags, tables):
            if len(allowed) != len(set(allowed)):
                raise ValueError('local records cannot be duplicated as capacity')
            if any(len(t) != len(b) or any(r not in self.domains[i] for i, r in zip(b, t))
                   for t in allowed):
                raise ValueError('invalid local record tuple')

    def count(self, tables):
        """Return root charge -> (exact count, actual witness), plus all messages."""
        self.validate_tables(tables)
        messages = [None] * len(self.bags)
        for v in reversed(range(len(self.bags))):
            out = {}
            for t in tables[v]:
                record = dict(zip(self.bags[v], t))
                sep = tuple(record[i] for i in self.separators[v])
                partial = {self.local_charge(v, record): (1, record)}
                for w in self.children[v]:
                    child_sep = tuple(record[i] for i in self.separators[w])
                    nxt = {}
                    for q, (count, witness) in partial.items():
                        for r, (other_count, other_witness) in messages[w].get(child_sep, {}).items():
                            z = self.group.add(q, r)
                            total = count * other_count
                            if z in nxt:
                                old, old_w = nxt[z]
                                nxt[z] = old + total, old_w
                            else:
                                nxt[z] = total, merge(witness, other_witness)
                    partial = nxt
                target = out.setdefault(sep, {})
                for q, (count, witness) in partial.items():
                    if q in target:
                        old, old_w = target[q]
                        target[q] = old + count, old_w
                    else:
                        target[q] = count, witness
            messages[v] = out
        return messages[0].get((), {}), messages

    def bridge(self, first_tables, second_tables):
        """Root (chargeA,chargeB) -> (minimum cost, full witnessA, full witnessB)."""
        self.validate_tables(first_tables)
        self.validate_tables(second_tables)
        messages = [None] * len(self.bags)
        for v in reversed(range(len(self.bags))):
            out = {}
            for ta, tb in product(first_tables[v], second_tables[v]):
                a, b = dict(zip(self.bags[v], ta)), dict(zip(self.bags[v], tb))
                sep = (tuple(a[i] for i in self.separators[v]),
                       tuple(b[i] for i in self.separators[v]))
                q = self.local_charge(v, a), self.local_charge(v, b)
                cost = sum(self.costs[i][a[i], b[i]] for i in self.owned[v])
                partial = {q: (cost, a, b)}
                for w in self.children[v]:
                    child_sep = (tuple(a[i] for i in self.separators[w]),
                                 tuple(b[i] for i in self.separators[w]))
                    nxt = {}
                    for (qa, qb), (cost, wa, wb) in partial.items():
                        for (ra, rb), (child_cost, ca, cb) in messages[w].get(child_sep, {}).items():
                            z = self.group.add(qa, ra), self.group.add(qb, rb)
                            value = cost + child_cost
                            if z not in nxt or value < nxt[z][0]:
                                nxt[z] = value, merge(wa, ca), merge(wb, cb)
                    partial = nxt
                target = out.setdefault(sep, {})
                for q, answer in partial.items():
                    if q not in target or answer[0] < target[q][0]:
                        target[q] = answer
            messages[v] = out
        return messages[0].get(((), ()), {}), messages

    def literal(self, tables):
        """Independent global enumeration, intended only for small controls."""
        fibres = {}
        for values in product(*(self.domains[i] for i in self.slots)):
            record = dict(zip(self.slots, values))
            if all(tuple(record[i] for i in b) in allowed
                   for b, allowed in zip(self.bags, tables)):
                q = self.group.zero
                for i in self.slots:
                    q = self.group.add(q, self.charges[i][record[i]])
                fibres.setdefault(q, []).append(record)
        return fibres
