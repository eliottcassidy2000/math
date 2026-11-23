"""Exact short walks between the odd neighbors of an even integer."""
from math import isqrt


def reflected_path(midpoint, half_roots):
    path = [midpoint - 1]
    for root in half_roots:
        path.append(4*root**2 - path[-1])
    assert path[-1] == midpoint + 1
    assert all(value > 0 and value % 2 for value in path)
    return path


def erase_loops(path):
    result = []
    for label in path:
        if label in result:
            result = result[:result.index(label)+1]
        else:
            result.append(label)
    assert result[0] == path[0] and result[-1] == path[-1]
    assert all(isqrt(x+y)**2 == x+y for x,y in zip(result, result[1:]))
    return result


def universal(midpoint):
    m = midpoint // 2
    assert m > 0 and midpoint == 2*m
    n = isqrt(m)
    if n*n < m:
        n += 1
    n += (n-m) % 2
    deficit = (n*n - m) // 2
    b = max(1, isqrt(deficit))
    a, c = divmod(deficit, b)
    roots = (n+a+1, n+a+b+1, n+b+1, n+c+1, n+c)
    return erase_loops(reflected_path(midpoint, roots))


def euclid(exponent):
    assert exponent >= 5 and exponent % 2
    n = 1 << (exponent-1)
    r = isqrt(n)
    center = 2*n*n - n
    roots = (n-r//2, n-r, n-r//2)
    path = reflected_path(center, roots)
    width = 4*n*r-3*n+1
    assert max(abs(v-center) for v in path) == width
    return center, width, path


def main():
    for midpoint in range(2, 20002, 2):
        path = universal(midpoint)
        assert len(path) in (2, 4, 6) and len(path) == len(set(path))
    print("five-edge construction through 20000; first paths", universal(2), universal(4))
    for exponent in (5, 7, 9, 11):
        center, width, winner = euclid(exponent)
        # Enumerate the two outer vertices in the entire proposed optimal window.
        left = {4*k*k-center+1 for k in range(isqrt((2*center+width)//4)+2)}
        right = {4*k*k-center-1 for k in range(isqrt((2*center+width)//4)+2)}
        left = {v for v in left if abs(v-center) <= width}
        right = {v for v in right if abs(v-center) <= width}
        found = [[center-1,u,v,center+1] for u in left for v in right if isqrt(u+v)**2 == u+v]
        assert found == [winner]
    center, width, path = euclid(61)
    print("E61 minimum three-edge width", width)
    print("E61 nearest path", *path)


if __name__ == "__main__":
    main()
