def hamilton(n):
    # Input n is the code 2**order * (2*adjacency+1) - 1.
    header = n + 1
    order = 0
    while header % 2 == 0:
        header = header // 2
        order += 1
    adjacency = header // 2
    paircount = order * (order - 1) // 2
    if adjacency >= 2 ** paircount:
        return 0
    total = 0
    candidate = 0
    bound = order ** order
    while candidate < bound:
        digits = candidate
        seen = 0
        good = 1
        position = 0
        previous = 0
        while position < order:
            vertex = digits % order
            digits = digits // order
            vertexbit = 2 ** vertex
            if (seen // vertexbit) % 2 == 0:
                seen += vertexbit
            else:
                good = 0
            if position > 0:
                if previous < vertex:
                    edgeindex = vertex * (vertex - 1) // 2 + previous
                    edge = (adjacency // (2 ** edgeindex)) % 2
                else:
                    if vertex < previous:
                        edgeindex = previous * (previous - 1) // 2 + vertex
                        edge = 1 - (adjacency // (2 ** edgeindex)) % 2
                    else:
                        edge = 0
                if edge == 0:
                    good = 0
            previous = vertex
            position += 1
        if good > 0:
            total += 1
        candidate += 1
    return total
