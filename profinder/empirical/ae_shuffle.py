"""Altschul-Erikson dinucleotide-preserving shuffle (Altschul & Erikson,
Mol. Biol. Evol. 1985).
"""

import random
from collections import defaultdict, Counter


def dinucleotide_counts(s):
    """Counter of the len(s)-1 adjacent dinucleotides in s."""
    return Counter(s[i:i + 2] for i in range(len(s) - 1))


def ae_shuffle(s, rng=None, max_tries=10000):
    """Return a dinucleotide-preserving shuffle of ``s``.
    """
    if rng is None:
        rng = random
    n = len(s)
    if n < 3:
        return s

    last = s[-1]
    graph = defaultdict(list)          # vertex -> list of successor symbols
    for i in range(n - 1):
        graph[s[i]].append(s[i + 1])
    vertices = list(graph.keys())

    def random_last_edges():
        """Pick, for each vertex != last, one out-edge to traverse last;
        return the mapping iff those edges form a tree rooted at ``last``
        (no cycle among non-last vertices), else None."""
        le = {v: rng.choice(graph[v]) for v in vertices if v != last}
        for v in le:
            seen = set()
            u = v
            while u != last:
                if u in seen or u not in le:
                    return None        # cycle / dead end -> not an arborescence
                seen.add(u)
                u = le[u]
        return le

    last_edge = None
    for _ in range(max_tries):
        le = random_last_edges()
        if le is not None:
            last_edge = le
            break
    if last_edge is None:
        return s                        # give up: return input (comp. preserved)

    # Per-vertex ordered out-edge list: permute all edges, then force one
    # instance of the designated last-edge target to the very end.
    order = {}
    for v in vertices:
        edges = list(graph[v])
        rng.shuffle(edges)
        if v != last:
            tgt = last_edge[v]
            edges.remove(tgt)
            edges.append(tgt)
        order[v] = edges

    # Walk the Euler path from the original first symbol.
    idx = defaultdict(int)
    out = [s[0]]
    cur = s[0]
    for _ in range(n - 1):
        nxt = order[cur][idx[cur]]
        idx[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)
