"""Sharp route-event bounds conditional on fixed, study-labelled DAG marginals.

The deterministic finite-state product retains accepting AND rejecting paths.
Projection equalities couple product arcs; this is an LP, not ordinary max-flow.
Only numpy/scipy are required. No complete-route enumeration is used here.
"""
from dataclasses import dataclass
from collections import deque
from typing import Callable, Hashable
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


@dataclass(frozen=True)
class Edge:
    tail: Hashable
    head: Hashable
    label: Hashable
    flow: float


@dataclass(frozen=True)
class Automaton:
    n_states: int
    initial: int
    accepting: frozenset
    transition: Callable[[int, int, Edge], int]
    name: str = "event"


def visit_groups(groups, require_all=True):
    """Track whether each label group is visited (same edge may hit both groups)."""
    groups = tuple(frozenset(g) for g in groups)
    if not groups:
        raise ValueError("Supply at least one group.")
    count = 1 << len(groups)
    def transition(q, i, e):
        for j, group in enumerate(groups):
            if e.label in group:
                q |= 1 << j
        return q
    accept = {count-1} if require_all else set(range(1, count))
    return Automaton(count, 0, frozenset(accept), transition,
                     "visits_all_groups" if require_all else "visits_any_group")


def ordered_groups(groups):
    """Visit successive label groups on distinct edges; intervening edges allowed."""
    groups = tuple(frozenset(g) for g in groups)
    if not groups:
        raise ValueError("Supply at least one group.")
    k = len(groups)
    return Automaton(k+1, 0, frozenset({k}),
                     lambda q, i, e: q+1 if q < k and e.label in groups[q] else q,
                     "ordered_groups")


def contains_segment(edge_indices):
    """Accept if the exact edge-index sequence occurs contiguously (KMP states)."""
    pattern = tuple(edge_indices)
    if not pattern:
        raise ValueError("The segment must be nonempty.")
    k = len(pattern)
    alphabet = set(pattern)
    cache = {}
    for q in range(k):
        for symbol in alphabet:
            seq = pattern[:q] + (symbol,)
            cache[q, symbol] = max(j for j in range(min(k, len(seq))+1)
                                   if j == 0 or seq[-j:] == pattern[:j])
    def transition(q, i, e):
        return k if q == k else cache.get((q, i), 0)
    return Automaton(k+1, 0, frozenset({k}), transition, "contains_segment")


def complement(event):
    return Automaton(event.n_states, event.initial,
                     frozenset(set(range(event.n_states))-set(event.accepting)),
                     event.transition, "not_" + event.name)


def validate_flow(edges, source, sink, tol=1e-9):
    if source == sink or not edges:
        raise ValueError("A nonempty source-sink flow is required.")
    nodes = list(dict.fromkeys([source, sink] + [v for e in edges for v in (e.tail, e.head)]))
    outgoing = {v: [] for v in nodes}
    incoming = {v: [] for v in nodes}
    for i, e in enumerate(edges):
        if not np.isfinite(e.flow) or e.flow <= 0 or e.tail == e.head:
            raise ValueError("Edges must have finite strictly positive flows and no loops.")
        outgoing[e.tail].append(i)
        incoming[e.head].append(i)
    if incoming[source] or outgoing[sink]:
        raise ValueError("The source has incoming edges or the sink has outgoing edges.")
    degrees = {v: len(incoming[v]) for v in nodes}
    queue = deque(v for v in nodes if degrees[v] == 0)
    order = []
    while queue:
        v = queue.popleft()
        order.append(v)
        for i in outgoing[v]:
            w = edges[i].head
            degrees[w] -= 1
            if degrees[w] == 0:
                queue.append(w)
    if len(order) != len(nodes):
        raise ValueError("The support must be acyclic.")
    reachable = {source}
    for v in order:
        if v in reachable:
            reachable.update(edges[i].head for i in outgoing[v])
    coreachable = {sink}
    for v in reversed(order):
        if any(edges[i].head in coreachable for i in outgoing[v]):
            coreachable.add(v)
    if any(v not in reachable or v not in coreachable for v in nodes):
        raise ValueError("Every support vertex must lie on a complete route.")
    total = sum(edges[i].flow for i in outgoing[source])
    balance_error = max(abs(sum(edges[i].flow for i in outgoing[v]) -
                            sum(edges[i].flow for i in incoming[v]) -
                            (total if v == source else -total if v == sink else 0)) for v in nodes)
    if balance_error > tol * max(1., total):
        raise ValueError(f"Flow is not conserved: {balance_error:g}.")
    return order, outgoing, total


def product_problem(edges, source, sink, event):
    """Sparse equalities and objective; forward reachable states, no acceptance pruning."""
    edges = tuple(edges)
    order, outgoing, total = validate_flow(edges, source, sink)
    if event.n_states < 1 or not 0 <= event.initial < event.n_states:
        raise ValueError("Invalid automaton initial state.")
    if any(not 0 <= q < event.n_states for q in event.accepting):
        raise ValueError("Invalid accepting state.")
    reachable = {v: set() for v in order}
    reachable[source].add(event.initial)
    arcs = []
    for v in order:
        for q in sorted(reachable[v]):
            for i in outgoing[v]:
                e = edges[i]
                qnext = event.transition(q, i, e)
                if not isinstance(qnext, (int, np.integer)) or not 0 <= qnext < event.n_states:
                    raise ValueError("Transitions must be total and deterministic on reachable pairs.")
                arcs.append((i, q, int(qnext)))
                reachable[e.head].add(int(qnext))
    product_nodes = [(v,q) for v in order for q in sorted(reachable[v])]
    conservation = {(v,q): len(edges)+j for j,(v,q) in
                    enumerate((v,q) for v,q in product_nodes if v != sink)}
    rhs = np.zeros(len(edges)+len(conservation))
    rhs[:len(edges)] = [e.flow for e in edges]
    rhs[conservation[source,event.initial]] = total
    rr, cc, vals = [], [], []
    cost = np.zeros(len(arcs))
    for j,(i,q,qnext) in enumerate(arcs):
        e = edges[i]
        rr.extend([i, conservation[e.tail,q]])
        cc.extend([j,j]); vals.extend([1.,1.])
        if e.head != sink:
            rr.append(conservation[e.head,qnext]); cc.append(j); vals.append(-1.)
        elif qnext in event.accepting:
            cost[j] = 1.
    matrix = coo_matrix((vals,(rr,cc)), shape=(len(rhs),len(arcs))).tocsr()
    return {"matrix":matrix,"rhs":rhs,"cost":cost,"arcs":arcs,"edges":edges,
            "source":source,"sink":sink,"event":event,"total":total,
            "n_product_nodes":len(product_nodes)}


def _certificate(problem, solution, cost):
    a, b, x = problem['matrix'], problem['rhs'], solution.x
    dual = solution.eqlin.marginals
    reduced = cost-a.T@dual
    return {"primal_residual":float(np.max(np.abs(a@x-b))),
            "nonnegative_violation":float(max(0.,-x.min())),
            "dual_violation":float(max(0.,-reduced.min())),
            "duality_gap":float(abs(cost@x-b@dual))}


def sharp_event_bounds(edges, source, sink, event, witnesses=False):
    problem = product_problem(edges, source, sink, event)
    result = {"event":event.name,"n_states":event.n_states,
              "variables":len(problem['arcs']),"equalities":len(problem['rhs']),
              "nonzeros":problem['matrix'].nnz,"product_nodes":problem['n_product_nodes']}
    for endpoint, sign in [('lower',1.),('upper',-1.)]:
        cost = sign*problem['cost']
        fit = linprog(cost, A_eq=problem['matrix'], b_eq=problem['rhs'],
                      bounds=(0,None), method='highs',
                      options={"primal_feasibility_tolerance":1e-9,
                               "dual_feasibility_tolerance":1e-9})
        if not fit.success:
            raise RuntimeError(f"Product LP {endpoint} failed: {fit.message}")
        result[endpoint] = float(sign*fit.fun)
        result[endpoint+'_certificate'] = _certificate(problem,fit,cost)
        if witnesses:
            result[endpoint+'_witness'] = decompose_product_solution(problem,fit.x)
    return result


def decompose_product_solution(problem, x, tol=1e-10):
    """Construct an attaining sparse route allocation, useful for independent checks."""
    residual = np.array(x,copy=True)
    adj = {}
    for j,(i,q,qn) in enumerate(problem['arcs']):
        e = problem['edges'][i]
        adj.setdefault((e.tail,q),[]).append(j)
    start = problem['source'],problem['event'].initial
    paths = []
    while any(residual[j] > tol for j in adj.get(start,[])):
        node = start; picked = []
        while node[0] != problem['sink']:
            candidates = [j for j in adj.get(node,[]) if residual[j] > tol]
            if not candidates:
                raise RuntimeError("Numerical conservation failure in witness decomposition.")
            j = candidates[0]; picked.append(j)
            i,q,qn = problem['arcs'][j]
            node = problem['edges'][i].head,qn
        mass = float(min(residual[j] for j in picked))
        for j in picked:
            residual[j] -= mass
        paths.append({"edges":[problem['arcs'][j][0] for j in picked], "mass":mass,
                      "accepted":node[1] in problem['event'].accepting})
    if residual.size and np.max(np.abs(residual)) > 100*tol:
        raise RuntimeError("Witness decomposition left unexplained flow.")
    return paths


def accepts(path, edges, event):
    state = event.initial
    for i in path:
        state = event.transition(state,i,edges[i])
    return state in event.accepting
