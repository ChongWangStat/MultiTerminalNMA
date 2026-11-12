import csv
from collections import defaultdict
from pathlib import Path
p=Path(__file__).with_name('COVID19_AE_Netlist.csv')
rows=list(csv.DictReader(p.open()))
assert len(rows)==11
flow=defaultdict(float)
for r in rows:
    cur=float(r['current'])
    flow[r['from']]-=cur
    flow[r['to']]+=cur
for node in ('A','B','C','D','E','RECOVERY','SOLIDARITY'):
    desired= -1 if node=='A' else 1 if node=='E' else 0
    assert abs(flow[node]-desired)<1e-10,(node,flow[node])
print('PASS: 11 resistive branches; 5 treatments + 2 multi-terminal studies; all reported terminal current balances exact.')
for node in sorted(flow): print(f'  {node:>12} signed inflow minus outflow = {flow[node]:+.6f}')
