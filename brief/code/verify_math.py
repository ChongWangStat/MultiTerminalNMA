#!/usr/bin/env python3
"""Self-contained verification of source-derived path family and graph identities.
Not a rerun of the archived depression or COVID analyses.
"""
import json
from pathlib import Path
import numpy as np
A=np.array([[1,1,0,0],[0,0,1,1],[1,0,1,0],[0,1,0,1]],float)
f=np.full(4,.5)
checks={}
for t in [0.0,0.13,0.25,0.5]:
 pi=np.array([t,.5-t,.5-t,t])
 checks[f'route_residual_t={t}']=float(np.max(np.abs(A@pi-f)))
checks['route_incidence_rank']=int(np.linalg.matrix_rank(A))
# Bipartite graph edges: treatment 0--study 3, treatment 1--study 3,
# treatment 1--study 4, treatment 2--study 4 (a simple chain)
node_count=5
edge_list=[(0,3,2.0),(1,3,1.0),(1,4,3.0),(2,4,4.0)]
B=np.zeros((len(edge_list),node_count))
g=np.zeros(len(edge_list))
for i,(u,v,conductance) in enumerate(edge_list):
 B[i,u]=1;B[i,v]=-1;g[i]=conductance
L=B.T@np.diag(g)@B
q=np.zeros(node_count);q[0]=1;q[2]=-1
h=np.linalg.pinv(L)@q
current=np.diag(g)@B@h
# Expected transient visits, for walk starting at treatment node 0 with sink 2
D=np.diag(np.diag(L));P=np.eye(node_count)-np.linalg.solve(D,L)
trans=[0,1,3,4]
visits=np.linalg.solve(np.eye(4)-P[np.ix_(trans,trans)].T,np.array([1.,0.,0.,0.]))
n=np.zeros(node_count);n[trans]=visits
d=np.diag(D)
net=np.array([g[i]*(n[u]/d[u]-n[v]/d[v]) for i,(u,v,_) in enumerate(edge_list)])
checks['random_walk_electrical_max_error']=float(np.max(np.abs(current-net)))
# Two-step transition and Kron-reduced walk
C=np.array([[2.,1.,0.],[0.,3.,4.]])
dk=C.sum(axis=1);dj=C.sum(axis=0)
Q=np.diag(1/dj)@C.T@np.diag(1/dk)@C
W=C.T@np.diag(1/dk)@C
np.fill_diagonal(W,0)
Q_noself=Q.copy();np.fill_diagonal(Q_noself,0)
Q_noself/=Q_noself.sum(axis=1)[:,None]
T=W/W.sum(axis=1)[:,None]
checks['two_step_kron_max_error']=float(np.max(np.abs(Q_noself-T)))
# Verify that the multi-terminal conductance is a basis-invariant
# external response, and that the independent-arm special case is
# equivalent to the usual covariance-weighted study Laplacian.
S=np.diag([.04,.09,.16,.25])
C4=np.array([[-1.,1.,0.,0.],[-1.,0.,1.,0.],[-1.,0.,0.,1.]])
Sigma=C4@S@C4.T
K=C4.T@np.linalg.solve(Sigma,C4)
H=np.array([[1.,.2,-.4],[0.,2.,.5],[0.,0.,1.7]])
Ctrans=H@C4
Sigmatrans=H@Sigma@H.T
Ktrans=Ctrans.T@np.linalg.solve(Sigmatrans,Ctrans)
w=1/np.diag(S)
Kpassive=np.diag(w)-np.outer(w,w)/sum(w)
checks['terminal_basis_invariance_max_error']=float(np.max(abs(K-Ktrans)))
checks['independent_arm_study_star_equivalence_max_error']=float(np.max(abs(K-Kpassive)))
checks['terminal_kirchhoff_conservation_error']=float(np.max(abs(np.ones(4)@K)))
q4=np.array([1.,-1.,0.,0.])
checks['two_arm_effective_resistance_error']=float(abs(q4@np.linalg.pinv(K)@q4-(S[0,0]+S[1,1])))
assert max(v for k,v in checks.items() if 'error' in k or 'residual' in k)<1e-10
assert checks['route_incidence_rank']==3
out=Path(__file__).resolve().parents[1]/'validation_results.json';out.write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
