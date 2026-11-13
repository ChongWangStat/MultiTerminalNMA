#!/usr/bin/env python3
"""Public compact check of v19.9 passive synthesis and NMA covariance sensitivity.

Synthetic multi-study data only; original psoriasis inputs are not needed.
Fixed/plug-in study covariances: no uncertainty from estimating them is modeled.
"""
import numpy as np
from itertools import permutations

# Four-arm correlated covariance: valid variance triangles, no positive-resistor mesh.
K=np.array([[1.95,0.05,-1,-1],[0.05,1.95,-1,-1],[-1,-1,3,-1],[-1,-1,-1,3.]])
C=np.column_stack((-np.ones(3),np.eye(3)))
S=C@np.linalg.pinv(K)@C.T
B=np.linalg.cholesky(np.linalg.inv(S)).T@C
assert min(np.linalg.eigvalsh(S))>0 and K[0,1]>0
np.testing.assert_allclose(B.T@B,K,atol=1e-12)
assert np.allclose(B@np.ones(4),0)
R=np.array([[(np.eye(4)[i]-np.eye(4)[j])@np.linalg.pinv(K)@(np.eye(4)[i]-np.eye(4)[j]) for j in range(4)] for i in range(4)])
assert all(R[i,j] <= R[i,l]+R[l,j]+1e-12 for i,j,l in permutations(range(4),3))
v=np.array([.2,-.8,1.1,.3]);np.testing.assert_allclose(v@K@v,np.linalg.norm(B@v)**2,atol=1e-12)

# Independent five-study NMA with correlated within-study contrasts.
rng=np.random.default_rng(20261113)
n=5;studies=[[0,1,2,3],[2,3,4],[0,4],[1,2,4],[0,3]]
X=[]; Sig=[];y=[]
for k,arms in enumerate(studies):
    m=len(arms);Ci=np.column_stack((-np.ones(m-1),np.eye(m-1))); Xi=Ci@np.eye(n)[arms]
    M=rng.normal(size=(m-1,m-1));Si=M@M.T+(.5+.2*k)*np.eye(m-1)
    X.append(Xi);Sig.append(Si);y.append(rng.normal(size=m-1))
q=np.eye(n)[0]-np.eye(n)[4];Z=np.eye(n)[:,:-1]
def fit(covs):
    L=sum(A.T@np.linalg.solve(Si,A) for A,Si in zip(X,covs))
    rhs=sum(A.T@np.linalg.solve(Si,yy) for A,Si,yy in zip(X,covs,y))
    theta=Z@np.linalg.solve(Z.T@L@Z,Z.T@rhs)
    phi=Z@np.linalg.solve(Z.T@L@Z,Z.T@q)
    return float(q@theta),float(q@phi),theta,phi,L
base,V,theta,phi,L=fit(Sig)
u=[np.linalg.solve(Si,A@phi) for Si,A in zip(Sig,X)]
energies=[float(ui@Si@ui) for ui,Si in zip(u,Sig)]
assert np.isclose(sum(energies),V,atol=1e-12)
rho=np.array(energies)/V;assert np.all(rho>=0) and np.isclose(sum(rho),1)
# 50 exact rank-one covariance perturbations: no GLS refit needed for formula.
max_err=0.;rnd=np.random.default_rng(919)
for k in range(5):
    Si=Sig[k];Ai=X[k];res=y[k]-Ai@theta
    for _ in range(10):
        h=rnd.normal(size=len(Si));z=np.linalg.solve(Si,h)
        g=Ai.T@z;gap=float(h@z-g@np.linalg.pinv(L)@g)
        a=float(u[k]@h);b=float(h@np.linalg.solve(Si,res))
        assert gap>=-1e-10
        t=float(rnd.uniform(.01,3))
        covs=list(Sig);covs[k]=Si+t*np.outer(h,h)
        newdelta,newV,*_=fit(covs)
        predV=V+t*a*a/(1+t*gap)
        preddelta=base-t*a*b/(1+t*gap)
        max_err=max(max_err,abs(newV-predV),abs(newdelta-preddelta))
        assert np.allclose([newV,newdelta],[predV,preddelta],atol=1e-10)
        assert newV>=V-1e-10
print('Four-arm passive factorization: PASS; positive resistor-only realization: FAIL as expected')
print('Comparison variance:',format(V,'.9f'))
print('GLS energy shares:',np.round(rho,7),'sum',format(sum(rho),'.12f'))
print('Exact covariance update checks:',50,'max error:',format(max_err,'.3e'))
print('ALL TESTS PASSED')
