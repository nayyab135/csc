import numpy as np
LAM=0.1; d0=LAM/2
def runlen(X):
    X=sorted(set(int(v) for v in X)); best=cur=1
    for i in range(1,len(X)):
        cur=cur+1 if X[i]==X[i-1]+1 else 1; best=max(best,cur)
    return best
def psl(pos):
    s=np.linspace(-1,1,40001)
    B=np.abs(np.exp(1j*2*np.pi*np.outer(s,pos)/LAM).sum(1))/len(pos)
    i0=np.argmax(B); i=i0
    while i+1<len(B) and B[i+1]<B[i]: i+=1
    j=i0
    while j-1>=0 and B[j-1]<B[j]: j-=1
    m=B.copy(); m[j:i+1]=0; return m.max()
def crb(pos,r=100.0):
    x=np.asarray(pos,float); A=np.vstack([np.ones_like(x),x,x*x]).T
    Mi=np.linalg.inv(A.T@A)
    return np.sqrt(Mi[2,2]/2)*LAM*r*r/np.pi
def rep(n,name):
    n=np.array(sorted(set(int(v) for v in n)))
    pos=(n-(n.max()+n.min())/2)*d0
    d={int(a-b) for a in n for b in n}; s={int(a+b) for a in n for b in n}
    ap=pos.max()-pos.min()
    print(f"{name:34s} {len(n):3d} {ap:7.2f} {2*ap*ap/LAM:9.0f} {psl(pos):7.3f} {runlen(d):8d} {runlen(s):8d} {crb(pos):8.2f}")
def nested(N):
    N1=N//2;N2=N-N1
    return sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)]))[:N]
def symnested(N):
    h=N//2;a1=h//2;a2=h-a1
    half=sorted(set(list(range(1,a1+1))+[(i+1)*(a1+1) for i in range(a2)]))[:h]
    return [-v for v in half[::-1]]+half
def sha(half): return [-v for v in half[::-1]]+list(half)
print(f"{'array':34s} {'N':>3s} {'aper[m]':>7s} {'Rayl[m]':>9s} {'PSL':>7s} {'uDOFdif':>8s} {'uDOFsum':>8s} {'sr[m]':>8s}")
rep(nested(24),"nested (NOT symmetric)")
rep(symnested(24),"sym-nested (near-field baseline)")
rep(sha([1,2,3,4,5,6,7,8,16,32,48,64]),"SHA greedy N=24  <-- PROPOSED")
print()
rep(nested(8),"nested N=8")
rep(symnested(8),"sym-nested N=8")
rep(sha([1,2,3,6]),"SHA greedy N=8")
