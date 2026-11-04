import numpy as np
LAM=0.1; d0=LAM/2; N=24; D=7.75
def centre(p):
    p=np.asarray(p,float); return p-(p.max()+p.min())/2
def rescale(p,D):
    p=centre(p); return p*(D/(p.max()-p.min()))
def nested(N):
    N1=N//2;N2=N-N1
    return np.array(sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)])))[:N]
def sumfirst(N):
    nc=N//2;Do=2*nc-1
    return np.array(sorted(set(list(range(nc))+[Do+i for i in range(N-nc)])))[:N]
def endload(N,fm,D):
    nm=int(round(N*fm)); ne=(N-nm)//2; nm=N-2*ne
    left=-D/2+d0*np.arange(ne); right=D/2-d0*np.arange(ne)[::-1]
    mid=d0*(np.arange(nm)-(nm-1)/2)
    return np.concatenate([left,mid,right])
def psl(pos):
    """peak sidelobe of the far-field beampattern (ambiguity risk)."""
    s=np.linspace(-1,1,20001)
    B=np.abs(np.exp(1j*2*np.pi*np.outer(s,pos)/LAM).sum(1))/len(pos)
    i0=np.argmax(B)
    # walk down from the main peak to its first null, mask it out
    i=i0
    while i+1<len(B) and B[i+1]<B[i]: i+=1
    j=i0
    while j-1>=0 and B[j-1]<B[j]: j-=1
    m=B.copy(); m[j:i+1]=0
    return B.max(), m.max()
print(f"N={N}, aperture D={D} m  -- peak sidelobe level = AMBIGUITY risk (lower is better)\n")
print(f"{'array':24s} {'PSL':>8s} {'PSL [dB]':>10s}  {'min gap [m]':>12s}")
for name,p in [("dense ULA",rescale(np.arange(N)*d0,D)),
               ("nested",rescale(nested(N)*d0,D)),
               ("sum-first (Idea B)",rescale(sumfirst(N)*d0,D)),
               ("end-loaded mid=0.50",endload(N,0.50,D)),
               ("end-loaded mid=0.34",endload(N,0.34,D)),
               ("end-loaded mid=0.25",endload(N,0.25,D))]:
    pk,sl=psl(p); g=np.diff(np.sort(p)); g=g[g>1e-9].min()
    print(f"{name:24s} {sl:8.3f} {20*np.log10(sl):10.1f}  {g:12.3f}")
