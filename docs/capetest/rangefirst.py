import numpy as np
LAM=0.1; d0=LAM/2
def fim_crb(pos):
    """phase_m = phi0 + w*x_m + psi*x_m^2 ; FIM = Gram of {1,x,x^2}. Returns CRB(w),CRB(psi)."""
    x=np.asarray(pos,float)
    A=np.vstack([np.ones_like(x),x,x*x]).T
    M=A.T@A
    try: Mi=np.linalg.inv(M)
    except np.linalg.LinAlgError: return np.inf,np.inf
    return Mi[1,1],Mi[2,2]
def sigma_r(crb_psi,r=100.0,snr=1.0):
    # sigma_psi = sqrt(crb_psi/(2*snr)) ; sigma_r = sigma_psi*lam*r^2/pi  (theta=0)
    return np.sqrt(crb_psi/(2*snr))*LAM*r*r/np.pi
def sigma_ang(crb_w,r=100.0,snr=1.0):
    # w = 2pi/lam*sin(th) -> sigma_sin = sigma_w*lam/(2pi)
    return np.sqrt(crb_w/(2*snr))*LAM/(2*np.pi)
def centre(p):
    p=np.asarray(p,float); return p-(p.max()+p.min())/2
# ---- geometries, N=24, all scaled to the SAME aperture D ----
N=24; D=7.75
def ula(N): return np.arange(N)
def nested(N):
    N1=N//2;N2=N-N1
    return np.array(sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)])))[:N]
def sumfirst(N):
    nc=N//2; Do=2*nc-1
    return np.array(sorted(set(list(range(nc))+[Do+i for i in range(N-nc)])))[:N]
def endload(N,frac_mid=0.5):
    """c-optimal-style: 1/4 at each end, rest in the middle cluster."""
    nm=int(round(N*frac_mid)); ne=(N-nm)//2; nm=N-2*ne
    left=np.arange(ne); mid=np.arange(nm); right=np.arange(ne)
    return np.concatenate([left, 1e6+mid, 2e6+right])   # placeholder, rescaled below
def build_endload(N,frac_mid,D):
    nm=int(round(N*frac_mid)); ne=(N-nm)//2; nm=N-2*ne
    g=d0
    left = -D/2 + g*np.arange(ne)
    right=  D/2 - g*np.arange(ne)[::-1]
    mid  = g*(np.arange(nm)-(nm-1)/2)
    return np.concatenate([left,mid,right])
def rescale(p,D):
    p=centre(p); s=D/(p.max()-p.min()); return p*s
rows=[]
for name,p in [("dense ULA",rescale(ula(N)*d0,D)),
               ("nested",rescale(nested(N)*d0,D)),
               ("sum-first (Idea B)",rescale(sumfirst(N)*d0,D))]:
    rows.append((name,p))
for fm in (0.5,0.34,0.25,0.0):
    rows.append((f"end-loaded mid={fm:.2f}",build_endload(N,fm,D)))
print(f"N={N}, all arrays scaled to aperture D={D} m, lambda={LAM} m, r=100 m, SNR=0 dB\n")
print(f"{'array':24s} {'CRB_ang(sin)':>13s} {'CRB_range[m]':>13s} {'range/angle':>12s}")
base=None
for name,p in rows:
    cw,cp=fim_crb(p)
    sa,sr=sigma_ang(cw),sigma_r(cp)
    if name=="nested": base=sr
    print(f"{name:24s} {sa:13.3e} {sr:13.3f} {('%.2fx'%(base/sr) if base else '-'):>12s}")
print("\n(range/angle column = range-RMSE improvement vs nested, higher is better)")
