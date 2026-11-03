import itertools
def D(S): return sorted({a-b for a in S for b in S})
def P(S): return sorted({a+b for a in S for b in S})
def runlen(X):
    X=sorted(set(X)); best=cur=1
    for i in range(1,len(X)):
        cur=cur+1 if X[i]==X[i-1]+1 else 1; best=max(best,cur)
    return best
def nested(N0):
    N1=N0//2; N2=N0-N1
    return sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)]))[:N0]
def twoblock(nc,delta):
    return list(range(nc))+[delta+i for i in range(nc)]

print("=== Idea B as written in the PDF: core {0,1,2,3}, offset 10 ===")
S=twoblock(4,10); print(" S =",S)
print(f"  sum  P: {P(S)[:5]}... holefree-run={runlen(P(S))}  (span {min(P(S))}..{max(P(S))})")
print(f"  diff D: holefree-run={runlen(D(S))}  totalDOF={len(D(S))}")

print("\n=== the offset that actually makes the sum hole-free (delta=7) ===")
S=twoblock(4,7); print(" S =",S)
print(f"  sum  P holefree-run={runlen(P(S))} span {min(P(S))}..{max(P(S))}  -> uDOFsum={runlen(P(S))}")
print(f"  diff D holefree-run={runlen(D(S))} totalDOF={len(D(S))} -> uDOFdiff={runlen(D(S))}")

print("\n=== cost vs nested at the SAME N=8 ===")
nb=nested(8); print(" nested S =",nb)
print(f"  nested uDOFdiff={runlen(D(nb))} totalDOF={len(D(nb))} ; uDOFsum={runlen(P(nb))}")

print("\n=== exhaustive search N=8: best uDOFsum, and its uDOFdiff cost ===")
best=[]
for S in itertools.combinations(range(0,22),8):
    if S[0]!=0: continue
    us,ud=runlen(P(S)),runlen(D(S))
    best.append((us,ud,S))
best.sort(key=lambda z:(-z[0],-z[1]))
print("  top 5 by uDOFsum:")
for us,ud,S in best[:5]: print(f"    uDOFsum={us:3d} uDOFdiff={ud:3d}  S={list(S)}")
bd=sorted(best,key=lambda z:(-z[1],-z[0]))
print("  top 5 by uDOFdiff:")
for us,ud,S in bd[:5]: print(f"    uDOFsum={us:3d} uDOFdiff={ud:3d}  S={list(S)}")
