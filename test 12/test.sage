import itertools

p=2
q=2*p+1
ch=3
F=GF(ch)
R=PolynomialRing(F,'y',p*q)

Mt=Matrix(R,q,p,R.gens())
pols=[]
for a in itertools.combinations(range(q),p):
	pols.append(Mt[list(a)].det())

Sts=[[] for i in range(p)]
idxs=[]
pols2=[]
for pol in pols:
	for j in range(Mt.nrows()):
		row=Mt[j]
		idxs.append(j)
		pols2.append(pol)
		for i in range(Mt.ncols()):
			Sts[i].append(row[i]*pol)

Mts=[Sequence(a).coefficients_monomials()[0] for a in Sts]

Fin=Mts[0]

for i in range(1,len(Mts)):
	Fin=Fin.augment(Mts[i])

print(Fin.nrows())
print(Fin.ncols())
print(Fin.rank())

print(Fin.left_kernel_matrix())

K=Fin.left_kernel_matrix()

# ~ print(K.nrows())
# ~ print(K.ncols())
# ~ for a in K:
	# ~ print(a.hamming_weight())

g=R.gens()
sol={g[i]:F.random_element() for i in range(len(g))}
rnd=Mt.subs(sol)
rnd=Matrix(rnd)
print(rnd)

newmat=-1
for u in K:
	st=[]
	cnt=0
	for i in range(len(pols)):
		V=rnd[0]-rnd[0]
		for j in range(Mt.nrows()):
			V+=u[cnt]*rnd[idxs[cnt]]
			cnt+=1
		st.append(V)
	
	st=Matrix(st)
	if newmat==-1:
		newmat=st
	else:
		newmat=newmat.augment(st)

print(newmat.nrows())
print(newmat.ncols())
print(newmat.rank())
