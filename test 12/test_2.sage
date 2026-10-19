import itertools
from tqdm import tqdm
# ~ Shape=[[1,0],[1,1]]
# ~ Shape=[[1,1]]
Shape=[[1,0],[1,1]]

tot=sum([sum(x) for x in Shape])
samp=4
ch=2
F=GF(ch)
R=PolynomialRing(F,'y',samp*tot)
deg=2
g=R.gens()



cnt=0
Mts=[]
for i in range(samp):
	Mt=[]
	for j in range(len(Shape)):
		row=[]
		for k in range(len(Shape[0])):
			if Shape[j][k]==0:
				row.append(0)
			else:
				row.append(g[cnt])
				cnt+=1
		Mt.append(row)
	Mts.append(Matrix(Mt))



# ~ print(Mts)

pols=[]

for a in itertools.combinations_with_replacement(range(len(g)),deg):
	pols.append(prod([g[x] for x in a]))

# ~ print(pols)

Sts=[[] for i in range(Mts[0].ncols())]

idxs=[]

SolPol=[]

for vpol in itertools.combinations_with_replacement(pols,Mts[0].nrows()):
	SolPol+=list(vpol)
	for j in range(len(Mts)):
		idxs.append(j)
		row=vector(vpol)*Mts[j]
		
	
		for i in range(len(row)):
			Sts[i].append(row[i])

SolPol=vector(SolPol)

Mts2=[Sequence(a).coefficients_monomials()[0] for a in Sts]

Fin=Mts2[0]

for i in range(1,len(Mts2)):
	Fin=Fin.augment(Mts2[i])

print("Rows",Fin.nrows())
print("Cols",Fin.ncols())
print("Rank",Fin.rank())

K=Fin.left_kernel_matrix()
K=Matrix(F,K)

print(K.nrows())
print(K.ncols())

Mt2=[]
for i in range(0,K.ncols(),samp):
	U=K[:,i:i+samp]
	Mt2.append(U.list())

Mt2=Matrix(Mt2)
print("Rows",Mt2.nrows())
print("Cols",Mt2.ncols())
print("Rank",Mt2.rank())

def sub_rnd_tensor(rnd,Test=False):
	newmat=-1
	for u in tqdm(list(K)):
		st=[]
		cnt=0
		assert(len(idxs)==len(u))
		for i in range(0,len(idxs),len(Mts)):
			M=rnd[0]-rnd[0]
			for j in range(len(rnd)):
				M+=rnd[j]*u[i+j]
			st.append([M])
		st=block_matrix(st)
		# ~ print("N",st.nrows())
		# ~ print("C",st.ncols())
		# ~ print(newmat.rank())
		if newmat==-1:
			newmat=st.transpose()
		else:
			newmat=newmat.stack(st.transpose())
		
		newmat=newmat.echelon_form()
		newmat=newmat[:newmat.rank()]
		# ~ print("R",newmat.rank())
		# ~ print(newmat)
	
	if Test:
		S=SolPol.subs(sol)
		print(newmat*S)
	
	print(newmat.nrows())
	print(newmat.ncols())
	print(newmat.rank())
				

# ~ def sub_rnd_mat(rnd):
	# ~ newmat=-1
	# ~ for u in K:
		# ~ st=[]
		# ~ cnt=0
		# ~ for i in range(len(pols)):
			# ~ V=rnd[0]-rnd[0]
			# ~ for j in range(Mt.nrows()):
				# ~ V+=u[cnt]*rnd[idxs[cnt]]
				# ~ cnt+=1
			# ~ st.append(V)
		
		# ~ st=Matrix(st)
		# ~ if newmat==-1:
			# ~ newmat=st
		# ~ else:
			# ~ newmat=newmat.augment(st)

	# ~ print(newmat.nrows())
	# ~ print(newmat.ncols())
	# ~ print(newmat.rank())

# ~ rnd=Matrix(F,0,Mt.ncols())

# ~ S=GL(len(Shape),F).random_element()
# ~ S=Matrix.identity(len(Shape))
# ~ T=GL(len(Shape[0]),F).random_element()
sol={g[i]:F.random_element() for i in range(len(g))}
rnd=[]
for i in range(samp):
	# ~ Mtr=[]
	# ~ for j in range(len(Shape)):
		# ~ row=[]
		# ~ for k in range(len(Shape[0])):
			# ~ if Shape[j][k]==0:
				# ~ row.append(0)
			# ~ else:
				# ~ row.append(F.random_element())
		# ~ Mtr.append(row)
	# ~ Mtr=Matrix(Mtr)
	
	# ~ Mtr=S*Mtr*T
	# ~ print(Mts[i],sol)
	rnd.append(Matrix(F,Mts[i].subs(sol)))

sub_rnd_tensor(rnd,True)

Mt2=[random_matrix(F,Mts[0].nrows(),Mts[0].ncols()) for i in range(samp)]


sub_rnd_tensor(Mt2)


