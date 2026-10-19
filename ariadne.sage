import secrets
from tqdm import tqdm
import itertools
q=5
n=15
F=GF(q)
Fe.<z>=GF(q^n)
R.<x>=PolynomialRing(Fe,sparse=True)

TT=0

D=q+2

Lm=int(log(D,2))+1

Pk=dict({})

DD=set({})
DD2=set({})

# ~ for i in range(Lm):
	# ~ if pow(q,i)<=D:
		# ~ Pk[tuple([i])]=Fe.random_element()

# ~ for i in range(Lm):
	# ~ for j in range(i,Lm):
		# ~ if pow(q,i)+pow(q,j)<=D:
			# ~ Pk[(i,j)]=Fe.random_element()

for i in range(Lm):
	for j in range(i,Lm):
		for k in range(j,Lm):
			if pow(q,i)+pow(q,j) + pow(q,k)<=D:
				print(i,j,k)
				Pk[(i,j,k)]=Fe.random_element()
				


Mtg=GL(n,F)

S=Mtg.random_element()
T=Mtg.random_element()

R=PolynomialRing(F,'x',n)
Bs=R.gens()

Re=PolynomialRing(Fe,'y',n)
Bse=Re.gens()


pert=[]
for _ in range(TT):
	Pldeg2=0
	for i in range(n):
		for j in range(n):
			Pldeg2+=Bse[i]*Bse[j]*Fe.random_element()

	Pldeg1=0
	for i in range(n):
		Pldeg1+=Bs[i]*F.random_element()
	pert.append((Fe.random_element(),Pldeg1,Pldeg2))

print(pert)

OrgPl=0

for exp in Pk:
	coef=Pk[exp]
	val=1
	for y in exp:
		val*=pow(x,pow(q,y))
	OrgPl+=val*coef

print(OrgPl)



def encrypt(x):
	xbar=S*vector(x)
	el=Fe(list(xbar))
	tot=0
	for exp in Pk:
		coef=Pk[exp]
		val=1
		for Y in exp:
			val*=pow(el,pow(q,Y))
		tot+=val*coef
	
	tot=Fe(tot)
	
	for (coef,Pldeg1,Pldeg2) in pert:
		tot+=coef*Pldeg1(*x)*Pldeg2(*x)

	# ~ print(Pldeg1(*x),Pldeg2(*x))
	
	
	tot=vector(list(tot))
	
	tot=T*tot
	
	return tot



Mns=[]

# ~ for lst in itertools.combinations_with_replacement(Bs,1):
	# ~ Mns.append(prod(lst))

# ~ for lst in itertools.combinations_with_replacement(Bs,2):
	# ~ Mns.append(prod(lst))

for lst in itertools.combinations_with_replacement(Bs,3):
	Mns.append(prod(lst))

# ~ print(Mns)
Mat=[]
Tar=[]

for U in range(3,4):
	for TT in tqdm(itertools.combinations_with_replacement(range(n),U)):
		arr=[0]*n
		for y in TT:
			arr[y]+=1
		
		row=[]
		for a in Mns:
			row.append(a(*arr))
		Mat.append(row)
		
		Tar.append(encrypt(vector(F,arr)))
			

Mat=Matrix(Mat)
print(Mat.nrows())
print(Mat.ncols())

print(Mat.rank())
X=Mat.solve_right(Matrix(Tar))

pols=[]
for a in X.transpose():
	pol=a*vector(Mns)
	pols.append(pol)

pols=vector(pols)

r=(F^n).random_element()

assert(pols(*r)==encrypt(r))

Mt=[]

for i in tqdm(range(n)):
	Base=[derivative(pol,Bs[i]) for pol in pols]
	for j in range(i,n):
		F2=[derivative(pol,Bs[j]) for pol in Base]
		F2=F2+list(Bs)
		Mat=Sequence(F2).coefficients_monomials()[0]
		Mat=Mat[:n]
		Mt.append(Mat.list())

Mt=Matrix(Mt)

print(Mt.nrows())
print(Mt.ncols())
print(Mt.rank())

K=Mt.right_kernel_matrix()

print("="*50)
print(K.nrows())
print(K.ncols())
print("="*50)
MTS=[]

for k in range(n*n):
	MT=[int(i==k) for i in range(n*n)]
	MT=Matrix(F,n,n,MT)
	R=[]
	for i,a in enumerate(Mt):
		X=Matrix(F,len(a)//n,n,a)*MT
		U=K*vector(X.list())
		R+=list(U)
	MTS.append(R)

MTS=Matrix(MTS)

print(MTS.nrows())
print(MTS.ncols())
print(MTS.rank())

field=MTS.left_kernel_matrix()

print(field.nrows())
print(field.ncols())



print(field)

# ~ for M in field:
	# ~ Mt=Matrix(F,n,n,list(M))
	# ~ print((T^(-1))*Mt*T)
	# ~ print(Mt.charpoly().factor())
	# ~ print()

# ~ for A in field:
	# ~ for B in field:
		# ~ MtA=Matrix(F,n,n,list(A))
		# ~ MtB=Matrix(F,n,n,list(B))
		# ~ if MtA*MtB!=MtB*MtA:
			# ~ print("ERROR")


# ~ Seq=[]
# ~ for i in tqdm(range(n)):
	# ~ for j in range(n):
		# ~ F=derivative(pols[j],Bs[i])
		# ~ Seq.append(F)

# ~ for i in tqdm(range(n)):
	# ~ alpha=[0]*n
	# ~ alpha[i]=1
	# ~ alpha=vector(alpha)
	# ~ v=vector(Bs)
	
	# ~ U=pols(*list(v+alpha))-pols(*Bs)-pols(*list(alpha))
	# ~ for pol in U:
		# ~ Seq.append(pol)

# ~ Seq=Sequence(Seq)

# ~ Mat,Mon=Seq.coefficients_monomials()

# ~ print(Mat.nrows())
# ~ print(Mat.ncols())
# ~ print(Mat.rank())
