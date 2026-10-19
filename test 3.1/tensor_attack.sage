import itertools
import random
pubk=load("public_key.sobj")

n=pubk.rank()
print("n =",n)
n3=pubk.ncols()

pubk=[[pubk[i][j] for j in range(pubk.ncols())] for i in range(n3)]
for i in range(n):
	for j in range(i,n):
		pubk[j][i]+=pubk[i][j]

pubk=Matrix(pubk)
Fe=pubk.base_ring()

tensor=[[pubk[i][j].list() for j in range(pubk.ncols())] for i in range(pubk.nrows())]
n1=len(tensor[0][0])
n2=n

Fq=tensor[0][0][0].base_ring()


# ~ for i in range(len(tensor)):
	# ~ for j in range(len(tensor[i])):
		# ~ assert(len(tensor)==n)

CorrectOrderTensor=[[[tensor[j][k][i] for k in range(n3)] for j in range(n2)] for i in range(n1)]

n1=n
n2=n

p=3
q=2*p+1

print(f"{p = }")
print(f"{q = }")

S1=itertools.combinations(range(q),p)
S2=itertools.combinations(range(q),p+1)

bmS1=[]
for a in S1:
	tot=0
	for i in range(q):
		if i in a:
			tot+=(1<<i)
	bmS1.append(tot)

bmS2=[]
for a in S2:
	tot=0
	for i in range(q):
		if i in a:
			tot+=(1<<i)
	bmS2.append(tot)



def build_matrix(tensor):
	
	pi=[i for i in range(n1)]
	Mat=[]

	for S in bmS1:
		for j in range(n2):
			row=[]
			
			for U in bmS2:
				for k in range(n3):
					tot=0
					sgn=1
					
					
					# ~ print(bin(U))
					for i in range(q):
						
						if U==(S|(1<<i)):
							tot+=sgn*tensor[pi[i]][j][k]
						if ((1<<i)&U)>0:
							sgn*=-1
					row.append(tot)
			Mat.append(row)

	Mat=Matrix(Mat)

	return Mat

def ToMatrix(VV):
	cnt=0
	Mt=[[0 for i in range(n2)] for j in range(p+1)]
	for S in itertools.combinations(range(q),p):
		for j in range(n2):
			for i in range(p+1):
				if set(list(S)+[i])==set(range(p+1)):
					Mt[i][j]=pow(Fq(-1),i-1)*VV[cnt]
			cnt+=1
	return Mt

M=build_matrix(CorrectOrderTensor)

print(M.nrows())
print(M.ncols())
print("Rank:",M.rank())

print("Factor:",binomial(q-1,p))
M2=build_matrix([[[1 for i in range(n3)] for j in range(n2)] for k in range(n1)])
V=VectorSpace(Fq,M.nrows())
V1=V.span(M.transpose())
Sp=[]
for i in range(M.nrows()):
	if M2[i][0]!=0:
		vec=[0]*M.nrows()
		vec[i]=1
		Sp.append(vector(Fq,vec))

Sp=Matrix(Sp)

V2=V.span(Sp)




Int=V1.intersection(V2)

print("I",Int)

Mts=[]
for v in Int.basis_matrix():
	Mt=Matrix(ToMatrix(v))
	Mts.append(Mt.transpose())

def solveMinRank(Mats):
	k=len(Mats)
	n=Mats[0].ncols()
	nguess=k//Mats[0].nrows()+1
	
	# ~ print("I",k,Mats[0].nrows(),Mats[0].ncols(),nguess)
	# ~ exit(0)
	TOT=0
	while True:
		TOT+=1
		guesses=[(Fq^n).random_element() for i in range(nguess)]
		SolMat=[]
		for Mt in Mats:
			row=[]
			for g in guesses:
				row+=list(Mt*g)
			row=vector(row)
			SolMat.append(row)
		SolMat=Matrix(SolMat)
		
		sol=SolMat.left_kernel_matrix()
		# ~ print(sol.nrows())
		if sol.nrows()>=1:
			sol=SolMat.left_kernel().random_element()
			S=Mats[0]-Mats[0]
			for (Mt,Ml) in zip(Mats,sol):
				S+=Mt*Ml
			
			if S.rank()<=1 and S!=0:
				return S
Mt=[]
for _ in range(20):
	M=solveMinRank(Mts)
	v=M.transpose().echelon_form()[0]
	if Mt==[]:
		Mt=Matrix([v])
	else:
		Mt=Mt.stack(Matrix([v]))

Mt=Mt.echelon_form()
Mt=Mt[:Mt.rank()]
print(Mt)
