

import itertools
import random
pubk=load("public_key.sobj")
Fq=pubk[0][0][0].base_ring()

# ~ n1=len(pubk[0])
# ~ n2=len(pubk[0][0])
# ~ n3=len(pubk)

# ~ CorrectOrderTensor=[[[pubk[k][i][j] for k in range(n3)] for j in range(n2)] for i in range(n1)]

n1=len(pubk)
n2=len(pubk[0])
n3=len(pubk[0][0])
CorrectOrderTensor=[[[pubk[i][j][k] for k in range(n3)] for j in range(n2)] for i in range(n1)]

p=4
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
	# ~ random.shuffle(pi)
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
print("CALCULATING RANK")
print("NR:",M.nrows())
print("NC:",M.ncols())
print("R:",M.rank())
print("F:",binomial(q-1,p))

print("R/F:",float(M.rank()/binomial(q-1,p)))
print("R/NC:",float(M.rank()/M.ncols()))

Ker=M.left_kernel_matrix()

print(Ker.nrows())
print(Ker.nrows()//n2)


V=VectorSpace(Fq,M.nrows())
V1=V.span(Ker)
print(V1)

for St in itertools.product([0,1],repeat=n2):
	if sum(St)>3:
		continue
	
	V2=[]
	for i in range(0,Ker.ncols()):
		if St[i%n2]==0:
			continue
		ei=[0]*Ker.ncols()
		ei[i]=1
		V2.append(ei)
	
	V2=V.span(V2)
	
	Int=V1.intersection(V2)
	print(Int,St)
	
	# ~ exit(0)
	

# ~ for i in range(0,Ker.ncols(),n2):
	# ~ print(Ker[:,i:i+n2])

# ~ M2=build_matrix([[[1 for i in range(n3)] for j in range(n2)] for k in range(n1)])
# ~ V=VectorSpace(Fq,M.nrows())
# ~ V1=V.span(M.transpose())
# ~ Sp=[]
# ~ for i in range(M.nrows()):
	# ~ if M2[i][0]!=0:
		# ~ vec=[0]*M.nrows()
		# ~ vec[i]=1
		# ~ Sp.append(vector(Fq,vec))

# ~ Sp=Matrix(Sp)

# ~ V2=V.span(Sp)


# ~ print(V1,V2)

# ~ Int=V1.intersection(V2)

# ~ print(Int)

# ~ print(Int.basis_matrix())
