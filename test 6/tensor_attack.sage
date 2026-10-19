

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
q=9



print(f"{p = }")
print(f"{q = }")

def build_matrix(tensor):
	
	pi=[i for i in range(n1)]
	random.shuffle(pi)
	Mat=[]

	for S in itertools.combinations(range(q),p):
		for j in range(n2):
			row=[]
			
			for U in itertools.combinations(range(q),p+1):
				for k in range(n3):
					tot=0
					sgn=1
					for i in range(q):
						
						if set(list(U))==set(list(S)+[i]):
							tot+=sgn*tensor[pi[i]][j][k]
						if i in U:
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

print("NR:",M.nrows())
print("NC:",M.ncols())
print("R:",M.rank())
print("F:",binomial(q-1,p))

print("R/F:",float(M.rank()/binomial(q-1,p)))
print("R/NC:",float(M.rank()/M.ncols()))

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


print(V1,V2)

Int=V1.intersection(V2)

print(Int)
