

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

p=5
q=11



print(f"{p = }")
print(f"{q = }")

print("NR:",binomial(q,p)*n2)
print("NC:",binomial(q,p+1)*n3)

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

# ~ print("NR:",M.nrows())
# ~ print("NC:",M.ncols())
print("R:",M.rank())
print("F:",binomial(q-1,p))

print("R/F:",float(M.rank()/binomial(q-1,p)))
print("R/NC:",float(M.rank()/M.ncols()))

print("DIFF",min(M.nrows(),M.ncols())-M.rank())
print("DIFF/F",float((min(M.nrows(),M.ncols())-M.rank())/binomial(q-1,p)))
exit(0)

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
