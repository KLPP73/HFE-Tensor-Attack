import itertools
import random

pubk=load("public_key.sobj")

m=len(pubk)
n=pubk[0].nrows()

# ~ CorrectOrderTensor=[[[pubk[k][i][j] for k in range(m)] for j in range(n)] for i in range(n)]
CorrectOrderTensor=[[[pubk[k][i][j] for k in range(m)] for j in range(n)] for i in range(n)]
# ~ CorrectOrderTensor=pubk
Fq=CorrectOrderTensor[0][0][0].base_ring()
n1=n
n2=n
n3=m

q=10
p=q*n3/(n2+n3)
p=int(round(p))

print(p,q)

print(binomial(q,p)*n2,binomial(q,p+1)*n3)


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

	print(Mat.nrows())
	print(Mat.ncols())
	print(Mat.rank())
	return Mat


M=build_matrix(CorrectOrderTensor)

# ~ print(M.nrows())
# ~ print(M.ncols())
# ~ print(M.rank())


# ~ exit(0)


# ~ print(M.rank(),binomial(q-1,p))

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
