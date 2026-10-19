import itertools

def dimensions(tensor):
	n1=len(tensor)
	n2=len(tensor[0])
	n3=len(tensor[0][0])
	return (n1,n2,n3)

def switch_dimensions(tensor,a,b):
	if a==b:
		return tensor
	n1,n2,n3=dimensions(tensor)
	if a==0 and b==1:
		newtensor=[[[tensor[i][j][k] for k in range(n3)] for i in range(n1)] for j in range(n2)]
		return newtensor
	
	if a==0 and b==2:
		newtensor=[[[tensor[i][j][k] for i in range(n1)] for j in range(n2)] for k in range(n3)]
		return newtensor
	
	if a==1 and b==2:
		newtensor=[[[tensor[i][j][k] for j in range(n2)] for k in range(n3)] for i in range(n1)]
		return newtensor

def construct_matrix(tensor,p,q):
	n1,n2,n3=dimensions(tensor)
	
	assert(n1>=q)
	

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
	
	
	
	
	
	Mat=[]
	# ~ print(n2*len(bmS1),n2)
	# ~ print(n3*len(bmS2),n3)
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
							tot+=sgn*tensor[i][j][k]
						if ((1<<i)&U)>0:
							sgn*=-1
					row.append(tot)
			Mat.append(row)

	Mat=Matrix(Mat)
	return Mat

def compute_stats(tensor,p,q):
	Mat=construct_matrix(tensor,p,q)
	return (Mat.nrows(),Mat.ncols(),Mat.rank())

def construct_attack(tensor,p,q):
	M=construct_matrix(tensor,p,q)
	n1,n2,n3=dimensions(tensor)
	M2=construct_matrix([[[1 for i in range(n3)] for j in range(n2)] for k in range(n1)],p,q)
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
	
	# ~ print(V2)
	
	Int=V1.intersection(V2)
	
	return Int

