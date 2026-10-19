import secrets
import random

q=65537
# ~ q=2
# ~ n1=10
# ~n2=10
n1=9
n2=9



F=GF(q)

V1=[(F**n1).random_element() for i in range(n1)]
V2=[(F**n2).random_element() for i in range(n2)]

# ~ V1=[identity_matrix(n1)[i] for i in range(n1)]
# ~ V2=[identity_matrix(n2)[i] for i in range(n2)]


# ~ u=3
# ~ St=[]
# ~ r=random.sample(range(n1),u)
# ~ r=range(u)
# ~ for i in range(n1):
	# ~ St+=[((i+r[j])%n1,i) for j in range(u)]

# ~ St=[(i,j) for j in range(n2) for i in range(2)]



St=[(i,(i+j)%n1) for i in range(n1) for j in range(4)]
# ~ St=[(i,(i+1)%n1) for i in range(n1)]+[(i,i) for i in range(n1)]+[(i,(i+2)%n1) for i in range(n1)]

print(St)

tensor=[]

for i in range(100):
	Mt=Matrix(F,n1,n2)
	
	for (j,k) in St:
		Mt+=Matrix(V1[j]).transpose()*Matrix(V2[k])*F.random_element()
	
	# ~ for j in range(n):
		# ~ for k in range(n):
			# ~ Mt+=Matrix(V1[j]).transpose()*Matrix(V2[k])*F.random_element()
	Mt=[[Mt[i][j] for j in range(n2)] for i in range(n1)]
	
	tensor.append(Mt)
	

# ~ tensor=[[[F.random_element() for k in range(n)] for j in range(n)] for i in range(2*n)]


save(tensor,"public_key")
