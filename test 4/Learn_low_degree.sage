import secrets
import random

q=2
n=10

F=GF(q)
Fe=GF(q^n)
S=GL(n,F).random_element()
T=GL(n,F).random_element()

R=PolynomialRing(Fe,'x')
x=R.gen()

def random_linear_function(S):
	arr=[0]*n
	for i in S:
		arr[i]=Fe([secrets.randbelow(q) for i in range(n)])
		# ~ arr.append(Fe.random_element()) #???
	return arr



def apply(f,x):
	U=T*vector(x)
	el=Fe(list(U))
	
	tot=0
	for i in range(n):
		if f[i]==0:
			continue
		tot+=f[i]*pow(el,pow(q,i))
	assert(len(tot.list())==n)
	return S*vector(tot.list())

def tomatrix(f):
	assert(len(f)==n)
	Mat=[]
	for i in range(n):
		ei=[0]*n
		ei[i]=1
		Mat.append(apply(f,ei))
	return Mat

# ~ St={0,2,7}
St=random.sample(range(n),k=5)
tensor=[tomatrix(random_linear_function(St)) for i in range(2*n)]
# ~ tensor=[[[F.random_element() for k in range(n)] for j in range(n)] for i in range(2*n)]


save(tensor,"public_key")
