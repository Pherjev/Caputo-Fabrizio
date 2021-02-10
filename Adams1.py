import numpy
import matplotlib.pyplot as plt

A = 0.01
b = 0.01
d = 0.01
p = 0.01
u = 0.01
m = 0.03
rr =0.1
alpha = 0.12
beta = 0.01
gamma = 0.01
delta = 0.001

# Parametros fraccionales

def M(x):
        return 1. #2./(2.-x)

#sigma = 0.6

#M1 = (2*(1-sigma))/float(M(sigma)*(2-sigma))

#M2 = (2*sigma)/float((2-sigma)*M(sigma))



# Funciones simples

def functionS(S,I,R):
	return A - (alpha*S*I)/float(1+beta*S+gamma*I**2) - d*S + p*R

def functionI(S,I,R):
	return (alpha*S*I)/float(1+beta*S+gamma*I**2) - (rr*u*I)/float(1+b*u*I) - d*I - delta*I - m*I

def functionR(S,I,R):
	return (rr*u*I)/float(1+b*u*I) - d*R - p*R + m*I


# Funciones fraccionarias



def xi1(S,I,R):
	return -alpha*I*(1+gamma*I**2)/float((1+beta*S+gamma*I**2)**2)-d

def xi2(S,I,R):
	ter1 = alpha*S*(1+beta*S-gamma*I**2)/float((1+beta*S+gamma*I**2)**2)
	ter2 = -(rr*u)/float((1+b*u*I)**2)
	return ter1 + ter2 -d - delta - m

def dfunctionS(S,I,R):
	return (M2*functionS(S,I,R))/float(1-M1*xi1(S,I,R))

def dfunctionI(S,I,R):
	return (M2*functionI(S,I,R))/float(1-M1*xi2(S,I,R))

def dfunctionR(S,I,R):
	return (M2*functionR(S,I,R))/(1+M1*(d+p))

h = 0.001



SS = []
II = []
RR = []

Sigma = [1, 0.9, 0.6, 0.3]

for ciclos in range(4):

	sigma = Sigma[ciclos]

	s = 0.9
	i = 0.09
	r = 0.01
	t = 0

	s1 = s +0
	i1 = i +0
	r1 = r +0

	s2 = s +0
	i2 = i +0
	r2 = r +0

	s3 = s +0
	i3 = i +0
	r3 = r +0


	St = [s]
	It = [i]
	Rt = [r]
	T  = [t]

	for j in range(200000):
	
		t += h


		if ciclos == 0:
			
		        s1 = s + h*functionS(s,i,r)
      			i1 = i + h*functionI(s,i,r)
        		r1 = r + h*functionR(s,i,r)

			s = s1 + 0
			i = i1 + 0
			r = r1 + 0

                        St.append(s)
                        It.append(i)
                        Rt.append(r)
                        T.append(t)


			
		
		else:
			ans= -(1./M(sigma))*(1-sigma+(4./3.)*h*sigma)*functionS(s2,i2,r2) + (5*h*sigma)/(12*M(sigma))*functionS(s3,i3,r3)
			ani= -(1./M(sigma))*(1-sigma+(4./3.)*h*sigma)*functionI(s2,i2,r2) + (5*h*sigma)/(12*M(sigma))*functionI(s3,i3,r3)
			anr= -(1./M(sigma))*(1-sigma+(4./3.)*h*sigma)*functionR(s2,i2,r2) + (5*h*sigma)/(12*M(sigma))*functionR(s3,i3,r3)

			s1 = s + (1./M(sigma))*(1-sigma + (23./12.)*sigma*h)*functionS(s,i,r) + ans
			i1 = i + (1./M(sigma))*(1-sigma + (23./12.)*sigma*h)*functionI(s,i,r) + ani
			r1 = r + (1./M(sigma))*(1-sigma + (23./12.)*sigma*h)*functionR(s,i,r) + anr

			s3 = s2 + 0
			i3 = i2 + 0
			r3 = r2 + 0
			

			s2 = s + 0
			i2 = i + 0
			r2 = r + 0

			s = s1 + 0			
			i = i1 + 0
			r = r1 + 0

			St.append(s)
			It.append(i)
			Rt.append(r)
			T.append(t)

	SS.append(St)
	II.append(It)
	RR.append(Rt)

fig,ax = plt.subplots()


ax.plot(T,SS[0], color = (0.4,0.7,0.3),label='1')
ax.plot(T,SS[1], ls='--', color = (0.4,0.7,0.3),label='0.9')
ax.plot(T,SS[2], ls='-.', color = (0.4,0.7,0.3),label='0.6')
ax.plot(T,SS[3], ls=':', color = (0.4,0.7,0.3),label='0.3')
#ax.plot(T,It)
#ax.plot(T,Rt)
#ax.set_title(u"SIRS Fractional Model")
ax.spines['left'].set_position(('outward',10))
ax.spines['bottom'].set_position(('outward',10))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_yticks(range(10),minor=True)
ax.set_ylim([0,1])
ax.legend(framealpha=1)
#ax.legend(bbox_to_anchor=(1.1, 1.05))
ax.legend()#loc="upper left", bbox_to_anchor=(0.8,0.2))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.savefig('AB-S1.png')
plt.show()

fig,ax = plt.subplots()

ax.plot(T,II[0], color = (0.7,0.3,0.3),label='1')
ax.plot(T,II[1],ls='--', color = (0.7,0.3,0.3),label='0.9')
ax.plot(T,II[2],ls='-.', color = (0.7,0.3,0.3),label='0.6')
ax.plot(T,II[3],ls=':', color = (0.7,0.3,0.3), label='0.3')
ax.spines['left'].set_position(('outward',10))
ax.spines['bottom'].set_position(('outward',10))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_yticks(range(10),minor=True)
ax.set_ylim([0,1])
ax.legend(framealpha=1)
#ax.legend(bbox_to_anchor=(1.1, 1.05))
ax.legend()#loc="upper left", bbox_to_anchor=(0.8,0.2))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.savefig('AB-I1.png')
plt.show()



fig,ax = plt.subplots()

ax.plot(T,RR[0], color = (0.3,0.4,0.7),label='1')
ax.plot(T,RR[1],ls='--', color = (0.3,0.4,0.7),label='0.9')
ax.plot(T,RR[2],ls='-.', color = (0.3,0.4,0.7),label='0.6')
ax.plot(T,RR[3],ls=':', color = (0.3,0.4,0.7), label='0.3')
ax.spines['left'].set_position(('outward',10))
ax.spines['bottom'].set_position(('outward',10))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_yticks(range(10),minor=True)
ax.set_ylim([0,1])
ax.legend(framealpha=1)
#ax.legend(bbox_to_anchor=(1.1, 1.05))
ax.legend()#loc="upper left", bbox_to_anchor=(0.8,0.2))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.savefig('AB-R1.png')
plt.show()
