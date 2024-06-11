from sympy import *




r=symbols("r",cls=Symbol,positive=True)

theta0A,theta0B,theta1A,theta1B=symbols("t0A,t0B,t1A,t1B",real=True)
mA,mB=symbols("mA,mB",cls=Symbol,positive=True)

theta0A=-mB/mA*theta0B-theta1A-mB/mA*theta1B


d0A0B=r*(theta0B-theta0A)

d0B1A=r*(theta1A-theta0B)

d1A1B=r*(theta1B-theta1A)

d1B0A=r*(2*pi+theta0A-theta1B)


c1,c2,a1,a2=symbols("c1,c2,a1,a2",cls=Symbol)




V=c1*(d0A0B-a1)**2+c2*(d0B1A-a2)**2\
    +c1*(d1A1B-a1)**2+c2*(d1B0A-a1)**2


V = simplify(V.expand())




f=c1*((1+mB/mA)*theta0B+theta1A+mB/mA*theta1B)**2\
    +c2*(theta1A-theta0B)**2\
    +c1*(theta1B-theta1A)**2\
    +c2*(2*pi+(-1-mB/mA)*theta1B-mB/mA*theta0B-theta1A)**2


g=2*c1*a1*((1+mB/mA)*theta0B+theta1A+mB/mA*theta1B)\
    +2*c2*a2*(theta1A-theta0B)\
    +2*c1*a1*(theta1B-theta1A)\
    +2*c2*a2*(2*pi+(-1-mB/mA)*theta1B-mB/mA*theta0B-theta1A)

rst=f*r**2-g*r+2*c1*a1**2+2*c2*a2**2

diff=V-rst

pprint(simplify(rst.expand()))