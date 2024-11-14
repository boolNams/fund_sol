import sys	

#================== var

E, v = var('E v')
#v = 0.2				# КОЭФ ПУАССОНА
#E = 1.0				# МОДУЛЬ ЮНГА

G = E / (2*(1+v))
Kd = (1/8)*pi*(1-v)*G

x1, x2, ksi1, ksi2 = var('x1 x2 ksi1 ksi2')

ksi1_ = -ksi1		# ksi ШТРИХ
ksi2_ = ksi2

x = x1
c = ksi1

R1 = x1 - ksi1_
R2 = x2 - ksi2_
r1 = x1 - ksi1
r2 = x2 - ksi2
R = sqrt(R1^2+R2^2)
teta = atan(R2/R1)

u11 = -(8*(1-v)^2 - (3-4*v))*ln(R) - ((3-4*v)*R1^2-2*c*x)/R^2 + 4*c*x*R1^2/R^4
u11 = Kd*u11

u12 = (3-4*v)*r1*r2/R^2 + 4*c*x*R1*r2/R^4 + 4*(1-v)*(1-2*v)*teta
u12 = Kd*u12

u21 = (3-4*v)*r1*r2/R^2 + 4*c*x*R1*r2/R^4 - 4*(1-v)*(1-2*v)*teta
u21 = Kd*u21

u22 = -(8*(1-v)^2 - (3-4*v))*ln(R) + ((3-4*v)*r2^2+2*c*x)/R^2 - 4*c*x*r2^2/R^4
u22 = Kd*u22

sol(x1,x2,ksi1,ksi2) = G*(diff(diff(u11,x1),x1) + diff(diff(u11,x2),x2)) + (G/(1-2*v)) * (diff(diff(u11,x1),x1) + diff(diff(u12,x2),x1))
print(sol.full_simplify())
'''
print(u11,'\n')
print(u12,'\n')
print(u21,'\n')
print(u22,'\n')
'''
a1, a2 = var('a1 a2')
r = sqrt(a1^2 + a2^2)
f = (1/2*pi)*ln(1/r)
sol = diff(diff(f,a1),a1) + diff(diff(f,a2),a2)
print(sol.full_simplify())