import sys

E, v = var('E v')
#v = 0.2				# КОЭФ ПУАССОНА
#E = 1.0				# МОДУЛЬ ЮНГА

G = E / (2*(1+v))

x1, x2, ksi1, ksi2 = var("x1 x2 ksi1 ksi2")
r = sqrt((x1 - ksi1)^2 + (x2 - ksi2)^2)

pd = -1/(8*pi*(1-v)*G)

u11 = pd * ((3-4*v)*ln(r)-diff(r,x1)*diff(r,x1))
u21 = pd * (-diff(r,x1)*diff(r,x2))

sol(x1,x2,ksi1,ksi2) = G*(diff(diff(u11,x1),x1) + diff(diff(u11,x2),x2)) + (G/(1-2*v)) * (diff(diff(u11,x1),x1) + diff(diff(u21,x2),x1))
if sol.full_simplify() == 0:
	print("1) test 1 point complete")


ksi1A, ksi2A = var("ksi1A ksi2A")
ksi1B, ksi2B = ksi1A, -ksi2A

u11A = u11(ksi1 = ksi1A, ksi2 = ksi2A)
u21A = u21(ksi1 = ksi1A, ksi2 = ksi2A)

u11B = u11(ksi1 = ksi1B, ksi2 = ksi2B)
u21B = u21(ksi1 = ksi1B, ksi2 = ksi2B)

u11(ksi1A, ksi2A, x1, x2) = u11A + u11B
u21(ksi1A, ksi2A, x1, x2) = u21A + u21B

sol(ksi1A,ksi2A, x1, x2) = G*(diff(diff(u11,x1),x1) + diff(diff(u11,x2),x2)) + (G/(1-2*v)) * (diff(diff(u11,x1),x1) + diff(diff(u21,x2),x1))
if sol.full_simplify() == 0:
	print("2) test 2 point complete")

sigma111A = 1
sigma121A = 1
sigma221A = 1

sigma111B = 1
sigma121B = 1
sigma221B = 1

sigma111 = sigma111A + sigma111B
sigma121 = sigma121A + sigma121B
sigma221 = sigma221A + sigma221B

eps111 = diff(u11,x1)
eps121 = (1/2)*(diff(u11,x2) + diff(u21,x1))
eps221 = diff(u21,x2)

sigma111 = 2*G*eps111 + ((2*G*v)/(1-2*v))*(eps111 + eps221)
sigma121 = 2*G*eps121
sigma221 = 2*G*eps221 + ((2*G*v)/(1-2*v))*(eps111 + eps221)

print('11',sigma111(ksi1A = 0, ksi2A = -1, x2 = 0,E = 1, v = 0.2).full_simplify())
print('12',sigma121(ksi1A = 0, ksi2A = -1, x2 = 0,E = 1, v = 0.2).full_simplify())
print('22',sigma221(ksi1A = 0, ksi2A = -1, x2 = 0,E = 1, v = 0.2).full_simplify())