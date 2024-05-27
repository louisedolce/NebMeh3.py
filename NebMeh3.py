import math as m
import numpy as np
e=0
epsilon = 23.43929111
epsilon=m.radians(epsilon)
k = 0.01720209895

file = open('nebo33.txt', 'r')
ok = file.read().split('\n')
a = float(ok[0]) #предполагаемая орбита
alf_1h = float(ok[1])
alf_1m = float(ok[2])
alf_1s = float(ok[3])
alf_2h = float(ok[4])
alf_2m = float(ok[5])
alf_2s = float(ok[6])
dlt_1d = float(ok[7])
dlt_1m = float(ok[8])
dlt_1s = float(ok[9])
dlt_2d = float(ok[10])
dlt_2m = float(ok[11])
dlt_2s = float(ok[12])
X_1 = float(ok[13])
X_2 = float(ok[14])
Y_1 = float(ok[15])
Y_2 = float(ok[16])
Z_1 = float(ok[17])
Z_2 = float(ok[18])
t_1 = float(ok[19])
t_2 = float(ok[20])
file.close()

def doli(x,y,z): #Функция перевода часов и градусов из минут и секунд в доли часа\градусов
    if x>0:
        c=x+y/60+z/3600
    else:
        c = x-y/60-z/3600
    return c
def h_g(x):
    c=m.radians(x*15)
    return c

alf = [h_g(doli(alf_1h,alf_1m,alf_1s)),h_g(doli(alf_2h,alf_2m,alf_2s))]
dlt = [m.radians(doli(dlt_1d,dlt_1m,dlt_1s)),m.radians(doli(dlt_2d,dlt_2m,dlt_2s))]
X = [X_1, X_2]
Y = [Y_1, Y_2]
Z = [Z_1, Z_2]
t = [t_1, t_2]

#Создание списков
A=[]
B=[]
C=[]
RcosTet=[]
R_2=[]
RsinTet_2=[]

# A,B,C,Rcos,Rsin,R
eps=0.005
for i in range(len(alf)):
    A.append(m.cos(alf[i])*m.cos(dlt[i]))
    B.append(m.sin(alf[i])*m.cos(dlt[i]))
    C.append(m.sin(dlt[i]))
    if abs(m.pow(A[i],2)+m.pow(B[i],2)+m.pow(C[i],2))-1<eps: #проверка
        print("Проверка А,В,С пройдена",str(i+1)+"-й набор")
    else:
        print("Проверка А,В,С не пройдена",str(i+1)+"-й набор")
        quit()
    RcosTet.append(-(A[i]*X[i]+B[i]*Y[i]+C[i]*Z[i]))
    R_2.append(X[i]**2+Y[i]**2+Z[i]**2)
    RsinTet_2.append(R_2[i]-RcosTet[i]**2)

# Создание доп списков
Ro=np.zeros(2) #строка из двух элементов
x=np.zeros(2)
y=np.zeros(2)
z=np.zeros(2)

def c(l):
    if l>=0:
        for i in range(2):
                Ro[i]=m.sqrt(l**2-RsinTet_2[i])-RcosTet[i]
                x[i]=A[i]*Ro[i]-X[i]
                y[i]=B[i]*Ro[i]-Y[i]
                z[i]=C[i]*Ro[i]-Z[i]
        sin_2fg=(1/4/l**2)*((x[0]-x[1])**2+(y[0]-y[1])**2+(z[0]-z[1])**2)
        fg=m.asin(m.sqrt(sin_2fg))
        fd=k*(t[1]-t[0])/(2*(l)**1.5)
        d=fg-fd
        return d
    else:
        print("Большая полуось принимает отрицательное значение")
        return quit()
a_2=0
a_1=a+0.1
while abs(a_2-a)!=0:
    d1=c(a)
    d2=c(a_1)
    a_2=a_1-d2*(a_1-a)/(d2-d1)
    a=a_1
    a_1=a_2

def c1(l):
    for i in range(2):
            Ro[i]=m.sqrt(m.pow(l,2)-RsinTet_2[i])-RcosTet[i]
            x[i]=A[i]*Ro[i]-X[i]
            y[i]=B[i]*Ro[i]-Y[i]
            z[i]=C[i]*Ro[i]-Z[i]
    sin_2fg=1/(4*m.pow(l,2))*(m.pow((x[0]-x[1]),2)+m.pow((y[0]-y[1]),2)+m.pow((z[0]-z[1]),2))
    fg=m.asin(m.sqrt(sin_2fg))
    fd=k*(t[1]-t[0])/(2*(l)**(3/2))
    f=(fg+fd)/2
    return f
f=c1(a)
Px=(x[0]+x[1])/(2*m.cos(f)*a)
Py=(y[0]+y[1])/(2*m.cos(f)*a)
Pz=(z[0]+z[1])/(2*m.cos(f)*a)
Qx=(x[1]-x[0])/(2*m.sin(f)*a)
Qy=(y[1]-y[0])/(2*m.sin(f)*a)
Qz=(z[1]-z[0])/(2*m.sin(f)*a)
Sum_P = m.pow(Px,2) + m.pow(Py,2) + m.pow(Pz,2)
Sum_Q=m.pow(Qx,2) + m.pow(Qy,2) + m.pow(Qz,2)
Mult_PQ=Px*Qx+Py*Qy+Pz*Qz
if abs(Sum_Q-1)<0.005 and abs(Sum_P-1)<0.005 and abs(Mult_PQ)<0.005:
    print("Проверка P и Q пройдена")
else:
    print("Проверка P и Q не пройдена")
    quit()

# Поиск элементов орбиты
def f(x, y):
    if y > 0 and x > 0:
        z = m.degrees(m.asin(y))
    elif y > 0 and x < 0:
        z = 180 - m.degrees(m.asin(y))
    elif y < 0 and x < 0:
        z = 180 + abs(m.degrees(m.asin(y)))
    else:
        z = 360 - abs(m.degrees(m.asin(y)))
    return z
sini=m.sqrt(m.pow((Qz*m.cos(epsilon)-Qy*m.sin(epsilon)),2)+m.pow((Pz*m.cos(epsilon)-Py*m.sin(epsilon)),2))
sin_omega=(Pz*m.cos(epsilon)-Py*m.sin(epsilon))/sini
cos_omega=(Qz*m.cos(epsilon)-Qy*m.sin(epsilon))/sini
omega=f(cos_omega,sin_omega)
sin_Bomega=(Py*m.cos(m.radians(omega))-Qy*m.sin(m.radians(omega)))/m.cos(epsilon)
cos_Bomega=Px*m.cos(m.radians(omega))-Qx*m.sin(m.radians(omega))
Bomega=f(cos_Bomega,sin_Bomega)
cosi=-(Px*m.sin(m.radians(omega))+Qx*m.cos(m.radians(omega)))/m.sin(m.radians(Bomega))
i=f(cosi,sini)
AvAnomaly=omega
t0=(t[1]+t[0])/2

# Ответ
rezultat3 = open('itog.txt', 'w')
rezultat3.write(str(e)+"\n")
rezultat3.write(str(a)+' а.е.'+"\n")
rezultat3.write(str(AvAnomaly)+'°'+"\n") #аругмент шир
rezultat3.write(str(Bomega)+'°'+"\n")
rezultat3.write(str(i)+'°'+"\n")
rezultat3.write('JD '+str(t0)+"\n")
rezultat3.close()