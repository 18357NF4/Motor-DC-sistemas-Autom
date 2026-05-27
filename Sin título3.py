import matplotlib.pyplot as plt
import numpy as np
#-----------------------------------DATOS------------------
J_motor = 1.4e-6
J_caja = 1.0e-6
J_carga = 50e-6
Nc = 10
J = J_motor + J_caja + (J_carga/(Nc**2))
B = 3e-6
k = 0.0208
R = 1.53

t_est = 0.5
zita = 1.2
wn = 4/(zita*t_est)

#p3 = (k**2+R*B)/(R*J)
p3=500

kenc = 636.62 * Nc
k_num = k /Nc

a1 = R * J
a0 = (R * B) + (k**2)

#---------------------------Constantes PID------------------------
ki = (a1 * p3 * wn**2) / k_num
kp = (a1 * (2*zita*wn*p3 + wn**2)) / k_num
kd = (a1 * (p3 + 2*zita*wn) - a0) / k_num
#-----------------------------Simulacion--------------------------
dt=0.0001
tiempo=np.arange(0,15,dt)
N=len(tiempo)

U_fisica = np.pi*np.ones(N)
carga=np.zeros(N)


theta=np.zeros(N)
omega=np.zeros(N)
voltaje=np.zeros(N)
integral=0
error_previo=0
carga=np.zeros(N)
for i in range(1,N):
    posicion_carga=theta[i-1]/Nc
    
    pulsos_recibidos = posicion_carga * kenc
    posicion_medida_micro = pulsos_recibidos / kenc
    error = U_fisica[i] - posicion_medida_micro
    proporcional=error*kp
    integral+=ki*error*dt
    derivativo=(error-error_previo)*kd/dt
    V=proporcional+integral+derivativo
    
    if V<-24:
        V=-24
    elif V>24:
        V=24
    
    voltaje[i]=V
    error_previo=error
    aceleracion=omega[i-1]*(-B/J-k/(R*J))+(k*V)/(R*Nc*J)
    omega[i]= omega[i-1]+aceleracion*dt
    theta[i]=theta[i-1]+omega[i-1]*dt
    carga[i]=theta[i]/Nc

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

ax1.plot(tiempo, carga, label='Respuesta del motor', color='blue')
ax1.plot(tiempo, U_fisica, label='Referencia', color='red', linestyle='--') 
ax1.set_title('Simulación PID Estándar')
ax1.set_ylabel('Radianes')
ax1.grid(True)
ax1.legend()

ax2.plot(tiempo, voltaje, label='Voltaje', color='green')
ax2.set_xlabel('Tiempo [segundos]')
ax2.set_ylabel('Voltaje [V]')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()
    
    
    
    
    
    