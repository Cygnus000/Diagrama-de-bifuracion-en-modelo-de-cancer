# codigo adaptado de https://stackoverflow.com/questions/62842261/in-python-how-to-make-a-bifurcation-diagram-of-the-lorenz-system-under-a-varyin
import numpy as np
import matplotlib.pyplot as plt
    # Definicion del sistema de cancer
    # aqui el parametro que vamos a variar es 'a' y los demas fijos
def cancer_system(x, y, z, a, b=2.5, c=0.6, d=1.5, e=4.5, f=1, g=0.2, h=0.5):
    x_dot = x * (1-x)- a * x * y - b * x * z
    y_dot = c * y * (1 - y) - d * x * y
    z_dot = (e * x * z) / (x + f) - g * x * z - h * z
    return x_dot, y_dot, z_dot

da = 0.0001#0.001  # tamaño de paso
a = np.arange(0.85, 1.05, da)  # rango de variacion max 1.04
dt = 0.001  # paso de tiempo
t = np.arange(0, 300, dt)# 100, dt)  # rango temporal

# arreglos de solucion
xs = np.empty(len(t) + 1)
ys = np.empty(len(t) + 1)
zs = np.empty(len(t) + 1)

# valores iniciales x0,y0,z0 para el sistema
xs[0], ys[0], zs[0] = (1,1,1)

# Guardamos las coordenadas de los puntos y graficamos con plt.plot
a_maxes = []
a_mins = []
y_maxes = []
y_mins = []

for A in a:
    print(f"{A=:.4f}") # Imprimimos el avance en el parametro a
    for i in range(len(t)):
        # Runge-kutta para resolver el sistema a lo largo de un tiempo t
        x_dot, y_dot, z_dot = cancer_system(xs[i], ys[i], zs[i], A)
        xs[i + 1] = xs[i] + (x_dot * dt)
        xs[i + 1] = 0.5 *  (xs[i] + (x_dot * dt) + xs[i + 1])
        ys[i + 1] = ys[i] + (y_dot * dt)
        ys[i + 1] = 0.5 *  (ys[i] + (y_dot * dt) + ys[i + 1])
        zs[i + 1] = zs[i] + (z_dot * dt)
        zs[i + 1] = 0.5 *  (zs[i] + (z_dot * dt) + zs[i + 1])
    # calculamos y guardamos los puntos de inflexion
    for i in range(1, len(ys) - 1):
        # guardamos los maximos locales
        if ys[i - 1] < ys[i] and ys[i] > ys[i + 1]:
            a_maxes.append(A)
            y_maxes.append(ys[i])
        # guardamos los minimos locales
        elif ys[i - 1] > ys[i] and ys[i] < ys[i + 1]:
            a_mins.append(A)
            y_mins.append(ys[i])

    # "usamos los valores finales para usarlos como valores iniciales
    # para permanecer cerca del atractor durante la siguiente iteracion"
    xs[0], ys[0], zs[0] = xs[i], ys[i], zs[i]

del a_maxes[:10]
del a_mins[:10]
del y_maxes[:10]
del y_mins[:10]

# Graficamos
    # tipo de letra latex
plt.rcParams.update({
"text.usetex" : True,
"font.family" : 'serif',
"font.serif"  : 'cm'
})
plt.scatter(a_maxes, y_maxes, color="black", s=0.8, alpha=0.2)
plt.scatter(a_mins, y_mins, color="red", s=0.8, alpha=0.2)
plt.legend(['$max(y(t))$','$min(y(t))$'], loc="center left") 
#plt.xlim(0.85, 1.1)
#plt.ylim(0, 1)
plt.title("Diagrama de bifurcación del modelo de Cáncer")
plt.xlabel("parámetro $a$")
plt.ylabel("$y(t)$")
plt.savefig('bifurcation_Y_vs_a.eps', format='eps')
plt.show()
