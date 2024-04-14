# codigo adaptado de https://stackoverflow.com/questions/62842261/in-python-how-to-make-a-bifurcation-diagram-of-the-lorenz-system-under-a-varyin
import numpy as np
import matplotlib.pyplot as plt
    # Definicion del sistema de cancer
    # aqui el parametro que vamos a variar es 'a' y los demas fijos
def cancer_system(x, y, z, c, a=1, b=2.5, d=1.5, e=4.5, f=1, g=0.2, h=0.5):
    x_dot = x * (1-x)- a * x * y - b * x * z
    y_dot = c * y * (1 - y) - d * x * y
    z_dot = (e * x * z) / (x + f) - g * x * z - h * z
    return x_dot, y_dot, z_dot

dc = 0.0001#0.001  # tamaño de paso
c = np.arange(0.3, 1.45, dc)  # rango de variacion  0.85 max 1.05
dt = 0.001  # paso de tiempo
t = np.arange(0, 300, dt)# 100, dt)  # rango temporal

# arreglos de solucion
xs = np.empty(len(t) + 1)
ys = np.empty(len(t) + 1)
zs = np.empty(len(t) + 1)

# valores iniciales x0,y0,z0 para el sistema
xs[0], ys[0], zs[0] = (1,1,1)

# Guardamos las coordenadas de los puntos y graficamos con plt.plot
cx_maxes = []
cx_mins = []
x_maxes = []
x_mins = []

cy_maxes = []
cy_mins = []
y_maxes = []
y_mins = []

cz_maxes = []
cz_mins = []
z_maxes = []
z_mins = []

for C in c:
    print(f"{C=:.4f}") # Imprimimos el avance en el parametro a
    for i in range(len(t)):
        # Runge-kutta para resolver el sistema a lo largo de un tiempo t
        x_dot, y_dot, z_dot = cancer_system(xs[i], ys[i], zs[i], C)
        xs[i + 1] = xs[i] + (x_dot * dt)
        xs[i + 1] = 0.5 *  (xs[i] + (x_dot * dt) + xs[i + 1])
        ys[i + 1] = ys[i] + (y_dot * dt)
        ys[i + 1] = 0.5 *  (ys[i] + (y_dot * dt) + ys[i + 1])
        zs[i + 1] = zs[i] + (z_dot * dt)
        zs[i + 1] = 0.5 *  (zs[i] + (z_dot * dt) + zs[i + 1])
    # calculamos y guardamos los maximos y minimos
    #X
    for i in range(1, len(xs) - 1):
        # guardamos los maximos locales
        if xs[i - 1] < xs[i] and xs[i] > xs[i + 1]:
            cx_maxes.append(C)
            x_maxes.append(xs[i])
        # guardamos los minimos locales
        elif xs[i - 1] > xs[i] and xs[i] < xs[i + 1]:
            cx_mins.append(C)
            x_mins.append(xs[i])
    #Y
    for i in range(1, len(ys) - 1):
        # guardamos los maximos locales
        if ys[i - 1] < ys[i] and ys[i] > ys[i + 1]:
            cy_maxes.append(C)
            y_maxes.append(ys[i])
        # guardamos los minimos locales
        elif ys[i - 1] > ys[i] and ys[i] < ys[i + 1]:
            cy_mins.append(C)
            y_mins.append(ys[i])
    #Z        
    for i in range(1, len(zs) - 1):
        # guardamos los maximos locales
        if zs[i - 1] < zs[i] and zs[i] > zs[i + 1]:
            cz_maxes.append(C)
            z_maxes.append(zs[i])
        # guardamos los minimos locales
        elif zs[i - 1] > zs[i] and zs[i] < zs[i + 1]:
            cz_mins.append(C)
            z_mins.append(zs[i])

    # "usamos los valores finales para usarlos como valores iniciales
    # para permanecer cerca del atractor durante la siguiente iteracion"
    xs[0], ys[0], zs[0] = xs[i], ys[i], zs[i]

# Eliminando los primeros puntos que contienen ruido
del cx_maxes[:10]
del cx_mins[:10]
del x_maxes[:10]
del x_mins[:10]
del cy_maxes[:10]
del cy_mins[:10]
del y_maxes[:10]
del y_mins[:10]
del cz_maxes[:10]
del cz_mins[:10]
del z_maxes[:10]
del z_mins[:10]
# salvando los datos en archivo
xpuntos = []
ypuntos = []
zpuntos = []
for valor in range(len(cx_maxes)):
    xpuntos.append((cx_maxes[valor],x_maxes[valor]))
for valor in range(len(cx_mins)):
    xpuntos.append((cx_mins[valor],x_mins[valor]))
np.savetxt('Bif_X_vs_c.dat', xpuntos)

for valor in range(len(cy_maxes)):
    ypuntos.append((cy_maxes[valor],y_maxes[valor]))
for valor in range(len(cy_mins)):
    ypuntos.append((cy_mins[valor],y_mins[valor]))
np.savetxt('Bif_Y_vs_c.dat', ypuntos)

for valor in range(len(cz_maxes)):
    zpuntos.append((cz_maxes[valor],z_maxes[valor]))
for valor in range(len(cz_mins)):
    zpuntos.append((cz_mins[valor],z_mins[valor]))
np.savetxt('Bif_Z_vs_c.dat', zpuntos)

# Graficamos
    # tipo de letra latex
plt.rcParams.update({
"text.usetex" : True,
"font.family" : 'serif',
"font.serif"  : 'cm'
})
# X
plt.scatter(cx_maxes, x_maxes, color="black", s=0.8, alpha=0.2)
plt.scatter(cx_mins, x_mins, color="red", s=0.8, alpha=0.2)
plt.legend(['$max(x(t))$','$min(x(t))$'], loc="best") 
#plt.xlim(0.85, 1.1)
#plt.ylim(0, 1)
plt.title("Diagrama de bifurcación del modelo de Cáncer")
plt.xlabel("parámetro $c$")
plt.ylabel("$x(t)$")
plt.savefig('bifurcation_X_vs_c.png', format='png')
plt.show()
plt.clf()
# Y
plt.scatter(cy_maxes, y_maxes, color="black", s=0.8, alpha=0.2)
plt.scatter(cy_mins, y_mins, color="red", s=0.8, alpha=0.2)
plt.legend(['$max(y(t))$','$min(y(t))$'], loc="best") 
#plt.xlim(0.85, 1.1)
#plt.ylim(0, 1)
plt.title("Diagrama de bifurcación del modelo de Cáncer")
plt.xlabel("parámetro $c$")
plt.ylabel("$y(t)$")
plt.savefig('bifurcation_Y_vs_c.png', format='png')
plt.show()
plt.clf()
# Z
plt.scatter(cz_maxes, z_maxes, color="black", s=0.8, alpha=0.2)
plt.scatter(cz_mins, z_mins, color="red", s=0.8, alpha=0.2)
plt.legend(['$max(z(t))$','$min(z(t))$'], loc="best") 
#plt.xlim(0.85, 1.1)
#plt.ylim(0, 1)
plt.title("Diagrama de bifurcación del modelo de Cáncer")
plt.xlabel("parámetro $c$")
plt.ylabel("$z(t)$")
plt.savefig('bifurcation_Z_vs_c.png', format='png')
plt.show()
