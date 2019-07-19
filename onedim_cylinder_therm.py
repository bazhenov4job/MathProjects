"""Скрипт для вычисления поля температур цилиндра путём численного решения
дифференциального уравнения теплопроводности для одномерного случая."""
"""Автор Баженов Вадим, email = bazhenov4job@gmail.com"""

from matricies import slau_gaus
from math import log1p


# Задаём теплофизические свойства материала
c = 400 # Дж/ кг*К
ro = 8930 # кг/ м3
lam = 400 # Вт/ м*К
t_init = 293 # К

# Задаём плотность тепловых потоков
q_l = 2 * 10 ** 4 # Вт/ м2
q_p = 4 * 10 ** 4

# Задаём начальный радиус, конечный радиус, число слоёв разбиения
r_init = 0.05
r_final = 0.1
n = 51
# находим приращение радиуса
del_r = (r_final - r_init) / (n-1)

# Задаём временной диапазон вычисления, сек
tau = 10
d_tau = 1

# Начальное распределение
temps = [t_init for i in range(n)]
temps[0] += 2 * q_l * r_init * d_tau / (ro * c * del_r)
temps[n - 1] -= 2 * q_p * r_final * d_tau / (ro * c * del_r)

# Вычисляем основную часть "критерия Фурье"
fo = d_tau * lam / (ro * c * del_r)

# Создаём матрицу левой части
free_coef = [[0 for i in range(n)] for i in range(n)]

free_coef[0][0] = (1 + 2 * fo / log1p(r_init / (r_init + del_r)))
free_coef[0][1] = (-2) * fo / log1p(r_init / (r_init + del_r))

for i in range(1, n-1):
    # Находим текущую величину радиуса n-го узла
    r_cur = r_init + del_r * i
    free_coef[i][i-1] = (-1) * fo / log1p(r_cur / (r_cur - del_r))
    free_coef[i][i] = (fo * (1 / log1p(r_cur / (r_cur - del_r)) +
                             1 / log1p(r_cur / (r_cur + del_r))) + 1)
    free_coef[i][i+1] = (-1) * fo / log1p(r_cur / (r_cur + del_r))

free_coef[n-1][n-2] = (-2) * fo / log1p(r_final / (r_final - del_r))
free_coef[n-1][n-1] = 2 * fo / log1p(r_final / (r_final - del_r)) + 1

# Вычисляем распределение на новом временном шаге
i = 0
while tau >= d_tau * i:
    new_temp = slau_gaus(free_coef, temps)
    temps = new_temp
    temps[0] += 2 * q_l * d_tau / (ro * c * l / (n-1))
    if temps[n-1] > 77.6:
        temps[n - 1] -= 2 * q_p * d_tau / (ro * c * l / (n-1))
    else:
        temps[n - 1] = 77.6
    if i % 10 == 0:
        result = open("result.txt", "a")
        string_list = []
        for temp in temps:
            string_list.append(str(temp))
        results = '\t'.join(string_list)
        results += '\n'
        result.write(results)
        result.close()
        print("Итерация № {} выполнена".format(i))
    i += 1
print(new_temp)



