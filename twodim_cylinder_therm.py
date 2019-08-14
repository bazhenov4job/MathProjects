"""Скрипт для вычисления поля температур цилиндрической двухмерной шайбы
путём численного решения дифференциального уравнения теплопроводности для
двухмерного случая.
|------------------------------------------------------------------|
|\  \  \  \  \  \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  \  |
| \  \  \  \  \  \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  \ |
|  \  \  \  \  \  \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  \|
|_________________ \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  |
|\/ \/ \/F\/ \/ \/| \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \ |
|/\ /\ /\F/\ /\ /\|  \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \|
 |↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓||                     |↑  ↑  ↑  ↑  ↑  ↑|          |
    Q_n           |                     |    Q_heat     |          |
                  |                     |               |          |
-------r_f------->|                     |               |          |
                                        |               |          |
-------------------r_in---------------->|               |          |
                                                        |          |
---------------------------r_out----------------------->|          |
                                                                   |
---------------------------r_final-------------------------------->|

Схема задачи приведена выше
Cu - медь
F - фторопласт
Q_n -  тепловой поток уходит в азот
Q_heat - тепловой поток приносится нагревателем
задача осесимметрична относительно центральной оси цилиндрического тела
"""
"""Автор Баженов Вадим, email: bazhenov4job@gmail.com"""

from matricies import slau_gaus
from math import log1p


# Задаём теплофизические свойства материалов
# Медь
c_Cu = 400  # Дж/ кг*К
ro_Cu = 8930  # кг/ м3
lam_Cu = 400  # Вт/ м*К
# Фторопласт
c_F = 970  # Дж/ кг*К
ro_F = 2200  # кг/ м3
lam_F = 0.25  # Вт/ м*К

t_init = 293  # К

# Задаём плотность тепловых потоков
q_h = 3.266 * 10 ** 3  # Вт/ м2

# Задаём геометрические характеристики, число слоёв разбиения
r_f = 0.025
r_in = 0.05
r_out = 0.08
r_final = 0.145
del_r = 0.005
del_x = 0.002
x_final = 0.01
# находим максимальный индекс узла на слое (узлов всего на слое: n + 1)
n = r_final / del_r
m = x_final / del_x

# Задаём временной диапазон вычисления, сек
tau = 10000
d_tau = 1

# Начальное распределение
temps = [t_init for i in range(n)]

"""
Continue from this point
"""

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
    temps[0] += 2 * q_l * r_init * d_tau / (ro * c * del_r)
    temps[n - 1] -= 2 * q_p * r_final * d_tau / (ro * c * del_r)
    new_temp = slau_gaus(free_coef, temps)
    temps = new_temp
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



