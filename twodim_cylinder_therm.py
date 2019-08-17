from matricies import slau_gaus
from math import log1p

"""Скрипт для вычисления поля температур цилиндрической двухмерной шайбы
путём численного решения дифференциального уравнения теплопроводности для
двухмерного случая.
|------------------------------------------------------------------|
|\  \  \  \  \  \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  \  |
| \  \  \  \  \  \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  \ |
|  \  \  \  \  \  \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  \|
|___\__\__\__\__\_ \  \  \  \  \Cu\  \  \  \  \  \  \  \  \  \  \  |
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



# Задаём теплофизические свойства материалов
# Медь
c_Cu = 400  # Дж/ кг*К
ro_Cu = 8930  # кг/ м3
lam_Cu = 400  # Вт/ м*К

a_Cu = lam_Cu / (ro_Cu * c_Cu)

# Фторопласт
c_F = 970  # Дж/ кг*К
ro_F = 2200  # кг/ м3
lam_F = 0.25  # Вт/ м*К

a_F = lam_F / (ro_F * c_F)

t_N2 = 76  # К - температура азота
alfa = 400  # Вт/ м2*к
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
# Всего узлов в модели
nodes = (m + 1) * (n + 1)
r_n = n + 1  # Переход на новый слой по х

# Задаём временной диапазон вычисления, сек
tau = 10000
d_tau = 1

# Начальное распределение
temps = [t_init for i in range(nodes)]

# Создаём матрицу левой части
free_coef = [[0 for i in range(nodes)] for i in range(nodes)]

free_coef[0][0] = (2 * d_tau * alfa / (ro_F * c_F * del_x)) + \
                  (2 * a_F * d_tau / del_x ** 2) + \
                  (8.0 / 5 * a_F * d_tau / (del_r ** 2 * log1p(0.5))) + \
                   1
free_coef[0][1] = (-1) * 8.0 / 5 * a_F * d_tau / (del_r ** 2 * log1p(0.5))
free_coef[0][r_n] = (-1) * 2 * a_F * d_tau / (del_x ** 2)

"""
1. Уравнения для узлов на фторопласте, который граничит с азотом
"""
for i in range(1, int((r_f / del_r)) - 1):
    # Находим текущую величину радиуса n-го узла
    r_cur = del_r * (i + 1)
    free_coef[i][i-1] = (-1) * a_F * d_tau / (log1p(r_cur / (r_cur - del_r)) * r_cur * del_r)
    free_coef[i][i] = (a_F * d_tau / (log1p(r_cur / (r_cur - del_r)) * r_f * del_r)) + \
                      (a_F * d_tau / log1p(r_cur / (r_cur + del_r))) + \
                      (a_F * 2 * d_tau / del_x ** 2) + \
                      (alfa * 2 * d_tau / (ro_F * c_F * del_x)) + 1
    free_coef[i][i+1] = (-1) * a_F * d_tau / (log1p(r_cur / (r_cur + del_r)) * r_cur * del_r)
    free_coef[i][i+r_n] = (-1) * 2 * a_F * d_tau / del_x ** 2

"""
2. Уравнение для узла на стыке фторопласт-медь на нижней поверхности цилиндра
"""
"""
Continue from here
"""

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



