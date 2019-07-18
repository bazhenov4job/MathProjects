"""Скрипт для вычисления поля температур одномерной пластины путём численного решения
дифференциального уравнения теплопроводности для одномерного случая."""
from matricies import slau_gaus


# Задаём теплофизические свойства материала
c = 400 # Дж/ кг*К
ro = 8930 # кг/ м3
lam = 400 # Вт/ м*К
t_init = 293 # К

# Задаём плотность тепловых потоков
q_l = 2 * 10 ** 4 # Вт/ м2
q_p = 4 * 10 ** 4

# Задаём длину пластины, m. Шаг разбиения
l = 0.1
n = 60

# Задаём временной диапазон вычисления, сек
tau = 1000
d_tau = 1

# Начальное распределение
temps = [t_init for i in range(n)]
temps[0] += 2 * q_l * d_tau / (ro * c * l / (n-1))
temps[n - 1] -= 2 * q_p * d_tau / (ro * c * l / (n-1))

# Вычисляем критерий Фурье
fo = d_tau * lam / (ro * c * (l / n) ** 2)

# Создаём матрицу левой части
free_coef = [[0 for i in range(n)] for i in range(n)]

free_coef[0][0] = (1 + 2 * fo)
free_coef[0][1] = (-2) * fo

for i in range(1, n-1):
    free_coef[i][i-1] = (-1) * fo
    free_coef[i][i] = (2 * fo + 1)
    free_coef[i][i+1] = (-1) * fo

free_coef[n-1][n-2] = (-2) * fo
free_coef[n-1][n-1] = 2 * fo + 1

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



