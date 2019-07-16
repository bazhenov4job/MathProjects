import copy


def transpose(matrix):
    # найдём число строк матрицы
    rows = len(matrix)
    # найдём число столбцов матрицы
    columns = len(matrix[0])
    # создадим пустой массив под болванку матрицы
    empty_matrix = [[0 for x in range(rows)] for y in range(columns)]
    for i in range(rows):
        for j in range(columns):
            if i == j:
                empty_matrix[i][j] = matrix[i][j]
            else:
                empty_matrix[j][i] = matrix[i][j]
    return empty_matrix


def determinant(matrix):
    rows = len(matrix)
    columns = len(matrix[0])
    if rows == columns:
        if rows == 1:
            det = matrix[0][0]
            return det
        elif rows == 2:
            det = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]
            return det
        else:
            # вычисление будет производится через миноры элементов
            det = 0
            for j in range(len(matrix[0])):
                # формируем минор элемента
                m_ij = [[matrix[n][k] for k in range(len(matrix[0])) if k != j] for n in range(len(matrix[0])) if n != 0]
                # вычисляем алгебраическое дополнение
                a_ij = matrix[0][j] * (-1) ** (0 + j)
                if a_ij == 0:
                    det = 0
                else:
                    # считаем детерминант, если минор элемента > 3 - привет, рекурсия
                    det += a_ij * determinant(m_ij)
                    # print(det)
            return det
    else:
        print("Матрица не является квадратной")
        return None


def rowsub(row_1, row_2):
    # вычитает из списка 1 список 2, возвращает результат вычитания
    row = []
    for i in range(len(row_1)):
        row.append(row_1[i] - row_2[i])
    return row


def rowmult(row, mult):
    # умножает список на число, возвращает умноженный список
    # приходится создавать глубокую копию строки, потому что
    # в противном случае она модифицирцется в функции
    new_row = copy.deepcopy(row)
    for i in range(len(new_row)):
        new_row[i] *= mult
    return new_row


def slau_gaus(left_part_matrix, right_part):
    # проверяем совместна ли система, вычисляя детерминант левой части
    #if determinant(left_part_matrix) != 0:
    # создаём дополненную матрицу
    matrix = []
    for i in range(len(left_part_matrix)):
        row = left_part_matrix[i]
        row.append(right_part[i])
        matrix.append(row)
    # начнаем прямую прогонку
    # итерируя по столбцам
    for j in range(len(matrix) - 1):
        # итерируем по строкам, не беря в рассмотрение последнюю
        for i in range(j, len(matrix) - 1):
            if matrix[i+1][j] == 0:
                continue
            else:
                multiplier = matrix[i+1][j] / matrix[j][j]
                matrix[i+1] = rowsub(matrix[i+1], rowmult(matrix[j], multiplier))

    # начинаем обратную прогонку
    # начинаем с последней строки матрицы
    for j in range(len(matrix) - 1, 0, -1):
        # начинаем с последнего столбца
        for i in range(j, 0, -1):
            if matrix[i-1][j] == 0:
                continue
            else:
                multiplier = matrix[i-1][j] / matrix[j][j]
                matrix[i-1] = rowsub(matrix[i-1], rowmult(matrix[j], multiplier))

    # тепер вычислим вектор - столбец неизвестных unknow
    unknow = []
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i == j:
                unknow.append(matrix[i][-1] / matrix[i][j])
    return unknow
    """else:
        print("Система уравнений несовместна")
        return None"""
