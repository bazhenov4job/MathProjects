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
                print(m_ij)
                # вычисляем алгебраическое дополнение
                a_ij = matrix[0][j] * (-1) ** (0 + j)
                print(a_ij)
                # считаем детерминант, если n > 3 - привет, рекурсия
                det += a_ij * determinant(m_ij)
            return det
    else:
        return None

