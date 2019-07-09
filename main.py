from matricies import transpose, determinant, slau_gaus, rowsub, rowmult

# print(transpose([[1, 2, 3], [4, 5, 6]]))
# print(determinant([[-2, -2, -3, -1],
#                   [4, 4, -1, 4],
#                   [5, 6, 0, 2],
#                   [1, 0, 1, 7]]))

print(slau_gaus([[-2, -2, -3, -1], [4, 3, -1, 4], [5, 6, 0, 2], [1, 0, 1, 7]], [1, 2, 3, 4]))
# print(rowsub([1, 2, 3], [0, 2, 3]))
# print(rowmult([1, 2, 3], 2))

