#include "MatrixSolver.h"
#include <cmath>

QString MatrixSolver::formatMatrix(const QVector<QVector<double>>& matrix,
                                   const QVector<double>& freeTerms) {
    if (matrix.isEmpty() || freeTerms.isEmpty() || matrix.size() != freeTerms.size()) {
        return "Invalid matrix format";
    }

    QString output;
    for (int i = 0; i < matrix.size(); ++i) {
        if (matrix[i].size() == 0) continue;

        output += "| ";
        for (int j = 0; j < matrix[i].size(); ++j) {
            output += QString::number(matrix[i][j], 'f', 4) + "\t";
        }
        output += "| \t" + QString::number(freeTerms[i], 'f', 4) + " |\n";
    }
    return output;
}

QString MatrixSolver::formatVector(const QVector<double>& vec) {
    if (vec.isEmpty()) return "[]";

    QString output = "[ ";
    for (int i = 0; i < vec.size(); ++i) {
        output += QString::number(vec[i], 'f', 4);
        if (i < vec.size() - 1) output += ", ";
    }
    output += " ]";
    return output;
}

bool MatrixSolver::solveSystem(Method method,
                               const QVector<QVector<double>>& inputMatrix,
                               const QVector<double>& inputFreeTerms,
                               QVector<double>& results,
                               QStringList* steps) {
    // Очищаем результаты и шаги перед началом
    results.clear();
    if (steps) steps->clear();

    // Проверка базовых условий
    if (inputMatrix.isEmpty() || inputFreeTerms.isEmpty() ||
        inputMatrix.size() != inputFreeTerms.size()) {
        if (steps) steps->append("Invalid input: matrix or free terms are empty or sizes mismatch");
        return false;
    }

    switch (method) {
    case Gauss:
        return solveGauss(inputMatrix, inputFreeTerms, results, steps);
    case JordanGauss:
        return solveJordanGauss(inputMatrix, inputFreeTerms, results, steps);
    case Cramer:
        return solveCramer(inputMatrix, inputFreeTerms, results, steps);
    case InverseMatrix:
        return solveInverseMatrix(inputMatrix, inputFreeTerms, results, steps);
    case LeastSquares:
        return solveLeastSquares(inputMatrix, inputFreeTerms, results, steps);
    default:
        if (steps) steps->append("Unknown solving method");
        return false;
    }
}

bool MatrixSolver::solveGauss(QVector<QVector<double>> matrix,
                              QVector<double> freeTerms,
                              QVector<double>& results,
                              QStringList* steps) {
    const int n = matrix.size();
    if (n == 0 || n != freeTerms.size()) return false;

    if (steps) {
        steps->append("Метод Гаусса:");
        steps->append("Начальная система:");
        steps->append(formatMatrix(matrix, freeTerms));
    }

    // Прямой ход
    for (int i = 0; i < n; ++i) {
        // Поиск опорного элемента
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // Перестановка строк
        if (maxRow != i) {
            std::swap(matrix[i], matrix[maxRow]);
            std::swap(freeTerms[i], freeTerms[maxRow]);
            if (steps) {
                steps->append(QString("Перестановка строк %1 и %2:").arg(i+1).arg(maxRow+1));
                steps->append(formatMatrix(matrix, freeTerms));
            }
        }

        // Проверка на вырожденность
        if (std::abs(matrix[i][i]) < 1e-10) {
            if (steps) steps->append("Система вырождена!");
            return false;
        }

        // Нормализация строки
        double divisor = matrix[i][i];
        if (steps) {
            steps->append(QString("Нормировка строки %1 (деление на %2):")
                              .arg(i+1).arg(divisor, 0, 'f', 4));
        }

        for (int j = i; j < n; ++j) {
            matrix[i][j] /= divisor;
        }
        freeTerms[i] /= divisor;

        if (steps) steps->append(formatMatrix(matrix, freeTerms));

        // Исключение элементов
        for (int k = i + 1; k < n; ++k) {
            double factor = matrix[k][i];
            if (steps) {
                steps->append(QString("Исключение в строке %1 (коэффициент %2 * строка %3):")
                                  .arg(k+1).arg(factor, 0, 'f', 4).arg(i+1));
            }

            for (int j = i; j < n; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
            }
            freeTerms[k] -= factor * freeTerms[i];

            if (steps) steps->append(formatMatrix(matrix, freeTerms));
        }
    }

    // Обратный ход
    if (steps) steps->append("Обратный ход:");
    results.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        results[i] = freeTerms[i];
        for (int j = i + 1; j < n; ++j) {
            results[i] -= matrix[i][j] * results[j];
        }

        if (steps) {
            steps->append(QString("x%1 = %2")
                              .arg(i+1)
                              .arg(results[i], 0, 'f', 4));
        }
    }

    return true;
}

bool MatrixSolver::solveJordanGauss(QVector<QVector<double>> matrix,
                                    QVector<double> freeTerms,
                                    QVector<double>& results,
                                    QStringList* steps) {
    const int n = matrix.size();
    if (n == 0 || n != freeTerms.size()) return false;

    if (steps) {
        steps->append("Метод Жордана-Гаусса:");
        steps->append("Начальная система:");
        steps->append(formatMatrix(matrix, freeTerms));
    }

    for (int i = 0; i < n; ++i) {
        // Поиск опорного элемента
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // Перестановка строк
        if (maxRow != i) {
            std::swap(matrix[i], matrix[maxRow]);
            std::swap(freeTerms[i], freeTerms[maxRow]);
            if (steps) {
                steps->append(QString("Перестановка строк %1 и %2:").arg(i+1).arg(maxRow+1));
                steps->append(formatMatrix(matrix, freeTerms));
            }
        }

        // Проверка на вырожденность
        if (std::abs(matrix[i][i]) < 1e-10) {
            if (steps) steps->append("Система вырождена!");
            return false;
        }

        // Нормализация строки
        double divisor = matrix[i][i];
        if (steps) {
            steps->append(QString("Нормировка строки %1 (деление на %2):")
                              .arg(i+1).arg(divisor, 0, 'f', 4));
        }

        for (int j = 0; j < n; ++j) {
            matrix[i][j] /= divisor;
        }
        freeTerms[i] /= divisor;

        if (steps) steps->append(formatMatrix(matrix, freeTerms));

        // Исключение элементов во всех строках
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;

            double factor = matrix[k][i];
            if (std::abs(factor) > 1e-10) {
                if (steps) {
                    steps->append(QString("Исключение в строке %1 (коэффициент %2 * строка %3):")
                                      .arg(k+1).arg(factor, 0, 'f', 4).arg(i+1));
                }

                for (int j = 0; j < n; ++j) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
                freeTerms[k] -= factor * freeTerms[i];

                if (steps) steps->append(formatMatrix(matrix, freeTerms));
            }
        }
    }

    // Результаты
    results = freeTerms;
    if (steps) {
        steps->append("Результат:");
        for (int i = 0; i < n; ++i) {
            steps->append(QString("x%1 = %2").arg(i+1).arg(results[i], 0, 'f', 4));
        }
    }

    return true;
}

double MatrixSolver::determinant(QVector<QVector<double>> matrix) {
    const int n = matrix.size();
    if (n == 0) return 0;

    double det = 1.0;

    for (int i = 0; i < n; ++i) {
        // Поиск опорного элемента
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }

        if (maxRow != i) {
            std::swap(matrix[i], matrix[maxRow]);
            det *= -1; // Смена знака при перестановке строк
        }

        // Проверка на вырожденность
        if (std::abs(matrix[i][i]) < 1e-10) {
            return 0;
        }

        // Исключение элементов
        for (int k = i + 1; k < n; ++k) {
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j < n; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }

        det *= matrix[i][i];
    }

    return det;
}

bool MatrixSolver::solveCramer(const QVector<QVector<double>>& matrix,
                               const QVector<double>& freeTerms,
                               QVector<double>& results,
                               QStringList* steps) {
    const int n = matrix.size();
    if (n == 0 || n != freeTerms.size()) return false;

    if (steps) steps->append("Метод Крамера:");

    // Вычисление определителя основной матрицы
    double mainDet = determinant(matrix);
    if (steps) steps->append(QString("Определитель основной матрицы: %1").arg(mainDet, 0, 'f', 4));

    if (std::abs(mainDet) < 1e-10) {
        if (steps) steps->append("Система вырождена!");
        return false;
    }

    results.resize(n);

    // Создаем копию матрицы для модификаций
    QVector<QVector<double>> modMatrix = matrix;

    for (int i = 0; i < n; ++i) {
        // Заменяем i-й столбец на свободные члены
        for (int j = 0; j < n; ++j) {
            modMatrix[j][i] = freeTerms[j];
        }

        if (steps) {
            steps->append(QString("Матрица для x%1:").arg(i+1));
            steps->append(formatMatrix(modMatrix, freeTerms));
        }

        // Вычисляем определитель модифицированной матрицы
        double modDet = determinant(modMatrix);
        if (steps) steps->append(QString("Определитель: %1").arg(modDet, 0, 'f', 4));

        // Вычисляем x_i
        results[i] = modDet / mainDet;

        if (steps) {
            steps->append(QString("x%1 = %2 / %3 = %4")
                              .arg(i+1)
                              .arg(modDet, 0, 'f', 4)
                              .arg(mainDet, 0, 'f', 4)
                              .arg(results[i], 0, 'f', 4));
        }

        // Восстанавливаем исходную матрицу
        for (int j = 0; j < n; ++j) {
            modMatrix[j][i] = matrix[j][i];
        }
    }

    return true;
}

QVector<QVector<double>> MatrixSolver::inverseMatrix(const QVector<QVector<double>>& inputMatrix) {
    const int n = inputMatrix.size();
    if (n == 0) return {};

    // Создаем копию для работы
    QVector<QVector<double>> matrix = inputMatrix;
    QVector<QVector<double>> inverse(n, QVector<double>(n, 0.0));

    // Создаем расширенную матрицу [A|I]
    for (int i = 0; i < n; ++i) {
        inverse[i][i] = 1.0;
    }

    // Приведение левой части к единичной матрице
    for (int i = 0; i < n; ++i) {
        // Поиск опорного элемента
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // Перестановка строк
        if (maxRow != i) {
            std::swap(matrix[i], matrix[maxRow]);
            std::swap(inverse[i], inverse[maxRow]);
        }

        // Проверка на вырожденность
        if (std::abs(matrix[i][i]) < 1e-10) {
            return {}; // Матрица вырождена
        }

        // Нормализация строки
        double divisor = matrix[i][i];
        for (int j = 0; j < n; ++j) {
            matrix[i][j] /= divisor;
            inverse[i][j] /= divisor;
        }

        // Исключение элементов
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;

            double factor = matrix[k][i];
            for (int j = 0; j < n; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
                inverse[k][j] -= factor * inverse[i][j];
            }
        }
    }

    return inverse;
}

bool MatrixSolver::solveInverseMatrix(const QVector<QVector<double>>& matrix,
                                      const QVector<double>& freeTerms,
                                      QVector<double>& results,
                                      QStringList* steps) {
    const int n = matrix.size();
    if (n == 0 || n != freeTerms.size()) return false;

    if (steps) steps->append("Метод обратной матрицы:");

    // Проверка определителя
    double det = determinant(matrix);
    if (steps) steps->append(QString("Определитель матрицы: %1").arg(det, 0, 'f', 4));

    if (std::abs(det) < 1e-10) {
        if (steps) steps->append("Матрица вырождена, обратной матрицы не существует!");
        return false;
    }

    // Вычисление обратной матрицы
    QVector<QVector<double>> inv = inverseMatrix(matrix);

    if (inv.isEmpty()) {
        if (steps) steps->append("Ошибка при вычислении обратной матрицы");
        return false;
    }

    if (steps) {
        steps->append("Обратная матрица:");
        for (int i = 0; i < n; ++i) {
            QString row;
            for (int j = 0; j < n; ++j) {
                row += QString::number(inv[i][j], 'f', 4) + "\t";
            }
            steps->append(row);
        }
    }

    // Умножение обратной матрицы на вектор свободных членов
    results.resize(n);
    for (int i = 0; i < n; ++i) {
        results[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            results[i] += inv[i][j] * freeTerms[j];
        }

        if (steps) {
            steps->append(QString("x%1 = %2")
                              .arg(i+1)
                              .arg(results[i], 0, 'f', 4));
        }
    }

    return true;
}

QVector<QVector<double>> MatrixSolver::transposeMatrix(const QVector<QVector<double>>& matrix) {
    if (matrix.isEmpty()) return {};

    const int rows = matrix.size();
    const int cols = matrix[0].size();

    QVector<QVector<double>> transposed(cols, QVector<double>(rows));

    for (int i = 0; i < rows; ++i) {
        if (matrix[i].size() != cols) continue; // Пропуск некорректных строк

        for (int j = 0; j < cols; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

QVector<QVector<double>> MatrixSolver::multiplyMatrices(const QVector<QVector<double>>& A,
                                                        const QVector<QVector<double>>& B) {
    if (A.isEmpty() || B.isEmpty()) return {};

    const int rowsA = A.size();
    const int colsA = A[0].size();
    const int rowsB = B.size();
    const int colsB = B[0].size();

    if (colsA != rowsB) return {};

    QVector<QVector<double>> result(rowsA, QVector<double>(colsB, 0.0));

    for (int i = 0; i < rowsA; ++i) {
        if (A[i].size() != colsA) continue; // Пропуск некорректных строк

        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                if (k >= rowsB || j >= B[k].size()) continue; // Проверка границ
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

bool MatrixSolver::solveLeastSquares(const QVector<QVector<double>>& matrix,
                                     const QVector<double>& freeTerms,
                                     QVector<double>& results,
                                     QStringList* steps) {
    const int m = matrix.size(); // Количество уравнений
    if (m == 0) return false;

    const int n = matrix[0].size(); // Количество переменных
    if (n == 0) return false;

    if (steps) steps->append("Метод наименьших квадратов:");

    // Транспонирование матрицы
    QVector<QVector<double>> A_T = transposeMatrix(matrix);
    if (A_T.isEmpty() || A_T.size() != n) {
        if (steps) steps->append("Ошибка транспонирования матрицы");
        return false;
    }

    if (steps) {
        steps->append("Транспонированная матрица A^T:");
        for (int i = 0; i < A_T.size(); ++i) {
            steps->append(formatVector(A_T[i]));
        }
    }

    // Вычисление A^T * A
    QVector<QVector<double>> ATA = multiplyMatrices(A_T, matrix);
    if (ATA.isEmpty() || ATA.size() != n || ATA[0].size() != n) {
        if (steps) steps->append("Ошибка при вычислении A^T * A");
        return false;
    }

    if (steps) {
        steps->append("Матрица A^T * A:");
        for (int i = 0; i < ATA.size(); ++i) {
            steps->append(formatVector(ATA[i]));
        }
    }

    // Вычисление A^T * b
    QVector<double> ATb(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (A_T[i].size() != m) continue;

        for (int j = 0; j < m; ++j) {
            ATb[i] += A_T[i][j] * freeTerms[j];
        }
    }

    if (steps) {
        steps->append("Вектор A^T * b:");
        steps->append(formatVector(ATb));
    }

    // Решаем систему A^T * A * x = A^T * b
    QVector<QVector<double>> ATA_matrix = ATA;
    return solveGauss(ATA_matrix, ATb, results, steps);
}
