#pragma once
#include <QVector>
#include <QStringList>
#include <cmath>
#include <algorithm>

class MatrixSolver {
public:
    enum Method {
        Gauss,
        JordanGauss,
        Cramer,
        InverseMatrix,
        LeastSquares
    };

    static bool solveSystem(Method method,
                            QVector<QVector<double>>& matrix,
                            QVector<double>& freeTerms,
                            QVector<double>& results,
                            QStringList* steps = nullptr);

private:
    static bool solveGauss(QVector<QVector<double>>& matrix,
                           QVector<double>& freeTerms,
                           QVector<double>& results,
                           QStringList* steps = nullptr);

    static bool solveJordanGauss(QVector<QVector<double>>& matrix,
                                 QVector<double>& freeTerms,
                                 QVector<double>& results,
                                 QStringList* steps = nullptr);

    static bool solveCramer(QVector<QVector<double>>& matrix,
                            QVector<double>& freeTerms,
                            QVector<double>& results,
                            QStringList* steps = nullptr);

    static bool solveInverseMatrix(QVector<QVector<double>>& matrix,
                                   QVector<double>& freeTerms,
                                   QVector<double>& results,
                                   QStringList* steps = nullptr);

    static bool solveLeastSquares(QVector<QVector<double>>& matrix,
                                  QVector<double>& freeTerms,
                                  QVector<double>& results,
                                  QStringList* steps = nullptr);

    static double determinant(QVector<QVector<double>> matrix);
    static QVector<QVector<double>> inverseMatrix(QVector<QVector<double>> matrix);
    static QVector<QVector<double>> transposeMatrix(const QVector<QVector<double>>& matrix);
    static QVector<QVector<double>> multiplyMatrices(const QVector<QVector<double>>& A,
                                                     const QVector<QVector<double>>& B);
    static QString formatMatrix(const QVector<QVector<double>>& matrix,
                                const QVector<double>& freeTerms);
    static QString formatVector(const QVector<double>& vec);
};
