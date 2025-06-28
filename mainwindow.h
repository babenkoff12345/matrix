#pragma once

#include <QMainWindow>
#include <QTableWidget>
#include <QVector>
#include "MatrixSolver.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void resizeMatrix();
    void solveSystem();

private:
    Ui::MainWindow *ui;
    void setupMatrixTable(int rows, int cols);
    bool getMatrixData(QVector<QVector<double>>& matrix, QVector<double>& freeTerms);
};
