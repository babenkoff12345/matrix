#include "MainWindow.h"
#include "MatrixSolver.h"
#include "ui_MainWindow.h"
#include <QHeaderView>
#include <QMessageBox>
#include <iomanip>
#include <sstream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    // Настройка таблицы
    ui->matrixTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->matrixTable->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    // Инициализация матрицы
    setupMatrixTable(3, 3);

    // Заполнение комбобокса методами
    ui->methodComboBox->addItem("Метод Гаусса", MatrixSolver::Gauss);
    ui->methodComboBox->addItem("Метод Жордана-Гаусса", MatrixSolver::JordanGauss);
    ui->methodComboBox->addItem("Метод Крамера", MatrixSolver::Cramer);
    ui->methodComboBox->addItem("Метод обратной матрицы", MatrixSolver::InverseMatrix);
    ui->methodComboBox->addItem("Метод наименьших квадратов", MatrixSolver::LeastSquares);

    // Соединение сигналов и слотов
    connect(ui->resizeButton, &QPushButton::clicked, this, &MainWindow::resizeMatrix);
    connect(ui->solveButton, &QPushButton::clicked, this, &MainWindow::solveSystem);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::resizeMatrix() {
    int rows = ui->equationsSpinBox->value();
    int cols = ui->variablesSpinBox->value();
    setupMatrixTable(rows, cols);
}

void MainWindow::setupMatrixTable(int rows, int cols) {
    ui->matrixTable->setRowCount(rows);
    ui->matrixTable->setColumnCount(cols + 1); // +1 для свободных членов

    // Установка заголовков
    QStringList headers;
    for (int i = 0; i < cols; ++i) {
        headers << QString("x%1").arg(i + 1);
    }
    headers << "=";
    ui->matrixTable->setHorizontalHeaderLabels(headers);

    // Заполнение нулями
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j <= cols; ++j) {
            QTableWidgetItem *item = new QTableWidgetItem("0");
            ui->matrixTable->setItem(i, j, item);
        }
    }
}

bool MainWindow::getMatrixData(
    QVector<QVector<double>>& matrix,
    QVector<double>& freeTerms)
{
    const int rows = ui->matrixTable->rowCount();
    const int cols = ui->matrixTable->columnCount() - 1; // Последний столбец - свободные члены

    matrix.resize(rows);
    freeTerms.resize(rows);

    for (int i = 0; i < rows; ++i) {
        matrix[i].resize(cols);
        for (int j = 0; j < cols; ++j) {
            QTableWidgetItem *item = ui->matrixTable->item(i, j);
            if (!item) return false;
            bool ok;
            double value = item->text().toDouble(&ok);
            if (!ok) return false;
            matrix[i][j] = value;
        }

        // Свободный член
        QTableWidgetItem *freeItem = ui->matrixTable->item(i, cols);
        if (!freeItem) return false;
        bool ok;
        double freeValue = freeItem->text().toDouble(&ok);
        if (!ok) return false;
        freeTerms[i] = freeValue;
    }
    return true;
}

void MainWindow::solveSystem() {
    QVector<QVector<double>> matrix;
    QVector<double> freeTerms;

    if (!getMatrixData(matrix, freeTerms)) {
        QMessageBox::critical(this, "Ошибка", "Недопустимые входные данные");
        return;
    }

    QVector<double> results;
    QStringList solutionSteps;

    // Получаем выбранный метод
    MatrixSolver::Method method = static_cast<MatrixSolver::Method>(
        ui->methodComboBox->currentData().toInt()
        );

    if (MatrixSolver::solveSystem(method, matrix, freeTerms, results, &solutionSteps)) {
        QString html = "<h2>Решение:</h2><table>";
        for (int i = 0; i < results.size(); ++i) {
            html += QString("<tr><td>x<sub>%1</sub> = %2</td></tr>")
            .arg(i+1)
                .arg(results[i], 0, 'f', 4);
        }
        html += "</table>";

        // Добавляем шаги решения
        html += "<h2>Ход решения:</h2><pre>" + solutionSteps.join("\n") + "</pre>";
        ui->resultsText->setHtml(html);
    } else {
        QString html = "<b>Решение не найдено!</b>";
        html += "<h2>Ход решения:</h2><pre>" + solutionSteps.join("\n") + "</pre>";
        ui->resultsText->setHtml(html);
    }
}
