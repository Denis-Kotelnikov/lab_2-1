#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <iomanip>

using namespace std;

// Функция для ввода матрицы
void inputMatrix(vector<vector<int>>& matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) matrix[i][j] = rand() / 10.;
    for (int i = 0; i < rows; ++i)
    {
        for (int i = 0; i < rows; ++i)
                cout << setw(8) << setprecision(5) << rand()/10;
            cout << endl;
        }
    }
}

// Функция для нахождения минимума матрицы
int findMatrixMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max();
    for (const auto& row : matrix) {
        for (int val : row) {
            if (val < minVal) {
                minVal = val;
            }
        }
    }
    return minVal;
}

// Функция для нахождения максимума матрицы
int findMatrixMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min();
    for (const auto& row : matrix) {
        for (int val : row) {
            if (val > maxVal) {
                maxVal = val;
            }
        }
    }
    return maxVal;
}

// Функция для нахождения максимума нижнетреугольной части
int findLowerTriangleMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min();
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j <= i; ++j) {
            if (matrix[i][j] > maxVal) {
                maxVal = matrix[i][j];
            }
        }
    }
    return maxVal;
}

// Функция для нахождения максимума верхнетреугольной части
int findUpperTriangleMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min();
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = i; j < matrix.size(); ++j) {
            if (matrix[i][j] > maxVal) {
                maxVal = matrix[i][j];
            }
        }
    }
    return maxVal;
}

// Функция для нахождения минимума нижнетреугольной части
int findLowerTriangleMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max();
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j <= i; ++j) {
            if (matrix[i][j] < minVal) {
                minVal = matrix[i][j];
            }
        }
    }
    return minVal;
}

// Функция для нахождения минимума верхнетреугольной части
int findUpperTriangleMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max();
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = i; j < matrix.size(); ++j) {
            if (matrix[i][j] < minVal) {
                minVal = matrix[i][j];
            }
        }
    }
    return minVal;
}

// Функция для нахождения минимума главной диагонали
int findMainDiagonalMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max();
    for (int i = 0; i < matrix.size(); ++i) {
        if (matrix[i][i] < minVal) {
            minVal = matrix[i][i];
        }
    }
    return minVal;
}

// Функция для нахождения максимума главной диагонали
int findMainDiagonalMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min();
    for (int i = 0; i < matrix.size(); ++i) {
        if (matrix[i][i] > maxVal) {
            maxVal = matrix[i][i];
        }
    }
    return maxVal;
}

// Функция для нахождения минимума второстепенной диагонали
int findSecondaryDiagonalMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max();
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        if (matrix[i][n - 1 - i] < minVal) {
            minVal = matrix[i][n - 1 - i];
        }
    }
    return minVal;
}

// Функция для нахождения максимума второстепенной диагонали
int findSecondaryDiagonalMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min();
    int n = matrix.size();

        for (int i = 0; i < n; ++i) {
            if (matrix[i][n - 1 - i] > maxVal) {
                maxVal = matrix[i][n - 1 - i];
            }
        }
    return maxVal;
}

// Функция для нахождения среднеарифметического значения элементов матрицы
double findMatrixAverage(const vector<vector<int>>& matrix) {
    double sum = 0.0;
    int count = 0;
    for (const auto& row : matrix) {
        for (int val : row) {
            sum += val;
            count++;
        }
    }
    return sum / count;
}

// Функция для нахождения среднеарифметического значения элементов нижнетреугольной части
double findLowerTriangleAverage(const vector<vector<int>>& matrix) {
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j <= i; ++j) {
            sum += matrix[i][j];
            count++;
        }
    }
    return sum / count;
}

// Функция для нахождения среднеарифметического значения элементов верхнетреугольной части
double findUpperTriangleAverage(const vector<vector<int>>& matrix) {
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = i; j < matrix.size(); ++j) {
            sum += matrix[i][j];
            count++;
        }
    }
    return sum / count;
}

// Функция для нахождения суммы строк матрицы
vector<int> findRowSums(const vector<vector<int>>& matrix) {
    vector<int> rowSums(matrix.size(), 0);
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            rowSums[i] += matrix[i][j];
        }
    }
    return rowSums;
}

// Функция для нахождения суммы столбцов матрицы
vector<int> findColumnSums(const vector<vector<int>>& matrix) {
    vector<int> colSums(matrix[0].size(), 0);
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            colSums[j] += matrix[i][j];
        }
    }
    return colSums;
}

// Функция для нахождения минимальных значений строк
vector<int> findRowMins(const vector<vector<int>>& matrix) {
    vector<int> rowMins(matrix.size(), numeric_limits<int>::max());
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            if (matrix[i][j] < rowMins[i]) {
                rowMins[i] = matrix[i][j];
            }
        }
    }
    return rowMins;
}

// Функция для нахождения минимальных значений столбцов
vector<int> findColumnMins(const vector<vector<int>>& matrix) {
    vector<int> colMins(matrix[0].size(), numeric_limits<int>::max());
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            if (matrix[i][j] < colMins[j]) {
                colMins[j] = matrix[i][j];
            }
        }
    }
    return colMins;
}

// Функция для нахождения максимальных значений строк
vector<int> findRowMaxs(const vector<vector<int>>& matrix) {
    vector<int> rowMaxs(matrix.size(), numeric_limits<int>::min());
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            if (matrix[i][j] > rowMaxs[i]) {
                rowMaxs[i] = matrix[i][j];
            }
        }
    }
    return rowMaxs;
}

// Функция для нахождения максимальных значений столбцов
vector<int> findColumnMaxs(const vector<vector<int>>& matrix) {
    vector<int> colMaxs(matrix[0].size(), numeric_limits<int>::min());
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            if (matrix[i][j] > colMaxs[j]) {
                colMaxs[j] = matrix[i][j];
            }
        }
    }
    return colMaxs;
}

// Функция для нахождения среднеарифметического значения строк
vector<double> findRowAverages(const vector<vector<int>>& matrix) {
    vector<double> rowAverages(matrix.size(), 0.0);

    for (size_t i = 0; i < matrix.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            sum += matrix[i][j];
        }
        rowAverages[i] = sum / matrix[i].size();

    }

    return rowAverages;
}

// Функция для нахождения среднеарифметического значения столбцов
vector<double> findColumnAverages(const vector<vector<int>>& matrix) {
    vector<double> colAverages(matrix[0].size(), 0.0);

    for (size_t j = 0; j < matrix[0].size(); ++j) {
        double sum = 0.0;
        for (size_t i = 0; i < matrix.size(); ++i) {
            sum += matrix[i][j];
        }
        colAverages[j] = sum / matrix.size();
    }

    return colAverages;
}

// Функция для нахождения суммы нижней и верхней треугольной части матрицы
pair<double, double> findTriangleSums(const vector<vector<int>>& matrix) {
    double lowerSum = 0.0, upperSum = 0.0;

    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j <= i; ++j) { // Нижняя треугольная часть
            lowerSum += matrix[i][j];
        }
        for (int j = i; j < matrix.size(); ++j) { // Верхняя треугольная часть
            upperSum += matrix[i][j];
        }
    }

    return { lowerSum, upperSum };
}

// Функция для нахождения элемента, наиболее близкого по значению к среднеарифметическому
int findClosestToAverage(const vector<vector<int>>& matrix, double average) {
    int closestValue = numeric_limits<int>::max();

    double minDiff = numeric_limits<double>::max();

    for (const auto& row : matrix) {
        for (int val : row) {
            double diff = abs(val - average);
            if (diff < minDiff ||
                (diff == minDiff && abs(val - closestValue) > diff)) { // Если разница такая же, выбираем меньшее по модулю значение
                closestValue = val;
                minDiff = diff;
            }
        }
    }

    return closestValue;
}

int main() {
    setlocale(LC_CTYPE, "Russian");
    const int rows = 5, cols = 5;
    vector<vector<int>> matrix(rows, vector<int>(cols));

    inputMatrix(matrix, rows, cols);

    cout << "Минимум матрицы: " << findMatrixMin(matrix) << endl;
    cout << "Максимум матрицы: " << findMatrixMax(matrix) << endl;
    cout << "Максимум нижнетреугольной части: " << findLowerTriangleMax(matrix) << endl;
    cout << "Максимум верхнетреугольной части: " << findUpperTriangleMax(matrix) << endl;
    cout << "Минимум нижнетреугольной части: " << findLowerTriangleMin(matrix) << endl;
    cout << "Минимум верхнетреугольной части: " << findUpperTriangleMin(matrix) << endl;
    cout << "Минимум главной диагонали: " << findMainDiagonalMin(matrix) << endl;
    cout << "Максимум главной диагонали: " << findMainDiagonalMax(matrix) << endl;
    cout << "Минимум второстепенной диагонали: " << findSecondaryDiagonalMin(matrix) << endl;
    cout << "Максимум второстепенной диагонали: " << findSecondaryDiagonalMax(matrix) << endl;

    double averageValue = findMatrixAverage(matrix);
    cout << "Среднеарифметическое значение элементов матрицы: " << averageValue << endl;

    cout << "Среднеарифметическое значение элементов нижнетреугольной части: "
        << findLowerTriangleAverage(matrix) << endl;

    cout << "Среднеарифметическое значение элементов верхнетреугольной части: "
        << findUpperTriangleAverage(matrix) << endl;

    auto rowSums = findRowSums(matrix);
    cout << "Суммы строк матрицы: ";
    for (int sum : rowSums)
        cout << sum << " ";
    cout << endl;

    auto colSums = findColumnSums(matrix);
    cout << "Суммы столбцов матрицы: ";
    for (int sum : colSums)
        cout << sum << " ";
    cout << endl;

    auto rowMins = findRowMins(matrix);
    cout << "Минимальные значения строк: ";
    for (int min : rowMins)
        cout << min << " ";
    cout << endl;

    auto colMins = findColumnMins(matrix);
    cout << "Минимальные значения столбцов: ";
    for (int min : colMins)
        cout << min << " ";
    cout << endl;

    auto rowMaxs = findRowMaxs(matrix);
    cout << "Максимальные значения строк: ";
    for (int max : rowMaxs)
        cout << max << " ";
    cout << endl;

    auto colMaxs = findColumnMaxs(matrix);
    cout << "Максимальные значения столбцов: ";
    for (int max : colMaxs)
        cout << max << " ";
    cout << endl;

    auto rowAverages = findRowAverages(matrix);

        cout << "Среднеарифметическое значение строк: ";
    for (double avg : rowAverages)
        cout << avg << " ";
    cout << endl;

    auto colAverages = findColumnAverages(matrix);
    cout << "Среднеарифметическое значение столбцов: ";
    for (double avg : colAverages)
        cout << avg << " ";
    cout << endl;

    auto triangleSums = findTriangleSums(matrix);
    cout << "Сумма нижней треугольной части: " << triangleSums.first
        << ", Сумма верхней треугольной части: " << triangleSums.second
        << endl;

    int closestValueToAvg = findClosestToAverage(matrix, averageValue);
    cout << "Элемент, наиболее близкий по значению к среднеарифметическому: "
        << closestValueToAvg << endl;

    return 0;
}