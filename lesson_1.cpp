#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <iomanip>
using namespace std;

// Функция для ввода матрицы
void inputMatrix(vector<vector<int>>& matrix, int rows, int cols) {
    // Заполнение матрицы случайными значениями
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = rand() / 10.; // Генерация случайного значения и присвоение его элементу матрицы
        }
    }
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cout << setw(8) << setprecision(5) << matrix[i][j]; // Форматированный вывод элемента матрицы
        }
        cout << endl; 
    }
}

// Функция для нахождения минимума матрицы
int findMatrixMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max(); // Инициализация minVal максимально возможным значением типа int
    for (const auto& row : matrix) { // Проход по каждой строке матрицы
        for (int val : row) { // Проход по каждому элементу строки
            if (val < minVal) { // Если текущий элемент меньше найденного минимума
                minVal = val; // Обновляем значение минимума
            }
        }
    }
    return minVal; // Возвращаем найденное минимальное значение
}

// Функция для нахождения максимума матрицы
int findMatrixMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min(); // Инициализация maxVal минимально возможным значением типа int
    for (const auto& row : matrix) { // Проход по каждой строке матрицы
        for (int val : row) { // Проход по каждому элементу строки
            if (val > maxVal) { // Если текущий элемент больше найденного максимума
                maxVal = val; // Обновляем значение максимума
            }
        }
    }
    return maxVal; // Возвращаем найденное максимальное значение
}

// Функция для нахождения максимума нижнетреугольной части
int findLowerTriangleMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min(); // Инициализация переменной maxVal минимально возможным значением для типа int
    for (int i = 0; i < matrix.size(); ++i) { // Проход по каждой строке матрицы
        for (int j = 0; j <= i; ++j) { // Проход по элементам нижнего треугольника (включая главную диагональ)
            if (matrix[i][j] > maxVal) { // Если текущий элемент больше текущего максимума
                maxVal = matrix[i][j]; // Обновляем значение максимума
            }
        }
    }
    return maxVal; // Возвращаем найденное максимальное значение
}



// Функция для нахождения максимума верхнетреугольной части
int findUpperTriangleMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min(); // Инициализация maxVal минимально возможным значением типа int
    for (int i = 0; i < matrix.size(); ++i) { // Проход по каждой строке матрицы
        for (int j = i; j < matrix.size(); ++j) { // Проход по элементам верхнего треугольника (включая главную диагональ)
            if (matrix[i][j] > maxVal) { // Если текущий элемент больше найденного максимума
                maxVal = matrix[i][j]; // Обновляем значение максимума
            }
        }
    }
    return maxVal; // Возвращаем найденное максимальное значение
}


// Функция для нахождения минимума нижнетреугольной части
int findLowerTriangleMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max(); // Инициализация minVal максимально возможным значением для типа int
    for (int i = 0; i < matrix.size(); ++i) { // Проход по каждой строке матрицы
        for (int j = 0; j <= i; ++j) { // Проход по элементам нижнего треугольника (включая главную диагональ)
            if (matrix[i][j] < minVal) { // Если текущий элемент меньше текущего минимума
                minVal = matrix[i][j]; // Обновляем значение минимума
            }
        }
    }
    return minVal; // Возвращаем найденное минимальное значение
}

// Функция для нахождения минимума верхнетреугольной части
int findUpperTriangleMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max(); // Инициализация minVal максимально возможным значением типа int
    for (int i = 0; i < matrix.size(); ++i) { // Проход по каждой строке матрицы
        for (int j = i; j < matrix.size(); ++j) { // Проход по элементам верхнего треугольника (включая главную диагональ)
            if (matrix[i][j] < minVal) { // Если текущий элемент меньше найденного минимума
                minVal = matrix[i][j]; // Обновляем значение минимума
            }
        }
    }
    return minVal; // Возвращаем найденное минимальное значение
}

// Функция для нахождения минимума главной диагонали
int findMainDiagonalMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max(); // Инициализация minVal максимально возможным значением для типа int
    for (int i = 0; i < matrix.size(); ++i) { // Проход по индексам главной диагонали
        if (matrix[i][i] < minVal) { // Если текущий элемент диагонали меньше найденного минимума
            minVal = matrix[i][i]; // Обновляем значение минимума
        }
    }
    return minVal; // Возвращаем найденное минимальное значение
}

// Функция для нахождения максимума главной диагонали
int findMainDiagonalMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min(); // Инициализация maxVal минимально возможным значением для типа int
    for (int i = 0; i < matrix.size(); ++i) { // Проход по индексам главной диагонали
        if (matrix[i][i] > maxVal) { // Если текущий элемент диагонали больше найденного максимума
            maxVal = matrix[i][i]; // Обновляем значение максимума
        }
    }
    return maxVal; // Возвращаем найденное максимальное значение
}


// Функция для нахождения минимума второстепенной диагонали
int findSecondaryDiagonalMin(const vector<vector<int>>& matrix) {
    int minVal = numeric_limits<int>::max(); // Инициализация minVal максимально возможным значением для типа int
    int n = matrix.size(); // Получаем размерность матрицы
    for (int i = 0; i < n; ++i) { // Проход по индексам второстепенной диагонали
        if (matrix[i][n - 1 - i] < minVal) { // Если текущий элемент диагонали меньше найденного минимума
            minVal = matrix[i][n - 1 - i]; // Обновляем значение минимума
        }
    }
    return minVal; // Возвращаем найденное минимальное значение
}

// Функция для нахождения максимума второстепенной диагонали
int findSecondaryDiagonalMax(const vector<vector<int>>& matrix) {
    int maxVal = numeric_limits<int>::min(); // Инициализация maxVal минимально возможным значением для типа int
    int n = matrix.size(); // Получаем размерность матрицы

    for (int i = 0; i < n; ++i) { // Проход по индексам второстепенной диагонали
        if (matrix[i][n - 1 - i] > maxVal) { // Если текущий элемент диагонали больше найденного максимума
            maxVal = matrix[i][n - 1 - i]; // Обновляем значение максимума
        }
    }
    return maxVal; // Возвращаем найденное максимальное значение
}

// Функция для нахождения среднеарифметического значения элементов матрицы
double findMatrixAverage(const vector<vector<int>>& matrix) {
    double sum = 0.0; // Инициализация суммы элементов матрицы
    int count = 0; // Счетчик количества элементов
    for (const auto& row : matrix) { // Проход по каждой строке матрицы
        for (int val : row) { // Проход по каждому элементу строки
            sum += val; // Добавляем значение к сумме
            count++; // Увеличиваем счетчик элементов
        }
    }
    return sum / count; // Возвращаем среднее арифметическое (сумма деленная на количество)
}

// Функция для нахождения среднеарифметического значения элементов нижнетреугольной части
double findLowerTriangleAverage(const vector<vector<int>>& matrix) {
    double sum = 0.0; // Инициализация суммы элементов нижнего треугольника
    int count = 0; // Счетчик количества элементов нижнего треугольника
    for (int i = 0; i < matrix.size(); ++i) { // Проход по строкам матрицы
        for (int j = 0; j <= i; ++j) { // Проход по элементам нижнего треугольника (включая главную диагональ)
            sum += matrix[i][j]; // Добавляем значение к сумме
            count++; // Увеличиваем счетчик элементов
        }
    }
    return sum / count; // Возвращаем среднее арифметическое для нижнего треугольника
}

// Функция для нахождения среднеарифметического значения элементов верхнетреугольной части
double findUpperTriangleAverage(const vector<vector<int>>& matrix) {
    double sum = 0.0; // Инициализация суммы элементов верхнего треугольника
    int count = 0; // Счетчик количества элементов верхнего треугольника
    for (int i = 0; i < matrix.size(); ++i) { // Проход по строкам матрицы
        for (int j = i; j < matrix.size(); ++j) { // Проход по элементам верхнего треугольника (включая главную диагональ)
            sum += matrix[i][j]; // Добавляем значение к сумме
            count++; // Увеличиваем счетчик элементов
        }
    }
    return sum / count; // Возвращаем среднее арифметическое для верхнего треугольника
}

// Функция для нахождения суммы строк матрицы
vector<int> findRowSums(const vector<vector<int>>& matrix) {
    vector<int> rowSums(matrix.size(), 0); // Инициализация вектора для хранения сумм строк с начальным значением 0
    for (size_t i = 0; i < matrix.size(); ++i) { // Проход по строкам матрицы
        for (size_t j = 0; j < matrix[i].size(); ++j) { // Проход по элементам каждой строки
            rowSums[i] += matrix[i][j]; // Суммируем элементы строки и сохраняем в соответствующем индексе вектора rowSums
        }
    }
    return rowSums; // Возвращаем вектор с суммами строк
}

// Функция для нахождения суммы столбцов матрицы
vector<int> findColumnSums(const vector<vector<int>>& matrix) {
    vector<int> colSums(matrix[0].size(), 0); // Инициализация вектора для хранения сумм столбцов с начальным значением 0
    for (size_t j = 0; j < matrix[0].size(); ++j) { // Проход по столбцам матрицы
        for (size_t i = 0; i < matrix.size(); ++i) { // Проход по элементам каждого столбца
            colSums[j] += matrix[i][j]; // Суммируем элементы столбца и сохраняем в соответствующем индексе вектора colSums
        }
    }
    return colSums; // Возвращаем вектор с суммами столбцов
}

// Функция для нахождения минимальных значений строк
vector<int> findRowMins(const vector<vector<int>>& matrix) {
    vector<int> rowMins(matrix.size(), numeric_limits<int>::max()); // Инициализация вектора для хранения минимальных значений строк с максимально возможным значением
    for (size_t i = 0; i < matrix.size(); ++i) { // Проход по строкам матрицы
        for (size_t j = 0; j < matrix[i].size(); ++j) { // Проход по элементам каждой строки
            if (matrix[i][j] < rowMins[i]) { // Если текущий элемент меньше найденного минимума в данной строке
                rowMins[i] = matrix[i][j]; // Обновляем значение минимума для этой строки
            }
        }
    }
    return rowMins; // Возвращаем вектор с минимальными значениями строк
}

// Функция для нахождения минимальных значений столбцов
vector<int> findColumnMins(const vector<vector<int>>& matrix) {
    // Инициализация вектора для хранения минимальных значений столбцов
    // Размер вектора равен количеству столбцов, начальное значение - максимально возможное
    vector<int> colMins(matrix[0].size(), numeric_limits<int>::max());

    // Проходим по каждому столбцу
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        // Проходим по каждой строке для текущего столбца
        for (size_t i = 0; i < matrix.size(); ++i) {
            // Если текущий элемент меньше найденного минимума для этого столбца
            if (matrix[i][j] < colMins[j]) {
                colMins[j] = matrix[i][j]; // Обновляем минимум
            }
        }
    }

    return colMins; // Возвращаем вектор с минимальными значениями столбцов
}

// Функция для нахождения максимальных значений строк
vector<int> findRowMaxs(const vector<vector<int>>& matrix) {
    // Инициализация вектора для хранения максимальных значений строк
    // Размер вектора равен количеству строк, начальное значение - минимально возможное
    vector<int> rowMaxs(matrix.size(), numeric_limits<int>::min());

    // Проходим по каждой строке матрицы
    for (size_t i = 0; i < matrix.size(); ++i) {
        // Проходим по каждому элементу строки
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            // Если текущий элемент больше найденного максимума для этой строки
            if (matrix[i][j] > rowMaxs[i]) {
                rowMaxs[i] = matrix[i][j]; // Обновляем максимум
            }
        }
    }

    return rowMaxs; // Возвращаем вектор с максимальными значениями строк
}

// Функция для нахождения максимальных значений столбцов
vector<int> findColumnMaxs(const vector<vector<int>>& matrix) {
    // Инициализация вектора для хранения максимальных значений столбцов
    // Размер вектора равен количеству столбцов, начальное значение - минимально возможное
    vector<int> colMaxs(matrix[0].size(), numeric_limits<int>::min());

    // Проходим по каждому столбцу
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        // Проходим по каждой строке для текущего столбца
        for (size_t i = 0; i < matrix.size(); ++i) {
            // Если текущий элемент больше найденного максимума для этого столбца
            if (matrix[i][j] > colMaxs[j]) {
                colMaxs[j] = matrix[i][j]; // Обновляем максимум
            }
        }
    }
    return colMaxs; // Возвращаем вектор с максимальными значениями столбцов
}

// Функция для нахождения среднеарифметического значения строк
vector<double> findRowAverages(const vector<vector<int>>& matrix) {
    // Инициализация вектора для хранения среднеарифметических значений строк
    vector<double> rowAverages(matrix.size(), 0.0);

    // Проходим по каждой строке матрицы
    for (size_t i = 0; i < matrix.size(); ++i) {
        double sum = 0.0; // Инициализация суммы элементов строки

        // Проходим по каждому элементу строки
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            sum += matrix[i][j]; // Суммируем элементы строки
        }

        // Вычисляем среднее арифметическое и сохраняем его в соответствующем индексе вектора
        rowAverages[i] = sum / matrix[i].size();
    }
    return rowAverages; // Возвращаем вектор со среднеарифметическими значениями строк
}


// Функция для нахождения среднеарифметического значения столбцов
vector<double> findColumnAverages(const vector<vector<int>>& matrix) {
    // Инициализация вектора для хранения среднеарифметических значений столбцов
    vector<double> colAverages(matrix[0].size(), 0.0);

    // Проходим по каждому столбцу матрицы
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        double sum = 0.0; // Инициализация суммы элементов столбца

        // Проходим по каждой строке для текущего столбца
        for (size_t i = 0; i < matrix.size(); ++i) {
            sum += matrix[i][j]; // Суммируем элементы столбца
        }

        // Вычисляем среднее арифметическое и сохраняем его в соответствующем индексе вектора
        colAverages[j] = sum / matrix.size();
    }

    return colAverages; // Возвращаем вектор со среднеарифметическими значениями столбцов
}


// Функция для нахождения суммы нижней и верхней треугольной части матрицы
pair<double, double> findTriangleSums(const vector<vector<int>>& matrix) {
    double lowerSum = 0.0, upperSum = 0.0; // Инициализация сумм для нижней и верхней треугольной части

    // Проходим по всем строкам матрицы
    for (int i = 0; i < matrix.size(); ++i) {
        // Суммируем элементы нижней треугольной части (включая главную диагональ)
        for (int j = 0; j <= i; ++j) {
            lowerSum += matrix[i][j];
        }
        // Суммируем элементы верхней треугольной части (включая главную диагональ)
        for (int j = i; j < matrix.size(); ++j) {
            upperSum += matrix[i][j];
        }
    }

    return { lowerSum, upperSum }; // Возвращаем пару с суммами нижней и верхней треугольной части
}


// Функция для нахождения элемента, наиболее близкого по значению к среднеарифметическому
int findClosestToAverage(const vector<vector<int>>& matrix, double average) {
    int closestValue = numeric_limits<int>::max(); // Инициализация переменной для хранения ближайшего значения

    double minDiff = numeric_limits<double>::max(); // Инициализация минимальной разницы

    // Проходим по всем элементам матрицы
    for (const auto& row : matrix) {
        for (int val : row) {
            double diff = abs(val - average); // Вычисляем разницу между элементом и средним

            // Если разница меньше найденной минимальной разницы или равна ей и элемент меньше по модулю, обновляем ближайшее значение
            if (diff < minDiff ||
                (diff == minDiff && abs(val - closestValue) > diff)) {
                closestValue = val;
                minDiff = diff;
            }
        }
    }

    return closestValue; // Возвращаем ближайшее значение к среднему арифметическому
}


int main() {
    setlocale(LC_CTYPE, "Russian"); // Для отображения русских символов в консоли.

    const int rows = 5, cols = 5; // Определяем константы для количества строк и столбцов матрицы.
    vector<vector<int>> matrix(rows, vector<int>(cols)); // Создаем двумерный матрицу размером 5x5, инициализируя все элементы нулями.

    inputMatrix(matrix, rows, cols); // Вызываем функцию для ввода значений в матрицу от пользователя.
    
    // Выводим
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

    double averageValue = findMatrixAverage(matrix); // Вычисляем среднеарифметическое значение всех элементов матрицы.
    cout << "Среднеарифметическое значение элементов матрицы: " << averageValue << endl; 
    cout << "Среднеарифметическое значение элементов нижнетреугольной части: "
        << findLowerTriangleAverage(matrix) << endl;
    cout << "Среднеарифметическое значение элементов верхнетреугольной части: "
        << findUpperTriangleAverage(matrix) << endl;

    auto rowSums = findRowSums(matrix); // Находим суммы всех строк матрицы и сохраняем в вектор rowSums.
    cout << "Суммы строк матрицы: ";

    for (int sum : rowSums) // Проходим по всем суммам строк.
        cout << sum << " ";
    cout << endl; 

    auto colSums = findColumnSums(matrix); // Находим суммы всех столбцов матрицы и сохраняем в вектор colSums.
    cout << "Суммы столбцов матрицы: "; 
    for (int sum : colSums) // Проходим по всем суммам столбцов.
        cout << sum << " "; 
    cout << endl; 

    auto rowMins = findRowMins(matrix); // Находим минимальные значения для каждой строки и сохраняем в вектор rowMins.
    cout << "Минимальные значения строк: "; 

    for (int min : rowMins) // Проходим по всем минимальным значениям строк.
        cout << min << " "; 
    cout << endl; 

    auto colMins = findColumnMins(matrix); // Находим минимальные значения для каждого столбца и сохраняем в вектор colMins.
    cout << "Минимальные значения столбцов: "; // Начинаем вывод минимальных значений столбцов.
    for (int min : colMins) // Проходим по всем минимальным значениям столбцов.
        cout << min << " "; // Выводим каждое минимальное значение через пробел.
    cout << endl; 

    auto rowMaxs = findRowMaxs(matrix); // Находим максимальные значения для каждой строки и сохраняем в вектор rowMaxs.
    cout << "Максимальные значения строк: "; 
    for (int max : rowMaxs) // Проходим по всем максимальным значениям строк.
        cout << max << " "; 
    cout << endl; 

    auto colMaxs = findColumnMaxs(matrix); // Находим максимальные значения для каждого столбца и сохраняем в вектор colMaxs.
    cout << "Максимальные значения столбцов: "; 
    for (int max : colMaxs) // Проходим по всем максимальным значениям столбцов.
        cout << max << " "; 
    cout << endl; 

    auto rowAverages = findRowAverages(matrix); // Находим среднеарифметические значения для каждой строки и сохраняем в вектор rowAverages.
    cout << "Среднеарифметическое значение строк: "; 
    for (double avg : rowAverages) // Проходим по всем среднеарифметическим значениям строк.
        cout << avg << " "; 
    cout << endl; 

    auto colAverages = findColumnAverages(matrix); // Находим среднеарифметические значения для каждого столбца и сохраняем в вектор colAverages.
    cout << "Среднеарифметическое значение столбцов: "; 
    for (double avg : colAverages) // Проходим по всем среднеарифметическим значениям столбцов.
        cout << avg << " "; 
    cout << endl; 

    auto triangleSums = findTriangleSums(matrix); // Находим суммы нижней и верхней треугольных частей матрицы и сохраняем их в пару triangleSums.
    cout << "Сумма нижней треугольной части: " << triangleSums.first
        << ", Сумма верхней треугольной части: " << triangleSums.second
        << endl; 

    int closestValueToAvg = findClosestToAverage(matrix, averageValue); // Находим элемент, наиболее близкий к среднеарифметическому значению всей матрицы.
    cout << "Элемент, наиболее близкий по значению к среднеарифметическому: "
        << closestValueToAvg << endl; 
    return 0; // Завершаем выполнение программы, возвращая 0.
}
