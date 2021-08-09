#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define false 0
#define true 1
using namespace std;

int scalar(int n, double *a, double *b) // Скалярное произведение
{
    int scalar = 0;
    for(int i = 0; i < n; i++)
    {
        scalar += a[i]*b[i];
    }
    return scalar;
}

void left_product(int n, double **a, double **b) // Произведение матрицы a на матрицу b слева
{
    double **result = (double**)malloc(n*sizeof(double*));
    for(int i = 0; i < n; i++)
    {
        result[i] = (double*)malloc(n*sizeof(double));
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < n; k++)
            {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            b[i][j] = result[i][j];
        }
    }
}

void minus2vec(int n, int gamma, double *v, double **T) // Вычитаем из матрицы T произведение векторов v и v транспонированного, делённое на гамму
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            T[i][j] -= v[i]*v[j]/gamma;
        }
    }
}

void outL(int n, double **l, double *y) // Решаем нижнетреугольную СЛАУ
{
    for(int i=0; i < n; i++)
    {
        y[i] = l[i][n]/l[i][i];
        for(int j = 0; j < i; j++)
        {
            y[i] = y[i] - ((l[i][j] * y[j]) / l[i][i]);
        }
    }
    for(int i=0; i<n; i++)
        {
            if(y[i] == y[i])
            {
                cout << "y" << i+1 << "=" << y[i] << endl;
            }
        }
    return;
}

void outR(int n, double **a) // Решаем верхнетреугольную СЛАУ
{
    double*x = (double*)malloc(n*sizeof(double));
    if (x == NULL)
    {
        return;
    }
    for(int i=n-1; i>=0; i--) //с конца
    {
        x[i]=a[i][n]/a[i][i];
        for(int j=n-1; j>i; j--)
        {
            x[i]=x[i] - a[i][j]*(x[j]/a[i][i]);
        }
    }
    for(int i=0; i<n; i++)
        {
            if((x[i] == x[i]))
            {
                cout << "x" << i+1 << "=" << x[i] << endl;
            }
        }
}

void output(int n, double **a, bool is_square) // Выводим матрицу
{
    cout << endl;
    for (int i=0; i<n; i++)  //переходим на следующую строку, i - номер строки
    {
        for (int j=0; j<n+(!is_square); j++)  //выводим строку, j - номер столбца - элемента в строке
        {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

bool symmetrical(int n, double **a) // Проверка на симметрию
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if (a[i][j] != a[j][i])
            {
                return 0;
            }
        }
    }
    return 1;
}

void check_rk(int n, double **a) // По т. Кронекера-Капелли проверяем ранг верхнетреугольной матрицы и выводим решение
{
    int rank1=0, rank2=0;
    double temp;
    for(int i=0; i<n; i++)
    {
        temp=0;
        for(int j=0; j < n; j++)
        {
            temp += a[i][j]*a[i][j];
        }
        if (temp != 0) // ранг треугольной матрицы равен количеству её ненулевых строк
        {
            rank1++;
        }
    }
    for(int i=0; i<n; i++)
    {
        temp=0;
        for(int j=0; j<n+1; j++)
        {
            temp += a[i][j]*a[i][j];
        }
        if (temp != 0)
        {
            rank2++;
        }
    }

    cout << "Rang(A)=" << rank1 << ", Rang(A|B)=" << rank2 << ", N = " << n << endl;
    if ((rank1==rank2) && (rank1==n))  // по теореме Кронекера-Капелли СЛАУ имеет решение тогда и только тогда когда ранг матрицы коэффициентов равен рангу расширенной матрицы системы
    {
        outR(n, a); // если есть решение - выводим
    }
    if (rank2 >= rank1 + 1)
    {
        cout << "System has no solution (Kronecker-Capelli th.)" << endl;
        return;
    }
    if ((rank1==rank2) && (rank1<n))
    {
        cout << "System has infinite number of solution (Kronecker-Capelli th.)" << endl;
    }
}

void Elimination(int n, double **a)   // Приводим к верхнетреугольному виду
{
    int p = 0; // p - номер опорного элемента, т.е. первого ненулевого в строке
    for (int q=0; q<n; q++)  // q - номер строки, которую мы собираемся вычитать из всех под ней
    {
        for (int i=1; q+i<n; i++) // i - номер строки под q-той
        {
            /* если опорный ненулевой, то всё хорошо и мы ничего не делаем */
            if(a[q][p] == 0) // если нулевой, то всё не очень хорошо и мы меняем его
            {
                cout << "Leading coefficient is zero. We need to change it:" << endl;
                int temp_sum = 0;
                for(int temp_col = q; temp_col < n; temp_col++) // считаем сумму квадратов чисел в столбце
                {
                    temp_sum += a[temp_col][p]*a[temp_col][p];
                }
                cout << "Sum of squares of elements in column = " << temp_sum << endl;

                if (temp_sum != 0) // если она не нулевая, то есть ненулевой элемент
                {
                    int temp_col = q;
                    while (a[temp_col][p] == 0) // нам надо его найти и мы его находим посредством простого пробегания по циклу
                    {
                        temp_col++;
                    }
                    cout << "Index of non-zero element in column is " << temp_col+1 << endl; // на всякий случай выводим его индекс в столбце
                    for (int j=p; j<n+1; j++) // меняем нашу строку с нулевым опорным на строку с ненулевым элементом из того же столбца
                    {
                        int c=a[q][j];
                        a[q][j]=a[temp_col][j];
                        a[temp_col][j]=c;
                    }
                    cout << "L" << q+1 << "<->" << "L" << temp_col+1 << endl;
                    output(n, a, false);
                }
                else // иначе столбец оказывается пустым и мы как будто его игнорируем, переходя на следующий и оставаясь в той же строке
                {
                    p++;
                    cout << "Changing the leading coefficient" << endl;
                }
            }
            cout << "Leading coefficient is a[q][p] = " << a[q][p] << ", q = " << q+1 << ", p = " << p+1 << endl;

            /* вычитаем строчку с ненулевым опорным из всех строчек ниже */
            double e, f;
            e=a[q+i][p];
            f=a[q][p];
            for (int j=0; j<n+1; j++)
            {
                a[q+i][j]=a[q+i][j]*f - a[q][j]*e;
            }
            /* теперь выводим то, что мы натворили */
            output(n, a, false);
            cout << "L" << q+i+1 << "*" << f << " - L" << q+1 << "*" << e << endl;
        }
        p++;
    }
}

void LRdecomposition(int n, double **a) // Разложение Холецкого
{
    double *y = (double*)malloc(n*sizeof(double));
    double**l = (double**)malloc(n*sizeof(double*)); // матрица разложения
    for (int i=0; i<n; i++)
    {
        l[i] = (double*)malloc((n+1)*sizeof(double));
    }
    int sum;
    if(!symmetrical(n, a)) // проверяем на симметричность
    {
        cout << "Not symmetrical" << endl;
        return;
    }
    cout << "Matrix is symmetrical" << endl;

    l[0][0] = sqrt(a[0][0]);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            sum = 0;
            if (j == i)
            {
                cout << "sum = ";
                for (int k = 0; k < j; k++)
                {
                    sum += l[j][k]*l[j][k];
                    cout << l[j][k] << "^2 + ";
                }
                cout << " 0 " << endl;
                l[j][j] = sqrt(a[j][j] - sum);
                cout << "l[" << j << "][" << j << "] = sqrt(" << a[j][j] << " - " << sum << ") = " << l[j][j] << endl;
            }
            else {

                for (int k = 0; k < j; k++)
                {
                    sum += (l[i][k] * l[j][k]);
                }
                l[i][j] = (a[i][j] - sum) / l[j][j];
                cout << "l[" << i << "][" << j << "] = (" << a[i][j] << " - " << sum << ") / " << l[j][j] << " = " << l[i][j] << endl;
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        l[i][n] = a[i][n];
    }

    output(n, l, false);
    outL(n, l, y);
    int tmp;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < i; j++)
        {
            tmp = l[i][j];
            l[i][j] = l[j][i];
            l[j][i] = tmp;
        }
    l[i][n] = y[i];
    }
    output(n, l, false);
    outR(n, l);
}

void QRrotation(int n, double **a) // Реализация метода вращений
{
    double**reserve = (double**)malloc(n*sizeof(double*)); // Резервная матрица
    for (int i=0; i<n; i++)
    reserve[i]=(double*)malloc((n+1)*sizeof(double));

    double **T = (double**)malloc(n*sizeof(double*)); // Матрица поворота
    for(int i = 0; i < n; i++)
    T[i] = (double*)malloc(n*sizeof(double));

    double **Q = (double**)malloc(n*sizeof(double*)); // Унитарная матрица
    for(int i = 0; i < n; i++)
    Q[i] = (double*)malloc(n*sizeof(double));

    double *Qb = (double*)malloc(n*sizeof(double)); // Будущий Q*b
    for(int i = 0; i < n; i++)
    Qb[i] = 0;

    for (int i = 0; i < n; i++) // Пусть они берут своё начало как единичные
    {
        for(int j = 0; j < n; j++)
        {
            reserve[i][j] = a[i][j];
            if (i == j)
            {
                T[i][i] = 1;
                Q[i][i] = 1;
            }
            else{
                T[i][j] = 0;
                Q[i][j] = 0;
            }
        }
        reserve[i][n] = a[i][n];
    }
    int iteration = 1;
    double k1 = a[0][n-1];
    double k2 = a[0][n-2];
    double sin, cos;

    for(int i = 0; i < n; i++)
    {
        for(int j = n-1; j > i; j--)
        {
            k1 = a[j][i]; // Считаем углы
            k2 = a[j-1][i];
            sin = -k1/sqrt(k1*k1 + k2*k2);
            cos = k2/sqrt(k1*k1 + k2*k2);
            T[j][j-1] = sin; // Матрица поворота
            T[j][j] = cos;
            T[j-1][j] = -sin;
            T[j-1][j-1] = cos;
            left_product(n, T, a); // Применяем T к A
            left_product(n, T, Q); // Применяем T к Q
            cout << "T" << iteration << ": " << endl;
            output(n, T, true);
            cout << "A" << iteration << ": " << endl;
            output(n, a, true);
            for (int i = 0; i < n; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        T[i][i] = 1;
                    }
                    else
                    {
                        T[i][j] = 0;
                    }
                }
            }
            iteration++;
        }
    }
    cout << "Q^(T):" << endl;
    output(n, Q, true);
    cout << "R = Q^T * A:" << endl;
    left_product(n, Q, reserve);

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(reserve[i][j] < 0.00001)
            {
                reserve[i][j] = 0;
            }
            Qb[i] += reserve[j][n] * Q[i][j];
        }
    }
    for(int i = 0; i < n; i++)
        reserve[i][n] = Qb[i];

    output(n, reserve, false);
    outR(n, reserve);
}

void QReflection(int n, double **a)
{
    double**reserve = (double**)malloc(n*sizeof(double*)); // Резервная матрица
    for (int i=0; i<n; i++)
    reserve[i]=(double*)malloc((n+1)*sizeof(double));

    double **U = (double**)malloc(n*sizeof(double*)); // Матрица U
    for(int i = 0; i < n; i++)
    U[i] = (double*)malloc(n*sizeof(double));

    double **P = (double**)malloc(n*sizeof(double*)); // Матрица отражения
    for(int i = 0; i < n; i++)
    P[i] = (double*)malloc(n*sizeof(double));

    double **Q = (double**)malloc(n*sizeof(double*)); // Унитарная матрица
    for(int i = 0; i < n; i++)
    Q[i] = (double*)malloc(n*sizeof(double));

    double *Qb = (double*)malloc(n*sizeof(double)); // Будущий Q*b
    for(int i = 0; i < n; i++)
    Qb[i] = 0;

    double gamma, sc;

    for (int i = 0; i < n; i++) // Пусть они берут своё начало как единичные
    {
        for(int j = 0; j < n; j++){
            reserve[i][j] = a[i][j];
            if (i == j){
                U[i][i] = 1;
                Q[i][i] = 1;
            }
            else{
                U[i][j] = 0;
                Q[i][j] = 0;
            }
        }
        reserve[i][n] = a[i][n];
    }

    for(int r =0; r < n-1; r++) // Начинаем цикл
    {
        double *u = (double*)malloc((n-r)*sizeof(double)); // векторы размерности n-r
        double *v = (double*)malloc((n-r)*sizeof(double));
        cout << "Vector s = {";
        for(int j = 0; j < n-r; j++)
        {
            u[j] = a[r][j]; // заполняем элементами из r-го столбца А
            cout << u[j];
            if (j == n-r-1) continue;
            cout << ", ";
        }
        cout << "} " << endl; // ПРАВКА СКАЛЯРНОГО ПРОИЗВЕДЕНИЯ
        sc = scalar(n-r, u, u); // считаем скалярное произведения вектора на себя
        cout << "Scalar product (s, s) = " << sc << endl;

        if (sc == 0) // если оно нулевое, то вектор v это e_r, gamma = 1/2
        {
            cout << "Vector v = e_" << r << "Y = 1/2";
            for(int i = 0; i < n-r; i++)
            {
                v[i] = 0;
                if(i == r)
                {
                    v[i] = 1;
                }
            }
            gamma = 0.5;
        }
        else
        {
            cout << "Vector u = {";
            for(int j = 0; j < n-r; j++)
            {
                u[j] = u[j]/sqrt(sc); // u = s/sqrt({s,s})
                cout << u[j];
                if (j == n-r-1) continue;
                cout << ", ";
            }
            cout << "} " << endl;
            /* Далее, v_i = 0, если i<r
             * v_i = u_{i-r+1}, если i>r
             * В случае i = r выполняется v_i = 1, если u_1 = 0, иначе v_r = (u_1/|u_1|)*(1 + |u_1|)
             * При этом gamma = 1 + |u_1| = |v_r|
             */
            cout << "Vector v = {";
            for(int i = 0; i < n-r; i++)
            {
                if(i < r)
                {
                    v[i] = 0;
                }
                if(i > r)
                {
                    v[i] = u[i-r];
                }
                else
                {
                    if(u[1] == 0)
                    {
                        v[i] = 1;
                    }
                    else
                    {
                        v[i] = (u[1]/abs(u[1])) * (1+abs(u[1]));
                    }
                }
                cout << v[i];
                if (i == n-r-1) continue;
                cout << ", ";
            }
            cout << "} " << endl;
            gamma = v[r];
        }

        /* После вычисления вектора v подстолбцы справа модифицируются по формуле a' = a - v*{a,v}/gamma */

        for (int i = r; i < n; i++) // столбец
        {
            double *x = (double*)malloc((n-r)*sizeof(int));
            for(int j = r; j < n; j++)
                x[j] = a[i][j];

            sc = scalar(n-r, x, v);
            for(int j = r; j < n; j++) // строка
            {
                a[i][j] -= v[j] * sc / gamma;
            }
        }

        minus2vec(n, gamma, v, U); // Считаем матрицу преобразования
        left_product(n, U, a); // Выполняем преобразование для A
        left_product(n, U, Q); // Выполняем преобразование для Q
        cout << "A" << r << ": ";
        output(n, a, true);
        cout << "U" << r << ": ";
        output(n, U, true);

        for (int i = 0; i < n; i++) // Снова делаем матрицу U единичной
        {
            for(int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    U[i][i] = 1;
                }
                else
                {
                    U[i][j] = 0;
                }
            }
        }
    }
    cout << "Q^(T):" << endl;
    output(n, Q, true);
    cout << "R = Q^T * A:" << endl;
    left_product(n, Q, reserve);

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(reserve[i][j] < 0.00001)
            {
                reserve[i][j] = 0;
            }
            Qb[i] += reserve[j][n] * Q[i][j];
        }
    }
    for(int i = 0; i < n; i++)
        reserve[i][n] = Qb[i];

    output(n, reserve, false);
    outR(n, reserve);
}

int main()
{
    string s;
    int n;
    cout << "Enter the number of equasions: " << endl;
    cin >> n;
    cout << "Enter matrix: " << endl;
    double**a = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++)
    a[i]=(double*)malloc((n+1)*sizeof(double));
    if (a==NULL)
    {
        cout << "error";
    }
    else
    {
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<n+1; j++)
            {
                cin >> a[i][j];
            }
        }
        cout << "Input matrix is: ";
        output(n, a, false);
        cout << "Enter Gauss for gaussian elimination" << endl;
        cout << "Enter lr for Cholesky's LR decomposition" << endl;
        cout << "Enter qr1 for QR decomposition by rotation method" << endl;
        cout << "Enter qr2 for QR decomposition by reflection method" << endl;
        while(s != "exit")
        {
            cin >> s;
            if(s == "Gauss")
            {
                Elimination(n, a);
                output(n, a, false);
                check_rk(n, a);
            }
            if(s == "lr")
            {
                LRdecomposition(n, a);
            }
            if(s == "qr1")
            {
                QRrotation(n, a);
            }
            if(s == "qr2")
            {
                QReflection(n, a);
            }
        }
        return 0;
    }
}
