#include <vector>
#include <functional>
#include <exception>
#include <fstream>
#include <string>

using Vec = std::vector<double>;

Vec operator*(double scalar, const Vec& array)
{
    Vec result(array.size());

    for (size_t i = 0; i < array.size(); ++i)
    {
        result[i] = scalar * array[i];
    }

    return result;
}

Vec operator+(const Vec& lhs, const Vec& rhs)
{
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("Размеры массивов не совпадают!");

    Vec result(lhs.size());

    for (size_t i = 0; i < result.size(); ++i)
    {
        result[i] = lhs[i] + rhs[i];
    }

    return result;
}

class OdeSolver
{
public:
    OdeSolver(std::function<Vec(double, Vec)> func
        , double xStart
        , double xEnd
        , double step
        , Vec y0
        )
        : _func(func)
        , _xStart(xStart)
        , _xEnd(xEnd)
        , _step(step)
        , _y0(y0)
    {
        _updateGrid();
    }

    Vec rk4()
    {
        std::vector<Vec> y(_grid.size());
        Vec result(_grid.size());
        Vec k1;
        Vec k2;
        Vec k3;
        Vec k4;

        y[0] = _y0;
        result[0] = y[0][0];

        for (size_t i = 1; i < _grid.size(); ++i)
        {
            k1 = _func(_grid[i - 1], y[i - 1]);
            k2 = _func(_grid[i - 1] + _step / 2, y[i - 1] + _step / 2 * k1);
            k3 = _func(_grid[i - 1] + _step / 2, y[i - 1] + _step / 2 * k2);
            k4 = _func(_grid[i - 1] + _step, y[i - 1] + _step * k3);

            y[i] = y[i - 1] + _step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

            result[i] = y[i][0];
        }

        return result;
    }

    Vec getGrid()
    {
        return _grid;
    }

private:
    void _updateGrid()
    {
        size_t gridSize = (_xEnd - _xStart) / _step + 1;

        _grid.resize(gridSize);

        _grid[0] = _xStart;

        for (size_t i = 1; i < gridSize; ++i)
        {
            _grid[i] = _grid[i - 1] + _step;
        }
    }

    std::function<Vec(double, Vec)> _func;
    double _xStart;
    double _xEnd;
    double _step;
    Vec _grid;
    Vec _y0;

};


Vec rightHandFunction(double x, const Vec& y)
{
    Vec result(5);
    result[0] = y[1];
    result[1] = y[2];
    result[2] = y[3];
    result[3] = y[4];
    result[4] = -(15 * y[4] + 90 * y[3] + 270 * y[2] + 405 * y[1] + 243 * y[0]);

    return result;
}

void saveToFile(const Vec& x, const Vec& y, const std::string& filename)
{
    std::ofstream outFile(filename);

    for (size_t i = 0; i < x.size(); ++i)
    {
        outFile << x[i] << " " << y[i] << "\n";
    }
}

int main()
{
    Vec y0{0, 3, -9, -8, 0};
    double xStart = 0;
    double xEnd = 5;
    double step = 0.001;

    OdeSolver odeSolver(rightHandFunction, xStart, xEnd, step, y0);

    Vec result = odeSolver.rk4();

    saveToFile(odeSolver.getGrid(), result, "cpp_results.txt");

    return 0;
}