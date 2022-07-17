#include "s21_matrix_plus.h"

S21Matrix::S21Matrix() : _rows(1), _cols(1) { Create_matrix(); }

S21Matrix::S21Matrix(int rows, int cols) : _rows(rows), _cols(cols) {
  if (_rows < 1 || _cols < 1) throw std::invalid_argument("Incorrect input");
  Create_matrix();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : _rows(other._rows), _cols(other._cols), _matrix(NULL) {
  Copy_matrix(other._matrix);
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : _rows(other._rows), _cols(other._cols), _matrix(other._matrix) {
  other._rows = 0;
  other._cols = 0;
  other._matrix = nullptr;
}

S21Matrix::~S21Matrix() {
  if (_matrix) {
    if (_rows > 0 && _cols > 0) Delete_matrix();
    _matrix = nullptr;
  }
}

bool S21Matrix::eq_matrix(const S21Matrix& other) {
  bool result = true;
  if (_rows == other._rows && _cols == other._cols) {
    for (int i = 0; i < _rows; i++) {
      for (int j = 0; j < _cols; j++) {
        if (fabs(_matrix[i][j] - other._matrix[i][j]) >= 1e-7) {
          result = false;
          return result;
        }
      }
    }
  } else {
    result = false;
  }
  return result;
}

void S21Matrix::sum_matrix(const S21Matrix& other) {
  if (!_matrix || !other._matrix || _rows != other._rows ||
      _cols != other._cols)
    throw std::invalid_argument("Incorrect size of matrix");
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] += other._matrix[i][j];
}

void S21Matrix::sub_matrix(const S21Matrix& other) {
  if (!_matrix || !other._matrix || _rows != other._rows ||
      _cols != other._cols)
    throw std::invalid_argument("Incorrect size of matrix");
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] -= other._matrix[i][j];
}

void S21Matrix::mul_number(const double num) {
  if (!_matrix) throw std::invalid_argument("Matrix is empty");
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] *= num;
}

void S21Matrix::mul_matrix(const S21Matrix& other) {
  if (!_matrix || !other._matrix || _cols != other._rows)
    throw std::invalid_argument("Incorrect size of matrix");
  S21Matrix result(_rows, other._cols);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      for (int k = 0; k < _cols; k++) {
        result._matrix[i][j] += _matrix[i][k] * other._matrix[k][j];
      }
    }
  }
  *this = result;
}

S21Matrix S21Matrix::transpose() {
  S21Matrix result(_cols, _rows);
  for (int i = 0; i < _cols; i++)
    for (int j = 0; j < _rows; j++) result._matrix[i][j] = _matrix[j][i];
  return result;
}

S21Matrix S21Matrix::calc_complements() {
  if (_cols != _rows)
    throw std::invalid_argument("Incorrect input! Matrix must be square");
  S21Matrix result(_rows, _cols);
  if (_rows == 1) {
    result._matrix[0][0] = _matrix[0][0];
  } else {
    for (int i = 0; i < _rows; i++) {
      for (int j = 0; j < _cols; j++) {
        S21Matrix tmp = minor(i, j);
        result._matrix[i][j] = pow(-1, i + j) * tmp.determinant();
      }
    }
  }
  return result;
}

double S21Matrix::determinant() {
  if (_cols != _rows)
    throw std::invalid_argument("Incorrect input! Matrix must be square");
  double result = 0.0;
  if (_rows == 1) {
    result = _matrix[0][0];
  } else if (_rows == 2) {
    result = _matrix[0][0] * _matrix[1][1] - _matrix[1][0] * _matrix[0][1];
  } else {
    for (int i = 0; i < _cols; i++) {
      S21Matrix tmp = minor(i, 0);
      result += _matrix[i][0] * tmp.determinant() * pow(-1, i);
    }
  }
  return result;
}

S21Matrix S21Matrix::inverse_matrix() {
  double deter = determinant();
  if (fabs(deter) < 1e-7)
    throw std::invalid_argument("Matrix determinant is 0");
  S21Matrix result = calc_complements();
  if (_rows == 1) {
    result._matrix[0][0] = 1.0 / _matrix[0][0];
  } else {
    result.mul_number(1 / deter);
    result = result.transpose();
  }
  return result;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    Delete_matrix();
    _rows = other._rows;
    _cols = other._cols;
    Create_matrix();
    for (int i = 0; i < _rows; i++) {
      for (int j = 0; j < _cols; j++) {
        _matrix[i][j] = other._matrix[i][j];
      }
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result = *this;
  result.sum_matrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result = *this;
  result.sub_matrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double number) {
  S21Matrix result = *this;
  result.mul_number(number);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result = *this;
  result.mul_matrix(other);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) { return eq_matrix(other); }

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  sum_matrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  sub_matrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double number) {
  mul_number(number);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  mul_matrix(other);
  return *this;
}

double& S21Matrix::operator()(int row, int col) {
  if (row >= _rows || col >= _cols || row < 0 || col < 0)
    throw std::invalid_argument("Incorrect index of matrix");
  return _matrix[row][col];
}

int S21Matrix::getRow() const { return _rows; }

int S21Matrix::getColumn() const { return _cols; }

double S21Matrix::getValue(int row, int col) {
  if (row >= _rows || col >= _cols || row < 0 || col < 0)
    throw std::invalid_argument("Wrong index");
  return _matrix[row][col];
}

void S21Matrix::setRows(int rows) {
  if (rows < 1) throw std::invalid_argument("Incorrect size of matrix");
  if (rows != _rows) {
    Delete_matrix();
    _rows = rows;
    Create_matrix();
  }
}

void S21Matrix::setColumns(int columns) {
  if (columns < 1) throw std::invalid_argument("Incorrect size of matrix");
  if (columns != _cols) {
    Delete_matrix();
    _cols = columns;
    Create_matrix();
  }
}

void S21Matrix::setValue(double value, int row, int col) {
  if (row >= _rows || col >= _cols || row < 0 || col < 0)
    throw std::invalid_argument("Wrong index");
  _matrix[row][col] = value;
}

//***************************************************//

void S21Matrix::Create_matrix() {
  _matrix = new double*[_rows];
  for (int i = 0; i < _rows; i++) _matrix[i] = new double[_cols]();
}

void S21Matrix::Delete_matrix() {
  for (int i = 0; i < _rows; i++) delete[] _matrix[i];
  delete[] _matrix;
}

S21Matrix S21Matrix::minor(int row, int col) {
  S21Matrix matrix(_rows - 1, _cols - 1);
  int m = 0;
  for (int i = 0; i < _rows; i++) {
    if (i != row) {
      int n = 0;
      for (int j = 0; j < _cols; j++) {
        if (j != col) {
          matrix._matrix[n][m] = _matrix[i][j];
          n++;
        }
      }
      m++;
    }
  }
  return matrix;
}

void S21Matrix::Copy_matrix(double** other) {
  Create_matrix();
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] = other[i][j];
}
