#ifndef SRC_S21_MATRIX_PLUS_H_
#define SRC_S21_MATRIX_PLUS_H_

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  int _rows, _cols;
  double** _matrix;

  void Create_matrix();
  void Delete_matrix();
  S21Matrix minor(int row, int col);
  void Copy_matrix(double** other);

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  bool eq_matrix(const S21Matrix& other);
  void sum_matrix(const S21Matrix& other);
  void mul_number(const double num);
  void sub_matrix(const S21Matrix& other);
  void mul_matrix(const S21Matrix& other);
  S21Matrix transpose();
  S21Matrix calc_complements();
  double determinant();
  S21Matrix inverse_matrix();

  bool operator==(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const double number);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const double number);
  S21Matrix& operator*=(const S21Matrix& other);
  double& operator()(int row, int col);

  int getRow() const;
  int getColumn() const;
  double getValue(int row, int col);
  void setRows(int rows);
  void setColumns(int columns);
  void setValue(double value, int row, int col);
};

#endif  // SRC_S21_MATRIX_PLUS_H_
