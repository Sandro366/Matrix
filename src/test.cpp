#include "gtest/gtest.h"
#include "s21_matrix_plus.h"

void Fill(S21Matrix& tmp, double value) {
  for (int i = 0; i < tmp.getRow(); i++)
    for (int j = 0; j < tmp.getColumn(); j++) tmp(i, j) = value++;
}

TEST(ConstructMatrix, test1) {
  double tmp[1][1] = {0};
  S21Matrix tmp_matrix = S21Matrix();
  EXPECT_EQ(tmp[0][0], tmp_matrix.getValue(0, 0));
}

TEST(ConstructMatrix, test2) {
  int row = 5;
  int col = 5;
  double tmp[row][col] = {{0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0}};
  S21Matrix tmp_matrix = S21Matrix(5, 5);
  EXPECT_EQ(row, tmp_matrix.getRow());
  EXPECT_EQ(col, tmp_matrix.getColumn());
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      EXPECT_EQ(tmp[i][j], tmp_matrix.getValue(i, j));
    }
  }
}

TEST(ConstructMatrix, test3) {
  S21Matrix A(3, 3);
  Fill(A, 1);
  S21Matrix B(A);
  bool check = A.eq_matrix(B);
  EXPECT_EQ(check, true);
  EXPECT_EQ(A.getRow(), B.getRow());
  EXPECT_EQ(A.getColumn(), B.getColumn());
}

TEST(ConstructMatrix, test4) {
  S21Matrix A(3, 3);
  Fill(A, 1);
  S21Matrix C(A);
  S21Matrix B = std::move(A);
  bool check = C.eq_matrix(B);
  EXPECT_EQ(check, true);
  EXPECT_EQ(C.getRow(), B.getRow());
  EXPECT_EQ(C.getColumn(), B.getColumn());
}

TEST(SumMatrix, test1) {
  S21Matrix A(3, 3), B(3, 3), tmp1(3, 3), tmp2(3, 3);
  Fill(A, 1);
  Fill(B, 2);
  S21Matrix C(3, 3);
  C(0, 0) = 3;
  C(0, 1) = 5;
  C(0, 2) = 7;
  C(1, 0) = 9;
  C(1, 1) = 11;
  C(1, 2) = 13;
  C(2, 0) = 15;
  C(2, 1) = 17;
  C(2, 2) = 19;
  tmp1 = A + B;
  tmp2 = A;
  tmp2 += B;
  bool check1 = tmp1.eq_matrix(C);
  bool check2 = tmp2.eq_matrix(C);
  EXPECT_EQ(check1, true);
  EXPECT_EQ(check2, true);
}

TEST(SubMatrix, test1) {
  S21Matrix A(3, 3), B(3, 3), tmp1(3, 3), tmp2(3, 3);
  Fill(A, 2);
  Fill(B, 1);
  S21Matrix C(3, 3);
  C(0, 0) = 1;
  C(0, 1) = 1;
  C(0, 2) = 1;
  C(1, 0) = 1;
  C(1, 1) = 1;
  C(1, 2) = 1;
  C(2, 0) = 1;
  C(2, 1) = 1;
  C(2, 2) = 1;
  tmp1 = A - B;
  tmp2 = A;
  tmp2 -= B;
  bool check1 = tmp1.eq_matrix(C);
  bool check2 = tmp2.eq_matrix(C);
  EXPECT_EQ(check1, true);
  EXPECT_EQ(check2, true);
}

TEST(EqMatrix, test1) {
  S21Matrix A(3, 3);
  Fill(A, 1);
  S21Matrix B = A;
  bool check = A == B;
  EXPECT_EQ(check, true);
}

TEST(MulNumber, test1) {
  S21Matrix A(3, 3), B(3, 3);
  Fill(A, 2.3);
  S21Matrix C(3, 3);
  C(0, 0) = 2.53;
  C(0, 1) = 3.63;
  C(0, 2) = 4.73;
  C(1, 0) = 5.83;
  C(1, 1) = 6.93;
  C(1, 2) = 8.03;
  C(2, 0) = 9.13;
  C(2, 1) = 10.23;
  C(2, 2) = 11.33;
  B = A;
  A = A * 1.1;
  B *= 1.1;
  bool check1 = A.eq_matrix(C);
  bool check2 = B.eq_matrix(C);
  EXPECT_EQ(check1, true);
  EXPECT_EQ(check2, true);
}

TEST(MulMatrix, test1) {
  S21Matrix A(4, 3), B(3, 4), D(4, 3);
  Fill(A, 1.1);
  Fill(B, -3.4);
  S21Matrix C(4, 4);
  C(0, 0) = 11.78;
  C(0, 1) = 18.08;
  C(0, 2) = 24.38;
  C(0, 3) = 30.68;
  C(1, 0) = 17.18;
  C(1, 1) = 32.48;
  C(1, 2) = 47.78;
  C(1, 3) = 63.08;
  C(2, 0) = 22.58;
  C(2, 1) = 46.88;
  C(2, 2) = 71.18;
  C(2, 3) = 95.48;
  C(3, 0) = 27.98;
  C(3, 1) = 61.28;
  C(3, 2) = 94.58;
  C(3, 3) = 127.88;
  D = A;
  A = A * B;
  D *= B;
  bool check1 = A.eq_matrix(C);
  bool check2 = D.eq_matrix(C);
  EXPECT_EQ(check1, true);
  EXPECT_EQ(check2, true);
}

TEST(TransposeMatrix, test1) {
  S21Matrix A(3, 3);
  Fill(A, 1);
  S21Matrix B = A.transpose();
  S21Matrix C(3, 3);
  C(0, 0) = 1;
  C(0, 1) = 4;
  C(0, 2) = 7;
  C(1, 0) = 2;
  C(1, 1) = 5;
  C(1, 2) = 8;
  C(2, 0) = 3;
  C(2, 1) = 6;
  C(2, 2) = 9;
  bool check = B.eq_matrix(C);
  EXPECT_EQ(check, true);
}

TEST(CalcCompl, test1) {
  S21Matrix A(3, 3);
  Fill(A, -10);
  A(0, 0) = 10;
  S21Matrix B = A.calc_complements();
  S21Matrix C(3, 3);
  C(0, 0) = -3;
  C(0, 1) = 6;
  C(0, 2) = -3;
  C(1, 0) = 6;
  C(1, 1) = -52;
  C(1, 2) = 66;
  C(2, 0) = -3;
  C(2, 1) = 106;
  C(2, 2) = -123;
  bool check = B.eq_matrix(C);
  EXPECT_EQ(check, true);
}

TEST(Determinant, test1) {
  S21Matrix A(3, 3);
  Fill(A, 5);
  A(2, 2) = 100;
  double det = A.determinant();
  EXPECT_EQ(det, -261);
}

TEST(Determinant, test2) {
  S21Matrix A(2, 2);
  A(0, 0) = 50;
  A(0, 1) = -10;
  A(1, 0) = 1;
  A(1, 1) = 5.5;
  double det = A.determinant();
  EXPECT_EQ(det, 285);
}

TEST(Determinant, test3) {
  S21Matrix A(1, 1);
  A(0, 0) = 5.55555;
  double det = A.determinant();
  EXPECT_EQ(det, 5.55555);
}

TEST(InverseMatrix, test1) {
  S21Matrix A(3, 3);
  Fill(A, -4);
  A(0, 0) = 16;
  S21Matrix B = A.inverse_matrix();
  S21Matrix C(3, 3);
  C(0, 0) = 1.0 / 20;
  C(0, 1) = -1.0 / 10;
  C(0, 2) = 1.0 / 20;
  C(1, 0) = -1.0 / 10;
  C(1, 1) = -17.0 / 15;
  C(1, 2) = 7.0 / 30;
  C(2, 0) = 1.0 / 20;
  C(2, 1) = 9.0 / 10;
  C(2, 2) = 1.0 / 20;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(C.getValue(i, j), B.getValue(i, j), 0.0000001);
    }
  }
}

TEST(AccessorsMutators, test1) {
  S21Matrix A;
  A.setColumns(3);
  A.setRows(4);
  EXPECT_EQ(A.getRow(), 4);
  EXPECT_EQ(A.getColumn(), 3);
}

TEST(AccessorsMutators, test2) {
  S21Matrix A(5, 10);
  A.setColumns(8);
  A.setRows(4);
  EXPECT_EQ(A.getRow(), 4);
  EXPECT_EQ(A.getColumn(), 8);
}

TEST(ExceptMatrix, test1) {
  S21Matrix A(4, 4);
  EXPECT_THROW(A.setColumns(-1), std::invalid_argument);
  EXPECT_THROW(A.setRows(0.5), std::invalid_argument);
}

TEST(ExceptMatrix, test2) {
  S21Matrix A(4, 4);
  EXPECT_THROW(A.getValue(-2, 3), std::invalid_argument);
  EXPECT_THROW(A.setValue(10, 5, 5), std::invalid_argument);
}

TEST(ExceptMatrix, test3) {
  S21Matrix A(4, 4), B(4, 5);
  Fill(A, 2);
  Fill(B, 1);
  EXPECT_THROW(A + B, std::invalid_argument);
  EXPECT_THROW(B - A, std::invalid_argument);
}

TEST(ExceptMatrix, test4) {
  S21Matrix A(4, 4), B(5, 4);
  Fill(A, 2);
  Fill(B, 1);
  EXPECT_THROW(A * B, std::invalid_argument);
}

TEST(ExceptMatrix, test5) {
  S21Matrix A(2, 4);
  Fill(A, 2);
  EXPECT_THROW(A.calc_complements(), std::invalid_argument);
  EXPECT_THROW(A.determinant(), std::invalid_argument);
}

TEST(ExceptMatrix, test6) {
  S21Matrix A(3, 3);
  Fill(A, 2);
  EXPECT_THROW(A.inverse_matrix(), std::invalid_argument);
}

int main() {
  testing::InitGoogleTest();
  return RUN_ALL_TESTS();
}
