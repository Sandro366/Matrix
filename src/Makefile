CC = g++
FLAGS = -Wall -Wextra -Werror -std=c++17
LIB = s21_matrix_plus
TESTFILE = test
GCOV_FLAG = -fprofile-arcs -ftest-coverage

all: clean $(LIB).a test gcov_report

clean:
	@echo "****** clean trash ******"
	rm -f *.a *.o $(TESTFILE) *.info *.gcno *.gcda 
	rm -rf report/
	rm -rf CPPLINT.cfg

s21_matrix_plus.a: clean
	$(CC) -c $(FLAGS) $(LIB).cpp -o $(LIB).o
	ar rcs $(LIB).a $(LIB).o
	ranlib $@

test: clean
	$(CC) $(GCOV_FLAG) $(LIB).cpp $(TESTFILE).cpp -o $(TESTFILE) -lgtest -lpthread
	./$(TESTFILE)

gcov_report: test
	@echo "****** gcov start ******"
	lcov -t "test" -o $(TESTFILE).info -c -d ./
	lcov --remove test.info "/usr/include/*" "/usr/local/*" "*asio/*" "*googletest/*" -o filtered_test.info
	genhtml -o report filtered_test.info
	open ./report/index.html

linter:
	@echo "****** linter start ******"
	cp -f ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extensions=cpp *.cpp *.h
	rm -rf CPPLINT.cfg

cppcheck:
	@echo "****** cppcheck start ******"
	cppcheck --enable=all --suppress=missingIncludeSystem --language=c++ $(LIB).cpp *.h

leaks: test
ifeq ($(shell uname), Linux)
	CK_FORK=no valgrind -s --tool=memcheck --leak-check=full --show-leak-kinds=all --show-reachable=yes ./$(TESTFILE)
else
	CK_FORK=no leaks --atExit -- ./$(TESTFILE)
endif
