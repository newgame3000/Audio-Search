all:
	g++ -g -pedantic -std=c++17 -O2 -Wall -Wextra main.cpp -o solution
clean:
	rm -rf solution
