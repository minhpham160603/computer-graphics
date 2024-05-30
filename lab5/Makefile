export LDFLAGS=-L/opt/homebrew/opt/libomp/lib
export CPPFLAGS=-I/opt/homebrew/opt/libomp/include

all:
	clang++ -o main main.cpp -std=c++17 -Xpreprocessor -fopenmp $(LDFLAGS) $(CPPFLAGS) -lomp 
clean:
	rm -f main