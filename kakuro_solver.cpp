#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>
#include <array>
//#include <bits/stdc++.h>

#include <chrono>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <stack>
#include <set>
#include <omp.h>

using namespace std;

enum direction { d_down, d_right, none };

#define COORD std::pair<int, int>

//#define DEBUG

int iter = 0;

///Auxiliary functions

void display_arr(int* arr, int n) {

	cout << "arr: ";

	for (int i = 0; i < n; i++) {
		cout << arr[i] << " ";
	}

	cout << endl;

}

void print_coords(COORD start, COORD end) {

	cout << "Start:" << start.first << "," << start.second << endl;
	cout << "End:" << end.first << "," << end.second << endl;

}

int find_length(COORD start, COORD end, direction dir) {

	if (dir == d_down)
		return end.first - start.first;
	if (dir == d_right)
		return end.second - start.second;

	return -1;
}

void convert_sol(int** mat, int**& sol_mat, int m, int n) {

	sol_mat = new int* [m]; //Rows
	for (int i = 0; i < m; i++) {
		sol_mat[i] = new int[n]; //Cols
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) { // ther was a bug here, it was m
			if (mat[i][j] == -2)
				sol_mat[i][j] = -2; //Empty value cell
			else
				sol_mat[i][j] = -1; //Hint or empty cell
		}
	}
}

void print_one_matrix(int** matrix, int m, int n) {
	std::cout << "Matrix: " << std::endl;
	for (int i = 0; i < m; i++) { //rows
		for (int j = 0; j < n; j++) { //cols
			std::cout << matrix[i][j] << "\t";
		}
		std::cout << "\n";
	}
}

void sol_to_file(int** mat, int** sol_mat, int m, int n, string fname) {

	//string fname = "visualize.kakuro";
	ofstream to_write(fname);

	to_write << m << " " << n << "\n";

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (mat[i][j] != -2)
				to_write << mat[i][j] << " ";
			else
				to_write << sol_mat[i][j] << " ";
		}
		to_write << "\n";
	}

	to_write.close();
}

void read_matrix(int**& matrix, std::ifstream& afile, int m, int n) {

	matrix = new int* [m]; //rows

	for (int i = 0; i < m; i++) {
		matrix[i] = new int[n]; //cols
	}

	int val;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			afile >> val;
			matrix[i][j] = val;
		}
	}
}

///Auxiliary functions

struct sum {
	COORD start;
	COORD end;

	int hint;
	int dir;
	int length;
	int* arr;

	void print_sum() {
		cout << "############################" << endl;
		cout << "Creating sum with: " << endl;
		print_coords(start, end);
		cout << "Hint: " << hint << endl;
		cout << "Direction: " << dir << endl;
		cout << "Length: " << length << endl;
		cout << "############################" << endl;
	}

	sum(COORD _start, COORD _end, int _hint, direction _dir) : start(_start), end(_end), hint(_hint), dir(_dir)
	{
		length = find_length(_start, _end, _dir);
		arr = new int[length];

#ifdef DEBUG
		cout << "############################" << endl;
		cout << "Creating sum with: " << endl;
		print_coords(start, end);
		cout << "Hint: " << hint << endl;
		cout << "Direction: " << dir << endl;
		cout << "Length: " << length << endl;
		cout << "############################" << endl;
#endif
	}

	//~sum(){
	//delete arr;
	//}
};


COORD find_end(int** matrix, int m, int n, int i, int j, direction dir) { //0 down 1 right

	if (dir == d_right) {
		for (int jj = j + 1; jj < n; jj++) {
			if (matrix[i][jj] != -2 || jj == n - 1) {
				if (matrix[i][jj] == -2 && jj == n - 1)
					jj++;
				COORD END = COORD(i, jj);
				return END;
			}
		}
	}

	if (dir == d_down) {
		for (int ii = i + 1; ii < m; ii++) {
			if (matrix[ii][j] != -2 || ii == m - 1) {
				if (matrix[ii][j] == -2 && ii == m - 1)
					ii++;
				COORD END = COORD(ii, j);
				return END;
			}
		}
	}

}


vector<sum> get_sums(int** matrix, int m, int n) {

	vector<sum> sums;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			int val = matrix[i][j];
			if (val != -1 && val != -2) {
				int hint = val;
				hint = hint / 10;

				if ((hint % 100) == 0) {
					hint = (int)(hint / 100);
					COORD START = COORD(i, j + 1);
					COORD END = find_end(matrix, m, n, i, j, d_right);
					sum _sum = sum(START, END, hint, d_right);
					sums.push_back(_sum);
				}

				else {
					int div = (int)(hint / 100);
					int rem = (int)(hint % 100);

					if (div == 0 && rem != 0) {
						COORD START = COORD(i + 1, j);
						COORD END = find_end(matrix, m, n, i, j, d_down);
						sum _sum = sum(START, END, rem, d_down);
						sums.push_back(_sum);
					}

					if (div != 0 && rem != 0) {
						COORD START1 = COORD(i + 1, j);
						COORD START2 = COORD(i, j + 1);
						COORD END1 = find_end(matrix, m, n, i, j, d_down);
						COORD END2 = find_end(matrix, m, n, i, j, d_right);
						sum _sum1 = sum(START1, END1, rem, d_down);
						sum _sum2 = sum(START2, END2, div, d_right);
						sums.push_back(_sum1);
						sums.push_back(_sum2);
					}
				}
			}


		}
	}
	return sums;
}

bool sum_is_a_smaller(const sum& a, const sum& b) {
	return a.hint > b.hint;
}

inline int set_sum_array_and_sum(const sum& s, int** sol_mat) {

	if (s.dir == d_down) {
		for (int i = 0; i < s.end.first - s.start.first; i++) {
			s.arr[i] = sol_mat[i + s.start.first][s.start.second];
		}
	}
	else {
		for (int i = 0; i < s.end.second - s.start.second; i++) {
			s.arr[i] = sol_mat[s.start.first][i + s.start.second];
		}
	}

	int ret = 0;
	for (short i = 0; i < s.length; i++) {
		ret += s.arr[i];
	}
	return ret;
}

inline bool sum_has_no_unfiled(const sum& S) {
	for (short i = 0; i < S.length; i++) {
		if (S.arr[i] == -2) {
			return false;
		}
	}
	return true;
}

// inputs are small so this is faster tahn a set or a hash table 
inline bool sum_has_no_duplicate(const sum& s) {
	bool has_duplicate = false;

	for (int i = 0; i < s.length - 1; ++i) {
		for (int j = i + 1; j < s.length; ++j) {
			if (s.arr[i] == s.arr[j]) {
				has_duplicate = true;
				break;
			}
		}
		if (has_duplicate) {
			break;
		}
	}
	return !has_duplicate;
}

inline bool is_sum_valid(const sum& s, int** sol_mat) {
	if (set_sum_array_and_sum(s, sol_mat) == s.hint)
		if (sum_has_no_unfiled(s))
			if (sum_has_no_duplicate(s))
				return true;
	return false;
}


inline bool is_all_sums_valid(const vector<sum>& sums, int** sol_mat) {
	
	for (int i = 0; i < sums.size(); i++) {
		if (!is_sum_valid(sums[i], sol_mat)) {
			return false;
		}
	}
	return true;
}

inline tuple<COORD, int > get_cor(const vector<sum>& sums, int** sol_mat) {
	for (int ss = 0; ss < sums.size(); ss++) {
		sum s = sums.at(ss);
		if (!is_sum_valid(s, sol_mat)) {
			for (short i = 0; i < s.length; i++) {
				if (s.arr[i] == -2) {
					if (s.dir == d_down) {
						return tuple < COORD, int>(COORD(s.start.first + i, s.start.second), ss);
					}
					else {
						return tuple < COORD, int>(COORD(s.start.first, s.start.second + i), ss);
					}

				}
			}
		}
	}
	return tuple < COORD, int>(COORD(-1, -1), -1);
}


sum  deep_copy(const sum& s) {
	sum ret = sum(s);
	ret.length = s.length; 
	ret.arr = new int[ret.length];
	for (int i = 0; i < ret.length; i++) {
		ret.arr[i] = s.arr[i];
	}
	return ret;
}

vector<sum> deep_copy(const vector<sum>& sums) {
	vector<sum> ret;
	for (int i = 0; i < sums.size(); i++) {
		ret.push_back(deep_copy(sums[i]));
	}
	return ret;
}


int** deep_copy(int**& mat, const int& m, const int& n) {
	int** ret = new int* [m];
	for (int i = 0; i < m; i++) {
		ret[i] = new int[n];
		for (int j = 0; j < n; j++) {
			ret[i][j] = mat[i][j];
		}
	}
	return ret;
}

void deep_delete(int**& mat, const int& m, const int& n) {
	for (int i = 0; i < m; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}

bool solution_single_thread_solver(int** mat, int** sol_mat, vector<sum>& sums, int& m, int& n) {
	stack<tuple<COORD, int, int>> callStack; // stack to store the state of each call
	
	// Initial state
	tuple<COORD, int> tmp = get_cor(sums, sol_mat);
	COORD next_cord = get<0>(tmp);
	int sums_index = get<1>(tmp);
	int i = 1;

	callStack.push(make_tuple(next_cord, sums_index, i));

	while (!callStack.empty()) {
		tie(next_cord, sums_index, i) = callStack.top();
		callStack.pop();

		if (sums_index == -1) {
			if (is_all_sums_valid(sums, sol_mat)) {
				return true;
			}
		}
		else {
			if (i < 10 ) {
				sol_mat[next_cord.first][next_cord.second] = i;
				i++;
				


				// Save the current state and push it back to the stack
				callStack.push(make_tuple(next_cord, sums_index, i));

				// Move on to the next state
				tmp = get_cor(sums, sol_mat);
				next_cord = get<0>(tmp);
				sums_index = get<1>(tmp);
				i = 1;
				callStack.push(make_tuple(next_cord, sums_index, i));
			}
			else {
				// Reset the cell to its initial value when backtracking
				sol_mat[next_cord.first][next_cord.second] = -2;
			}
		}
	}

	return false;
}


// should work ? 
bool solution_multi_thread_solver(int** mat, int** sol_mat, vector<sum>& sums, int& m, int& n) {
	bool found = false;
	int num_of_remaining_threads = omp_get_max_threads();
	
	// max 100 threds can ben inpruved
	#pragma omp parallel for num_threads(num_of_remaining_threads) shared(mat, sol_mat, sums, m, n, found) collapse(2)
	for (int i = 1; i < 10; i++) {
		for (int j = 1; j < 10; j++) {
			if (!found) {
				int** local_sol_mat = new int* [m];
				for (int k = 0; k < m; k++) {
					local_sol_mat[k] = new int[n];
					memcpy(local_sol_mat[k], sol_mat[k], n * sizeof(int));
				}
				vector<sum> copy_of_susm = deep_copy(sums);

				// this is worng 
				COORD startign_chagne = copy_of_susm[0].start;
				local_sol_mat[startign_chagne.first][startign_chagne.second] = i;
				if (copy_of_susm[0].dir == d_down) {
					local_sol_mat[startign_chagne.first + 1][startign_chagne.second] = j;
				}
				else {
					local_sol_mat[startign_chagne.first][startign_chagne.second + 1] = j;
				}

				if (solution_single_thread_solver(mat, local_sol_mat, copy_of_susm, m, n)) {
				#pragma omp critical
					{
						found = true;
						for (int k = 0; k < m; k++) {
							memcpy(sol_mat[k], local_sol_mat[k], n * sizeof(int));
						}						
					}
				}

				for (int k = 0; k < m; k++) {
					delete[] local_sol_mat[k];
				}
				delete[] local_sol_mat;
			}
		}
	}
	
	return found;
}




bool solution(int** mat, int** sol_mat, vector<sum> sums, int m, int n) {

	//TO DO: Write the solution
	//You can use any algorithm and data type
	//Write your solution to file in main function using sol_to_mat() after solving it	
	sort(sums.begin(), sums.end(), sum_is_a_smaller);
	/*
	sol_mat[1][1] = 6; //6313
	sol_mat[1][2] = 3;
	sol_mat[2][1] = 1;
	sol_mat[2][2] = 3;
	cout << "is valid :" << is_all_sums_valid(sums, sol_mat) << endl;
	return false;
	*/
	// Have to deep copy the sums as well 
	cout << "#########################" << endl;
	print_one_matrix(sol_mat, m, n);
	cout << "#########################" << endl;
	print_one_matrix(mat, m, n);
	cout << "##########################" << endl;




	int** sol_copy = deep_copy(sol_mat, m, n);
	vector<sum> sums_copy = deep_copy(sums);
	auto start = chrono::high_resolution_clock::now();
	bool got = solution_multi_thread_solver(mat, sol_copy, sums_copy, m, n);	
	auto end = chrono::high_resolution_clock::now();
	cout << "solution_multi_thread_solver execution time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "micro seconds" << endl;
	cout << "solution got: " << got << endl;
	print_one_matrix(sol_copy,m,n);
	/*
	print_one_matrix(sol_mat, m, n );
	start = chrono::high_resolution_clock::now();
	got = solution_single_thread_solver(mat, sol_mat, sums, m, n) ;
	end = chrono::high_resolution_clock::now();
	cout << "solution_single_thread_solver execution time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "micro seconds" << endl;
	cout << "solution got: " << got << endl;
	*/
	
	return got;

}


int main(int argc, char** argv) {

	std::string filename(argv[1]);
	std::ifstream file;
	file.open(filename.c_str());

	int m, n;
	file >> m;
	file >> n;

	int** mat;
	read_matrix(mat, file, m, n);
	print_one_matrix(mat, m, n);

	int** sol_mat;
	convert_sol(mat, sol_mat, m, n);
	print_one_matrix(sol_mat, m, n);

	vector<sum> sums = get_sums(mat, m, n);

	// get very persice time 
	solution(mat, sol_mat, sums, m, n);

	print_one_matrix(sol_mat, m, n);
	sol_to_file(mat, sol_mat, m, n, "solution.kakuro");



	delete mat;
	delete sol_mat;

	return 0;
}
