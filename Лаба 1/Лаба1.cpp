#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

typedef double real;

class FDM {
public:
	int nx, ny;
	real kx, ky;
	real x0, xm, xe;
	real y0, ym, ye;
	real hx, hy;

	int bottom, top, left, right;

	real lambda = 1;
	real gamma = 2;
	real w = 0.9;
	real e = 10e-20;
	int maxiter = 1000;

	vector<vector<real>> grid;

	vector<real> x;
	vector<real> xk;
	vector<real> f;
	vector<real> d1, d2, d3, d4, d5;

	void input();
	void create_grid();
	void create_matrix();
	void create_vector();

	void first_condition_bottom();
	void first_condition_top();
	void first_condition_left();
	void first_condition_right();

	void second_condition_bottom();
	void second_condition_top();
	void second_condition_left();
	void second_condition_right();

	void null_nodes();

	void output();

	real func_f(real x, real y);
	real der_ux(real x, real y);
	real der_uy(real x, real y);
	real u(real x, real y);

	real Ax_before(int i, vector<real> x);
	real Ax_after(int i, vector<real> x);
	real discrepancy();
	void iteration_step(vector<real> y, vector<real> yk);
	void Jacobi_method();
	void Gauss_Seidel_method();

	real error_count();
};


real FDM::u(real x, real y) {
	//return 5;
	//return x + y;
	return x * x + y * y;
	//return x * x * x + y * y * y;
}


real FDM::der_ux(real x, real y) {
	//return 0;
	//return 1;
	return 2 * x;
	//return -lambda * 6 + gamma * 3 * x * x;
}


real FDM::der_uy(real x, real y) {
	//return 0;
	//return 1;
	return 2 * y;
	//return 3 * y * y;
}


real FDM::func_f(real x, real y) {
	//return 5 * gamma;
	//return gamma * (x + y);
	return -4 * lambda + gamma * (x * x + y * y);
	//return -lambda * (6 * x + 6 * y) + gamma * (x * x * x + y * y * y);
}


void FDM::input() {
	ifstream input("input.txt");
	input >> x0 >> xm >> xe;
	input >> y0 >> ym >> ye;
	input >> nx >> kx;
	input >> ny >> ky;
	input.close();

	ifstream boundary_conditions("boundary_conditions.txt");
	boundary_conditions >> bottom >> top >> left >> right;
	boundary_conditions.close();
}


void FDM::create_grid() {
	if (kx == 1) hx = (xe - x0) / nx;
	else hx = (xe - x0) * (kx - 1) / (pow(kx, nx) - 1);
	
	if (ky == 1) hy = (ye - y0) / ny;
	else hy = (ye - y0) * (ky - 1) / (pow(ky, ny) - 1);

	grid.resize((nx + 1) * (ny + 1));

	for (int i = 0; i < grid.size(); i++) grid[i].resize(2);

	real x = x0;
	real y = y0;

	for (int i = 0; i <= ny; i++) {
		x = x0;
		grid[(nx + 1) * i][0] = x;
		grid[(nx + 1) * i][1] = y;
		for (int j = 1; j <= nx; j++) {
			x += pow(kx, j - 1) * hx;
			grid[(nx + 1) * i + j][0] = x;
			grid[(nx + 1) * i + j][1] = y;
		}
		y += pow(ky, i) * hy;
	}
}


void FDM::create_matrix() {
	d1.resize((nx + 1) * (ny + 1), 0);
	d5.resize((nx + 1) * (ny + 1), 0);
	d2.resize((nx + 1) * (ny + 1), 0);
	d4.resize((nx + 1) * (ny + 1), 0);
	d3.resize((nx + 1) * (ny + 1), 0);

	int indx = 0;
	while (grid[indx][0] < xm) indx++;

	int indy = 0, k = 0;
	while (grid[k][1] < ym) {
		k += nx + 1;
		indy++;
	}

	for (int i = nx + 2; i <= indy * (nx + 1) + indx; i++) {
		if (!(i % (nx + 1) == 0) && !(i % (nx + 1) == nx)) {
			d1[i] = -2 * lambda / (grid[i][1] - grid[i - nx - 1][1]) / (grid[i + nx + 1][1] - grid[i - nx - 1][1]);
			d2[i] = -2 * lambda / (grid[i][0] - grid[i - 1][0]) / (grid[i + 1][0] - grid[i - 1][0]);
			d3[i] = 2 * lambda * ((1 / (grid[i][0] - grid[i - 1][0]) / (grid[i + 1][0] - grid[i][0])) + (1 / (grid[i][1] - grid[i - nx - 1][1]) / (grid[i + nx + 1][1] - grid[i][1]))) + gamma;
			d4[i + 1] = -2 * lambda / (grid[i + 1][0] - grid[i][0]) / (grid[i + 1][0] - grid[i - 1][0]);
			d5[i + nx + 1] = -2 * lambda / (grid[i + nx + 1][1] - grid[i][1]) / (grid[i + nx + 1][1] - grid[i - nx - 1][1]);
		}
	}

	for (int i = indy + 1; i < ny; i++) {
		for (int j = i * (nx + 1) + 1; j < i * (nx + 1) + indx; j++) {
			d1[j] = -2 * lambda / (grid[j][1] - grid[j - nx - 1][1]) / (grid[j + nx + 1][1] - grid[j - nx - 1][1]);
			d2[j] = -2 * lambda / (grid[j][0] - grid[j - 1][0]) / (grid[j + 1][0] - grid[j - 1][0]);
			d3[j] = 2 * lambda * ((1 / (grid[j][0] - grid[j - 1][0]) / (grid[j + 1][0] - grid[j][0])) + (1 / (grid[j][1] - grid[j - nx - 1][1]) / (grid[j + nx + 1][1] - grid[j][1]))) + gamma;
			d4[j + 1] = -2 * lambda / (grid[j + 1][0] - grid[j][0]) / (grid[j + 1][0] - grid[j - 1][0]);
			d5[j + nx + 1] = -2 * lambda / (grid[j + nx + 1][1] - grid[j][1]) / (grid[j + nx + 1][1] - grid[j - nx - 1][1]);
			
		}
	}

	if(bottom == 1) first_condition_bottom();
	else second_condition_bottom();

	if (top == 1) first_condition_top();
	else second_condition_top();

	if (right == 1) first_condition_right();
	else second_condition_right();

	if (left == 1) first_condition_left();
	else second_condition_left();
}


void FDM::create_vector() {
	f.resize((nx + 1) * (ny + 1));
	x.resize((nx + 1) * (ny + 1));
	xk.resize((nx + 1) * (ny + 1));
	for (int i = 0; i < grid.size(); i++) {
		f[i] = func_f(grid[i][0], grid[i][1]);
		x[i] = 1;
		xk[i] = 1;
	}
}


void FDM::first_condition_bottom() {
	for (int i = 0; i <= nx; i++) { 
		d3[i] = 1;
		f[i] = u(grid[i][0], grid[i][1]);
	}
}


void FDM::first_condition_top() {
	for (int i = ny * (nx + 1); i < grid.size(); i++) {
		if (grid[i][0] <= xm) {
			d3[i] = 1;
			f[i] = u(grid[i][0], grid[i][1]);
		}
	}

	int indx = 0;
	while (grid[indx][0] < xm) indx++;

	int indy = 0, k = 0;
	while (grid[k][1] < ym) {
		k += nx + 1;
		indy++;
	}

	for(int i = indy * (nx + 1) + indx + 1; i < (indy + 1) * (nx + 1); i++)
	{
		d3[i] = 1;
		f[i] = u(grid[i][0], grid[i][1]);
	}
}


void FDM::first_condition_left() {
	for (int i = 0; i <= ny; i++) {
		d3[i * (nx + 1)] = 1;
		f[i * (nx + 1)] = u(grid[i * (nx + 1)][0], grid[i * (nx + 1)][1]);
	}
}


void FDM::first_condition_right() {
	for (int i = 0; i <= ny; i++) {
		if (grid[i * (nx + 1) + nx][1] <= ym) {
			d3[i * (nx + 1) + nx] = 1;
			f[i * (nx + 1) + nx] = u(grid[i * (nx + 1) + nx][0], grid[i * (nx + 1) + nx][1]);
		}
	}

	int indx = 0;
	while (grid[indx][0] < xm) indx++;

	int indy = 0, k = 0;
	while (grid[k][1] < ym) { 
		k += nx + 1; 
		indy++;
	}

	for (int i = (indy + 1) * (nx + 1) + indx; i < grid.size(); i += nx + 1) {
		d3[i] = 1;
		f[i] = u(grid[i][0], grid[i][1]);
	}
}


void FDM::second_condition_bottom() {
	if (left == 2) {
		f[0] = -lambda * der_uy(grid[0][0], grid[0][1]);
		d5[nx + 1] = -lambda / hy;
		d3[0] = lambda / hy;
	}

	for (int i = 1; i < nx; i++) {
		f[i] = -lambda * der_uy(grid[i][0], grid[i][1]);
		d5[i + nx + 1] = -lambda / hy;
		d3[i] = lambda / hy;
	}

	if (right == 2) {
		f[nx] = -lambda * der_uy(grid[nx][0], grid[nx][1]);
		d5[2 * nx + 1] = -lambda / hy;
		d3[nx] = lambda / hy;
	}
}


void FDM::second_condition_top() {
	if (left == 2) {
		f[(nx + 1) * ny] = lambda * der_uy(grid[(nx + 1) * ny][0], grid[(nx + 1) * ny][1]);
		d1[(nx + 1) * ny] = -lambda / (grid[(nx + 1) * ny][1] - grid[(nx + 1) * (ny - 1)][1]);
		d3[(nx + 1) * ny] = lambda / (grid[(nx + 1) * ny][1] - grid[(nx + 1) * (ny - 1)][1]);
	}

	int indx = 0;
	while (grid[indx][0] < xm) indx++;

	int indy = 0, k = 0;
	while (grid[k][1] < ym) {
		k += nx + 1;
		indy++;
	}

	for (int i = (nx + 1) * ny + 1; i <= (nx + 1) * ny + indx - 1; i++) {
		f[i] = lambda * der_uy(grid[i][0], grid[i][1]);
		d1[i] = -lambda / (grid[i][1] - grid[i - nx - 1][1]);
		d3[i] = lambda / (grid[i][1] - grid[i - nx - 1][1]);
	}

	if (right == 2) {
		f[(nx + 1) * ny + indx] = lambda * der_uy(grid[(nx + 1) * ny + indx][0], grid[(nx + 1) * ny + indx][1]);
		d1[(nx + 1) * ny + indx] = -lambda / (grid[(nx + 1) * ny + indx][1] - grid[(nx + 1) * (ny - 1) + indx][1]);
		d3[(nx + 1) * ny + indx] = lambda / (grid[(nx + 1) * ny + indx][1] - grid[(nx + 1) * (ny - 1) + indx][1]);

		f[(nx + 1) * (indy + 1) - 1] = lambda * der_uy(grid[(nx + 1) * (indy + 1) - 1][0], grid[(nx + 1) * (indy + 1) - 1][1]);
		d1[(nx + 1) * (indy + 1) - 1] = -lambda / (grid[(nx + 1) * (indy + 1) - 1][1] - grid[(nx + 1) * indy - 1][1]);
		d3[(nx + 1) * (indy + 1) - 1] = lambda / (grid[(nx + 1) * (indy + 1) - 1][1] - grid[(nx + 1) * indy - 1][1]);
	}

	for (int i = indy * (nx + 1) + indx + 1; i < (indy + 1) * (nx + 1) - 1; i++) {
		f[i] = lambda * der_uy(grid[i][0], grid[i][1]);
		d1[i] = -lambda / (grid[i][1] - grid[i - nx - 1][1]);
		d3[i] = lambda / (grid[i][1] - grid[i - nx - 1][1]);
	}
}


void FDM::second_condition_left() {
	/*if (bottom == 2) {
		f[0] = -lambda * der_fx(grid[0][0], grid[0][1]);
		d4[1] = -lambda / hx;
		d3[0] = lambda / hx;
	}*/

	for (int i = nx + 1; i < (nx + 1) * ny; i += nx + 1) {
		f[i] = -lambda * der_ux(grid[i][0], grid[i][1]);
		d4[i + 1] = -lambda / hx;
		d3[i] = lambda / hx;
	}

	/*if (top == 2) {
		f[(nx + 1) * ny] = -lambda * der_fx(grid[(nx + 1) * ny][0], grid[(nx + 1) * ny][1]);
		d4[(nx + 1) * ny + 1] = -lambda / hx;
		d3[(nx + 1) * ny] = lambda / hx;
	}*/
}


void FDM::second_condition_right() {
	/*if (bottom == 2) {
		f[nx] = lambda * der_fx(grid[nx][0], grid[nx][1]);
		d2[nx] = -lambda / (grid[nx][0] - grid[nx - 1][0]);
		d3[nx] = lambda / (grid[nx][0] - grid[nx - 1][0]);
	}*/

	int indx = 0;
	while (grid[indx][0] < xm) indx++;

	int indy = 0, k = 0;
	while (grid[k][1] < ym) {
		k += nx + 1;
		indy++;
	}

	/*if (top == 2) {
		f[(indy + 1) * (nx + 1) - 1] = lambda * der_fx(grid[(indy + 1) * (nx + 1) - 1][0], grid[(indy + 1) * (nx + 1) - 1][1]);
		d2[(indy + 1) * (nx + 1) - 1] = -lambda / (grid[(indy + 1) * (nx + 1) - 1][0] - grid[(indy + 1) * (nx + 1) - 2][0]);
		d3[(indy + 1) * (nx + 1) - 1] = lambda / (grid[(indy + 1) * (nx + 1) - 1][0] - grid[(indy + 1) * (nx + 1) - 2][0]);

		f[ny * (nx + 1) + indx] = lambda * der_fx(grid[ny * (nx + 1) + indx][0], grid[ny * (nx + 1) + indx][1]);
		d2[ny * (nx + 1) + indx] = -lambda / (grid[ny * (nx + 1) + indx][0] - grid[ny * (nx + 1) + indx - 1][0]);
		d3[ny * (nx + 1) + indx] = lambda / (grid[ny * (nx + 1) + indx][0] - grid[ny * (nx + 1) + indx - 1][0]);
	}*/

	for (int i = 2 * nx + 1; i < (indy + 1) * (nx + 1) - 1; i += nx + 1) {
		f[i] = lambda * der_ux(grid[i][0], grid[i][1]);
		d2[i] = -lambda / (grid[i][0] - grid[i - 1][0]);
		d3[i] = lambda / (grid[i][0] - grid[i - 1][0]);
	}

	for (int i = (indy + 1) * (nx + 1) + indx; i < ny * (nx + 1) + indx; i += nx + 1) {
		f[i] = lambda * der_ux(grid[i][0], grid[i][1]);
		d2[i] = -lambda / (grid[i][0] - grid[i - 1][0]);
		d3[i] = lambda / (grid[i][0] - grid[i - 1][0]);
	}
}


void FDM::null_nodes() {
	for (int i = 0; i < grid.size(); i++) {
		if (grid[i][0] > xm && grid[i][0] <= xe && grid[i][1] > ym && grid[i][1] <= ye) {
			f[i] = 0;
			d3[i] = 1;
		}
	}
}


real FDM::Ax_before(int i, vector<real> x) {

	real s = 0;
	int m = nx - 1;
	int n = grid.size();

	if (i == 0) {
		return s;
	}

	if (i == n - 1) {
		s += d1[i] * x[i - 2 - m];
		s += d2[i] * x[i - 1];
		return s;
	}

	if (0 < i < n - 1) {
		if (i - m - 2 >= 0) s += d1[i] * x[i - m - 2];
		s += d2[i] * x[i - 1];
		return s;
	}
}


real FDM::Ax_after(int i, vector<real> x) {

	real s = 0;
	int m = nx - 1;
	int n = grid.size();

	if (i == 0) {
		s += d3[i] * x[i];
		s += d4[i + 1] * x[i + 1];
		s += d5[i + m + 2] * x[i + m + 2];
		return s;
	}

	if (i == n - 1) {
		s += d3[i] * x[i];
		return s;
	}

	if (0 < i < n - 1) {
		s += d3[i] * x[i];
		s += d4[i + 1] * x[i + 1];
		if (i + m + 2 < n) s += d5[i + m + 2] * x[i + m + 2];
		return s;
	}
}


real FDM::discrepancy() {
	int n = grid.size();
	real norm_f = 0;
	for (int i = 0; i < n; i++) norm_f += f[i] * f[i];

	real norm_nev = 0;
	for (int i = 0; i < n; i++) norm_nev += (f[i] - Ax_before(i, x) - Ax_after(i, x)) * (f[i] - Ax_before(i, x) - Ax_after(i, x));

	return sqrt(norm_nev) / sqrt(norm_f);
}


void FDM::iteration_step(vector<real> yk, vector<real> y) {
	int n = grid.size();
	for (int i = 0; i < n; i++) {
		real sum = w * (f[i] - Ax_before(i, yk) - Ax_after(i, y)) / d3[i];
		xk[i] = x[i] + sum;
	}
}


void FDM::Jacobi_method() {
	int iteration = 0;
	int n = grid.size();
	while (iteration <= maxiter && discrepancy() >= e) {
		iteration_step(x, x);

		for (int i = 0; i < n; i++) x[i] = xk[i];
		cout << "Iteration number: " << iteration << endl;
		cout << "Discrepancy value: " << discrepancy() << endl << endl;
		iteration++;
	}
}


void FDM::Gauss_Seidel_method() {
	int iteration = 0;
	int n = grid.size();
	while (iteration <= maxiter && discrepancy() >= e) {
		iteration_step(xk, x);

		for (int i = 0; i < n; i++) x[i] = xk[i];
		/*cout << "Iteration number: " << iteration << endl;
		cout << "Discrepancy value: " << discrepancy() << endl << endl;*/
		iteration++;
	}
}


real FDM::error_count() {
	real sum = 0;  
	for (int i = 0; i < grid.size(); i++) {
		if (!(grid[i][0] > xm && grid[i][0] <= xe && grid[i][1] > ym && grid[i][1] <= ye))
			sum += (x[i] - u(grid[i][0], grid[i][1])) * (x[i] - u(grid[i][0], grid[i][1]));
	}
	return sqrt(sum);
}


void FDM::output() {
	for (int i = 0; i < grid.size(); i++) {
		cout << i << " : " << grid[i][0] << ", " << grid[i][1];
		cout << endl;
	}

	cout << endl;

	cout << "f: " << endl;
	for (int i = 0; i < f.size(); i++) {
		cout << i << " : " << f[i];
		cout << endl;
	}

	cout << endl;
	cout << "d1: " << endl;
	for (int i = 0; i < d1.size(); i++) {
		cout << i << " : " << d1[i];
		cout << endl;
	}

	cout << endl;

	cout << "d2: " << endl;
	for (int i = 0; i < d2.size(); i++) {
		cout << i << " : " << d2[i];
		cout << endl;
	}

	cout << endl;

	cout << "d3: " << endl;
	for (int i = 0; i < d3.size(); i++) {
		cout << i << " : " << d3[i];
		cout << endl;
	}

	cout << endl;

	cout << "d4: " << endl;
	for (int i = 0; i < d3.size(); i++) {
		cout << i << " : " << d4[i];
		cout << endl;
	}

	cout << endl;

	cout << "d5: " << endl;
	for (int i = 0; i < d3.size(); i++) {
		cout << i << " : " << d5[i];
		cout << endl;
	}


	ofstream output("output.txt");
	for (int i = 0; i < grid.size(); i++) {
		if (!(grid[i][0] > xm && grid[i][0] <= xe && grid[i][1] > ym && grid[i][1] <= ye)) {
			cout << fixed;
			cout << setprecision(14) << "[" << i << "]: " << x[i] << "		" << u(grid[i][0], grid[i][1]) << "		" << abs(x[i] - u(grid[i][0], grid[i][1])) << endl;
			output << fixed;
			output << setprecision(14) << x[i] << endl;
		}
	}
	output.close();

	cout << "Error value = " << setprecision(14) << error_count();
}


int main() {
	FDM m;
	m.input();
	m.create_grid();
	m.create_vector();
	m.create_matrix();
	m.null_nodes();
	m.Gauss_Seidel_method();
	m.output();
}