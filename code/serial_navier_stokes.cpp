#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

const int Nx = 201;
const int Ny = 101;

const double Lx = 0.1, Ly = 0.05;
const double rho = 1000, nu = 1e-6;
const double P_max = 0.5;
const double t_end = 50.0;
const double dt_min = 1.e-3;
const double courant = 0.01;
const double dt_out = 0.5;

vector<vector<double>> P, P_old, u, u_old, v, v_old, PPrhs;
double dx, dy, dt, t;

void grids_to_file(int out)
{
	//Write the output for a single time step to file
	stringstream fname;
	fstream f1;
	fname << "./out/P" << "_" << out << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
			f1 << P[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
	fname.str("");
	fname << "./out/u" << "_" << out << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
			f1 << v[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
	fname.str("");
	fname << "./out/u" << "_" << out << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
			f1 << v[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
}

void setup(void)
{
	P.resize(Nx, vector<double>(Ny, 0.0));
	P_old.resize(Nx, vector<double>(Ny, 0.0));
	u.resize(Nx, vector<double>(Ny, 0.0));
	u_old.resize(Nx, vector<double>(Ny, 0.0));
	v.resize(Nx, vector<double>(Ny, 0.0));
	v_old.resize(Nx, vector<double>(Ny, 0.0));
	PPrhs.resize(Nx, vector<double>(Ny, 0.0));

	dx = Lx / (Nx - 1);
	dy = Ly / (Ny - 1);

	for (int j = 0; j < Ny; j++)
		P[0][j] = P_max;

	P_old = P;

	t = 0.0;
}

void calculate_ppm_RHS_central(void)
{
	for (int i = 1; i < Nx - 1; i++)
		for (int j = 1; j < Ny - 1; j++)
		{
			PPrhs[i][j] = rho / dt * ((u[i + 1][j] - u[i - 1][j]) / (2. * dx) + (v[i][j + 1] - v[i][j - 1]) / (2. * dy));
		}
}

void set_pressure_BCs(void)
{
	for (int i = 0; i < Nx; i++)
	{
		P[i][0] = P[i][1];
		P[i][Ny - 1] = P[i][Ny - 2];
	}

	for (int j = Ny / 2; j < Ny; j++)
		P[Nx - 1][j] = P[Nx - 2][j];
}

int pressure_poisson_jacobi(double rtol = 1.e-5)
{
	double tol = 10. * rtol;
	int it = 0;

	//swap P and P_old without copying the data - quick
	swap(P, P_old);

	while (tol > rtol)
	{
		double sum_val = 0.0;
		tol = 0.0;
		it++;

		//Jacobi iteration
		for (int i = 1; i < Nx - 1; i++)
			for (int j = 1; j < Ny - 1; j++)
			{
				P[i][j] = 1.0 / (2.0 + 2.0 * (dx * dx) / (dy * dy)) * (P_old[i + 1][j] + P_old[i - 1][j] + 
					(P_old[i][j + 1] + P_old[i][j - 1]) * (dx * dx) / (dy * dy)
					- (dx * dx) * PPrhs[i][j]);

				sum_val += fabs(P[i][j]);
				tol += fabs(P[i][j] - P_old[i][j]);
			}

		set_pressure_BCs();

		tol = tol / max(1.e-10, sum_val);

		swap(P, P_old);
	}

	return it;
}

void calculate_intermediate_velocity(void)
{
	for (int i = 1; i < Nx - 1; i++)
		for (int j = 1; j < Ny - 1; j++)
		{
			//viscous diffusion
			u[i][j] = u_old[i][j] + dt * nu * ((u_old[i + 1][j] + u_old[i - 1][j] - 2.0 * u_old[i][j]) / (dx * dx) + (u_old[i][j + 1] + u_old[i][j - 1] - 2.0 * u_old[i][j]) / (dy * dy));
			v[i][j] = v_old[i][j] + dt * nu * ((v_old[i + 1][j] + v_old[i - 1][j] - 2.0 * v_old[i][j]) / (dx * dx) + (v_old[i][j + 1] + v_old[i][j - 1] - 2.0 * v_old[i][j]) / (dy * dy));
			//advection - upwinding
			if (u[i][j] > 0.0)
			{
				u[i][j] -= dt * u_old[i][j] * (u_old[i][j] - u_old[i - 1][j]) / dx;
				v[i][j] -= dt * u_old[i][j] * (v_old[i][j] - v_old[i - 1][j]) / dx;
			}
			else
			{
				u[i][j] -= dt * u_old[i][j] * (u_old[i + 1][j] - u_old[i][j]) / dx;
				v[i][j] -= dt * u_old[i][j] * (v_old[i + 1][j] - v_old[i][j]) / dx;
			}
			if (v[i][j] > 0.0)
			{
				u[i][j] -= dt * v_old[i][j] * (u_old[i][j] - u_old[i][j - 1]) / dy;
				v[i][j] -= dt * v_old[i][j] * (v_old[i][j] - v_old[i][j - 1]) / dy;
			}
			else
			{
				u[i][j] -= dt * v_old[i][j] * (u_old[i][j + 1] - u_old[i][j]) / dy;
				v[i][j] -= dt * v_old[i][j] * (v_old[i][j + 1] - v_old[i][j]) / dy;
			}
		}
}

void set_velocity_BCs(void)
{
	for (int j = 0; j < Ny; j++)
		u[0][j] = u[1][j];

	for (int j = 0; j < Ny / 2; j++)
		u[Nx - 1][j] = u[Nx - 2][j];
}

double project_velocity(void)
{
	double vmax = 0.0;
	for (int i = 1; i < Nx - 1; i++)
		for (int j = 1; j < Ny - 1; j++)
		{
			u[i][j] = u[i][j] - dt * (1. / rho) * (P[i + 1][j] - P[i - 1][j]) / (2. * dx);
			v[i][j] = v[i][j] - dt * (1. / rho) * (P[i][j + 1] - P[i][j - 1]) / (2. * dy);

			double vel = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);

			vmax = max(vmax, vel);
		}

	set_velocity_BCs();

	return vmax;
}


void solve_NS(void)
{
	double vel_max = 0.0;
	int time_it = 0;
	int its;
	int out_it = 0;
	double t_out = dt_out;

	grids_to_file(out_it);

	while (t < t_end)
	{
		if (vel_max > 0.0)
		{
			dt = min(courant * min(dx, dy) / vel_max, dt_min);
		}
		else dt = dt_min;

		t += dt;
		time_it++;
		swap(u, u_old);
		swap(v, v_old);

		calculate_intermediate_velocity();
		calculate_ppm_RHS_central();
		its = pressure_poisson_jacobi(1.e-5);
		vel_max = project_velocity();

		if (t >= t_out)
		{
			out_it++;
			t_out += dt_out;
			cout << time_it << ": " << t << " Jacobi iterations: " << its << " vel_max: " << vel_max << endl;
			grids_to_file(out_it);
		}
	}
}

int main()
{
	setup();
	solve_NS();
	return 0;
}
