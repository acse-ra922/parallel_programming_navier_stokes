#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include "mpi.h"

using namespace std;

//defining the domain, boundary and initial conditions
const int Nx = 201;
const int Ny = 101;

const double Lx = 0.1, Ly = 0.05;
const double rho = 1000, nu = 1e-6;
const double P_max = 0.5;
const double t_end = 50.0;
const double dt_min = 1.e-3;
const double courant = 0.01;
const double dt_out = 0.5;


class Cmatrix
{
public:
	double* mat_1D;
	double** mat_2D;
	int n, m;

	Cmatrix()
	{
		mat_1D = nullptr;
		mat_2D = nullptr;
	}
	Cmatrix(int imax, int jmax)
	{
		n = imax;
		m = jmax;
		mat_1D = new double[n * m];
		mat_2D = new double* [n];
		for (int i = 0; i < n; i++)
			mat_2D[i] = &mat_1D[i * m];

		for (int i = 0; i < n * m; i++)
			mat_1D[i] = 0.0;
	}

	~Cmatrix()
	{
		delete[] mat_1D;
		delete[] mat_2D;
	}
};

class Cvel
{
public:
	double u, v;

	Cvel()
	{
		u = 0.0;
		v = 0.0;
	}

	static void buildMPIType();

	static MPI_Datatype MPI_Type;
};

MPI_Datatype Cvel::MPI_Type;

void Cvel::buildMPIType()
{
	int block_lengths[2];
	MPI_Aint offsets[2];
	MPI_Aint addresses[2], add_start;
	MPI_Datatype typelist[2];

	Cvel temp;

	typelist[0] = MPI_DOUBLE;
	block_lengths[0] = 1;
	MPI_Get_address(&temp.u, &addresses[0]);

	typelist[1] = MPI_DOUBLE;
	block_lengths[1] = 1;
	MPI_Get_address(&temp.v, &addresses[1]);

	MPI_Get_address(&temp, &add_start);
	for (int i = 0; i < 2; i++) offsets[i] = addresses[i] - add_start;

	MPI_Type_create_struct(2, block_lengths, offsets, typelist, &MPI_Type);
	MPI_Type_commit(&MPI_Type);
}

class Cmatrix_vel
{
public:
	Cvel* mat_1D;
	Cvel** mat_2D;
	int n, m;

	Cmatrix_vel()
	{
		mat_1D = nullptr;
		mat_2D = nullptr;
	}
	Cmatrix_vel(int imax, int jmax)
	{
		n = imax;
		m = jmax;
		mat_1D = new Cvel[n * m];
		mat_2D = new Cvel * [n];
		for (int i = 0; i < n; i++)
			mat_2D[i] = &mat_1D[i * m];
	}

	~Cmatrix_vel()
	{
		delete[] mat_1D;
		delete[] mat_2D;
	}
};

Cmatrix* P, * P_old, * PPrhs;
Cmatrix_vel* vel, * vel_old;
double dx, dy, dt, t;

int id, p;
int tag_num;

vector<int> num_i, start_i;

void grids_to_file(int out)
{
	//Write the output for a single time step to file
	stringstream fname;
	fstream f1;
	fname << "./out/P" << "_" << out << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 1; i < num_i[id] + 1; i++)
	{
		for (int j = 0; j < Ny; j++)
			f1 << P->mat_2D[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
	fname.str("");
	fname << "./out/u" << "_" << out << "_" << id << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 1; i < num_i[id] + 1; i++)
	{
		for (int j = 0; j < Ny; j++)
			f1 << vel->mat_2D[i][j].u << "\t";
		f1 << endl;
	}
	f1.close();
	fname.str("");
	fname << "./out/v" << "_" << out << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 1; i < num_i[id] + 1; i++)
	{
		for (int j = 0; j < Ny; j++)
			f1 << vel->mat_2D[i][j].v << "\t";
		f1 << endl;
	}
	f1.close();
}

void setup(void)
{
	//Set up communications
	tag_num = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	Cvel::buildMPIType();

	int i_rem = Nx;
	int i_s = 0;
	num_i.resize(p);
	start_i.resize(p);

	for (int n = 0; n < p; n++)
	{
		start_i[n] = i_s;
		num_i[n] = i_rem / (p - n);
		i_rem -= num_i[n];
		i_s += num_i[n];
	}

	P = new Cmatrix(num_i[id] + 2, Ny);
	P_old = new Cmatrix(num_i[id] + 2, Ny);
	PPrhs = new Cmatrix(num_i[id] + 2, Ny);

	vel = new Cmatrix_vel(num_i[id] + 2, Ny);
	vel_old = new Cmatrix_vel(num_i[id] + 2, Ny);

	dx = Lx / (Nx - 1);
	dy = Ly / (Ny - 1);

	if (id == 0)
		for (int j = 0; j < Ny; j++)
			P->mat_2D[1][j] = P_max;

	for (int i = 0; i < num_i[id] + 2 * Ny; i++)
		P_old->mat_1D[i] = P->mat_1D[i];

	t = 0.0;
}

void calculate_ppm_RHS_central(void)
{
	//
	for (int i = (id == 0 ? 2 : 1); i < (id == p-1 ? num_i[id] : num_i[id] + 1); i++)
		for (int j = 1; j < Ny - 1; j++)
		{
			PPrhs->mat_2D[i][j] = rho / dt * ((vel->mat_2D[i + 1][j].u - vel->mat_2D[i - 1][j].u) / (2. * dx)
				+ (vel->mat_2D[i][j + 1].v - vel->mat_2D[i][j - 1].v) / (2. * dy));
		}
}

void set_pressure_BCs(void)
{
	//This function sets pressure boundary conditions
	//Impose Neumann Boundary Conditions on relevant boundaries if they are present in the domain
	for (int i = 1; i < num_i[id] + 1; i++)
	{
		P->mat_2D[i][0] = P->mat_2D[i][1];
		P->mat_2D[i][Ny - 1] = P->mat_2D[i][Ny - 2];
	}

	if (id == p-1)
		for (int j = Ny / 2; j < Ny; j++)
			P->mat_2D[num_i[id]][j] = P->mat_2D[num_i[id] - 1][j];
}

int pressure_poisson_jacobi(double rtol = 1.e-5)
{
	//
	double resid_vals[2];
	int it = 0;
	MPI_Request requests[4];
	int req_cnt;

	resid_vals[0] = rtol * 10;

	while (resid_vals[0] > rtol)
	{
		swap(P, P_old);
		resid_vals[0] = 0.0;
		resid_vals[1] = 0.0;
		it++;

		//Jacobi iteration
		for (int i = (id == 0 ? 2 : 1); i < (id == p - 1 ? num_i[id] : num_i[id] + 1); i++)
		{
			for (int j = 1; j < Ny - 1; j++)
			{
				P->mat_2D[i][j] = 1.0 / (2.0 + 2.0 * (dx * dx) / (dy * dy)) * (P_old->mat_2D[i + 1][j] + P_old->mat_2D[i - 1][j] +
					(P_old->mat_2D[i][j + 1] + P_old->mat_2D[i][j - 1]) * (dx * dx) / (dy * dy)
					- (dx * dx) * PPrhs->mat_2D[i][j]);

				resid_vals[1] += fabs(P->mat_2D[i][j]);
				resid_vals[0] += fabs(P->mat_2D[i][j] - P_old->mat_2D[i][j]);
			}
		}

		set_pressure_BCs();

		req_cnt = 0;
		//LHS
		if (id != 0)
		{
			MPI_Irecv(&P->mat_2D[0][0], Ny, MPI_DOUBLE, id - 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
			req_cnt++;
			MPI_Isend(&P->mat_2D[1][0], Ny, MPI_DOUBLE, id - 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
			req_cnt++;
		}
		//RHS
		if (id != p - 1)
		{
			MPI_Irecv(&P->mat_2D[num_i[id] + 1][0], Ny, MPI_DOUBLE, id + 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
			req_cnt++;
			MPI_Isend(&P->mat_2D[num_i[id]][0], Ny, MPI_DOUBLE, id + 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
			req_cnt++;
		}

		MPI_Waitall(req_cnt, requests, MPI_STATUSES_IGNORE);

		MPI_Allreduce(MPI_IN_PLACE, resid_vals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		resid_vals[0] = resid_vals[0] / max(1.e-10, resid_vals[1]);
	}

	return it;
}

void calculate_intermediate_velocity(void)
{
	MPI_Request requests[4];
	int req_cnt;

	for (int i = (id == 0 ? 2 : 1); i < (id == p - 1 ? num_i[id] : num_i[id] + 1); i++)
		for (int j = 1; j < Ny - 1; j++)
		{
			//viscous diffusion
			vel->mat_2D[i][j].u = vel_old->mat_2D[i][j].u + dt * nu * ((vel_old->mat_2D[i + 1][j].u + vel_old->mat_2D[i - 1][j].u
				- 2.0 * vel_old->mat_2D[i][j].u) / (dx * dx) + (vel_old->mat_2D[i][j + 1].u + vel_old->mat_2D[i][j - 1].u - 2.0 * vel_old->mat_2D[i][j].u) / (dy * dy));
			vel->mat_2D[i][j].v = vel_old->mat_2D[i][j].v + dt * nu * ((vel_old->mat_2D[i + 1][j].v + vel_old->mat_2D[i - 1][j].v
				- 2.0 * vel_old->mat_2D[i][j].v) / (dx * dx) + (vel_old->mat_2D[i][j + 1].v + vel_old->mat_2D[i][j - 1].v - 2.0 * vel_old->mat_2D[i][j].v) / (dy * dy));
			//advection - upwinding
			if (vel_old->mat_2D[i][j].u > 0.0)
			{
				vel->mat_2D[i][j].u -= dt * vel_old->mat_2D[i][j].u * (vel_old->mat_2D[i][j].u - vel_old->mat_2D[i - 1][j].u) / dx;
				vel->mat_2D[i][j].v -= dt * vel_old->mat_2D[i][j].u * (vel_old->mat_2D[i][j].v - vel_old->mat_2D[i - 1][j].v) / dx;
			}
			else
			{
				vel->mat_2D[i][j].u -= dt * vel_old->mat_2D[i][j].u * (vel_old->mat_2D[i+1][j].u - vel_old->mat_2D[i][j].u) / dx;
				vel->mat_2D[i][j].v -= dt * vel_old->mat_2D[i][j].u * (vel_old->mat_2D[i+1][j].v - vel_old->mat_2D[i][j].v) / dx;
			}

			if (vel_old->mat_2D[i][j].v > 0.0)
			{
				vel->mat_2D[i][j].u -= dt * vel_old->mat_2D[i][j].v * (vel_old->mat_2D[i][j].u - vel_old->mat_2D[i][j-1].u) / dx;
				vel->mat_2D[i][j].v -= dt * vel_old->mat_2D[i][j].v * (vel_old->mat_2D[i][j].v - vel_old->mat_2D[i][j-1].v) / dx;
			}
			else
			{
				vel->mat_2D[i][j].u -= dt * vel_old->mat_2D[i][j].v * (vel_old->mat_2D[i][j+1].u - vel_old->mat_2D[i][j].u) / dx;
				vel->mat_2D[i][j].v -= dt * vel_old->mat_2D[i][j].v * (vel_old->mat_2D[i][j+1].v - vel_old->mat_2D[i][j].v) / dx;
			}
		}

	req_cnt = 0;
	//LHS
	if (id != 0)
	{
		MPI_Irecv(&vel->mat_2D[0][0], Ny, Cvel::MPI_Type, id - 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
		MPI_Isend(&vel->mat_2D[1][0], Ny, Cvel::MPI_Type, id - 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
	}
	//RHS
	if (id != p - 1)
	{
		MPI_Irecv(&vel->mat_2D[num_i[id] + 1][0], Ny, Cvel::MPI_Type, id + 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
		MPI_Isend(&vel->mat_2D[num_i[id]][0], Ny, Cvel::MPI_Type, id + 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
	}

	MPI_Waitall(req_cnt, requests, MPI_STATUSES_IGNORE);
}

void set_velocity_BCs(void)
{
	//This function sets velocity boundary conditions
	//Impose Neumann Boundary Conditions on relevant boundaries if they are present in the domain
	if (id == 0)
		for (int j = 0; j < Ny; j++)
			vel->mat_2D[1][j].u = vel->mat_2D[2][j].u;

	if (id == p - 1)
		for (int j = 0; j < Ny / 2; j++)
			vel->mat_2D[num_i[id]][j].u = vel->mat_2D[num_i[id] - 1][j].u;
}

double project_velocity(void)
{
	MPI_Request requests[4];
	int req_cnt;
	double vmax = 0.0;

	for (int i = (id == 0 ? 2 : 1); i < (id == p - 1 ? num_i[id] : num_i[id] + 1); i++)
		for (int j = 1; j < Ny - 1; j++)
		{
			vel->mat_2D[i][j].u = vel->mat_2D[i][j].u - dt * (1. / rho) * (P->mat_2D[i + 1][j] - P->mat_2D[i - 1][j]) / (2. * dx);
			vel->mat_2D[i][j].v = vel->mat_2D[i][j].v - dt * (1. / rho) * (P->mat_2D[i][j + 1] - P->mat_2D[i][j - 1]) / (2. * dy);

			double vel_mag = sqrt(vel->mat_2D[i][j].u * vel->mat_2D[i][j].u + vel->mat_2D[i][j].v * vel->mat_2D[i][j].v);

			vmax = max(vmax, vel_mag);
		}

	set_velocity_BCs();

	req_cnt = 0;
	//LHS
	if (id != 0)
	{
		MPI_Irecv(&vel->mat_2D[0][0], Ny, Cvel::MPI_Type, id - 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
		MPI_Isend(&vel->mat_2D[1][0], Ny, Cvel::MPI_Type, id - 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
	}
	//RHS
	if (id != p - 1)
	{
		MPI_Irecv(&vel->mat_2D[num_i[id] + 1][0], Ny, Cvel::MPI_Type, id + 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
		MPI_Isend(&vel->mat_2D[num_i[id]][0], Ny, Cvel::MPI_Type, id + 1, tag_num, MPI_COMM_WORLD, &requests[req_cnt]);
		req_cnt++;
	}

	MPI_Waitall(req_cnt, requests, MPI_STATUSES_IGNORE);

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

		MPI_Allreduce(MPI_IN_PLACE, &dt, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		t += dt;
		time_it++;
		swap(vel, vel_old);

		calculate_intermediate_velocity();
		calculate_ppm_RHS_central();
		its = pressure_poisson_jacobi(1.e-5);
		vel_max = project_velocity();

		if (t >= t_out)
		{
			out_it++;
			t_out += dt_out;
			if (id == 0)
			{
				cout << time_it << ": " << t << " Jacobi iterations: " << its << " vel_max: " << vel_max << endl;
				cout.flush();
			}
			grids_to_file(out_it);
		}
	}
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	setup();
	solve_NS();


	MPI_Type_free(&Cvel::MPI_Type);
	MPI_Finalize();

	return 0;
}

