#include "Fluid3D.h"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#define MIN2(x, y) ( (x) < (y) ? (x) : (y) )
#define MAX2(x, y) ( (x) > (y) ? (x) : (y) )
#define MAX3(x, y, z) ( (x) > MAX2( (y), (z) ) ? (x) :  MAX2( (y), (z) ) )

inline void CLAMP (float &x, float a, float b)
{
	if ( x < a ) x = a;
	if ( x > b)  x = b;

}


Fluid3D::Fluid3D() : Fluid3D(32, 32, 32)
{
	m_numberOfJacobiIters = 40;
	m_beta  = 0.4f;
}


Fluid3D::Fluid3D(unsigned int nx, unsigned int ny, unsigned int nz)
{
	m_Nx = nx;
	m_Ny = ny;
	m_Nz = nz;
	m_NMax = MAX3(m_Nx, m_Ny, m_Nz);

	// voxel sizes(step size)
	// the solver assumes uniform cells, so the longest dimension is mapped to 1
	m_Dx = m_Dy = m_Dz = 1.0f / (float)m_NMax;

	// allocate 2 extra grid layer to take account the boundary conditions
	m_numCells = (m_Nx + 2) * (m_Ny + 2) * (m_Nz + 2);

	// value-initialization zeroes out all grids
	m_U0 = new float[m_numCells]();
	m_V0 = new float[m_numCells]();
	m_W0 = new float[m_numCells]();
	m_U1 = new float[m_numCells]();
	m_V1 = new float[m_numCells]();
	m_W1 = new float[m_numCells]();
	m_Fx = new float[m_numCells]();
	m_Fy = new float[m_numCells]();
	m_Fz = new float[m_numCells]();
	m_P0 = new float[m_numCells]();
	m_P1 = new float[m_numCells]();
	m_T0 = new float[m_numCells]();
	m_T1 = new float[m_numCells]();
	m_S0 = new float[m_numCells]();
	m_S1 = new float[m_numCells]();
	m_Div= new float[m_numCells]();

	m_dt    = 0.1f;
	m_numberOfJacobiIters = 20;
	m_alpha = 0.1f;
	m_beta  = 0.6f;
	m_Tamb  = 0.0f;
	m_visc  = 0.0001f;
	m_diff  = 0.0001f;

}


Fluid3D:: ~Fluid3D()
{
	// deallocate grids
	delete [] m_U0;
	delete [] m_V0;
	delete [] m_W0;
	delete [] m_U1;
	delete [] m_V1;
	delete [] m_W1;

	delete [] m_P0;
	delete [] m_P1;

	delete [] m_T0;
	delete [] m_T1;

	delete [] m_S0;
	delete [] m_S1;

	delete [] m_Div;

	delete [] m_Fx;
	delete [] m_Fy;
	delete [] m_Fz;

}

float Fluid3D::getCellValueAt(float *grid, int x, int y, int z)
{
	float value = -0.999999999f;

	// bounds check
	if (grid)
	{
		value = grid[IX(x, y, z)];
	}

	return value;
}

void Fluid3D::setCellValueAt(float *grid, float value, int x, int y, int z)
{
	// bounds check
	if (grid)
	{
		grid[IX(x, y, z)] = value;
	}

}


unsigned Fluid3D::getNumCellinX() const
{
	return m_Nx;
}


unsigned Fluid3D::getNumCellinY() const
{
	return m_Ny;
}


unsigned Fluid3D::getNumCellinZ() const
{
	return m_Nz;
}


unsigned Fluid3D::getGridTotalCells() const
{
	return m_numCells;
}


const float* Fluid3D::getVelocityX() const
{
	return m_U1;
}


const float* Fluid3D::getVelocityY() const
{
	return m_V1;
}

const float* Fluid3D::getVelocityZ() const
{
	return m_W1;
}

const float* Fluid3D::getVelocityXOld() const
{
	return m_U0;
}


const float* Fluid3D::getVelocityYOld() const
{
	return m_V0;
}


const float* Fluid3D::getVelocityZOld() const
{
	return m_W0;
}


const float* Fluid3D::getPressure() const
{
	return m_P1;
}


const float* Fluid3D::getPressureOld() const
{
	return m_P0;
}



float* Fluid3D::getTemperature()
{
	return m_T1;
}

float* Fluid3D::getTemperatureOld()
{
	return m_T0;
}


float* Fluid3D::getDensity()
{
	return m_S1;
}


float* Fluid3D::getDensityOld()
{
	return m_S0;
}


void Fluid3D::clearGrid(float *grid)
{
	if (!grid) return;

	for(int i=0; i < m_numCells; ++i)
	{
		grid[i] = 0.0;
	}

}

// add force vector F<Fx, Fy ,Fz> to the grid
void Fluid3D::addForce(float fx, float fy, float fz, float dt)
{
	if (!m_U1 || !m_V1 || !m_W1) return;

	for(int i=0; i < m_numCells; ++i)
	{
		m_U1[i] += fx*dt;
		m_V1[i] += fy*dt;
		m_W1[i] += fz*dt;
	}
}


void Fluid3D::addForce(float *u, float *v, float *w, float *fx, float *fy, float *fz, float dt)
{
	if (!fx || !fy || !fz) return;

	for(int i=0; i < m_numCells; ++i)
	{
		u[i] += fx[i]*dt;
		v[i] += fy[i]*dt;
		w[i] += fz[i]*dt;
	}
}


void Fluid3D::addBuoyancyForce()
{
	// ensure they are all positive
	if ( m_alpha < 0) m_alpha = 0.1f;
	if ( m_beta < 0 ) m_beta  = 0.4f;

	// add fbuoy = ( -alpha*rho + beta(T-Tamb) )*z where z is a vector in an upward direction i.e. <0,1,0>
	// current density and temperature are stored in S1 and T1
	for(int i = 1; i <= m_Nx; ++i)
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int k=1; k <= m_Nz; ++k)
			{
				float fb = -m_alpha* getCellValueAt(m_S1, i, j, k) + m_beta*(getCellValueAt(m_T1, i, j, k) - m_Tamb);
				setCellValueAt(m_Fy, fb, i, j, k);
			}
		}
	}

}


// apply a semi-lagrangian advection to solve PDE:
// dq/dt + u * grad(q) = 0
void Fluid3D::transport(float *s1, float *s0, float *u, float *v, float *w, float dt, bool isScalar)
{
	float x0, y0, z0, snew;

	// velocities are expressed in unit cube space, convert the time step to grid space
	const float dt0 = dt / m_Dx;

	for(int k=1; k <= m_Nz; ++k)
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int i=1; i <= m_Nx; ++i)
			{
				// trace a path backwards starting at the cell center through the velocity field
				x0 = i - dt0 * getCellValueAt(u, i, j, k);
				y0 = j - dt0 * getCellValueAt(v, i, j, k);
				z0 = k - dt0 * getCellValueAt(w, i, j, k);

				// clamp to boundaries if traced path exceeds grid boundary
				CLAMP(x0, 0.5f, m_Nx + 0.5f);
				CLAMP(y0, 0.5f, m_Ny + 0.5f);
				CLAMP(z0, 0.5f, m_Nz + 0.5f);

				// apply trilinear interpolation to find the final value of scalar S1
				snew = triLerp(s0, x0, y0, z0);
				setCellValueAt(s1, snew, i, j, k);
			}
		}
	}


	// apply boundary conditions
	applyBoundaryConditions(s1, isScalar);
}


// tri-linear interpolator
float Fluid3D::triLerp(float *s, float x, float y, float z)
{
	// locate neighbors to interpolate
	const int x0 = (int)x;
	const int x1 = x0 + 1;
	const int y0 = (int)y;
	const int y1 = y0 + 1;
	const int z0 = (int)z;
	const int z1 = z0 + 1;

	// get interpolation weights
	const float s1 = x - (float)x0;
	const float s0 = 1.0f - s1;
	const float t1 = y - (float)y0;
	const float t0 = 1.0f - t1;
	const float u1 = z - (float)z0;
	const float u0 = 1.0f - u1;

	const int i000 = IX(x0, y0, z0);
	const int i010 = IX(x0, y1, z0);
	const int i100 = IX(x1, y0, z0);
	const int i110 = IX(x1, y1, z0);
	const int i001 = IX(x0, y0, z1);
	const int i011 = IX(x0, y1, z1);
	const int i101 = IX(x1, y0, z1);
	const int i111 = IX(x1, y1, z1);

	// interpolate
	return ( u0 * (s0 * (t0 * s[i000] + t1 * s[i010]) + s1 * (t0 * s[i100] + t1 * s[i110])) +
					u1 * (s0 * (t0 * s[i001] + t1 * s[i011]) + s1 * (t0 * s[i101] + t1 * s[i111])) );

}


void Fluid3D::simulate(int numFrames, float dt)
{
	char strPBRTFileName[64]; // output filename

	for(int i=1; i <= numFrames; ++i)
	{
		// increment file id
		snprintf(strPBRTFileName, sizeof(strPBRTFileName), "Out/Frame_%d.pbrt", i);

		VelocityStep(m_U1, m_V1, m_W1, m_U0, m_V0, m_W0, dt);
		ScalarStep(m_T1, m_T0, m_U1, m_V1, m_W1, dt);
		ScalarStep(m_S1, m_S0, m_U1, m_V1, m_W1, dt);

		// export density to pbrt volume grid format
		if (exportGridToPBRTFormat(m_S1, m_Nx, m_Ny, m_Nz, strPBRTFileName) < 0)
			return;

		std::cout << "Frame-" << i << " has been exported to " << strPBRTFileName << "\n";

	}

}


void Fluid3D::ScalarStep(float *s1, float *s0, float *u, float *v, float *w, float dt)
{
	// add heat and density sources
	addSourcesToSimulation(s1);

	// diffuse s1 into s0
	diffuse(s0, s1, m_diff, dt, true);

	// advect scalar field s0 throughout the velocity field u
	// and write the output value to s1
	transport(s1, s0, u, v, w, dt, true);

}


void Fluid3D::VelocityStep(float *u1, float *v1, float *w1, float *u0, float *v0, float *w0, float dt)
{
	// add forces
	//addForce(0.0f, -9.81f, 0.0f, dt); // g: gravity acceleration
	addBuoyancyForce();
	addForce(u1, v1, w1, m_Fx, m_Fy, m_Fz, dt);

	// diffuse velocity into u0
	diffuse(u0, u1, m_visc, dt, false);
	diffuse(v0, v1, m_visc, dt, false);
	diffuse(w0, w1, m_visc, dt, false);

	project(u0, v0, w0);

	// transport each component of the velocity seperately back into u1
	transport(u1, u0, u0, v0, w0, dt, false);
	transport(v1, v0, u0, v0, w0, dt, false);
	transport(w1, w0, u0, v0, w0, dt, false);

	project(u1, v1, w1);

}


// diffuse u0 by solving (I - visc*dt*laplacian) u = u0
void Fluid3D::diffuse(float *u, float *u0, float visc, float dt, bool isScalar)
{
	if (visc <= 0.0f)
	{
		for(int i=0; i < m_numCells; ++i) u[i] = u0[i];
		return;
	}

	float alpha =  (m_Dx*m_Dx)/(visc*dt);

	// use u0 as the initial guess
	for(int i=0; i < m_numCells; ++i) u[i] = u0[i];

	jacobi(u, u0, alpha, 6.0f + alpha, isScalar);

}


// projection
void Fluid3D::project(float *u, float *v, float *w)
{

	// take the divergence of the velocity i.e. del dot u
	for(int i=1; i <= m_Nx; ++i)
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int k=1; k <= m_Nz; ++k)
			{
				float val = -0.5f*m_Dx * ( getCellValueAt(u,i+1,j,k) - getCellValueAt(u,i-1,j,k) +
							               getCellValueAt(v,i,j+1,k) - getCellValueAt(v,i,j-1,k) +
							               getCellValueAt(w,i,j,k+1) - getCellValueAt(w,i,j,k-1) )  ;

				setCellValueAt(m_Div, val, i, j, k);

				// clear out pressure values
				// since we are using an initial quess of p as 0
				setCellValueAt(m_P0, 0.0f, i, j, k);
			}
		}
	}
	applyBoundaryConditions(m_Div, true);
	applyBoundaryConditions(m_P0, true);

	// solve for the pressure poisson equation
	jacobi(m_P0, m_Div, 1.0f, 6.0f, true);

	// subtract gradient of pressure from the intermediate velocity u
	for(int i=1; i <= m_Nx; ++i)
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int k=1; k <= m_Nz; ++k)
			{
				float wx = getCellValueAt(u, i, j, k) - 0.5f/m_Dx * ( getCellValueAt(m_P0, i+1, j, k) - getCellValueAt(m_P0, i-1, j, k) );
				float wy = getCellValueAt(v, i, j, k) - 0.5f/m_Dy * ( getCellValueAt(m_P0, i, j+1, k) - getCellValueAt(m_P0, i, j-1, k) );
				float wz = getCellValueAt(w, i, j, k) - 0.5f/m_Dz * ( getCellValueAt(m_P0, i, j, k+1) - getCellValueAt(m_P0, i, j, k-1) );

				setCellValueAt(u, wx, i, j, k);
				setCellValueAt(v, wy, i, j, k);
				setCellValueAt(w, wz, i, j, k);
			}
		}
	}

	// apply boundary conditions
	applyBoundaryConditions(u, false);
	applyBoundaryConditions(v, false);
	applyBoundaryConditions(w, false);

}


// iterative solver for 3D Poisson equation
// Ax = b
// note: x is updated in place, which makes this a Gauss-Seidel iteration
void Fluid3D::jacobi(float *x, float *x0, float alpha, float beta, bool isScalar)
{
	float xnext=0.0f;

	for(int n=0; n < m_numberOfJacobiIters; ++n)
	{
		for(int i=1; i <= m_Nx; ++i)
		{
			for(int j=1; j <= m_Ny; ++j)
			{
				for(int k=1; k <= m_Nz; ++k)
				{
					// 3D Jacobi iteration update based on the following formula
					// (x[i][j][k])^n = (x[i-1][j][k] + x[i+1][j][k] + x[i][j-1][k] + x[i][j+1][k] + x[i][j][k-1]+
					//				     x[i][j][k+1] + alpha*x0[i][j][k]) * 1.0/beta
					xnext = ( alpha*getCellValueAt(x0, i, j, k)  +  getCellValueAt(x, i-1, j, k) +  getCellValueAt(x, i+1, j, k) +
						      getCellValueAt(x, i, j-1, k) + getCellValueAt(x, i, j+1, k) + getCellValueAt(x, i, j, k-1) + getCellValueAt(x, i, j, k+1))/beta;

					// set output value x,  after nth iteration (x[i][j][k])^(n+1)
					setCellValueAt(x, xnext , i, j, k);
				}
			}
		}

		// apply correct boundary conditions after the end of each iteration
		applyBoundaryConditions(x, isScalar);

	}
}


// apply correct boundary conditions
// we are only considering here two boundary conditions:
// No-Slip condition and pure neumann boundary condition
void Fluid3D::applyBoundaryConditions(float *u, bool isScalar)
{
	float sign = isScalar ? 1.0f : -1.0f;// switch btw pure neumann and no-slip conditions flag

	// Left & Right Walls
	for(int z=0; z < (m_Nz + 2); ++z)
	{
		for(int y=0; y < (m_Ny + 2); ++y)
		{
			// left side
			setCellValueAt(u, sign*getCellValueAt(u, 1, y, z),  0, y, z);

			// right side
			setCellValueAt(u, sign*getCellValueAt(u, m_Nx, y, z),  m_Nx+1, y, z);

		}
	}

	// Top & Bottom Walls
	for(int z=0; z < (m_Nz + 2); ++z)
	{
		for(int x=0; x < (m_Nx + 2); ++x)
		{
			// bottom side
			setCellValueAt(u, sign*getCellValueAt(u, x, 1, z),  x, 0, z);

			// top side
			setCellValueAt(u, sign*getCellValueAt(u, x, m_Ny, z), x , m_Ny+1, z);

		}
	}

	// Front & Back Walls
	for(int y=0; y < (m_Ny + 2); ++y)
	{
		for(int x=0; x < (m_Nx + 2); ++x)
		{
			// front side
			setCellValueAt(u, sign*getCellValueAt(u, x, y, 1),  x, y, 0);

			// back side
			setCellValueAt(u, sign*getCellValueAt(u, x, y, m_Nz),  x, y, m_Nz+1);

		}
	}

}

void Fluid3D::setBoundaryToZero(float *u)
{
	// Left & Right Walls
	for(int z=0; z < (m_Nz + 2); ++z)
	{
		for(int y=0; y < (m_Ny + 2); ++y)
		{
			// left side
			setCellValueAt(u, 0.0f,  0, y, z);

			// right side
			setCellValueAt(u, 0.0f,  m_Nx+1, y, z);

		}
	}

	// Top & Bottom Walls
	for(int z=0; z < (m_Nz + 2); ++z)
	{
		for(int x=0; x < (m_Nx + 2); ++x)
		{
			// bottom side
			setCellValueAt(u, 0.0f,  x, 0, z);

			// top side
			setCellValueAt(u, 0.0f, x , m_Ny+1, z);

		}
	}

	// Front & Back Walls
	for(int y=0; y < (m_Ny + 2); ++y)
	{
		for(int x=0; x < (m_Nx + 2); ++x)
		{
			// front side
			setCellValueAt(u, 0.0f,  x, y, 0);

			// back side
			setCellValueAt(u, 0.0f,  x, y, m_Nz+1);

		}
	}
}

void Fluid3D::addSourcesToSimulation(float* field)
{

	/*====== 3D Gaussian splat ===========*/
	float radius = m_Nx / 4.0f;
	float point[3] = { 0.5f*m_Nx, 0.0f, 0.5f*m_Nz};

	for (int z = 1; z <= m_Nz; z++)
	{
		for (int y =1; y <= m_Ny; y++)
		{
			 for (int x = 1; x <= m_Nx; x++)
			 {

				float dst = sqrtf ( (point[0]-x) * (point[0]-x) + (point[1]-y) * (point[1]-y) + (point[2]-z) * (point[2]-z) );

				if (dst < radius)
				{
					float a = (radius-dst)*0.5f;
					a = MIN2(a, 1.0f);
					setCellValueAt(field, a, x, y, z);
				}
			}
		 }
	  }

}

int Fluid3D::exportGridToDF3Format(const float *grid, unsigned short nx, unsigned short ny, unsigned short nz, const char *filename)
{
	// add .df3 extension to filename
	std::string path = std::string(filename) + ".df3";

	// try to open file for writing
	FILE *fp = fopen(path.c_str(), "wb");

	if( fp == NULL)
	{
		std::cerr << "Error(!): Unable to open the file " << path << "\n";
		return -1;
	}

	// DF3 stores everything in big endian order
	auto toBigEndian = [](unsigned short v) -> unsigned short
	{
		const unsigned char bytes[2] = { (unsigned char)(v >> 8), (unsigned char)(v & 0xFF) };
		unsigned short out;
		memcpy(&out, bytes, sizeof(out));
		return out;
	};

	// dimensions
	const int byteSize = 2;
	const unsigned short int nx1 = toBigEndian(nx), ny1 = toBigEndian(ny), nz1 = toBigEndian(nz);

	fwrite((const void*)&nx1, sizeof(unsigned short int), 1, fp);
	fwrite((const void*)&ny1, sizeof(unsigned short int), 1, fp);
	fwrite((const void*)&nz1, sizeof(unsigned short int), 1, fp);

	const int nitems = nx*ny*nz;
	const float mul = (float)( (1<<(8*byteSize))-1);

	std::vector<unsigned short int> buf(nitems);

	// skip the boundary layer of the padded grid
	for (int k = 0; k < nz; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++)
			{
				float val = grid[(i+1) + (j+1)*(nx+2) + (k+1)*(nx+2)*(ny+2)];
				CLAMP(val, 0.0f, 1.0f);
				buf[k*(nx*ny)+j*nx+i] = toBigEndian((unsigned short int)(val*mul));
			}
		}
	}

	fwrite((const void*)buf.data(), sizeof(unsigned short int), nitems, fp);

	fclose(fp);

	return 1;

}


int Fluid3D::exportGridToPBRTFormat(const float* grid, int nx, int ny, int nz, const char *filename)
{
	// try to open file for writing
	FILE *fp = fopen(filename, "w");

	if( fp == NULL)
	{
		std::cerr << "Error(!): Unable to open the file " << filename << "\n";
		return -1;
	}

	const int size = nx*ny*nz;

	// copy the interior of the padded grid into a buffer
	std::vector<float> field(size);

	for (int k = 0; k < nz; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++)
			{
				field[k*(nx*ny) + j*nx + i] = std::fabs(grid[(i+1) + (j+1)*(nx+2) + (k+1)*(nx+2)*(ny+2)]);
			}
		}
	}

	// normalize values
	float maxDensVal = 0.0f;
	const float targetNorm = 0.5f;

	for (int i = 0; i < size; i++)
	{
		if (field[i] > maxDensVal) maxDensVal = field[i];
	}

	if (maxDensVal > 0.0f)
	{
		for (int i = 0; i < size; i++)
		{
			field[i] = field[i] / maxDensVal * targetNorm;
		}
	}

	int maxRes = MAX3(nx, ny, nz);

	const float xSize = 1.0f / (float)maxRes * (float)nx;
	const float ySize = 1.0f / (float)maxRes * (float)ny;
	const float zSize = 1.0f / (float)maxRes * (float)nz;

	// pbrt-v4 grid medium named "smoke", Include this file from a scene and attach it
	// to a shape with MediumInterface "smoke" ""
	// sigma_a, sigma_s and scale control how dark and thick the smoke looks
	fprintf(fp, "MakeNamedMedium \"smoke\" \"string type\" \"uniformgrid\"\n");
	fprintf(fp, " \"rgb sigma_a\" [ 1 1 1 ] \"rgb sigma_s\" [ 4 4 4 ] \"float scale\" 4\n");
	fprintf(fp, " \"integer nx\" %i\n", nx);
	fprintf(fp, " \"integer ny\" %i\n", ny);
	fprintf(fp, " \"integer nz\" %i\n", nz);
	fprintf(fp, " \"point3 p0\" [ 0.0 0.0 0.0 ] \"point3 p1\" [%f %f %f ] \n", xSize, ySize, zSize);
	fprintf(fp, " \"float density\" [ \n");

	for(int i = 0; i < size; i++) {
		fprintf(fp, "%f ", field[i]);
	}

	// end of volume density
	fprintf(fp, "] \n \n");

	fclose(fp);

	return 1;

}
