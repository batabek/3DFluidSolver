#include "Fluid3D.h"
#include <iostream>
#include <cmath>
#include <cstring>

#define MIN2(x, y) ( (x) < (y) ? x : y )
#define MAX2(x, y) ( (x) > (y) ? x : y )
#define MAX3(x, y, z) ( (x) > MAX2( (y), (z) ) ? (x) :  MAX2( (y), (z) ) ) 

inline void CLAMP (float &x, float a, float b)
{
	if ( x < a ) x = a;
	if ( x > b)  x = b;

}

inline void SWAP(float* &grid1, float* &grid2)
{
	float *temp = grid1;
	grid1 = grid2;
	grid2 = temp;
}


Fluid3D::Fluid3D()
{
	// default grid size
	m_Nx = m_Ny = m_Nz = 32;
	m_NMax = MAX3(m_Nx, m_Ny, m_Nz); 

	// voxel sizes(step size)
	// Here we are discretizing inside a unit cube 
	m_Dx = m_Dy = m_Dz = 1.0f/(float)m_NMax;
	
	// allocate 2 extra grid layer to accomodate the boundary conditions
	m_numCells = (m_Nx + 2)  * (m_Ny + 2) * (m_Nz + 2);

	m_U0 = new float[m_numCells];
	m_V0 = new float[m_numCells];
	m_W0 = new float[m_numCells];
	m_U1 = new float[m_numCells];
	m_V1 = new float[m_numCells];
	m_W1 = new float[m_numCells];
	m_Fx = new float[m_numCells];
	m_Fy = new float[m_numCells];
	m_Fz = new float[m_numCells]; 
	m_P0 = new float[m_numCells];
	m_P1 = new float[m_numCells];
	m_T0 = new float[m_numCells];
	m_T1 = new float[m_numCells];
    m_S0 = new float[m_numCells];
	m_S1 = new float[m_numCells];  
    m_Div = new float[m_numCells];

	// zero out all grids
	for(int i=0; i < m_numCells; ++i)
	{
		m_U0[i] = m_V0[i] = m_W0[i] = m_U1[i]= m_V1[i] = m_W1[i] = 0.0f;
		m_Fx[i] = m_Fy[i] = m_Fz[i] = 0.0f;
		m_P0[i] = m_P1[i] = 0.0f;
		m_T0[i] = m_T1[i] = 0.0f;
		m_S0[i] = m_S1[i] = 0.0f;
		m_Div[i]= 0.0f;
		
	}
	    
    m_dt    = 0.1f;
    m_numberOfJacobiIters = 40;
	m_alpha = 0.1f;
	m_beta  = 0.4f;
	m_Tamb  = 0.0f;

}


Fluid3D::Fluid3D(unsigned int nx, unsigned int ny, unsigned int nz)
{
	// allocate 2 extra grid layer to take account the boundary conditions
	m_Nx = nx;
	m_Ny = ny;
	m_Nz = nz;
	m_NMax = MAX3(m_Nx, m_Ny, m_Nz); 
	
	// voxel sizes
	m_Dx = 1.0f / (float)m_Nx;
	m_Dy = 1.0f / (float)m_Ny;
	m_Dz = 1.0f / (float)m_Nz;
	
	m_numCells = (m_Nx + 2) * (m_Ny + 2) * (m_Nz + 2);
		
	m_U0 = new float[m_numCells];
	m_V0 = new float[m_numCells];
	m_W0 = new float[m_numCells];
	m_U1 = new float[m_numCells];
	m_V1 = new float[m_numCells];
	m_W1 = new float[m_numCells];
	m_Fx = new float[m_numCells];
	m_Fy = new float[m_numCells];
	m_Fz = new float[m_numCells]; 
	m_P0 = new float[m_numCells];
	m_P1 = new float[m_numCells];
	m_T0 = new float[m_numCells];
	m_T1 = new float[m_numCells];
    m_S0 = new float[m_numCells];
	m_S1 = new float[m_numCells];  
    m_Div= new float[m_numCells];

	// zero out all grids
	for(int i=0; i < m_numCells; ++i)
	{
		m_U0[i] = m_V0[i] = m_W0[i] = m_U1[i]= m_V1[i] = m_W1[i] = 0.0f;
		m_Fx[i] = m_Fy[i] = m_Fz[i] = 0.0f;
		m_P0[i] = m_P1[i] = 0.0f;
		m_T0[i] = m_T1[i] = 0.0f;
		m_S0[i] = m_S1[i] = 0.0f;
		m_Div[i]= 0.0f;
		
	}
	    
    m_dt    = 0.1f;
    m_numberOfJacobiIters = 20;
	m_alpha = 0.1f;
	m_beta  = 0.6f;
	m_Tamb  = 0.0f;
	
}


Fluid3D:: ~Fluid3D()
{
	// deallocate grids
	if (m_U0) delete [] m_U0;
	if (m_V0) delete [] m_V0;
	if (m_W0) delete [] m_W0;
	if (m_U1) delete [] m_U1;
	if (m_V1) delete [] m_V1;
	if (m_W1) delete [] m_W1;


	if (m_P0)	delete	[] m_P0;
	if (m_P1)	delete	[] m_P1;
	

	if (m_T0)	delete	[] m_T0;
	if (m_T1)	delete	[] m_T1;

	if (m_S0)	delete	[] m_S0;
	if (m_S1)	delete	[] m_S1;

	if(m_Div) delete [] m_Div;

	if (m_Fx) delete [] m_Fx;
	if (m_Fy) delete [] m_Fy;
	if (m_Fz) delete [] m_Fz;
	

}
              
float Fluid3D::getCellValueAt(float *grid, int x, int y, int z)
{
	float value = -0.999999999f;
	
	// bounds check
	if (grid)
	{
		value = grid[x + y*m_Nx + m_Nx*m_Ny*z];
	}

	return value;
}

void Fluid3D::setCellValueAt(float *grid, float value, int x, int y, int z)
{
	// bounds check
	if (grid)
	{
		grid[x + y*m_Nx + m_Nx*m_Ny*z] = value;
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
	if (!m_U0 || !m_V0 || !m_W0) return;

	for(int i=0; i < m_numCells; ++i)
	{
		m_U0[i] += fx*dt;
		m_V0[i] += fy*dt;
		m_W0[i] += fz*dt;
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
	for(int i = 1; i <= m_Nx; ++i) 
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int k=1; k <= m_Nz; ++k)
			{
				float fb = -m_alpha* getCellValueAt(m_S0, i, j, k) + m_beta*(getCellValueAt(m_T0, i, j, k) - m_Tamb);
				fb /= (float)m_NMax;
				setCellValueAt(m_Fy, fb, i, j, k);
			}
		}
	}
	
}


// apply a semi-lagrangian advection to solve PDE:
// dq/dt + u * grad(q) = 0
void Fluid3D::transport(float *s1, float *s0, float *u, float *v, float *w, float dt)
{
	float x0, y0, z0, snew;
		
	for(int k=1; k <= m_Nz; ++k) 
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int i=1; i <= m_Nx; ++i)
			{
				// sample the grid at the cell centers
				//x = (i + 0.5f) * m_Dx; y = (j + 0.5f) * m_Dy; z = (k + 0.5f) * m_Dz;
				
				// trace a path starting at x,y,z through the velocity field U0
				x0 = i - dt * getCellValueAt(u, i, j, k);
				y0 = j - dt * getCellValueAt(v, i, j, k);
				z0 = k - dt * getCellValueAt(w, i, j, k);
				
				// clamp to boundaries if traced path exceeds grid boundary 
				if (x0 < 0.5f) x0 = 0.5f; if (x0 > m_Nx - 1.5f) x0 = m_Nx - 1.5f; 
				if (y0 < 0.5f) y0 = 0.5f; if (y0 > m_Ny - 1.5f) y0 = m_Ny - 1.5f; 
				if (z0 < 0.5f) z0 = 0.5f; if (z0 > m_Nz - 1.5f) z0 = m_Nz - 1.5f;
				
				// apply trilinear interpolation to find the final value of scalar S1
				snew = triLerp(s0, x0, y0, z0);
				setCellValueAt(s1, snew, i, j, k);
			}
		}
	}
	

	// apply boundary conditions
	applyBoundaryConditions(s1, false);
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

	const int slabSize = m_Nx*m_Ny;
	const int i000 = x0 + y0 * m_Nx + z0 * slabSize;
	const int i010 = x0 + y1 * m_Nx + z0 * slabSize;
	const int i100 = x1 + y0 * m_Nx + z0 * slabSize;
	const int i110 = x1 + y1 * m_Nx + z0 * slabSize;
	const int i001 = x0 + y0 * m_Nx + z1 * slabSize;
	const int i011 = x0 + y1 * m_Nx + z1 * slabSize;
	const int i101 = x1 + y0 * m_Nx + z1 * slabSize;
	const int i111 = x1 + y1 * m_Nx + z1 * slabSize;

	// interpolate (indices could be computed once)
	return ( u0 * (s0 * (t0 * s[i000] + t1 * s[i010]) + s1 * (t0 * s[i100] + t1 * s[i110])) +
					u1 * (s0 * (t0 * s[i001] + t1 * s[i011]) + s1 * (t0 * s[i101] + t1 * s[i111])) );

}


void Fluid3D::simulate(int numFrames, float dt)
{
	char strPBRTFileName[50]; // output filename

	for(int i=1; i <= numFrames; ++i)
	{
		// increment file id
		sprintf(strPBRTFileName, ".//Out//Frame_%d", i);
						
		VelocityStep(m_U1, m_V1, m_W1, m_U0, m_V0, m_W0, dt);
		ScalarStep(m_T1, m_T0, m_U1, m_V1, m_W1, 100.0f, dt);
		ScalarStep(m_S1, m_S0, m_U1, m_V1, m_W1, 100.0f, dt);
			
		// export density to pbrt volume grid format
		exportGridToPBRTFormat(m_S0, m_Nx, m_Ny, m_Nz, strPBRTFileName);
		
		//g_Fluid3d.exportGridToDF3Format(g_Fluid3d.getScalar(), NX, NY, NZ, filename); 
		std::cout << "Frame-" << i << " has been exported to " << "Out\\Frame_" << i << "\n";
			 
	}
		
}


void Fluid3D::ScalarStep(float *s1, float *s0, float *u, float *v, float *w, float source, float dt)
{
	// advect scalar field s0 throughout the velocity field u 
	// and write the output value to s1
	// add heat and density sources
	addSourcesToSimulation(s1, m_Nx, m_Ny, m_Nz);
	SWAP(s0, s1); diffuse(s1, s0, 0.5f, dt);
	SWAP(s0, s1); transport(s1, s0, u, v, w, dt);

}


void Fluid3D::VelocityStep(float *u1, float *v1, float *w1, float *u0, float *v0, float *w0, float dt)
{
	// add forces
	//addForce(0.0f, 0.99f, 0.0f, dt); // g: gravity acceleration 
	addBuoyancyForce();
	addForce(u1, v1, w1, m_Fx, m_Fy, m_Fz, dt);
	
	diffuse(u0, u1, 0.5f, dt);
	diffuse(v0, v1, 0.5f, dt);
	
	SWAP(u0, u1);SWAP(v0, v1);SWAP(w0, w1);
	project(u1, v1, w1, u0, v0, w0, dt);
	
	// transport each component of the velocity seperately
	transport(u0, u1, u0, v0, w0, dt);
	transport(v0, v1, u0, v0, w0, dt);
	transport(w0, w1, u0, v0, w0, dt);
	
	//project(u1, v1, w1, u0, v0, w0, dt);

}


// diffuse
void Fluid3D::diffuse(float *u, float *u0, float visc, float dt)
{
	// not required for the smoke
	float alpha =  (m_Dx*m_Dx*m_Dx)/(visc*dt);
	
	jacobi(u, u0, alpha, 6.0f + alpha, false);
	
}


// projection 
void Fluid3D::project(float *u1, float *v1, float *w1, float *u0, float *v0, float *w0, float dt)
{
	
	// take the divergence of the velocity i.e. del dot u
	for(int i=1; i <= m_Nx; ++i)
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int k=1; k <= m_Nz; ++k)
			{
				float val = -0.5f*m_Dx * ( getCellValueAt(u1,i+1,j,k) - getCellValueAt(u1,i-1,j,k) +
							               getCellValueAt(v1,i,j+1,k) - getCellValueAt(v1,i,j-1,k) + 
							               getCellValueAt(w1,i,j,k+1) - getCellValueAt(w1,i,j,k-1) )  ;
						
					setCellValueAt(m_Div, val, i, j, k);
					
					// clear out pressure values
					// since we are using an initial quess of p as 0
					setCellValueAt(m_P0, 0.0f, i, j, k);
			}
		}
	}
	// are these necessary for the divergence? 
	applyBoundaryConditions(m_Div, true);
	applyBoundaryConditions(m_P0, true);
	
	// solve for the pressure poisson equation
	jacobi(m_P0, m_Div, -m_Dx*m_Dx*m_Dx, 6.0f, true);
	
	// subtract gradient of pressure from the intermediate velocity u
	for(int i=1; i <= m_Nx; ++i)
	{
		for(int j=1; j <= m_Ny; ++j)
		{
			for(int k=1; k <= m_Nz; ++k)
			{
				float wx = getCellValueAt(u1, i, j, k) - 0.5f/(float)m_Dx * ( getCellValueAt(m_P0, i+1, j, k) - getCellValueAt(m_P0, i-1, j, k) );
				float wy = getCellValueAt(v1, i, j, k) - 0.5f/(float)m_Dx * ( getCellValueAt(m_P0, i, j+1, k) - getCellValueAt(m_P0, i, j-1, k) );
				float wz = getCellValueAt(w1, i, j, k) - 0.5f/(float)m_Dx * ( getCellValueAt(m_P0, i, j, k+1) - getCellValueAt(m_P0, i, j, k-1) );

				setCellValueAt(u1, wx, i, j, k);
				setCellValueAt(v1, wy, i, j, k);
				setCellValueAt(w1, wz, i, j, k);
			}
		}
	}

	// apply boundary conditions
	applyBoundaryConditions(u1, false);
	applyBoundaryConditions(v1, false);
	applyBoundaryConditions(w1, false);

}


// jacobi iterative solver for 3D Poisson equation
// Ax = b
void Fluid3D::jacobi(float *x, float *x0, float alpha, float beta, bool isPressure)
{
	int n = 0 ;
	float xnext=0.0f;
	
	while( ++n < m_numberOfJacobiIters )
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
		applyBoundaryConditions(x, isPressure);
		
	}
}


// apply correct boundary conditions 
// we are only considering here two boundary conditions:
// No-Slip condition and pure neumann boundary condition
void Fluid3D::applyBoundaryConditions(float *u, bool isPressure)
{
	float sign = isPressure ? 1.0f : -1.0f;// switch btw pure neumann and no-slip conditions flag

	// Left & Right Walls
	for(int z=0; z < (m_Nz + 2); ++z)
	{
		for(int y=0; y < (m_Ny + 2); ++y)
		{
			// left side
			setCellValueAt(u, sign*getCellValueAt(u, 1, y, z),  0, y, z);
			
			// right side
			setCellValueAt(u, sign*getCellValueAt(u, m_Nx-1, y, z),  m_Nx+1, y, z);
			
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
			setCellValueAt(u, sign*getCellValueAt(u, x, m_Ny-1, z), x , m_Ny+1, z);
			
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
			setCellValueAt(u, sign*getCellValueAt(u, x, y, m_Nz-1),  x, y, m_Nz+1);
			
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

void Fluid3D::addSourcesToSimulation(float* field, int xres, int yres, int zres)
{
	
	/*====== 3D Gaussian splat ===========*/
	float radius = m_Nx / 4.0f;
	float point[3] = { 0.5f*m_Nx, 0.0f, 0.5f*m_Nz};

	for (int z = 1; z <= zres; z++)
	{
		for (int y =1; y <= yres; y++)
		{
			 for (int x = 1; x <= xres; x++)
			 {
				
				float dst = sqrtf ( (point[0]-x) * (point[0]-x) + (point[1]-y) * (point[1]-y) + (point[2]-z) * (point[2]-z) ); 
				
				if (dst < radius)
				{
					float a = (radius-dst)*0.5f;
					a = MIN2(a, 1.0f);
					int index = x + y * xres + z * xres*yres;
					field[index] = a;
				}
			}
		 }
	  }
	
}

int Fluid3D::exportGridToDF3Format(const float *grid, unsigned short nx, unsigned short ny, unsigned short nz, char *filename)
{
	// add .df3 extension to filename
	char *ext = ".df3";
	strcat(filename, ext);
	
	// try to open file for writing
	FILE *fp = fopen(filename, "w");

	if( fp == NULL)
	{
		std::cerr << "Error(!): Unable to open the file " << filename << "\n";
		return -1;
	}

	// dimensions
    const int byteSize = 2;
    const unsigned short int onx=nx,ony=ny,onz=nz;
    unsigned short int nx1,ny1,nz1; 
    
	nx1 = onx >> 8;
    ny1 = ony >> 8;
    nz1 = onz >> 8;
    nx1 += (onx << 8);
    ny1 += (ony << 8);
    nz1 += (onz << 8);
    
	fwrite((void*)&nx1, sizeof(short), 1, fp);
    fwrite((void*)&ny1, sizeof(short), 1, fp);
    fwrite((void*)&nz1, sizeof(short), 1, fp);

    const int nitems = onx*ony*onz;
    const float mul = (float)( (1<<(8*byteSize))-1); 

    unsigned short int *buf = new unsigned short int[nitems];
    
	for (int k = 0; k < onz; k++) 
      for (int j = 0; j < ony; j++) 
        for (int i = 0; i < onx; i++)
		{
			float val = grid[k*(onx*ony)+j*onx+i] ;
			CLAMP(val, 0.0, 1.0);
			buf[k*(onx*ony)+j*onx+i] = (short int)(val*mul);
        }
    
	fwrite((void*)buf, sizeof(unsigned short int)* nitems, 1, fp);

    fclose(fp);
    delete[] buf;

	return 1;

}


int Fluid3D::exportGridToPBRTFormat(float* grid, int nx, int ny, int nz, char *filename)
{
    // add .pbrt extension to filename
	char *ext = ".pbrt";
	strcat(filename, ext);
	
	// try to open file for writing
	FILE *fp = fopen(filename, "wb");

	if( fp == NULL)
	{
		std::cerr << "Error(!): Unable to open the file " << filename << "\n";
		return -1;
	}
    
	int size = nx*ny*nz;
	
	// create a buffer for the field
    float *field = new float[size];
    
	// normalize values
    float maxDensVal = abs(grid[0]);
    float targetNorm = 0.5f;

	
    for (int i = 0; i < size; i++)
	{
      //float c = grid[i];
	  if( abs(grid[i])>maxDensVal) maxDensVal = abs(grid[i]);
      field[i] = 0.;
    }

    if(maxDensVal>0.0)
	{
      for (int i = 0; i < size; i++)
	  {
		field[i] = abs(grid[i]) / maxDensVal * targetNorm;
      }
    }
    
    
    int maxRes = MAX3(nx, ny, nz);
    
    const float xSize = 1.0f/ (float)maxRes * (float)nx;
    const float ySize = 1.0f / (float)maxRes * (float)ny;
    const float zSize = 1.0f / (float)maxRes * (float)nz;

    // dimensions
    fprintf(fp, "Volume \"volumegrid\" \n");
    fprintf(fp, " \"integer nx\" %i\n", nx);
    fprintf(fp, " \"integer ny\" %i\n", ny);
    fprintf(fp, " \"integer nz\" %i\n", nz);
    fprintf(fp, " \"point p0\" [ 0.0 0.0 0.0 ] \"point p1\" [%f %f %f ] \n", xSize, ySize, zSize);
    fprintf(fp, " \"float density\" [ \n");
	
	for(int i = 0; i < size; i++) {
		fprintf(fp, "%f ", abs(grid[i]));
	}
	
	
	// end of volume density 
	fprintf(fp, "] \n \n");

    fclose(fp);
    delete[] field;
	
	return 1;

  }
