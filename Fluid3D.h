/**
 * @file Fluid3D.cpp
 * @author Berk Atabek (atabek @domain.com)
 * @brief This is an implementation of a 3D Fluid Solver 
 *       based on Jos Stams' "Stable Fluids" paper.
 *       It uses a colocated grid instead of a staggered grid 
 *		 since this is straightforward to implement.
 */
#ifndef FLUID3D_H
#define FLUID3D_H

class Fluid3D 
{
	public:
		// constructors
		Fluid3D();
		Fluid3D(unsigned nx, unsigned ny, unsigned nz);
		~Fluid3D();

		// accessors
		unsigned getNumCellinX()const;
		unsigned getNumCellinY()const;
		unsigned getNumCellinZ()const;
		unsigned getGridTotalCells() const;

		const float* getVelocityX() const;
		const float* getVelocityY() const;
		const float* getVelocityZ() const;
		const float* getVelocityXOld() const;
		const float* getVelocityYOld() const;
		const float* getVelocityZOld() const;

		const float* getPressure() const;
		const float* getPressureOld() const;
		
		float* getTemperature();
		float* getTemperatureOld();
		
		float* getDensity(); 
		float* getDensityOld();
			
		// setters
		

		// main solver interface
		void simulate(int numFrames, float dt);	
		void ScalarStep(float *s1, float *s0, float *u, float *v, float *w, float source, float dt);
		void VelocityStep(float *u1, float *v1, float *w1, float *u0, float *v0, float *w0, float dt);

		// add a constant force to velocity grid
		void addForce(float fx, float fy, float fz, float dt);
		
		// add a force field 
		void addForce(float *u, float *v, float *w, float *fx, float *fy, float *fz, float dt);
		
		// add buoyant force to force grid
		void addBuoyancyForce();
		
		// transport (semi lagrangian advection)
		// advect scalar or vector quantity s0 through fluid's velocity field <u,v,w>,
		// output value is stored in s1
		void transport(float *s1, float *s0, float *u, float *v, float *w, float dt);
		
		// diffuse
		void diffuse(float *u, float *u0, float visc, float dt);

		// projection 
		void project(float *u1, float *v1, float *w1, float *u0, float *v0, float *w0, float dt);

		// 3D Jacobi iterative solver
		void jacobi(float *x, float *x0, float alpha, float beta, bool isPressure);

		// tri-linear interpolator
		float triLerp(float *s, float x, float y, float z);

		// apply boundary conditions
		void applyBoundaryConditions(float *u0, bool isPressure);
		void setBoundaryToZero(float *u);

		// add concentration of particles and heat sources
		void addSourcesToSimulation(float* field, int xres, int yres, int zres);

		// creates a DF3 density file to visualize in POV-Ray
		int exportGridToDF3Format(const float *grid, unsigned short int nx, unsigned short int ny, unsigned short int nz, char *filename);

		// creates a volume for PBRT
		int exportGridToPBRTFormat(float* grid, int nx, int ny, int nz, char *filename);
		
	private:

		// helpers for retrieving grid value based on the 3D index value
		float getCellValueAt(float *grid, int x, int y, int z);
		void setCellValueAt(float *grid, float value, int x, int y, int z);

		// zero out the values in the grid
		void clearGrid(float *grid);

		// number of cells in each coordinate
		int m_Nx, m_Ny, m_Nz;
		int m_NMax; // for making uniform grid

		// size of each voxel for each coordinate
		// we are asuming a grid size of 1 for each dimensions, so that 
		// it maps unit cube when rendering density  
		float m_Dx, m_Dy, m_Dz;

		// number of total grid cells
		int m_numCells;
				
		float *m_U0, *m_V0, *m_W0;  // 3D velocity grid <u,v,w>
		float *m_U1, *m_V1, *m_W1; 
		float *m_Fx, *m_Fy, *m_Fz;  // 3D force grid(f)	         
		float *m_P0, *m_P1;			// pressure grid
		float *m_T0, *m_T1;			// temperature grid
		float *m_S0, *m_S1;			// scalar grid(density)
		float *m_Div;				// divergence i.e. del dot u
						
		// buoyancy force parameters
		float m_alpha, m_beta;
		float m_Tamb; // ambient temperature

		// time step
		float m_dt;

		// number of iterations required for the Jacobi solver
		int m_numberOfJacobiIters;

};


#endif