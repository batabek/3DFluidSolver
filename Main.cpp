/**
 * @file Main.cpp
 * @author Berk Atabek (atabek dot berk at gmail dot com)
 * @brief Driver program for 3D Fluid Solver
 * @version 0.1
 * @date 2020-04-08
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "Fluid3D.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

// grid sizes
const int NX = 1<<4;
const int NY = 1<<4;
const int NZ = 1<<4;

// simulation parameters
const int numFrames = 100;
const float dt = 0.1f;

int main(int argc, char **argv)
{
	
	// create output directories
	system("rmdir /q /s Out");
	system("mkdir Out");
	//system("cd Out");

	// create a fluid container
	Fluid3D  g_Fluid3d(NX, NY, NZ);
	
	// start simulating 
	std::cout << "Beginning Simulation...\n";
	g_Fluid3d.simulate(numFrames, dt);
	std::cout << "Simulation has finished...\n";

	return (int)0;	

}