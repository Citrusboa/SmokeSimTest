#ifndef GRID_DATA_MATRIX_H
#define GRID_DATA_MATRIX_H

#include "grid_data.h"

class GridDataMatrix {
public:
	GridData diag;
	GridData plusI;
	GridData plusJ;
	GridData plusK;

	GridDataMatrix() {
		diag.initialize();
		plusI.initialize();
		plusJ.initialize();
		plusK.initialize();
	}
};

#endif