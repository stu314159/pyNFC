/*
 * LBM_DataHandler.h
 *
 *  Created on: Jul 10, 2017
 *      Author: stu
 */

#ifndef LBM_DATAHANDLER_H_
#define LBM_DATAHANDLER_H_
/*
 * Data structure designed to hold individual lattice state data
 * for a full time-step.  Separate individual lattice point state data
 * from the lattice-structure class that says what to do with the data.
 * Ultimately, many threads may want to use the same lattice class instantiation
 * but need to processes separate lattice points.  This is intended to
 * facilitate that behavior.
 *
 */


class LBM_DataHandler
{
public:
	LBM_DataHandler(const int numSpd);
	~LBM_DataHandler();
	float ux;
	float uy;
	float uz;
	float rho;
	float u_bc;
	float rho_bc;
	int nodeType;
	float omega;
	float * f; // dangerous, but publicly available
	float * fEq;
	float * fOut;
	float piFlat[9];


private:
	const int numSpd;
};



#endif /* LBM_DATAHANDLER_H_ */
