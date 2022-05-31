#include "PhysParameters.h"


std::vector<TMatrixD> MakeMatricesA(const std::vector<TVector2> v1, const std::vector<TVector2> v2)
{	
	std::vector<TMatrixD> v;
	for (int i = 0; i < v1.size(); i++)
	{	
		std::vector<double> temp{ v1[i].X(), v1[i].Y(), v2[i].X(), v2[i].Y()};
		TMatrixD m(2, 2, &temp[0]);
		v.push_back(m);
	}
	return v;
}

std::vector<TMatrixD> MakeMatricesB(const std::vector<double> v1, const std::vector<double> v2)
{	
	std::vector<TMatrixD> v;
	for (int i = 0; i < v1.size(); ++i)
	{	
		std::vector<double> temp{ v1[i], v2[i] };
		TMatrixD m(2, 1, &temp[0]);
		v.push_back(m);
	}
	return v;
}

std::vector<TMatrixD> SovleMatricesEq(const std::vector<TMatrixD> As, const std::vector<TMatrixD> Bs)
{	
	std::vector<TMatrixD> Xs;
	for (int i = 0; i < As.size(); ++i)
	{
		TMatrixD A = As[i];
		TMatrixD B = Bs[i];
		cout<<"Rap: "<<Rap[i]<<endl;
		cout<<"[	"<<A(0,0)<<"	"<<A(0,1)<<"	]"<<"	"<<"["<<B(0,0)<<"]"<<endl;
		cout<<"[	"<<A(1,0)<<"	"<<A(1,1)<<"	]"<<"	"<<"["<<B(1,0)<<"]"<<endl;
		TMatrixD X = A.Invert() * B;
		cout<<"X.T: ["<<X(0,0)<<" "<<X(1,0)<<"]"<<endl;
		cout<<endl;
		Xs.push_back(X);
	}
	return Xs;
}

std::vector<std::vector<double>> MaxtricesToVecs(std::vector<TMatrixD> m)
{	
	std::vector<std::vector<double>> vs;
	for (int i = 0; i < m.size(); ++i)
	{
		std::vector<double> v;
		v.push_back(m[i](0,0)); v.push_back(m[i](1,0));
		vs.push_back(v);
		cout<<"Xsect: "<<m[i](0,0)<<" "<<m[i](1,0)<<endl;
	}
	return vs;
}

void photonEnergyContrib()
{

	cout<< "photonEnergyContrib ------> Processing 0nXnSum_XnXn Case"<<endl;
	// Now Solving matrix equation A * X = B
	std::vector<TMatrixD> As_0nXnSum_XnXn = MakeMatricesA	(flux_0nXnSum, 	flux_XnXn);
	std::vector<TMatrixD> Bs_0nXnSum_XnXn = MakeMatricesB	(Xsec_0nXnSum,		Xsec_XnXn);
	std::vector<TMatrixD> Xs_0nXnSum_XnXn = SovleMatricesEq	(As_0nXnSum_XnXn, Bs_0nXnSum_XnXn);

	// cout<< "photonEnergyContrib ------> Processing 0nXn_XnXn Case"<<endl;
	// // // Now Solving matrix equation A * X = B
	// std::vector<TMatrixD> As_0nXn_XnXn = MakeMatricesA	(flux_0nXnSum, 	flux_XnXn);
	// std::vector<TMatrixD> Bs_0nXn_XnXn = MakeMatricesB	(Xsec_0nXn,		Xsec_XnXn);
	// std::vector<TMatrixD> Xs_0nXn_XnXn = SovleMatricesEq	(As_0nXn_XnXn, Bs_0nXn_XnXn);

	// cout<< "photonEnergyContrib ------> Processing Xn0n_XnXn Case"<<endl;
	// // // Now Solving matrix equation A * X = B
	// std::vector<TMatrixD> As_Xn0n_XnXn = MakeMatricesA	(flux_0nXnSum, 	flux_XnXn);
	// std::vector<TMatrixD> Bs_Xn0n_XnXn = MakeMatricesB	(Xsec_Xn0n,		Xsec_XnXn);
	// std::vector<TMatrixD> Xs_Xn0n_XnXn = SovleMatricesEq	(As_Xn0n_XnXn, Bs_Xn0n_XnXn);

	// cout<< "photonEnergyContrib ------> Processing 0n0n_XnXn Case"<<endl;
	// std::vector<double> matrix_elements_0n0n_XnXn{flux_0n0n};
	// matrix_elements_0n0n_XnXn.insert(matrix_elements_0n0n_XnXn.end(), flux_XnXn.begin(), flux_XnXn.end());
	// // Now Solving matrix equation A * X = B
	// std::vector<TMatrixD> As_0n0n_XnXn = MakeMatricesA	(flux_0n0n, 	flux_XnXn);
	// std::vector<TMatrixD> Bs_0n0n_XnXn = MakeMatricesB	(Xsec_0n0n,		Xsec_XnXn);
	// std::vector<TMatrixD> Xs_0n0n_XnXn = SovleMatricesEq	(As_0n0n_XnXn, Bs_0n0n_XnXn);

	// cout<< "photonEnergyContrib ------> Processing 0n0n_0nXnSum Case"<<endl;
	// std::vector<double> matrix_elements_0n0n_0nXnSum{flux_0n0n};
	// matrix_elements_0n0n_0nXnSum.insert(matrix_elements_0n0n_0nXnSum.end(), flux_0nXnSum.begin(), flux_0nXnSum.end());
	// // Now Solving matrix equation A * X = B
	// std::vector<TMatrixD> As_0n0n_0nXnSum = MakeMatricesA	(flux_0n0n, 	flux_0nXnSum);
	// std::vector<TMatrixD> Bs_0n0n_0nXnSum = MakeMatricesB	(Xsec_0n0n,		Xsec_0nXnSum);
	// std::vector<TMatrixD> Xs_0n0n_0nXnSum = SovleMatricesEq	(As_0n0n_0nXnSum, Bs_0n0n_0nXnSum);
}
