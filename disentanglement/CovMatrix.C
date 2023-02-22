// #include "../common/Matrix.C"
#include "../common/SortByKeys.C"
#include "UncerAnalyzer.C"

class CalculateCovMatrix
{
public:
	TMatrixD cov_matrix;
	TMatrixD corr_matrix;

	AnalysisData *AnaData = nullptr;
	TString name, param;
	bool isSigned = false;	// Fully correlated uncertainty use unsigned covariance
	std::vector<double> sort_keys;

	// Constructor 
	CalculateCovMatrix(const TString name_, TMatrixD matrix_) : name{name_}, cov_matrix{matrix_} {}
	CalculateCovMatrix(const TString name_, AnalysisData *AnaData_, const TString param_) : name{name_}, AnaData{AnaData_}, param{param_} {Init();}
	CalculateCovMatrix(const TString name_, AnalysisData *AnaData_, const TString param_, std::vector<double> sort_keys_) : name{name_}, AnaData{AnaData_}, param{param_}, sort_keys{sort_keys_} {Init();}
	CalculateCovMatrix(const TString name_, AnalysisData *AnaData_, const TString param_, std::vector<double> sort_keys_, bool isSigned_) : name{name_}, AnaData{AnaData_}, param{param_}, sort_keys{sort_keys_},isSigned{isSigned_} {Init();}

	void Init()
	{
		int nrows = AnaData->GetSize(param);
		int ncols = AnaData->GetSize(param);

		cov_matrix.ResizeTo(nrows, ncols);
		corr_matrix.ResizeTo(nrows, ncols);

		cov_matrix = MatrixFactory(nrows, ncols, 1);
		corr_matrix = MatrixFactory(nrows, ncols, 1);
	}

	void Print()	const
	{
		cov_matrix.Print();
		corr_matrix.Print();
	}

	// void SetParam(const TString param_)
	// {
	// 	this->param = param_;
	// 	CalculateCov();
	// }

	// void SetSigned(bool isSigned_)
	// {
	// 	this->isSigned = isSigned_;
	// 	CalculateCov();
	// }

	void CalculateCov(std::vector<double> Norm_v)
	{
		if (AnaData == nullptr)		std::runtime_error("CalculateCovMatrix --> AnaData is not set!");
		if (sort_keys.empty())		std::runtime_error("CalculateCovMatrix --> sort_keys is not set!");

		std::vector<double> vec = AnaData->Get(param);
		if (isSigned) vec = AnaData->Get(param + "_Signed");

		//sort the vector based on the sort_keys
		vec = SortByKeys(vec, sort_keys);

		for (int i = 0; i < cov_matrix.GetNrows(); ++i)
		{
			for (int j = 0; j < cov_matrix.GetNcols(); ++j)
			{
				double cov = (vec[i]/100.00)*Norm_v[i] * (vec[j]/100.00) * Norm_v[j];
				cov_matrix(i, j) = cov;
			}
		}
		cout << "Covariance matrix from " << name << " is calculated for " << param << endl;
		Print();

		//update the correlation matrix
		corr_matrix = CovToCorr(cov_matrix);
	}

	void SetSortKeys(std::vector<double> sort_keys_)
	{
		this->sort_keys = sort_keys_;
	}

	static TMatrixD CovToCorr(const TMatrixD& cov)
	{
		TMatrixD corr = cov;
		int n = cov.GetNrows();

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				corr(i,j) = cov(i,j) / (TMath::Sqrt(cov(i,i)) * TMath::Sqrt(cov(j,j)));
			}
		}

		corr.Print();

		return corr;
	}

	//create a matrix with all elements equal to value
	static TMatrixD MatrixFactory(const int nrows, const int ncols, const int value)
	{
		TMatrixD matrix(nrows, ncols);
		for (int i = 0; i < nrows; ++i)
		{
			for (int j = 0; j < ncols; ++j)
			{
				matrix(i, j) = value;
			}
		}
		return matrix;
	}

	// Method to read a matrix from a text file
	static TMatrixD Read(const TString filename)
	{
		std::ifstream infile(filename.Data());
		if (infile.is_open())
		{
			// Read the dimensions of the matrix from the first line of the file
			int nrows, ncols;
			infile >> nrows >> ncols;

			// Create a new matrix with the given dimensions
			TMatrixD result = CalculateCovMatrix::MatrixFactory(nrows, ncols, 0);

			// Read the elements of the matrix from the file
			for (int i = 0; i < nrows; ++i)
			{
				for (int j = 0; j < ncols; ++j)
				{
					double d;
					infile >> d;
					result(i, j) = d;
				}
			}

			infile.close();
			return result;
		}
		else
		{
			throw std::runtime_error( "Error: could not open file " + filename + " for reading.");
		}
	}

	// Method to write the matrix to a text file
	static void WriteToFile(const TMatrixD matrix, const TString filename)
	{
		std::ofstream outfile(filename.Data());
		if (outfile.is_open())
		{
			outfile << matrix.GetNrows() << " " << matrix.GetNcols() << std::endl;
			for (int i = 0; i < matrix.GetNrows(); ++i)
			{
				for (int j = 0; j < matrix.GetNcols(); ++j)
				{
					outfile << matrix(i, j) << " ";
				}
				outfile << std::endl;
			}
			outfile.close();
		} 
		else
		{
			throw std::runtime_error( "Error: could not open file " + filename + " for writing.");
		}
	}

	static void DrawCov(const TMatrixD &matrix)
	{
		TCanvas* c = new TCanvas("c", "c", 1000, 1000);
		c->SetRightMargin(0.15);

		TH2D* hist = new TH2D("hist", "Heat Map", matrix.GetNrows(), 0, matrix.GetNrows(), matrix.GetNcols(), 0, matrix.GetNcols());
		for (int i = 0; i < matrix.GetNrows(); ++i)
		{
			for (int j = 0; j < matrix.GetNcols(); ++j)
			{
				hist->SetBinContent((i+1), (j+1), matrix[i][j]);
			}
		}
		hist->Draw("colz");

		c->SaveAs("outFigures/CovMatrix.pdf");
		
		delete c;
		delete hist;
	}	

	static void DrawCorr(const TMatrixD &matrix)
	{
		TCanvas* c = new TCanvas("c", "c", 1000, 1000);
		c->SetRightMargin(0.15);

		TH2D* hist = new TH2D("hist", "Heat Map", matrix.GetNrows(), 0, matrix.GetNrows(), matrix.GetNcols(), 0, matrix.GetNcols());
		for (int i = 0; i < matrix.GetNrows(); ++i)
		{
			for (int j = 0; j < matrix.GetNcols(); ++j)
			{
				hist->SetBinContent((i+1), (j+1), matrix[i][j]);
			}
		}
		hist->GetZaxis()->SetRangeUser(-1.00001, 1.00001);
		hist->Draw("colz");
		c->SaveAs("outFigures/CorrMatrix.pdf");

		delete c;
		delete hist;
	}
};

void CovMatrix()
{
	std::vector<TString> constCorrUncer = {"Lumi", "BR"};
	std::vector<TString> fullCorrUncer = {"Signal Ext", "HFveto", "TnP"};
	std::vector<TString> partCorrUncer = {"PhotonFlux", "n-PileUp"};

	std::vector<CalculateCovMatrix> fullCovMatrices;
	std::vector<CalculateCovMatrix> partCovMatrices;

	struct UncerAnalyzer UncerAna;
	UncerAna.Handle();

	auto sort_keys = UncerAna.Default_AnaData.Get("W");

	//calculate the covariance matrix for each uncertainty
	//constCorrUncer is in the vector of fullly correlated matrix
	for (auto& name : constCorrUncer)
	{
		class CalculateCovMatrix cov(name, &UncerAna.AnaData_Breakdown.at(name), name + "_SysUncer", sort_keys);
		cov.CalculateCov(UncerAna.Default_AnaData.Get("Sigma"));
		fullCovMatrices.push_back(cov);
	}

	for (auto& name : fullCorrUncer)
	{
		class CalculateCovMatrix cov(name, &UncerAna.AnaData_Breakdown.at(name), "Sigma_SysUncer", sort_keys);
		cov.CalculateCov(UncerAna.Default_AnaData.Get("Sigma"));
		fullCovMatrices.push_back(cov);
	}

	for (auto& name : partCorrUncer)
	{
		class CalculateCovMatrix cov(name, &UncerAna.AnaData_Breakdown.at(name), "Sigma_SysUncer", sort_keys, true);
		cov.CalculateCov(UncerAna.Default_AnaData.Get("Sigma"));
		partCovMatrices.push_back(cov);
	}

	std::vector<CalculateCovMatrix> CovMatrices = fullCovMatrices;
	CovMatrices.insert(CovMatrices.end(), partCovMatrices.begin(), partCovMatrices.end());

	//sum the covariance matrix
	TMatrixD cov_total = CalculateCovMatrix::MatrixFactory(CovMatrices[0].cov_matrix.GetNrows(), CovMatrices[0].cov_matrix.GetNcols(), 0.0);
	for (auto& cov : CovMatrices)
	{
		cov_total += cov.cov_matrix;
	}

	cov_total.Print();
	auto corr_total = CalculateCovMatrix::CovToCorr(cov_total);

	auto cov_sigEx = fullCovMatrices[2].cov_matrix;
	auto corr_sigEx = CalculateCovMatrix::CovToCorr(cov_sigEx);
	auto cov_HFveto = fullCovMatrices[3].cov_matrix;
	auto corr_HFveto = CalculateCovMatrix::CovToCorr(cov_HFveto);
	auto cov_TnP = fullCovMatrices[4].cov_matrix;
	auto corr_TnP = CalculateCovMatrix::CovToCorr(cov_TnP);

	auto cov_flux = partCovMatrices[0].cov_matrix;
	auto corr_flux = CalculateCovMatrix::CovToCorr(cov_flux);
	auto cov_neutron = partCovMatrices[1].cov_matrix;
	auto corr_neutron = CalculateCovMatrix::CovToCorr(cov_neutron);

	// CalculateCovMatrix::DrawCorr(corr_sigEx);
	// CalculateCovMatrix::DrawCov(cov_sigEx);
	// CalculateCovMatrix::DrawCorr(corr_HFveto);
	// CalculateCovMatrix::DrawCov(cov_HFveto);
	// CalculateCovMatrix::DrawCorr(corr_TnP);
	// CalculateCovMatrix::DrawCov(cov_TnP);

	CalculateCovMatrix::DrawCorr(corr_total);
	CalculateCovMatrix::DrawCov(cov_total);
	// CalculateCovMatrix::DrawCorr(corr_flux);
	// CalculateCovMatrix::DrawCov(cov_flux);
	// CalculateCovMatrix::DrawCorr(corr_neutron);
	// CalculateCovMatrix::DrawCov(cov_neutron);

	// auto flux_cov = partCovMatrices[0].cov_matrix;
	// auto flux_corr = CalculateCovMatrix::CovToCorr(flux_cov);
	// CalculateCovMatrix::DrawCorr(flux_corr);
	// CalculateCovMatrix::DrawCov(flux_cov);
	// flux_corr.Print();

}