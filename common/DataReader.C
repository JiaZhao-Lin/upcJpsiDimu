/*
	A data structure for reading Columnized data from files
	Do:
		SingleDataReader(FileDir, ColumnNames)
	JiaZhao Lin
	Oct. 2022
*/

struct ColumnPrototype
{
	std::vector<TString> ColumnNames;

	ColumnPrototype( std::vector<TString> ColumnNames_) : ColumnNames{ColumnNames_} {}
	TString getHeader()
	{
		TString header{""};
		for(auto name: ColumnNames){	header += name + ":";	}
		return header.Remove(header.Length()-1);
	}
	int size()	{return ColumnNames.size();}
};

struct DataReader
{
	virtual void Read() 	= 0;
};

struct SingleDataReader	:	DataReader
{
	/*
		In put the directory of the file and the column name (as a vector of TString)
	*/
	TString FileDir;
	ColumnPrototype cp;
	std::map<TString, std::vector<double>> Map;

	SingleDataReader(TString FileDir_,	ColumnPrototype cp_) 				: FileDir{FileDir_}, cp{cp_} 			{ 	Read(); }
	SingleDataReader(TString FileDir_,	std::vector<TString> ColumnNames_) 	: FileDir{FileDir_}, cp{ColumnNames_} 	{	Read(); }
	void Read() override
	{
		Init_Map();

		TTree *tree = new TTree();
		tree->ReadFile(	FileDir, cp.getHeader().Data() );
		int N = tree->GetEntries();
		std::vector<float> temp(cp.size());
		for (int i = 0; i < cp.size(); ++i){ tree->SetBranchAddress(cp.ColumnNames[i],	&temp[i]); }
		for (int i = 0; i < N; ++i)
		{
			tree->GetEntry(i);
			for (int j = 0; j < cp.size(); ++j)
			{
				Map.at(cp.ColumnNames[j]).push_back(double(temp[j]));
			}
		}
	}
	std::vector<double> GetVec(TString s) const
	{
		if (Map.size() == 0)	std::runtime_error("SingleDataReader: No MAP!?");
		return Map.at(s);
	}
	
private:
	void Init_Map(){ if(Map.size()!=0) cout << Form( "WARNING: Overwriting: %s", FileDir.Data() ) << endl; for (auto name : cp.ColumnNames) { Map[name] = {}; } }
};

struct MultiDataReader : DataReader
{
	std::vector<TString> FileDirs;
	std::vector<TString> ColumnNames;
	std::vector<SingleDataReader> SDRs;

	MultiDataReader(std::vector<TString> FileDirs_, std::vector<TString> ColumnNames_) : FileDirs{FileDirs_}, ColumnNames{ColumnNames_} { Read(); }

	void Read() override { for (auto FileDir : FileDirs) {	SDRs.push_back( {FileDir, ColumnNames} ); }	}
	std::vector<double> GetVec(int i, TString s) const { return SDRs[i].GetVec(s); }
};


void DataReader(){
	//Testing
	SingleDataReader a("../physicsFigures/inputfiles/LTA_Jpsi_weak_shadowing.dat", {"y", "AnAn", "0n0n", "0nXnSum", "XnXn"});
	cout<<a.cp.getHeader().Data()<<endl;
	MultiDataReader b({"../physicsFigures/inputfiles/LTA_Jpsi_weak_shadowing.dat", "../physicsFigures/inputfiles/LTA_Jpsi_strong_shadowing.dat"}, {"y", "AnAn", "0n0n", "0nXnSum", "XnXn"});
	cout<<a.GetVec("y")[2]<<endl;
	cout<<b.GetVec(0,"y")[7]<<endl;
}