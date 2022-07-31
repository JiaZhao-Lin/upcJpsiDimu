
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
	virtual void Read() = 0;
};

struct SingleDataReader	:	DataReader
{
	TString FileName;
	ColumnPrototype cp;
	std::map<TString, std::vector<double>> Map;

	SingleDataReader(TString FileName_,	ColumnPrototype cp_) : FileName{FileName_}, cp{cp_} { Init_Map(); }
	SingleDataReader(TString FileName_,	std::vector<TString> ColumnNames_) : FileName{FileName_}, cp{ColumnNames_} { Init_Map(); }
	void Read() override
	{
		TTree *tree = new TTree();
		tree->ReadFile(	FileName, cp.getHeader().Data() );
		int N = tree->GetEntries();
		std::vector<float> temp(cp.size());
		for (int i = 0; i < cp.size(); ++i){ tree->SetBranchAddress(cp.ColumnNames[i],	&temp[i]); }
		for (int i = 0; i < N; ++i)
		{
			tree->GetEntry(i);
			for (int j = 0; j < cp.size(); ++j)
			{
				Map.at(cp.ColumnNames[j]).push_back(double(temp[i]));
			}
		}
	}
private:
	void Init_Map(){	for(auto name : cp.ColumnNames) { Map[name] = {}; }	}
};

void DataReader(){
	SingleDataReader a("TEST",	{"a","b","c"});
	cout<<a.cp.getHeader().Data()<<endl;
}