void saveMap(std::map<TString, std::vector<double>> map, TString fileName = "rootfiles/Results_Map.root")
{
	cout << Form("saveMap-------->Saving Map To The RootFile------->%s",	fileName.Data()) << endl;
	TFile *file = TFile::Open(fileName.Data(), "RECREATE");
	std::vector<TString> keys;

	for (auto it = map.begin(); it != map.end(); ++it)
	{
		TString key 	=	it->first;
		std::cout << "Saving the key:	" << key << endl;
		keys.push_back(key);
  		std::vector<double> 	value	=	it->second;
		file->WriteObject(&value, key.Data()); // I store the vector in the TFile
	}
	file->WriteObject(&keys, "keys");
	file->Close();
	cout << Form("saveMap-------->DONE SAVING The RootFile------->%s",	fileName.Data()) << endl;
}

std::map<TString, std::vector<double>> readMap(TString fileName = "rootfiles/Results_Map.root")
{
	cout << Form("readMap-------->Reading Map From The RootFile------->%s",	fileName.Data()) << endl;
	TFile *file = TFile::Open(fileName, "READ");
	std::map<TString, std::vector<double>> map;

	std::vector<TString> *keys;
	file->GetObject("keys", keys); // I try to retrieve the vector
	for(auto it = keys->begin(); it != keys->end(); ++it) {
		std::vector<Double_t> *temp;

		file->GetObject(*it, temp);
		std::cout << "Retrieving the key:	" << *it << endl;
		map[*it] = *temp;
	}
	cout << Form("readMap-------->DONE READING The RootFile------->%s",	fileName.Data()) << endl;
	return map;
}