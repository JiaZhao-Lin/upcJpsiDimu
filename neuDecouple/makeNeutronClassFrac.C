
void makeNeutronClassFrac()
{
	TH2D* hNeutronClassFrac_inZB =  new TH2D("hNeutronClassFrac_inZB", "", 2,0,2, 2,0,2);

	hNeutronClassFrac_inZB->Fill(1,1, 0.0889779);
	hNeutronClassFrac_inZB->Fill(1,2, 0.05084);
	hNeutronClassFrac_inZB->Fill(2,1, 0.05306);
	hNeutronClassFrac_inZB->Fill(2,2, 0.00631);

}
