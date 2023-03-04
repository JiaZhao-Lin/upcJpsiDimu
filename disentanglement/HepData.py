# Generate HepData for sumbission

from hepdata_lib import Submission, Table, Variable, Uncertainty
import ROOT

# a class to create a submission

class HepDataCreator:
    
	def __init__(self, name, infile, infile_sys):
		self.name = name
		self.infile = infile
		self.infile_sys = infile_sys
		self.data_dict,	self.data_sys_dict = self.read()

	def read(self):
        # Open the ROOT file and get the TTree
		root_file = ROOT.TFile(self.infile, "READ")
		keys = root_file.Get("keys")
		keys = [str(key) for key in keys]
		data_dict = {	key : list(root_file.Get(key)) for key in keys	}

		root_file_sys = ROOT.TFile(self.infile_sys, "READ")
		keys_sys = root_file_sys.Get("keys")
		keys_sys = [str(key) for key in keys_sys]
		data_dict_sys = {	key : list(root_file_sys.Get(key)) for key in keys_sys	}

		# sort_keys = list(data_dict["W"])
		# data_dict_sorted = dict()

		#sort the data_dict_sorted by the values of W
		# for [key, value] in data_dict.items():
		# 	if len(value) != 6:
		# 		data_dict_sorted[key] = value
		# 		continue
		# 	if key == "W":
		# 		data_dict_sorted[key] = sorted(value)
		# 		continue
		# 	zipped = zip(value,	sort_keys)
		# 	data_dict_sorted[key] = [x for x,y in sorted(zipped, key=lambda x: x[1])]

		# print(data_dict_sorted["Sigma"])
		# print(data_dict_sorted["W"])
		# print(data_dict_sorted["Dy"])
		# print(data_dict_sorted["DSigmaDy_AnAn"])
		# print(data_dict_sorted)

		return data_dict, data_dict_sys

	def create_variables(self, param, independent=False, binned=False, units="", sys_breakdown=False):
		var = Variable(param, is_independent=independent, is_binned=binned, units=units)
		var.values = self.data_dict[param]
		
		unc = Uncertainty(param+"_stat_unc", is_symmetric=True)
		unc.values = self.data_dict[param+"_Err"]
		
		var.add_uncertainty(unc)
		if sys_breakdown:
			sys_thoery_unc = Uncertainty(param+"_theory_unc", is_symmetric=True)
			sys_thoery_unc.values = self.data_sys_dict[param+"_TheorySysErr"]
			var.add_uncertainty(sys_thoery_unc)

			sys_Experiment_unc = Uncertainty(param+"_experiment_unc", is_symmetric=True)
			sys_Experiment_unc.values = self.data_sys_dict[param+"_ExperiSysErr"]
			var.add_uncertainty(sys_Experiment_unc)

		else:
			sys_unc = Uncertainty(param+"_sys_unc", is_symmetric=True)
			sys_unc.values = self.data_sys_dict[param+"_SysErr"]
			var.add_uncertainty(sys_unc)

		return var


	def create_table(self):

		Dy_var = Variable("Dy", is_independent=True, is_binned=False, units="")
		Dy_var.values = self.data_dict["Dy"]
		# Dy_unc = Uncertainty("Dy_stat_unc", is_symmetric=True)
		# Dy_unc.values = self.data_dict["Dy_Err"]
		# print(self.data_dict["Dy_Err"])
		# Dy_var.add_uncertainty(Dy_unc)

		DSigmaDy_AnAn_var = self.create_variables("DSigmaDy_AnAn", units="mb")
		DSigmaDy_0n0n_var = self.create_variables("DSigmaDy_0n0n", units="mb")
		DSigmaDy_0nXnSum_var = self.create_variables("DSigmaDy_0nXnSum", units="mb")
		DSigmaDy_XnXn_var = self.create_variables("DSigmaDy_XnXn", units="mb")

		DSigmaDy_table = Table("DSigmaDy")
		DSigmaDy_table.add_variable(Dy_var)
		DSigmaDy_table.add_variable(DSigmaDy_AnAn_var)
		DSigmaDy_table.add_variable(DSigmaDy_0n0n_var)
		DSigmaDy_table.add_variable(DSigmaDy_0nXnSum_var)
		DSigmaDy_table.add_variable(DSigmaDy_XnXn_var)


		W_var = Variable("W", is_independent=True, is_binned=False, units="GeV")
		W_var.values = self.data_dict["W"]

		Sigma_var = self.create_variables("Sigma", units="mb", sys_breakdown=True)

		Sigma_table = Table("Sigma")
		Sigma_table.add_variable(W_var)
		Sigma_table.add_variable(Sigma_var)


		X_var = Variable("X", is_independent=True, is_binned=False, units="")
		X_var.values = self.data_dict["X"]

		R_var = self.create_variables("R", sys_breakdown=True)

		R_table = Table("R")
		R_table.add_variable(X_var)
		R_table.add_variable(R_var)

		return [DSigmaDy_table, Sigma_table, R_table]


class HepDataCreatorMatrix:
	def __init__(self, name, infile):
		self.name = name
		self.infile = infile
		self.matrix = self.read()
		self.nrows = len(self.matrix)
		self.ncols = len(self.matrix[0])

	def read(self):
		#read a matrix from a text file
		#the first line is the header, showing the number of rows and columns
		matrix = []
		with open(self.infile, "r") as f:
			header = f.readline()
			header = header.split()
			nrows = int(header[0])
			ncols = int(header[1])

			#read the data
			for i in range(nrows):
				line = f.readline()
				line = line.split()
				line = [float(x) for x in line]
				matrix.append(list(line))
			
		return matrix

	
	def create_table(self, table_name="Total covariance matrix"):
		# Create the table object and add the variables
		x = []
		y = []
		for i in range(self.nrows):
			for j in range(self.ncols):
				x.append(i+1)
				y.append(j+1)

		#flatten the matrix
		z = [j for sub in self.matrix for j in sub]

		x_var = Variable("x", is_independent=True, is_binned=False, units="")
		x_var.values = x

		y_var = Variable("y", is_independent=True, is_binned=False, units="")
		y_var.values = y

		z_var = Variable("z", is_independent=False, is_binned=False, units="")
		z_var.values = z

		table = Table(table_name)
		for var in [x_var,y_var,z_var]:
			table.add_variable(var)

		return [table]

def submit(tables_):
	# Submit the table to the HEPData database
	tables = list(tables_)
	submission = Submission()
	for table in tables:
		submission.add_table(table)
	
	submission.create_files("./outFiles/HepData", remove_old=True)

if (__name__ == "__main__"):
	print("HEPData.py")
	hep = HepDataCreator("test", "./outFiles/Result_CMS.root", 	"./outFiles/Result_CMS_SysUncer.root")
	hep_matrix = HepDataCreatorMatrix("test", "./outFiles/CovMatrix.txt")
	hep_matrix_experi = HepDataCreatorMatrix("test", "./outFiles/CovMatrixExperi.txt")
	hep_matrix_flux = HepDataCreatorMatrix("test", "./outFiles/CovMatrixFlux.txt")

	data_tables = hep.create_table()
	matrix_table = hep_matrix.create_table()
	matrix_table_experi = hep_matrix_experi.create_table("Experimental covariance matrix")
	matrix_table_flux = hep_matrix_flux.create_table("Flux covariance matrix")

	tables = list(data_tables + matrix_table + matrix_table_experi + matrix_table_flux)
	submit(tables)
