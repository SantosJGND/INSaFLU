'''
Created on 03/05/2020

@author: mmp
'''
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
from utils.lock_atomic_transaction import LockedAtomicTransaction

class DefaultSoftware(object):
	'''
	classdocs
	'''
	software_names = SoftwareNames()
	
	def test_all_defaults(self, user):
		### test all defaults
		self.test_default_trimmomatic_db(user)
		
	def test_default_trimmomatic_db(self, user):
		"""
		test if exist, if not persist in database
		"""
		try:
			Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user)
		except Software.DoesNotExist:
			### create a default one for this user
			with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
				vect_parameters = self._get_trimmomatic_default(user)
				self._persist_parameters(vect_parameters)
				
	def get_trimmomatic_parameters(self, user):
		"""
		get trimmomatic parameters
		"""
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user)
		except Software.DoesNotExist:
			return ""

		## get parameters for a specific user
		parameters = Parameter.objects.filter(software=software)
		
		### parse them
		dict_out = {}
		vect_order_ouput = []
		for parameter in parameters:
			### don't set the not set parameters
			if (not parameter.not_set_value is None and parameter.parameter == parameter.not_set_value): continue
			
			### create a dict with parameters
			if (parameter.name in dict_out): 
				dict_out[parameter.name][1].append(parameter.parameter)
				dict_out[parameter.name][0].append(parameter.union_char)
			else:
				vect_order_ouput.append(parameter.name) 
				dict_out[parameter.name] = [[parameter.union_char], [parameter.parameter]]
			
		return_parameter = ""
		for par_name in vect_order_ouput:
			if (len(dict_out[par_name][1]) == 1 and len(dict_out[par_name][1][0]) == 0):
				return_parameter += " {}".format(par_name)
			else:
				return_parameter += " {}".format(par_name)
				for _ in range(len(dict_out[par_name][0])):
					return_parameter += "{}{}".format(dict_out[par_name][0][_], dict_out[par_name][1][_])
		return return_parameter.strip()

	def set_default_software(self, software, user):
		vect_parameters = self._get_trimmomatic_default(user)
		parameters = Parameter.objects.filter(software=software)
		for parameter in parameters:
			if parameter.can_change:
				for parameter_to_set_default in vect_parameters:
					if (parameter_to_set_default.sequence_out == parameter.sequence_out):
						parameter.parameter = parameter_to_set_default.parameter
						parameter.save()
						break
		
	def _persist_parameters(self, vect_parameters):
		"""
		presist a specific software by default 
		"""
		software = None
		dt_out_sequential = {}
		for parameter in vect_parameters:
			assert parameter.sequence_out not in dt_out_sequential
			if software is None:
				software = parameter.software 
				software.save()
			parameter.software = software
			parameter.save()
			
			## set sequential number
			dt_out_sequential[parameter.sequence_out] = 1

	def get_parameters(self, software_name, user):
		"""
		"""
		if (software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name): return self.get_trimmomatic_parameters(user)
		return ""


	def get_all_software(self):
		"""
		get all softwares available by this class
		"""
		vect_software = []
		vect_software.append(self.software_names.get_trimmomatic_name())
		return vect_software


	def _get_trimmomatic_default(self, user):
		"""
		SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_name
		software.version = SoftwareNames.SOFTWARE_TRIMMOMATIC_VERSION
		software.owner = user
		
		vect_parameters =  []
		parameter = Parameter()
		parameter.name = "SLIDINGWINDOW"
		parameter.parameter = "5"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[3:50]"
		parameter.range_max = "50"
		parameter.range_min = "3"
		parameter.description = "SLIDINGWINDOW:<windowSize> specifies the number of bases to average across"
		vect_parameters.append(parameter)
		
		
		parameter = Parameter()
		parameter.name = "SLIDINGWINDOW"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[10:100]"
		parameter.range_max = "100"
		parameter.range_min = "10"
		parameter.description = "SLIDINGWINDOW:<requiredQuality> specifies the average quality required"
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "LEADING"
		parameter.parameter = "3"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "LEADING:<quality> Remove low quality bases from the beginning."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "TRAILING"
		parameter.parameter = "3"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 4
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "TRAILING:<quality> Remove low quality bases from the end."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "MINLEN"
		parameter.parameter = "35"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[5:500]"
		parameter.range_max = "500"
		parameter.range_min = "5"
		parameter.description = "SMINLEN:<length> This module removes reads that fall below the specified minimal length."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "CROP"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 6
		parameter.range_available = "[0:400]"
		parameter.range_max = "400"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "CROP:<length> Cut the read to a specified length."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "HEADCROP"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 7
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "HEADCROP:<length> Cut the specified number of bases from the start of the read."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "AVGQUAL"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 8
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "AVGQUAL:<length> Drop the read if the average quality is below the specified level."
		vect_parameters.append(parameter)
		
		
		parameter = Parameter()
		parameter.name = "TOPHRED33"
		parameter.parameter = ""
		parameter.type_data = Parameter.PARAMETER_null
		parameter.software = software
		parameter.union_char = ""
		parameter.can_change = False
		parameter.sequence_out = 9
		parameter.description = "This (re)encodes the quality part of the FASTQ file to base 33."
		vect_parameters.append(parameter)
		return vect_parameters





