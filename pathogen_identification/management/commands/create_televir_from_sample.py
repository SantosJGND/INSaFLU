'''
Created on January 30, 2023

@author: daniel.sobral
'''
import os
import logging
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from django.db import transaction
from managing_files.models import Sample
from pathogen_identification.models import Projects, PIProject_Sample
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.constants_settings import ConstantsSettings
from utils.process_SGE import ProcessSGE
from constants.constants import Constants

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Create a TELEVIR project and add a sample to it. Returns project id."

	# logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)

	def add_arguments(self, parser):
		parser.add_argument('--project_name', nargs='?', type=str,
							required=True, help='Project Name')
		parser.add_argument('--sample_name', nargs='?', type=str,
							required=True, help='Sample Name')
		parser.add_argument('--user_login', nargs='?', type=str, required=True,
							help='User login of the project and sample owner')		
					

	# A command must define handle()
	def handle(self, *args, **options):

		project_name = options['project_name']
		sample_name = options['sample_name']
		account = options['user_login']

		try:
			
			process_SGE = ProcessSGE()
			user = User.objects.get(username=account)
			sample = Sample.objects.get(name=sample_name, owner=user)

			project_count = Projects.objects.filter(
				name__iexact=project_name,
                is_deleted=False,
                owner__username=user.username,
            ).count()
			
			if(project_count > 0):

				self.stdout.write("Project '{}' already exists, reusing...".format(project_name))

				project = Projects.objects.filter(
				    name__iexact=project_name,
                    is_deleted=False,
                    owner__username=user.username,
                )[0]

				with transaction.atomic():
					project_sample = PIProject_Sample()
					project_sample.project = project
					project_sample.sample = sample
					project_sample.name = sample.name
					project_sample_input = sample.file_name_1
					if sample.is_valid_2:
						project_sample_input += ";" + sample.file_name_2                    
					project_sample.input = project_sample_input
					sample_technology = "ONT"
					if(sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina):
						sample_technology = "Illumina/IonTorrent"			
					if(project.technology != sample_technology):
						self.stdout.write("Project has different technology {} from sample technology {}...".format(project.technology, sample_technology))
					project_sample.technology = sample.type_of_fastq
					project_sample.report = "report"
					project_sample.save()

				utils = Utils_Manager()
				runs_to_deploy = utils.check_runs_to_deploy(user, project)

				if runs_to_deploy:
					taskID = process_SGE.set_submit_televir_job(
                    	user=user,
                    	project_pk=project.pk,
					)
					self.stdout.write("Project submitted as task {}.".format(project.id, taskID))   

			else:

				sample = Sample.objects.get(name=sample_name, owner=user)

				with transaction.atomic():
					project = Projects()
					project.name = project_name
					project.owner = user
					project.owner_id = user.id
					# TODO Check where these constants are, or define them somewhere...
					technology = "ONT"
					if(sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina):
						technology = "Illumina/IonTorrent"
					project.technology = technology
					project.save()
					project_sample_input = sample.file_name_1
					if sample.is_valid_2:
						project_sample_input += ";" + sample.file_name_2
					project_sample = PIProject_Sample()
					project_sample.project = project
					project_sample.sample = sample
					project_sample.name = sample.name
					project_sample.input = project_sample_input
					project_sample.technology = sample.type_of_fastq
					project_sample.report = "report"
					project_sample.save()
					self.stdout.write("Project created with id {}.".format(project.id))   

				utils = Utils_Manager()
				runs_to_deploy = utils.check_runs_to_deploy(user, project)

				if runs_to_deploy:
					taskID = process_SGE.set_submit_televir_job(
                    	user=user,
                    	project_pk=project.pk,
                	)
					self.stdout.write("Project submitted as task {}.".format(project.id, taskID))   

		except User.DoesNotExist as e:
			self.stdout.write("Error: User '{}' does not exist.".format(account))                
		except Exception as e:
			self.stdout.write("Error: {}.".format(e))
           
