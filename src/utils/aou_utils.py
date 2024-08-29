"""
Utilities for running workflows
"""

import os
import shutil
import subprocess
import sys
import tempfile

def ZipWDL(wdl_dependencies_file, file_list):
	"""
	Put all WDL dependencies into a zip file

	Arguments
	---------
	wdl_dependencies_fie : str
	    Zip file to put other wdls in
	"""
	dirname = tempfile.mkdtemp()
	for f in file_list:
		shutil.copyfile(f, dirname+"/"+os.path.basename(f))
	shutil.make_archive(os.path.splitext(wdl_dependencies_file)[0], "zip", root_dir=dirname)

def RunWorkflow(wdl_file, json_file, json_options_file, wdl_dependencies_file="", dryrun=False):
	"""
	Run workflow on AoU

	Arguments
	---------
	wdl_file : str
	    Main WDL file to run
	json_file : str
	    JSON file path with input arguments
	json_options_file : str
	    JSON with additional options for cromshell
	wdl_dependencies_file : str (optional)
	    Zip file with other WDLs that are imported
	dryrun : bool
	    Just print the command, don't actually run cromshell
	"""
	cmd = "cromshell submit {wdl} {json} -op {options}".format(wdl=wdl_file, json=json_file, \
		options=json_options_file)
	if wdl_dependencies_file.strip() != "":
		cmd += " -d {otherwdl}".format(otherwdl=wdl_dependencies_file)
	if dryrun:
		sys.stderr.write("Run: %s\n"%cmd)
		return
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))

def UploadGS(local_path, gcp_path):
	"""
	Upload a local file to GCP

	Arguments
	---------
	local_path : str
	   Local path
	gcp_path : str
	   GCP path to upload to
	"""
	cmd = "gsutil cp {src} {dest}".format(src=local_path, dest=gcp_path)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))